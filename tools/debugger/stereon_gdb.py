# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""GDB pretty-printers for Stereon (ADR-0003, ADR-0004).

Load them in a session or in ~/.gdbinit:

    source /path/to/Stereon/tools/debugger/stereon_gdb.py

Values print as Stereon's own std::format output does, so a debugger
session, a log line and a test failure read alike:

    vector          (1, 2.5, -3)
    interval        [1, 2]          (empty: [])
    box             [(0, 0, 0), (1, 2, 3)]
    vertex          Vertex#12       (invalid: Vertex#-)
    builder handle  Face#0@0
    coedge          Coedge#42 (reversed)
    curve           Curve(Circle) {Center = (1, 2, 3), Radius = 2.5, ...}
    arena           PagedArena of 6 entries in 2 pages, 1 shared

Two details of the GCC 16 module build shape the code. Type names carry
the module they are attached to ("Stereon::Interval@Stereon.Core"), which
the matchers strip. And a type-erased Curve or Surface holds its value in
a model whose dynamic type GDB cannot resolve on its own, because the
vtable names template arguments as "3ull" where the debug information says
"3"; the printer reads the vtable symbol, takes the wrapped type's name
from it and looks that type up directly.

The formatting rules live in this file only; tools/debugger/stereon_lldb.py
mirrors them for LLDB, and tests/debugger/check_printers.py runs both
against the same expected output.
"""

from __future__ import annotations

import math
import re

import gdb

MODULE_SUFFIX = re.compile(r"@[A-Za-z_][\w.]*")
ENTITY_KINDS = ["Vertex", "Edge", "Coedge", "Loop", "Face", "Shell", "Solid"]
INVALID_INDEX = 0xFFFFFFFF
SENSE_BIT = 0x80000000


def plain_name(gdb_type) -> str:
    """The type's name without typedefs, const or module decorations."""
    t = gdb_type.strip_typedefs().unqualified()
    return MODULE_SUFFIX.sub("", str(t))


def number(x: float) -> str:
    """A double as std::format("{}") writes it: shortest round-trip form."""
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if x == int(x) and abs(x) < 1e16:
        return str(int(x))
    return repr(x)


def tuple_text(value, axes) -> str:
    return "(" + ", ".join(number(float(value[a])) for a in axes) + ")"


def interval_text(value) -> str:
    lower = float(value["m_Lower"])
    upper = float(value["m_Upper"])
    if lower > upper:
        return "[]"
    return f"[{number(lower)}, {number(upper)}]"


def box_text(value, axes) -> str:
    intervals = [value[a] for a in axes]
    if any(float(i["m_Lower"]) > float(i["m_Upper"]) for i in intervals):
        return "[]"
    low = "(" + ", ".join(number(float(i["m_Lower"])) for i in intervals) + ")"
    high = "(" + ", ".join(number(float(i["m_Upper"])) for i in intervals) + ")"
    return f"[{low}, {high}]"


def kind_of(name: str) -> str:
    match = re.search(r"\(Stereon::EntityKind\)(\d+)", name)
    if match and int(match.group(1)) < len(ENTITY_KINDS):
        return ENTITY_KINDS[int(match.group(1))]
    match = re.search(r"Stereon::EntityKind::(\w+)", name)
    return match.group(1) if match else "Entity"


def handle_text(kind: str, index: int) -> str:
    return f"{kind}#-" if index == INVALID_INDEX else f"{kind}#{index}"


def short_geometry_name(full: str) -> str:
    """Line, Circle2d, Cylinder, ... for built-in types; the name otherwise."""
    name = MODULE_SUFFIX.sub("", full)
    builtin = re.fullmatch(r"Stereon::Basic(Line|Circle|Ellipse|BSplineCurve|Hyperbola|Parabola)<double, (\d)[uUlL]*>", name)
    if builtin:
        return builtin.group(1) + ("2d" if builtin.group(2) == "2" else "")
    builtin = re.fullmatch(r"Stereon::Basic(Plane|Cylinder|Cone|Sphere|Torus)<double>", name)
    if builtin:
        return builtin.group(1)
    return name


def value_text(value) -> str:
    """Any value: a Stereon summary when one applies, else a field list."""
    text = summary(value)
    if text is not None:
        return text
    code = value.type.strip_typedefs().code
    if code == gdb.TYPE_CODE_FLT:
        return number(float(value))
    if code == gdb.TYPE_CODE_BOOL:
        return "true" if bool(value) else "false"
    if code == gdb.TYPE_CODE_STRUCT:
        return struct_text(value)
    return str(value)


def struct_text(value) -> str:
    fields = [f for f in value.type.strip_typedefs().fields() if not f.is_base_class and not f.artificial]
    return "{" + ", ".join(f"{f.name} = {value_text(value[f.name])}" for f in fields) + "}"


# --- Type-erased wrappers --------------------------------------------------

VOID_PP = None


def split_template_arguments(text: str) -> list[str]:
    """The top-level template arguments of "Name<a, b<c, d>, e>"."""
    start = text.index("<")
    depth, current, parts = 0, "", []
    for ch in text[start + 1 : -1]:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(current.strip())
            current = ""
        else:
            current += ch
    parts.append(current.strip())
    return parts


def without_integer_suffixes(text: str) -> str:
    """3ull or 3ul (Windows and Linux spell size_t differently) becomes 3."""
    return re.sub(r"\b(\d+)[uUlL]+\b", r"\1", text)


FOUND_TYPES = {}


def find_type(name: str):
    """The gdb.Type named name, ignoring module decorations.

    A demangled vtable may name a type without the module it is attached to
    (Linux) or with it (MinGW), while the debug information always has it, so
    a direct lookup is tried first and then the type list is searched.
    """
    if name in FOUND_TYPES:
        return FOUND_TYPES[name]
    found = None
    try:
        found = gdb.lookup_type(name)
    except gdb.error:
        base = re.match(r"^([^<]*)", name).group(1).split("::")[-1]
        wanted = MODULE_SUFFIX.sub("", name).replace(" ", "")
        listing = gdb.execute(f"info types {base}", to_string=True)
        for line in listing.splitlines():
            candidate = re.sub(r"^\s*\d*:?\s*(?:struct|class|union|typedef)?\s*", "", line)
            candidate = candidate.rstrip(";").strip()
            if MODULE_SUFFIX.sub("", candidate).replace(" ", "") == wanted:
                try:
                    found = gdb.lookup_type(candidate)
                    break
                except gdb.error:
                    continue
    FOUND_TYPES[name] = found
    return found


def wrapped_value(model_pointer):
    """(type name, gdb value or None) of what a Curve or Surface model holds."""
    global VOID_PP
    if VOID_PP is None:
        VOID_PP = gdb.lookup_type("void").pointer().pointer()
    vptr = int(model_pointer.cast(VOID_PP).dereference())
    symbol = gdb.execute(f"info symbol {vptr}", to_string=True)
    match = re.search(r"vtable for (.*?) \+ \d+", symbol)
    if not match:
        return "?", None
    model = match.group(1)
    arguments = split_template_arguments(model)
    # CurveModel<D, T, Size> or SurfaceModel<T, Size>.
    type_name = arguments[-2]
    value_type = find_type(without_integer_suffixes(type_name))
    if value_type is None:
        return type_name, None
    size = int(without_integer_suffixes(arguments[-1]))
    pointer_size = VOID_PP.sizeof
    offset = max(pointer_size, value_type.alignof)
    address = int(model_pointer) + offset
    if value_type.sizeof + pointer_size <= size and value_type.alignof <= 16:
        return type_name, gdb.Value(address).cast(value_type.pointer()).dereference()
    # SharedValue: a shared_ptr<const T>, whose first word is the element.
    element = gdb.Value(address).cast(VOID_PP).dereference()
    return type_name + ", shared", gdb.Value(int(element)).cast(value_type.pointer()).dereference()


def erased_text(label: str, value) -> str:
    model = value["m_Storage"]["m_Model"]
    if int(model) == 0:
        return f"{label} (empty)"
    type_name, held = wrapped_value(model)
    shared = type_name.endswith(", shared")
    short = short_geometry_name(type_name.removesuffix(", shared"))
    head = f"{label}({short}{', shared' if shared else ''})"
    return head if held is None else f"{head} {value_text(held)}"


# --- Containers --------------------------------------------------------------


def vector_elements(vector):
    impl = vector["_M_impl"]
    start, finish = impl["_M_start"], impl["_M_finish"]
    return [(start + i).dereference() for i in range(int(finish - start))]


def arena_text(value) -> str:
    pages = vector_elements(value["m_Pages"])
    size = int(value["m_Size"])
    shared = 0
    for page in pages:
        counts = page["_M_refcount"]["_M_pi"]
        if int(counts) != 0 and int(counts.dereference()["_M_use_count"]) > 1:
            shared += 1
    editing = ", editing" if bool(value["m_Editing"]) else ""
    return f"PagedArena of {size} entries in {len(pages)} pages, {shared} shared{editing}"


def arena_entries(value):
    """("[i]", entry) for every live entry, read page by page."""
    pages = vector_elements(value["m_Pages"])
    if not pages:
        return []
    size = int(value["m_Size"])
    # The page is a std::array<T, N>; N is the extent of its element array.
    elements_type = pages[0]["_M_ptr"].dereference()["_M_elems"].type.strip_typedefs()
    per_page = elements_type.range()[1] + 1
    entries = []
    for index in range(size):
        page = pages[index // per_page]["_M_ptr"].dereference()["_M_elems"]
        entries.append((f"[{index}]", page[index % per_page]))
    return entries


# --- Dispatch ----------------------------------------------------------------

MATCHERS = [
    (re.compile(r"^Stereon::Basic(Vec|Point)2<.*>$"), lambda v: tuple_text(v, ["X", "Y"])),
    (re.compile(r"^Stereon::Basic(Vec|Point)3<.*>$"), lambda v: tuple_text(v, ["X", "Y", "Z"])),
    (re.compile(r"^Stereon::Interval$"), interval_text),
    (re.compile(r"^Stereon::Box2$"), lambda v: box_text(v, ["X", "Y"])),
    (re.compile(r"^Stereon::Box3$"), lambda v: box_text(v, ["X", "Y", "Z"])),
    (re.compile(r"^Stereon::Id<.*>$"),
     lambda v: handle_text(kind_of(str(v.type.strip_typedefs())), int(v["m_Index"]))),
    (re.compile(r"^Stereon::BuilderId<.*>$"), lambda v: (
        handle_text(kind_of(str(v.type.strip_typedefs())), int(v["m_Index"]))
        + ("" if int(v["m_Index"]) == INVALID_INDEX else f"@{int(v['m_Generation'])}"))),
    (re.compile(r"^Stereon::Oriented<.*>$"), lambda v: oriented_text(v)),
    (re.compile(r"^Stereon::BasicErasedCurve<3>$"), lambda v: erased_text("Curve", v)),
    (re.compile(r"^Stereon::BasicErasedCurve<2>$"), lambda v: erased_text("Curve2d", v)),
    (re.compile(r"^Stereon::Surface$"), lambda v: erased_text("Surface", v)),
    (re.compile(r"^Stereon::PagedArena<.*>$"), arena_text),
    (re.compile(r"^Stereon::JournalId$"),
     lambda v: "J#-" if int(v["Value"]) == 0 else f"J#{int(v['Value'])}"),
    (re.compile(r"^Stereon::Sign$"), lambda v: {-1: "Negative", 0: "Zero", 1: "Positive"}.get(int(v), str(int(v)))),
    (re.compile(r"^Stereon::Expansion<.*>$"), lambda v: expansion_text(v)),
]


def oriented_text(value) -> str:
    kind = kind_of(str(value.type.strip_typedefs()))
    bits = int(value["m_Bits"])
    if bits == INVALID_INDEX:
        return f"{kind}#-"
    text = f"{kind}#{bits & ~SENSE_BIT}"
    return text + " (reversed)" if bits & SENSE_BIT else text


def expansion_text(value) -> str:
    length = int(value["m_Length"])
    components = value["m_Components"]["_M_elems"]
    return "Expansion [" + ", ".join(number(float(components[i])) for i in range(length)) + "]"


def summary(value):
    """The Stereon summary of value, or None for other types."""
    try:
        name = plain_name(value.type)
    except gdb.error:
        return None
    for pattern, render in MATCHERS:
        if pattern.match(name):
            return render(value)
    return None


class SummaryPrinter:
    def __init__(self, value, text):
        self.value = value
        self.text = text

    def to_string(self):
        return self.text


class ArenaPrinter(SummaryPrinter):
    def children(self):
        return arena_entries(self.value)

    def display_hint(self):
        return "array"


def lookup(value):
    try:
        name = plain_name(value.type)
    except gdb.error:
        return None
    if not name.startswith("Stereon::"):
        return None
    text = summary(value)
    if text is None:
        return None
    if name.startswith("Stereon::PagedArena<"):
        return ArenaPrinter(value, text)
    return SummaryPrinter(value, text)


def register(objfile=None):
    target = objfile if objfile is not None else gdb
    if lookup not in target.pretty_printers:
        target.pretty_printers.append(lookup)


register()
