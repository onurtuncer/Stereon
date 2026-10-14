# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""LLDB summaries for Stereon (ADR-0003, ADR-0004).

Load them in a session or in ~/.lldbinit:

    command script import /path/to/Stereon/tools/debugger/stereon_lldb.py

The summaries are the ones tools/debugger/stereon_gdb.py prints, which are
Stereon's own std::format output; tests/debugger/check_printers.py runs both
debuggers against the same expected text.

LLDB-specific details: type names carry the module they are attached to
("Stereon::Interval@Stereon.Core"), which the matchers allow and strip.
LLDB's demangler does not understand GCC's module mangling, so the vtable
symbol of a type-erased model arrives mangled; it is demangled with
binutils' c++filt, which uses the same libiberty demangler as GDB and ships
with GCC. Without c++filt a Curve still prints, as "Curve(?)".
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess

import lldb

MODULE_SUFFIX = re.compile(r"@[A-Za-z_][\w.]*")
ENTITY_KINDS = ["Vertex", "Edge", "Coedge", "Loop", "Face", "Shell", "Solid"]
INVALID_INDEX = 0xFFFFFFFF
SENSE_BIT = 0x80000000
CATEGORY = "stereon"
# An optional module decoration after an identifier. LLDB's regexes are POSIX
# extended: no (?:...) groups and no \w.
M = r"(@[A-Za-z_][A-Za-z0-9_.]*)?"


def plain_name(value) -> str:
    return MODULE_SUFFIX.sub("", value.GetType().GetCanonicalType().GetUnqualifiedType().GetName())


def number(x: float) -> str:
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if x == int(x) and abs(x) < 1e16:
        return str(int(x))
    return repr(x)


def double_of(value) -> float:
    error = lldb.SBError()
    x = value.GetData().GetDouble(error, 0)
    return x if error.Success() else float("nan")


def member(value, name):
    return value.GetNonSyntheticValue().GetChildMemberWithName(name)


def tuple_text(value, axes) -> str:
    return "(" + ", ".join(number(double_of(member(value, a))) for a in axes) + ")"


def interval_text(value) -> str:
    lower, upper = double_of(member(value, "m_Lower")), double_of(member(value, "m_Upper"))
    if lower > upper:
        return "[]"
    return f"[{number(lower)}, {number(upper)}]"


def box_text(value, axes) -> str:
    intervals = [member(value, a) for a in axes]
    bounds = [(double_of(member(i, "m_Lower")), double_of(member(i, "m_Upper"))) for i in intervals]
    if any(lower > upper for lower, upper in bounds):
        return "[]"
    low = "(" + ", ".join(number(b[0]) for b in bounds) + ")"
    high = "(" + ", ".join(number(b[1]) for b in bounds) + ")"
    return f"[{low}, {high}]"


def kind_of(name: str) -> str:
    match = re.search(r"\(Stereon::EntityKind\)(\d+)", name)
    if match and int(match.group(1)) < len(ENTITY_KINDS):
        return ENTITY_KINDS[int(match.group(1))]
    match = re.search(r"Stereon::EntityKind::(\w+)", name)
    return match.group(1) if match else "Entity"


def handle_text(kind: str, index: int) -> str:
    return f"{kind}#-" if index == INVALID_INDEX else f"{kind}#{index}"


def id_text(value) -> str:
    return handle_text(kind_of(plain_name(value)), member(value, "m_Index").GetValueAsUnsigned())


def builder_id_text(value) -> str:
    index = member(value, "m_Index").GetValueAsUnsigned()
    text = handle_text(kind_of(plain_name(value)), index)
    if index != INVALID_INDEX:
        text += f"@{member(value, 'm_Generation').GetValueAsUnsigned()}"
    return text


def oriented_text(value) -> str:
    kind = kind_of(plain_name(value))
    bits = member(value, "m_Bits").GetValueAsUnsigned()
    if bits == INVALID_INDEX:
        return f"{kind}#-"
    text = f"{kind}#{bits & ~SENSE_BIT}"
    return text + " (reversed)" if bits & SENSE_BIT else text


def short_geometry_name(full: str) -> str:
    name = MODULE_SUFFIX.sub("", full)
    builtin = re.fullmatch(r"Stereon::Basic(Line|Circle|Ellipse)<double, (\d)[uUlL]*>", name)
    if builtin:
        return builtin.group(1) + ("2d" if builtin.group(2) == "2" else "")
    builtin = re.fullmatch(r"Stereon::Basic(Plane|Cylinder|Cone|Sphere|Torus)<double>", name)
    if builtin:
        return builtin.group(1)
    return name


def split_template_arguments(text: str) -> list[str]:
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


DEMANGLED = {}


def demangle(symbol: str) -> str:
    """The demangled vtable name, using c++filt for GCC module mangling."""
    if symbol.startswith("vtable for "):
        return symbol
    mangled = symbol.split("$")[-1]
    if mangled in DEMANGLED:
        return DEMANGLED[mangled]
    tool = shutil.which("c++filt")
    text = ""
    if tool is not None:
        text = subprocess.run([tool, mangled], capture_output=True, text=True).stdout.strip()
    DEMANGLED[mangled] = text
    return text


def without_integer_suffixes(text: str) -> str:
    """3ull or 3ul (Windows and Linux spell size_t differently) becomes 3."""
    return re.sub(r"\b(\d+)[uUlL]+\b", r"\1", text)


TYPE_INDEX = {}


def type_index(target):
    """Every class and struct type of the target, by its name without module
    decorations and spaces. Built once per target: LLDB's FindTypes does not
    match a base name whose stored name carries a decoration."""
    key = target.GetExecutable().fullpath
    if key not in TYPE_INDEX:
        index = {}
        for m in range(target.GetNumModules()):
            types = target.GetModuleAtIndex(m).GetTypes(lldb.eTypeClassClass | lldb.eTypeClassStruct)
            for i in range(types.GetSize()):
                t = types.GetTypeAtIndex(i)
                index.setdefault(MODULE_SUFFIX.sub("", t.GetName()).replace(" ", ""), t)
        TYPE_INDEX[key] = index
    return TYPE_INDEX[key]


def find_type(target, name: str):
    """The SBType named name, ignoring module decorations; see find_type in
    stereon_gdb.py."""
    found = target.FindFirstType(name)
    if found.IsValid():
        return found
    return type_index(target).get(MODULE_SUFFIX.sub("", name).replace(" ", ""))


def wrapped_value(value, model_address: int):
    target = value.GetTarget()
    process = value.GetProcess()
    error = lldb.SBError()
    vptr = process.ReadPointerFromMemory(model_address, error)
    symbol = target.ResolveLoadAddress(vptr).GetSymbol().GetName() or ""
    match = re.search(r"vtable for (.*)$", demangle(symbol))
    if not match:
        return "?", None
    arguments = split_template_arguments(match.group(1))
    type_name = arguments[-2]
    value_type = find_type(target, without_integer_suffixes(type_name))
    if value_type is None:
        return type_name, None
    size = int(without_integer_suffixes(arguments[-1]))
    pointer_size = process.GetAddressByteSize()
    alignment = max(pointer_size, 8)
    address = model_address + alignment
    if value_type.GetByteSize() + pointer_size <= size:
        held = target.CreateValueFromAddress("value", lldb.SBAddress(address, target), value_type)
        return type_name, held
    element = process.ReadPointerFromMemory(address, error)
    held = target.CreateValueFromAddress("value", lldb.SBAddress(element, target), value_type)
    return type_name + ", shared", held


def erased_text(label: str, value) -> str:
    model = member(member(value, "m_Storage"), "m_Model").GetValueAsUnsigned()
    if model == 0:
        return f"{label} (empty)"
    type_name, held = wrapped_value(value, model)
    shared = type_name.endswith(", shared")
    short = short_geometry_name(type_name.removesuffix(", shared"))
    head = f"{label}({short}{', shared' if shared else ''})"
    return head if held is None else f"{head} {value_text(held)}"


def vector_elements(vector):
    impl = member(vector, "_M_impl")
    start = member(impl, "_M_start")
    finish = member(impl, "_M_finish")
    element_type = start.GetType().GetPointeeType()
    size = element_type.GetByteSize()
    count = (finish.GetValueAsUnsigned() - start.GetValueAsUnsigned()) // size if size else 0
    target = vector.GetTarget()
    return [
        target.CreateValueFromAddress(
            f"[{i}]", lldb.SBAddress(start.GetValueAsUnsigned() + i * size, target), element_type)
        for i in range(count)
    ]


def arena_text(value) -> str:
    pages = vector_elements(member(value, "m_Pages"))
    size = member(value, "m_Size").GetValueAsUnsigned()
    shared = 0
    for page in pages:
        counts = member(member(page, "_M_refcount"), "_M_pi")
        if counts.GetValueAsUnsigned() != 0:
            if member(counts.Dereference(), "_M_use_count").GetValueAsSigned() > 1:
                shared += 1
    editing = ", editing" if member(value, "m_Editing").GetValueAsUnsigned() else ""
    return f"PagedArena of {size} entries in {len(pages)} pages, {shared} shared{editing}"


def expansion_text(value) -> str:
    length = member(value, "m_Length").GetValueAsUnsigned()
    components = member(member(value, "m_Components"), "_M_elems")
    return "Expansion [" + ", ".join(
        number(double_of(components.GetChildAtIndex(i))) for i in range(length)) + "]"


MATCHERS = [
    (re.compile(r"^Stereon::Basic(Vec|Point)2<.*>$"), lambda v: tuple_text(v, ["X", "Y"])),
    (re.compile(r"^Stereon::Basic(Vec|Point)3<.*>$"), lambda v: tuple_text(v, ["X", "Y", "Z"])),
    (re.compile(r"^Stereon::Interval$"), interval_text),
    (re.compile(r"^Stereon::Box2$"), lambda v: box_text(v, ["X", "Y"])),
    (re.compile(r"^Stereon::Box3$"), lambda v: box_text(v, ["X", "Y", "Z"])),
    (re.compile(r"^Stereon::Id<.*>$"), id_text),
    (re.compile(r"^Stereon::BuilderId<.*>$"), builder_id_text),
    (re.compile(r"^Stereon::Oriented<.*>$"), oriented_text),
    (re.compile(r"^Stereon::BasicErasedCurve<3>$"), lambda v: erased_text("Curve", v)),
    (re.compile(r"^Stereon::BasicErasedCurve<2>$"), lambda v: erased_text("Curve2d", v)),
    (re.compile(r"^Stereon::Surface$"), lambda v: erased_text("Surface", v)),
    (re.compile(r"^Stereon::PagedArena<.*>$"), arena_text),
    (re.compile(r"^Stereon::JournalId$"), lambda v: (
        "J#-" if member(v, "Value").GetValueAsUnsigned() == 0
        else f"J#{member(v, 'Value').GetValueAsUnsigned()}")),
    (re.compile(r"^Stereon::Sign$"),
     lambda v: {-1: "Negative", 0: "Zero", 1: "Positive"}.get(v.GetValueAsSigned(), str(v.GetValueAsSigned()))),
    (re.compile(r"^Stereon::Expansion<.*>$"), expansion_text),
]


def summary_of(value):
    name = plain_name(value)
    for pattern, render in MATCHERS:
        if pattern.match(name):
            return render(value)
    return None


def value_text(value) -> str:
    text = summary_of(value)
    if text is not None:
        return text
    basic = value.GetType().GetCanonicalType().GetBasicType()
    if basic in (lldb.eBasicTypeDouble, lldb.eBasicTypeFloat):
        return number(double_of(value))
    if basic == lldb.eBasicTypeBool:
        return "true" if value.GetValueAsUnsigned() else "false"
    raw = value.GetNonSyntheticValue()
    if raw.GetNumChildren() > 0 and value.GetType().GetCanonicalType().GetTypeClass() in (
            lldb.eTypeClassStruct, lldb.eTypeClassClass):
        parts = [f"{raw.GetChildAtIndex(i).GetName()} = {value_text(raw.GetChildAtIndex(i))}"
                 for i in range(raw.GetNumChildren())]
        return "{" + ", ".join(parts) + "}"
    return raw.GetValue() or ""


def summary(value, internal_dict):
    text = summary_of(value)
    return text if text is not None else ""


# The regexes LLDB matches against type names, module decorations allowed.
REGEXES = [
    rf"^Stereon::Basic(Vec|Point)[23]{M}<.+>$",
    rf"^Stereon::Interval{M}$",
    rf"^Stereon::Box[23]{M}$",
    rf"^Stereon::Id{M}<.+>$",
    rf"^Stereon::BuilderId{M}<.+>$",
    rf"^Stereon::Oriented{M}<.+>$",
    rf"^Stereon::BasicErasedCurve{M}<[23]>$",
    rf"^Stereon::Surface{M}$",
    rf"^Stereon::PagedArena{M}<.+>$",
    rf"^Stereon::JournalId{M}$",
    rf"^Stereon::Sign{M}$",
    rf"^Stereon::Expansion{M}<.+>$",
]


def __lldb_init_module(debugger, internal_dict):
    for regex in REGEXES:
        debugger.HandleCommand(
            f'type summary add -w {CATEGORY} -x "{regex}" -F {__name__}.summary')
    debugger.HandleCommand(f"type category enable {CATEGORY}")
