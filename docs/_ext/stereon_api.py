# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""The API reference of the Stereon manual, read from the sources and Doxygen.

Stereon's public API is flat: every library exports its names into the
``Stereon`` namespace, most of it as free functions overloaded per type
(``Eval``, ``Project``, ``ToBSpline``), and the types relate through concepts
and type erasure rather than base classes (ADR-0003). A namespace dump, or
Doxygen's inheritance graphs, show none of that. This extension organises
the reference the way the code is organised instead, by library and by
module partition, and draws the relations the code actually has.

Directives:

``.. stereon-library-graph::``
    The libraries and their dependencies, from ``stereon_add_library(...
    DEPENDS ...)`` in ``src/<library>/CMakeLists.txt``. ``:highlight:``
    marks one library.

``.. stereon-partition-graph:: <library>``
    The partitions of the library's module and the partitions and modules
    each imports, from the ``import`` lines of its ``.cppm`` files.

``.. stereon-type-map:: <library>``
    How the library's types fit together: for curves and surfaces the
    type-erased wrappers, the built-in types they dispatch over (the
    ``TypeList`` in ``Visit.cppm``), the class templates behind them and the
    concepts a user type models; for topology, the entity records and the
    handles that link them (``Records.cppm``).

``.. stereon-api:: <library>``
    The reference: an index of the operations, then one section per
    partition with its description, a summary table, every type rendered by
    Breathe with its collaboration graph, and the free functions grouped
    into overload sets. Without Doxygen XML it renders a note instead.

The graphs need only the sources and Graphviz; the reference needs the
Doxygen XML that ``conf.py`` locates (``breathe_projects``).
"""

from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from functools import lru_cache

from docutils import nodes
from docutils.parsers.rst import directives
from docutils.statemachine import StringList
from sphinx.util import logging
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles

logger = logging.getLogger(__name__)

# The libraries in layering order (ADR-0015; tools/check_layering.py).
LIBRARIES = ("core", "robust", "curve", "surface", "intersect", "topo",
             "build", "boolean", "blend", "mesh", "io")

# The wrappers of ADR-0003: name, the TypeList of built-in types it
# dispatches over, and the concept a user type models to be wrapped. The
# extension warns when a name here is missing from the sources.
TYPE_MAPS = {
    "curve": (("Curve", "BuiltinCurveTypes", "CurveGeometry"),
              ("Curve2d", "BuiltinCurve2dTypes", "Curve2dGeometry")),
    "surface": (("Surface", "BuiltinSurfaceTypes", "SurfaceGeometry"),),
}

# Namespaces that are not API (Doxyfile.in excludes them too).
HIDDEN_SCOPES = ("Stereon::Detail", "Stereon::Testing")


# -- The sources ----------------------------------------------------------------


@dataclass
class Partition:
    name: str            # "Line", or "" for the primary module interface unit
    path: str            # relative to src/, forward slashes
    summary: str         # the partition's leading /// comment, as one paragraph
    imports: list[str]   # partitions of the same module it imports
    modules: list[str]   # other modules it imports ("Stereon.Robust")


@dataclass
class Library:
    name: str
    module: str          # "Stereon.Curve"
    depends: list[str]
    summary: str = ""
    partitions: list[Partition] = field(default_factory=list)


_MODULE_RE = re.compile(r"^export\s+module\s+(Stereon\.\w+)(?::(\w+))?\s*;", re.M)
_IMPORT_RE = re.compile(r"^(?:export\s+)?import\s+(:\w+|Stereon\.\w+)\s*;", re.M)


def _leading_comment(text: str) -> str:
    """The file's own /// comment: the first /// block that is followed by a blank
    line or by the ``export module`` line, rather than by a declaration."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].lstrip().startswith("///"):
            block = []
            while i < len(lines) and lines[i].lstrip().startswith("///"):
                block.append(lines[i].lstrip()[3:].strip())
                i += 1
            following = lines[i].strip() if i < len(lines) else ""
            if following == "" or following.startswith("export module"):
                return " ".join(part for part in block if part)
            continue
        i += 1
    return ""


def _library_spec(src: str, name: str) -> tuple[list[str], list[str]]:
    """MODULES and DEPENDS of stereon_add_library(<name> ...)."""
    with open(os.path.join(src, name, "CMakeLists.txt"), encoding="utf-8") as f:
        text = re.sub(r"#[^\n]*", "", f.read())
    match = re.search(r"stereon_add_library\s*\(\s*" + name + r"\b(.*?)\)", text, re.S)
    if not match:
        return [], []
    words = match.group(1).split()
    groups: dict[str, list[str]] = {}
    key = None
    for word in words:
        if word.isupper():
            key = word
            groups.setdefault(key, [])
        elif key:
            groups[key].append(word)
    return groups.get("MODULES", []), groups.get("DEPENDS", [])


@lru_cache(maxsize=None)
def read_libraries(src: str) -> dict[str, Library]:
    libraries: dict[str, Library] = {}
    for name in LIBRARIES:
        modules, depends = _library_spec(src, name)
        library = Library(name, "", depends)
        units: dict[str, Partition] = {}
        order: list[str] = []
        for unit in modules:
            path = os.path.join(src, name, unit)
            if not os.path.isfile(path):
                continue
            with open(path, encoding="utf-8") as f:
                text = f.read()
            declared = _MODULE_RE.search(text)
            if not declared:
                continue
            library.module = declared.group(1)
            partition = declared.group(2) or ""
            imports = _IMPORT_RE.findall(text)
            unit_record = Partition(
                partition, f"{name}/{unit}", _leading_comment(text),
                [i[1:] for i in imports if i.startswith(":")],
                [i for i in imports if not i.startswith(":")])
            units[partition] = unit_record
            if partition == "":
                library.summary = unit_record.summary
                # The primary interface's export-import order is the module's own order.
                order = [i[1:] for i in imports if i.startswith(":")]
        library.partitions = [units[p] for p in order if p in units]
        library.partitions += [u for p, u in sorted(units.items()) if p and p not in order]
        if "" in units:
            library.partitions.insert(0, units[""])
        libraries[name] = library
    return libraries


def _type_lists(src: str, library: str) -> tuple[dict[str, list[str]], dict[str, str]]:
    """The library's ``using X = TypeList<...>`` lists, and its aliases of class
    templates (``using Line = BasicLine<double, 3>`` maps Line to BasicLine)."""
    lists: dict[str, list[str]] = {}
    aliases: dict[str, str] = {}
    folder = os.path.join(src, library)
    for unit in sorted(os.listdir(folder)):
        if not unit.endswith(".cppm"):
            continue
        with open(os.path.join(folder, unit), encoding="utf-8") as f:
            text = f.read()
        for name, body in re.findall(r"using\s+(\w+)\s*=\s*TypeList<([^>]*)>\s*;", text):
            lists[name] = [t.strip() for t in body.split(",") if t.strip()]
        for name, template in re.findall(r"using\s+(\w+)\s*=\s*(Basic\w+)\s*<", text):
            aliases[name] = template
    return lists, aliases


def _records(src: str) -> dict[str, list[tuple[str, str]]]:
    """The topology records and their handle fields: {"Edge": [("Start", "Vertex"), ...]}."""
    with open(os.path.join(src, "topo", "Records.cppm"), encoding="utf-8") as f:
        text = f.read()
    records: dict[str, list[tuple[str, str]]] = {}
    for kind, body in re.findall(r"struct\s+(\w+)Record\s*\{(.*?)\n\};", text, re.S):
        records[kind] = re.findall(r"^\s*(\w+)Id\s+(\w+)\s*;", body, re.M)
        records[kind] = [(member, target) for target, member in records[kind]]
    return records


# -- Graphs ---------------------------------------------------------------------

_GRAPH_STYLE = (
    'graph [fontname="Helvetica", fontsize=10, nodesep=0.25, ranksep=0.35, bgcolor="transparent"];',
    'node [fontname="Helvetica", fontsize=10, shape=box, style="rounded,filled", '
    'fillcolor="#f2f5fa", color="#5b6b82", penwidth=1];',
    'edge [color="#5b6b82", arrowsize=0.7, fontname="Helvetica", fontsize=9];',
)


def _quote(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _first_sentence(text: str) -> str:
    match = re.match(r"(.+?[.:])(\s|$)", text)
    return (match.group(1) if match else text).rstrip(":")


def _transitive_reduction(edges: set[tuple[str, str]]) -> set[tuple[str, str]]:
    """Drop every edge implied by a longer path, so the graph shows only direct structure."""
    succ: dict[str, set[str]] = {}
    for a, b in edges:
        succ.setdefault(a, set()).add(b)

    def reachable(start: str, skip: tuple[str, str]) -> set[str]:
        seen: set[str] = set()
        stack = [n for n in succ.get(start, ()) if (start, n) != skip]
        while stack:
            n = stack.pop()
            if n not in seen:
                seen.add(n)
                stack.extend(succ.get(n, ()))
        return seen

    return {(a, b) for a, b in edges if b not in reachable(a, (a, b))}


def _graphviz_block(dot: str, caption: str) -> list[str]:
    lines = [".. graphviz::", "   :align: center", f"   :caption: {caption}", ""]
    lines += ["   " + line for line in dot.splitlines()]
    return lines + [""]


def library_graph(libraries: dict[str, Library], highlight: str | None, link: str,
                  neighbourhood: bool = False) -> str:
    """Every library, or with ``neighbourhood`` only the highlighted one, the libraries
    it uses directly and the ones that use it directly."""
    if neighbourhood and highlight in libraries:
        keep = {highlight, *libraries[highlight].depends}
        keep |= {n for n, lib in libraries.items() if highlight in lib.depends}
        libraries = {n: lib for n, lib in libraries.items() if n in keep}
        libraries = {n: Library(lib.name, lib.module, [d for d in lib.depends if d in keep], lib.summary)
                     for n, lib in libraries.items()}
    out = ["digraph libraries {", *_GRAPH_STYLE, "rankdir=BT;"]
    for name, library in libraries.items():
        attrs = [f"label={_quote(library.module or name)}",
                 f"tooltip={_quote(_first_sentence(library.summary.split(': ', 1)[-1]) or name)}",
                 f"URL={_quote(link.format(name=name))}", 'target="_top"']
        if name == highlight:
            attrs += ['fillcolor="#2f5d9a"', 'fontcolor="white"', 'color="#2f5d9a"']
        out.append(f"{_quote(name)} [{', '.join(attrs)}];")
    for name, library in libraries.items():
        for dep in library.depends:
            out.append(f"{_quote(name)} -> {_quote(dep)};")
    out.append("}")
    return "\n".join(out)


def partition_graph(library: Library) -> str:
    out = ["digraph partitions {", *_GRAPH_STYLE, "rankdir=LR;"]
    names = {p.name for p in library.partitions if p.name}
    edges: set[tuple[str, str]] = set()
    external: set[str] = set()
    for unit in library.partitions:
        if not unit.name:
            continue
        for dep in unit.imports:
            if dep in names:
                edges.add((unit.name, dep))
        for module in unit.modules:
            edges.add((unit.name, module))
            external.add(module)
    edges = _transitive_reduction(edges)
    for unit in library.partitions:
        if unit.name:
            tooltip = _first_sentence(unit.summary) or unit.name
            out.append(f"{_quote(unit.name)} [label={_quote(':' + unit.name)}, "
                       f"tooltip={_quote(tooltip)}, URL={_quote('#api-' + library.name + '-' + unit.name.lower())}, "
                       'target="_top"];')
    for module in sorted(external):
        out.append(f'{_quote(module)} [shape=box, style="dashed,rounded", fillcolor="white"];')
    for a, b in sorted(edges):
        out.append(f"{_quote(a)} -> {_quote(b)};")
    out.append("}")
    return "\n".join(out)


def type_map(src: str, library: str) -> tuple[str, list[str]]:
    """The type-map graph of a library and the names it could not find."""
    missing: list[str] = []
    out = ["digraph types {", *_GRAPH_STYLE, "rankdir=LR;", "newrank=true;"]
    if library == "topo":
        records = _records(src)
        for kind in records:
            out.append(f"{_quote(kind)} [label={_quote(kind + 'Record')}];")
        for kind, links in records.items():
            for member, target in links:
                if target in records:
                    out.append(f"{_quote(kind)} -> {_quote(target)} [label={_quote(member)}];")
        if not records:
            missing.append("topo/Records.cppm records")
        out.append("}")
        return "\n".join(out), missing

    # Left to right: wrapper, its built-in types, the class templates behind them.
    # A second wrapper (Curve2d) is mirrored on the right, its edges drawn
    # backwards, so both share the column of templates in the middle.
    lists, aliases = _type_lists(src, library)
    templates: list[str] = []
    for side, (wrapper, list_name, concept) in enumerate(TYPE_MAPS.get(library, ())):
        builtins = lists.get(list_name)
        if builtins is None:
            missing.append(list_name)
            continue

        def edge(a: str, b: str, attrs: str = "", mirrored: bool = side == 1) -> str:
            if mirrored:
                return f"{_quote(b)} -> {_quote(a)} [dir=back{', ' + attrs if attrs else ''}];"
            return f"{_quote(a)} -> {_quote(b)}{' [' + attrs + ']' if attrs else ''};"

        out.append(f'{_quote(wrapper)} [fillcolor="#2f5d9a", fontcolor="white", color="#2f5d9a"];')
        out.append(f'{_quote(concept)} [shape=ellipse, style="filled", fillcolor="#fff6e0", color="#a8842c"];')
        out.append(edge(wrapper, concept, "style=dashed"))
        for builtin in builtins:
            out.append(edge(wrapper, builtin))
            template = aliases.get(builtin)
            if template:
                if template not in templates:
                    templates.append(template)
                out.append(edge(builtin, template, "style=dotted, arrowhead=empty, arrowtail=empty"))
    for template in templates:
        out.append(f'{_quote(template)} [style="rounded", fillcolor="white"];')
    out.append("}")
    return "\n".join(out), missing


# -- Doxygen XML ----------------------------------------------------------------


@dataclass
class Compound:
    kind: str            # class, struct, union, concept
    name: str            # Stereon::BasicLine
    path: str
    line: int
    brief: ET.Element | None
    has_graph: bool


@dataclass
class Member:
    kind: str            # function, typedef, variable, enum
    name: str
    path: str
    line: int
    element: ET.Element


@dataclass
class ApiModel:
    compounds: list[Compound]
    members: list[Member]
    type_names: set[str]


def _location(element: ET.Element) -> tuple[str, int]:
    """Where a declaration is, preferring the declaration over a definition in a .cpp."""
    location = element.find("location")
    if location is None:
        return "", 0
    path = location.get("declfile") or location.get("file") or ""
    line = location.get("declline") or location.get("line") or "0"
    return path.replace("\\", "/"), int(line)


def _graph_nodes(element: ET.Element, tag: str) -> int:
    graph = element.find(tag)
    return 0 if graph is None else len(graph.findall("node"))


@lru_cache(maxsize=4)
def read_model(xml_dir: str, stamp: float) -> ApiModel:
    del stamp  # part of the cache key only: a new Doxygen run gives a new model
    index = ET.parse(os.path.join(xml_dir, "index.xml")).getroot()
    compounds: list[Compound] = []
    members: list[Member] = []
    for entry in index.findall("compound"):
        kind = entry.get("kind")
        name = entry.findtext("name") or ""
        if kind in ("class", "struct", "union", "concept"):
            if name.count("::") != 1 or not name.startswith("Stereon::"):
                continue  # nested types are rendered inside their parent
            root = ET.parse(os.path.join(xml_dir, entry.get("refid") + ".xml")).getroot()
            definition = root.find("compounddef")
            path, line = _location(definition)
            compounds.append(Compound(kind, name, path, line, definition.find("briefdescription"),
                                      _graph_nodes(definition, "collaborationgraph") > 1))
        elif kind == "namespace" and name == "Stereon":
            root = ET.parse(os.path.join(xml_dir, entry.get("refid") + ".xml")).getroot()
            for section in root.find("compounddef").findall("sectiondef"):
                for member in section.findall("memberdef"):
                    if member.get("kind") not in ("function", "typedef", "variable", "enum"):
                        continue
                    if re.search(r"[:<]", member.findtext("name") or ""):
                        continue  # an explicit instantiation Doxygen files under the namespace
                    path, line = _location(member)
                    members.append(Member(member.get("kind"), member.findtext("name") or "",
                                          path, line, member))
    type_names = {c.name.split("::", 1)[1] for c in compounds}
    type_names |= {m.name for m in members if m.kind in ("typedef", "enum")}
    return ApiModel(compounds, members, type_names)


# -- Doxygen descriptions as reStructuredText -----------------------------------

_BEFORE_OK = set(" \t\n-:/'\"<([{")
_AFTER_OK = set(" \t\n-.,:;!?\\/'\")]}>")


class Describer:
    """Doxygen description XML to reStructuredText lines.

    Doxygen nests block elements (lists, code, notes) inside paragraphs, so a
    paragraph is flushed at every block child and resumed after it. Inline
    markup is separated from adjacent text by an escaped space wherever
    reStructuredText would otherwise not recognise it.
    """

    BLOCKS = {"itemizedlist", "orderedlist", "programlisting", "verbatim", "simplesect",
              "parameterlist", "table", "blockquote", "heading", "sect1", "sect2", "sect3",
              "hruler", "image", "xrefsect", "para", "title", "preformatted"}

    def __init__(self, type_names: set[str]):
        self.type_names = type_names

    # Inline -------------------------------------------------------------------

    @staticmethod
    def escape(text: str) -> str:
        text = text.replace("\\", "\\\\")
        for ch in "*`|_":
            text = text.replace(ch, "\\" + ch)
        return text

    @staticmethod
    def plain(element: ET.Element) -> str:
        parts: list[str] = []

        def walk(e: ET.Element) -> None:
            if e.tag == "sp":
                parts.append(" ")
            elif e.text:
                parts.append(e.text)
            for child in e:
                walk(child)
                if child.tail:
                    parts.append(child.tail)

        walk(element)
        return "".join(parts)

    def markup(self, child: ET.Element) -> tuple[str, bool]:
        """(text, is_markup) for one inline element."""
        tag = child.tag
        if tag in ("sp", "linebreak"):
            return " ", False
        if tag == "anchor":
            return "", False
        text = " ".join(self.plain(child).split())
        if not text:
            return "", False
        if tag == "computeroutput":
            if "``" in text:
                return self.escape(text), False
            return f"``{text}``", True
        if tag == "ref":
            target = re.sub(r"\(\)$", "", text)
            if re.fullmatch(r"[A-Za-z_]\w*", target) and (child.get("kindref") == "member"
                                                           or target in self.type_names):
                # Qualified: a short name does not resolve inside the generated content.
                return f":cpp:any:`{text} <Stereon::{target}>`", True
            return self.escape(text), False
        if tag in ("bold", "emphasis"):
            mark = "**" if tag == "bold" else "*"
            return f"{mark}{self.escape(text)}{mark}", True
        if tag == "ulink":
            url = child.get("url", "")
            label = self.escape(text).replace("<", "\\<")
            return (f"`{label} <{url}>`__", True) if url else (self.escape(text), False)
        if tag == "formula":
            body = text.strip().strip("$")
            body = re.sub(r"^\\\[|\\\]$", "", body).strip()
            return f":math:`{body}`", True
        return self.escape(text), False

    @staticmethod
    def join(tokens: list[tuple[str, bool]]) -> str:
        out = ""
        previous_markup = False
        for text, is_markup in tokens:
            if not text:
                continue
            if is_markup and out and (out[-1] not in _BEFORE_OK or previous_markup):
                out += "\\ "
            elif previous_markup and not is_markup and text[0] not in _AFTER_OK:
                out += "\\ "
            out += text
            previous_markup = is_markup
        out = " ".join(out.split())
        if re.match(r"^([-*+]|\d+[.)]|#\.)\s", out):
            out = "\\" + out
        return out

    # Blocks -------------------------------------------------------------------

    def blocks(self, container: ET.Element | None) -> list[str]:
        """The lines of a description container (briefdescription, para, listitem...)."""
        if container is None:
            return []
        lines: list[str] = []
        tokens: list[tuple[str, bool]] = []

        def separate() -> None:
            if lines and lines[-1] != "":
                lines.append("")

        def flush() -> None:
            paragraph = self.join(tokens)
            tokens.clear()
            if paragraph.strip():
                separate()
                lines.extend([paragraph, ""])

        if container.text and container.text.strip():
            tokens.append((self.escape(container.text), False))
        for child in container:
            if child.tag in self.BLOCKS:
                flush()
                block = self.block(child)
                if block:
                    separate()
                    lines.extend(block + [""])
            else:
                tokens.append(self.markup(child))
            if child.tail and child.tail.strip():
                tokens.append((self.escape(child.tail), False))
            elif child.tail and tokens:
                tokens.append((" ", False))
        flush()
        while lines and lines[-1] == "":
            lines.pop()
        return lines

    @staticmethod
    def indent(lines: list[str], by: str) -> list[str]:
        return [by + line if line else "" for line in lines]

    def block(self, child: ET.Element) -> list[str]:
        tag = child.tag
        if tag == "para":
            return self.blocks(child)
        if tag in ("itemizedlist", "orderedlist"):
            bullet = "- " if tag == "itemizedlist" else "#. "
            out: list[str] = []
            for item in child.findall("listitem"):
                body = self.blocks(item) or [""]
                out.append(bullet + body[0])
                out.extend(self.indent(body[1:], " " * len(bullet)))
            return out
        if tag == "programlisting":
            code = [self.plain(line) for line in child.findall("codeline")]
            return [".. code-block:: cpp", ""] + ["   " + c for c in code]
        if tag in ("verbatim", "preformatted"):
            return ["::", ""] + ["   " + c for c in (child.text or "").splitlines()]
        if tag == "simplesect":
            return self.simplesect(child)
        if tag == "parameterlist":
            field_name = {"param": "param", "templateparam": "tparam", "exception": "throws",
                          "retval": "retval"}.get(child.get("kind", "param"), "param")
            out = []
            for item in child.findall("parameteritem"):
                names = [" ".join(self.plain(n).split()) for n in item.iter("parametername")]
                body = self.blocks(item.find("parameterdescription")) or [""]
                out.append(f":{field_name} {', '.join(names)}: {body[0]}")
                out.extend(self.indent(body[1:], "   "))
            return out
        if tag == "table":
            return self.table(child)
        if tag == "blockquote":
            return self.indent(self.blocks(child), "   ")
        if tag in ("heading", "title"):
            text = " ".join(self.plain(child).split())
            return [f".. rubric:: {text}"] if text else []
        if tag in ("sect1", "sect2", "sect3"):
            title = child.findtext("title") or ""
            body: list[str] = []
            for part in child:
                if part.tag != "title":
                    body += self.block(part) + [""]
            return ([f".. rubric:: {title}", ""] if title else []) + body
        if tag == "xrefsect":
            title = child.findtext("xreftitle") or "Note"
            body = self.blocks(child.find("xrefdescription"))
            return ([".. admonition:: " + title, ""] + self.indent(body, "   ")) if body else []
        return []

    def simplesect(self, child: ET.Element) -> list[str]:
        kind = child.get("kind")
        body = self.blocks(child)
        if not body:
            return []
        if kind in ("note", "warning", "attention", "remark"):
            name = "warning" if kind == "warning" else "note"
            return [f".. {name}::", ""] + self.indent(body, "   ")
        if kind == "see":
            return [".. seealso::", ""] + self.indent(body, "   ")
        if kind == "return":
            return [":returns: " + body[0]] + self.indent(body[1:], "   ")
        labels = {"pre": "Precondition", "post": "Postcondition", "invariant": "Invariant",
                  "since": "Since", "deprecated": "Deprecated"}
        if kind in labels:
            return [f"**{labels[kind]}:** " + body[0]] + body[1:]
        title = child.findtext("title")
        return ([f".. rubric:: {title}", ""] if title else []) + body

    def table(self, child: ET.Element) -> list[str]:
        rows = child.findall("row")
        if not rows:
            return []
        width = max(len(r.findall("entry")) for r in rows)
        out = [".. list-table::", ""]
        for row in rows:
            entries = row.findall("entry")
            for i in range(width):
                body = (self.blocks(entries[i]) if i < len(entries) else []) or [""]
                out.append(("   * - " if i == 0 else "     - ") + body[0])
                out.extend(self.indent(body[1:], "       "))
        return out

    def description(self, element: ET.Element) -> list[str]:
        brief = self.blocks(element.find("briefdescription"))
        detail = self.blocks(element.find("detaileddescription"))
        if brief and detail:
            return brief + [""] + detail
        return brief or detail

    def brief_line(self, element: ET.Element | None) -> str:
        """The brief description on one line, for summary tables."""
        if element is None:
            return ""
        lines = self.blocks(element)
        return " ".join(line.strip() for line in lines if line.strip() and not line.startswith(".."))


# -- Declarations ---------------------------------------------------------------


def _template_prefix(element: ET.Element) -> str:
    params = element.find("templateparamlist")
    if params is None:
        return ""
    parts = []
    for param in params.findall("param"):
        kind = " ".join(Describer.plain(param.find("type")).split()) if param.find("type") is not None else ""
        name = param.findtext("declname") or ""
        text = kind if (not name or kind.endswith(" " + name) or kind == name) else f"{kind} {name}"
        default = param.find("defval")
        if default is not None:
            text += " = " + " ".join(Describer.plain(default).split())
        parts.append(text)
    prefix = f"template<{', '.join(parts)}>"
    requires = element.find("requiresclause")
    if requires is not None:
        prefix += " requires " + " ".join(Describer.plain(requires).split())
    return prefix + " "


def _text(element: ET.Element, tag: str) -> str:
    found = element.find(tag)
    return " ".join(Describer.plain(found).split()) if found is not None else ""


def declaration(member: Member) -> tuple[str, str]:
    """(directive, signature) for a namespace member, without the Stereon:: scope."""
    e = member.element
    template = _template_prefix(e)
    if member.kind == "function":
        return "function", f"{template}{_text(e, 'type')} {member.name}{_text(e, 'argsstring')}".strip()
    if member.kind == "typedef":
        aliased = _text(e, "type")
        return "type", f"{template}{member.name} = {aliased}".strip()
    if member.kind == "variable":
        # A short initializer is part of the declaration; a long one (a whole
        # catalogue) is left to the source.
        initializer = _text(e, "initializer")
        initializer = " " + initializer if initializer and len(initializer) <= 60 else ""
        return "var", f"{template}{_text(e, 'type')} {member.name}{initializer}".strip()
    strong = "enum-class" if e.get("strong") == "yes" else "enum"
    base = _text(e, "type")
    return strong, f"{member.name}{' : ' + base if base else ''}"


@lru_cache(maxsize=None)
def _parses(directive: str, signature: str) -> bool:
    """Whether Sphinx's C++ domain can parse the signature; a declaration it cannot
    is shown as a code block rather than producing a warning."""
    from sphinx.domains.cpp._parser import DefinitionParser
    from sphinx.util.cfamily import DefinitionError

    class _Config:
        cpp_id_attributes: list[str] = []
        cpp_paren_attributes: list[str] = []

    objtype = {"function": "function", "type": "type", "var": "member",
               "enum": "enum", "enum-class": "enum"}[directive]
    try:
        parser = DefinitionParser(signature, location=None, config=_Config())
        parser.parse_declaration(objtype, directive.replace("enum-class", "enum"))
        parser.skip_ws()
        parser.assert_end()
    except (DefinitionError, Exception):  # noqa: BLE001  (any failure means: show it verbatim)
        return False
    return True


def comment_rst(text: str) -> str:
    """A source comment as one line of reStructuredText: `code` spans become literals,
    everything else is escaped."""
    tokens = []
    for i, part in enumerate(re.split(r"`([^`]+)`", text)):
        tokens.append((f"``{part}``", True) if i % 2 else (Describer.escape(part), False))
    return Describer.join(tokens)


_OPERATOR_WORDS = {"+": "plus", "-": "minus", "*": "times", "/": "div", "%": "mod", "<": "lt",
                   ">": "gt", "=": "eq", "!": "not", "&": "and", "|": "or", "^": "xor",
                   "~": "compl", "[": "index", "(": "call"}


def _label(*parts: str) -> str:
    """A reference label; operator names keep their symbols as words (operator+ and
    operator- must not both become "operator")."""
    words = []
    for part in parts:
        if part.startswith("operator"):
            part = "operator" + "".join(_OPERATOR_WORDS.get(ch, ch) for ch in part[len("operator"):])
        words.append(re.sub(r"[^a-z0-9]+", "-", part.lower()).strip("-"))
    return "-".join(w for w in words if w)


# -- The reference --------------------------------------------------------------


@dataclass
class Group:
    """A partition's slice of the reference."""
    partition: Partition
    types: list[Compound] = field(default_factory=list)
    aliases: list[Member] = field(default_factory=list)
    functions: dict[str, list[Member]] = field(default_factory=dict)


def _groups(library: Library, model: ApiModel) -> list[Group]:
    by_path: dict[str, Group] = {}
    groups: list[Group] = []
    for unit in library.partitions:
        group = Group(unit)
        by_path[unit.path] = group
        groups.append(group)

    def group_for(path: str) -> Group | None:
        if not path.startswith(library.name + "/"):
            return None
        if path in by_path:
            return by_path[path]
        stem = os.path.splitext(path)[0]
        if stem + ".cppm" in by_path:
            return by_path[stem + ".cppm"]
        # Not API: a partition the module does not export (curve's :Homogeneous),
        # or a declaration made only in an implementation file.
        if "/include/" not in path:
            return None
        # A public header of the library (core's Exec.hpp, say): a section of its own.
        group = Group(Partition(os.path.basename(path), path, "", [], []))
        by_path[path] = group
        groups.append(group)
        return group

    for compound in sorted(model.compounds, key=lambda c: (c.path, c.line)):
        group = group_for(compound.path)
        if group is not None:
            group.types.append(compound)
    for member in sorted(model.members, key=lambda m: (m.path, m.line)):
        group = group_for(member.path)
        if group is None:
            continue
        if member.kind == "function":
            group.functions.setdefault(member.name, []).append(member)
        else:
            group.aliases.append(member)
    return [g for g in groups if g.types or g.aliases or g.functions]


class Reference:
    """The reStructuredText of one library's reference."""

    def __init__(self, library: Library, model: ApiModel):
        self.library = library
        self.model = model
        self.describe = Describer(model.type_names)
        self.groups = _groups(library, model)
        self.unparsed = 0

    def section_label(self, group: Group) -> str:
        return _label("api", self.library.name, group.partition.name or "module")

    def function_label(self, group: Group, name: str) -> str:
        return _label("api", self.library.name, group.partition.name or "module", "fn", name)

    def title(self, group: Group) -> str:
        unit = group.partition
        if not unit.name:
            return f"The module interface, ``{self.library.module}``"
        if unit.name.endswith((".hpp", ".h")):
            return f"``<{unit.path.split('/include/', 1)[-1]}>``"
        return unit.name

    def lines(self) -> list[str]:
        out: list[str] = [".. cpp:namespace:: Stereon", ""]
        out += self.operations()
        for group in self.groups:
            out += self.group(group)
        return out

    def operations(self) -> list[str]:
        """An index of the free functions: what you can do, and for which types."""
        names: dict[str, list[Group]] = {}
        briefs: dict[str, str] = {}
        for group in self.groups:
            for name, overloads in group.functions.items():
                names.setdefault(name, []).append(group)
                if name not in briefs:
                    for member in overloads:
                        brief = self.describe.brief_line(member.element.find("briefdescription"))
                        if brief:
                            briefs[name] = brief
                            break
        if not names:
            return []
        out = ["Operations at a glance", "~~~~~~~~~~~~~~~~~~~~~~", "",
               "Every free function of the library, by name, with the sections that declare it. "
               "Most are overloaded per type: the same call works on every type listed.", "",
               ".. list-table::", "   :header-rows: 1", "   :widths: 18 32 50", "   :class: stereon-ops", "",
               "   * - Function", "     - Declared in", "     - What it does"]
        for name in sorted(names, key=str.lower):
            where = ", ".join(f":ref:`{g.partition.name or 'module'} <{self.function_label(g, name)}>`"
                              for g in names[name])
            out += [f"   * - ``{name}``", f"     - {where}", f"     - {briefs.get(name, '')}"]
        return out + [""]

    def group(self, group: Group) -> list[str]:
        title = self.title(group)
        out = [f".. _{self.section_label(group)}:", "", title, "~" * max(len(title), 4), ""]
        unit = group.partition
        if unit.summary:
            out += [comment_rst(unit.summary), ""]
        if unit.path.endswith(".cppm") and unit.name:
            out += [f"*Partition* ``{self.library.module}:{unit.name}``, in ``src/{unit.path}``.", ""]

        summary = self.summary(group)
        if summary:
            out += summary
        for compound in group.types:
            kind = {"class": "doxygenclass", "struct": "doxygenstruct", "union": "doxygenunion",
                    "concept": "doxygenconcept"}[compound.kind]
            # Breathe qualifies the name itself, so its directives run at global scope
            # (inside "Stereon" they would declare Stereon::Stereon::X).
            out += [".. cpp:namespace:: nullptr", ""]
            out += [f".. {kind}:: {compound.name}", "   :project: Stereon"]
            if compound.kind != "concept":
                out += ["   :members:", "   :undoc-members:"]
                if compound.has_graph:
                    out += ["   :allow-dot-graphs:"]
            out += ["", ".. cpp:namespace:: Stereon", ""]
        for member in group.aliases:
            out += self.member(member)
        if group.functions:
            out += [".. rubric:: Functions", ""]
            for name, overloads in group.functions.items():
                out += [f".. _{self.function_label(group, name)}:", ""]
                for member in overloads:
                    out += self.member(member)
        return out

    def summary(self, group: Group) -> list[str]:
        rows = []
        for compound in group.types:
            short = compound.name.split("::", 1)[1]
            role = "cpp:concept" if compound.kind == "concept" else "cpp:class"
            rows.append((f":{role}:`{short} <{compound.name}>`", compound.kind,
                         self.describe.brief_line(compound.brief)))
        for member in group.aliases:
            kind = {"typedef": "alias", "variable": "constant", "enum": "enum"}[member.kind]
            brief = self.describe.brief_line(member.element.find("briefdescription"))
            if not brief and member.kind == "typedef":
                brief = f"Alias of ``{' '.join(_text(member.element, 'type').split())}``."
            rows.append((f":cpp:any:`{member.name} <Stereon::{member.name}>`", kind, brief))
        for name, overloads in group.functions.items():
            brief = next((b for b in (self.describe.brief_line(m.element.find("briefdescription"))
                                      for m in overloads) if b), "")
            count = f" ({len(overloads)} overloads)" if len(overloads) > 1 else ""
            rows.append((f":ref:`{name} <{self.function_label(group, name)}>`", "function" + count, brief))
        if len(rows) < 2:
            return []
        out = [".. list-table::", "   :header-rows: 1", "   :widths: 25 15 60", "   :class: stereon-summary", "",
               "   * - Name", "     - Kind", "     - Summary"]
        for name, kind, brief in rows:
            out += [f"   * - {name}", f"     - {kind}", f"     - {brief}"]
        return out + [""]

    def member(self, member: Member) -> list[str]:
        directive, signature = declaration(member)
        body = self.describe.description(member.element)
        if member.kind == "enum":
            body = body + [""] if body else []
            for value in member.element.findall("enumvalue"):
                body += [f".. cpp:enumerator:: {value.findtext('name')}", ""]
                value_body = self.describe.description(value)
                body += self.describe.indent(value_body, "   ") + ([""] if value_body else [])
        if _parses(directive, signature):
            out = [f".. cpp:{directive}:: {signature}", ""]
            out += self.describe.indent(body, "   ")
            return out + [""]
        self.unparsed += 1
        out = [".. code-block:: cpp", "   :class: stereon-signature", "", "   " + signature, ""]
        if body:
            out += [".. container:: stereon-signature-body", ""] + self.describe.indent(body, "   ") + [""]
        return out


# -- Directives -----------------------------------------------------------------


def _src(env) -> str:
    return os.path.join(os.path.dirname(env.srcdir), "src")


def _xml(env) -> str | None:
    project = env.config.breathe_default_project
    xml_dir = env.config.breathe_projects.get(project)
    if xml_dir and os.path.isfile(os.path.join(xml_dir, "index.xml")):
        return xml_dir
    return None


class _ParsedDirective(SphinxDirective):
    def parse(self, lines: list[str]) -> list[nodes.Node]:
        container = nodes.section() if any(set(line) == {"~"} for line in lines) else nodes.container()
        nested_parse_with_titles(self.state, StringList(lines, source=self.get_source_info()[0]), container)
        return container.children


class LibraryGraph(_ParsedDirective):
    has_content = False
    option_spec = {"highlight": directives.unchanged, "caption": directives.unchanged}

    def run(self) -> list[nodes.Node]:
        libraries = read_libraries(_src(self.env))
        highlight = self.options.get("highlight")
        # From an api/<library> page the others are siblings; from api.rst they are under api/.
        link = "{name}.html" if self.env.docname.startswith("api/") else "api/{name}.html"
        # A library's page shows its neighbourhood; the overview shows them all.
        neighbourhood = highlight is not None
        caption = self.options.get(
            "caption", "The libraries this one uses directly (above) and those that use it directly "
                       "(below); click one to open it. The overview shows them all."
            if neighbourhood else "The libraries and what each depends on (click to open).")
        return self.parse(_graphviz_block(library_graph(libraries, highlight, link, neighbourhood), caption))


class PartitionGraph(_ParsedDirective):
    required_arguments = 1
    option_spec = {"caption": directives.unchanged}

    def run(self) -> list[nodes.Node]:
        library = read_libraries(_src(self.env))[self.arguments[0]]
        if len([p for p in library.partitions if p.name]) < 2:
            return []
        caption = self.options.get(
            "caption", f"The partitions of ``{library.module}``: an arrow is an import; implied "
                       "imports are left out, and dashed boxes are other modules.")
        return self.parse(_graphviz_block(partition_graph(library), caption))


class TypeMap(_ParsedDirective):
    required_arguments = 1
    option_spec = {"caption": directives.unchanged}

    def run(self) -> list[nodes.Node]:
        dot, missing = type_map(_src(self.env), self.arguments[0])
        for name in missing:
            logger.warning("stereon-type-map: %s not found in the sources", name, location=self.get_location())
        return self.parse(_graphviz_block(dot, self.options.get("caption", "How the types fit together.")))


class ApiReference(_ParsedDirective):
    required_arguments = 1

    def run(self) -> list[nodes.Node]:
        library = read_libraries(_src(self.env))[self.arguments[0]]
        xml_dir = _xml(self.env)
        if xml_dir is None:
            return self.parse([".. note::", "",
                               "   This build has no Doxygen XML, so the reference below the overview is",
                               "   empty. Configure with ``-DSTEREON_BUILD_DOCS=ON`` and Doxygen installed,",
                               "   or set ``STEREON_RUN_DOXYGEN=1`` when running ``sphinx-build``.", ""])
        self.env.note_dependency(os.path.join(xml_dir, "index.xml"))
        for unit in library.partitions:
            self.env.note_dependency(os.path.join(_src(self.env), unit.path))
        model = read_model(xml_dir, os.path.getmtime(os.path.join(xml_dir, "index.xml")))
        reference = Reference(library, model)
        lines = reference.lines()
        if reference.unparsed:
            logger.info("stereon-api %s: %d declarations shown verbatim (not parsable by the C++ domain)",
                        library.name, reference.unparsed)
        return self.parse(lines)


def setup(app):
    app.add_directive("stereon-library-graph", LibraryGraph)
    app.add_directive("stereon-partition-graph", PartitionGraph)
    app.add_directive("stereon-type-map", TypeMap)
    app.add_directive("stereon-api", ApiReference)
    return {"version": "1", "parallel_read_safe": True, "parallel_write_safe": True}
