# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
"""Render implementation notes beside the method family that explains their use."""

import ast
from pathlib import Path
import re
import subprocess
from urllib.parse import quote

from docutils import nodes
from sphinx.util.docutils import SphinxDirective


def notes(path):
    """Yield source line, declaration and prose from C++ documentation or Python docstrings."""
    source = path.read_text(encoding="utf-8")
    if path.suffix == ".py":
        for item in ast.walk(ast.parse(source)):
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                description = ast.get_docstring(item)
                if description:
                    yield item.lineno, item.name, description
        return
    if path.suffix not in {".cpp", ".cppm", ".cc", ".cxx", ".ixx", ".h", ".hpp", ".hxx", ".inl", ".ipp", ".c"}:
        return
    lines = source.splitlines()
    i = 0
    while i < len(lines):
        stripped = lines[i].lstrip()
        internal = stripped.startswith("// ") and lines[i][:1].isspace()
        if not stripped.startswith("///") and not internal:
            i += 1
            continue
        start = i
        block = []
        prefix = "// " if internal else "///"
        while i < len(lines) and lines[i].lstrip().startswith(prefix):
            block.append(lines[i].lstrip()[3:].lstrip())
            i += 1
        declaration = lines[i].strip() if i < len(lines) else ""
        if internal and len(block) == 1:
            continue
        # A template prefix alone does not identify the documented declaration.
        if declaration.startswith("template") and i + 1 < len(lines):
            declaration += " " + lines[i + 1].strip()
        yield start + 1, declaration, "\n".join(block)


class AlgorithmSources(SphinxDirective):
    has_content = True

    def run(self):
        """Render linked declarations and source notes for the named method family."""
        root = Path(self.env.srcdir).parent
        index = nodes.bullet_list()
        result = [index]
        for raw in self.content:
            name = raw.strip()
            if not name:
                continue
            path = (root / name).resolve()
            if not path.is_relative_to(root.resolve()):
                raise self.error(f"Source outside repository: {name}")
            self.env.note_dependency(str(path))
            uri = f"https://github.com/onurtuncer/Stereon/blob/{self.env.config.algorithm_revision}/{quote(name)}"
            link = nodes.reference("", name, refuri=uri)
            index += nodes.list_item("", nodes.paragraph("", "", link))
            entries = list(notes(path))
            if entries:
                result.append(nodes.rubric("", name))
                listing = nodes.definition_list()
                for line, decl, description in entries:
                    target = nodes.reference("", decl or f"Line {line}", refuri=f"{uri}#L{line}")
                    definition = nodes.definition()
                    for paragraph in re.split(r"\n\s*\n", description.strip()):
                        definition += nodes.paragraph("", " ".join(paragraph.splitlines()))
                    listing += nodes.definition_list_item("", nodes.term("", "", target), definition)
                result.append(listing)
        return result


def setup(app):
    """Register the directive and pin source links to the documentation revision."""
    root = Path(app.confdir).parent
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        revision = "main"
    app.add_config_value("algorithm_revision", revision, "env")
    app.add_directive("algorithm-sources", AlgorithmSources)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
