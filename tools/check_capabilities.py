#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Every geometric operation declares its derivative capability (ADR-0021).

ADR-0021's foundation is complete when "a new geometric operation cannot be
added without declaring its derivative capability, dependencies,
correspondence and validity assumptions (enforced by a registration check in
tests)". This is that check.

It reads the public interface of the geometry and topology libraries
(src/curve, src/surface, src/intersect, src/topo, src/build): every function declared
at namespace scope in an `export namespace Stereon { ... }` block of a module
interface unit. Each function name must be either

- declared in a capability catalogue, a `Declare("Name", ...)` entry of a
  module's Capabilities partition (src/<library>/Capabilities.cppm), which
  states what the operation's result depends on, how its samples
  correspond, under which assumptions a derivative would be valid, and
  which products exist (none yet: NotSupported, never a silent zero); or
- listed below as not an operation, with the reason: an evaluator or
  query of existing geometry, an accessor, a traversal, a validator or a
  utility that makes no new geometry.

and every catalogue entry must name an exported function. A new public
function therefore fails this check until someone decides which it is.
Overloads share a name and a declaration.

This is a script rather than a compile-time check because no compiler the
project uses can enumerate a module's exports (ADR-0016's reflection does
not reach namespaces yet), as the other source checks are.

Usage: check_capabilities.py [--root <repository root>] [--list] | --self-test
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

LIBRARIES = ("curve", "surface", "intersect", "topo", "build")

_EVALUATOR = "an evaluator of existing geometry: it makes nothing, and its derivative is the geometry's own"
_QUERY = "a query of existing geometry or of its representation"
_BASIS = "a B-spline basis utility over knots, which are structural (ADR-0021, section 4)"
_ACCESSOR = "a topology accessor or traversal: it reads a shape and makes nothing"
_CATALOGUE = "the capability catalogue itself"

# Public functions that are not geometric operations, by the reason they
# make no new geometry or shape whose derivative could be asked for. Adding a
# name here is a review decision, as adding a declaration is.
NOT_OPERATIONS: dict[str, str] = {
    # Evaluators and differential properties.
    **dict.fromkeys(
        (
            "Eval", "Deriv", "Derivatives", "TotalDerivatives", "Curvature", "UnitTangent",
            "FundamentalForms", "GaussianCurvature", "MeanCurvature", "PrincipalCurvatures",
            "UnitNormal", "ImplicitResidual", "RadialTurned", "Mapped",
        ),
        _EVALUATOR,
    ),
    # Queries of a geometry's representation.
    **dict.fromkeys(
        (
            "BoundingBox", "Domain", "Period", "Periods", "IsRational", "PoleCountU", "PoleCountV",
            "Apex", "Axis", "GeometryNameOf", "BuiltinIndexOf", "VisitBuiltin",
        ),
        _QUERY,
    ),
    **dict.fromkeys(
        (
            "BasisDerivatives", "BasisFunctions", "BasisFunctionsOfAllDegrees", "ExpandKnots",
            "FindSpan", "KnotMultiplicity",
        ),
        _BASIS,
    ),
    **dict.fromkeys(
        (
            "AllEntities", "CoedgesOf", "CurveOf", "CachedCurveOf", "CacheCertificateOf",
            "CacheBoundOf", "EdgeOf", "EdgesOf", "EndOf", "FaceOf", "FacesOf", "IsReversed",
            "KindOf", "LoopOf", "LoopsOf", "PCurveOf", "PositionOf", "PrimaryCoedgeOf",
            "RadialCoedgesOf", "RangeOf", "ShellOf", "ShellsOf", "SolidOf", "StartOf",
            "StoredEdgeUseSense", "SurfaceOf", "ToleranceOf", "HasCachedCurve",
        ),
        _ACCESSOR,
    ),
    "Describe": "a diagnostic text of a shape",
    "ToString": "the name of an enumerator",
    "ValidateTopology": "a validator: it reports on a shape and makes nothing (ADR-0018)",
    "RuleId": "the identifier of a validation rule",
    "RuleName": "the name of a validation rule",
    **dict.fromkeys(
        ("CurveCapabilities", "SurfaceCapabilities", "IntersectCapabilities", "TopoCapabilities",
         "KernelCapabilities", "BuildCapabilities"),
        _CATALOGUE,
    ),
}

DECLARE = re.compile(r'\bDeclare\(\s*"(?P<name>[A-Za-z_]\w*)"')


def blank_comments_and_strings(text: str) -> str:
    """Replace comments, string and character literals with spaces, keeping line breaks."""

    def blank(match: re.Match[str]) -> str:
        return re.sub(r"[^\n]", " ", match.group(0))

    pattern = re.compile(
        r'//[^\n]*|/\*.*?\*/|R"(?P<delim>[^(\s]*)\(.*?\)(?P=delim)"|"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'',
        re.DOTALL,
    )
    return pattern.sub(blank, text)


def export_blocks(code: str) -> list[str]:
    """The bodies of the `export namespace Stereon { ... }` blocks (not Stereon::Detail)."""
    blocks = []
    for match in re.finditer(r"\bexport\s+namespace\s+Stereon\s*\{", code):
        depth = 1
        start = match.end()
        i = start
        while i < len(code) and depth:
            if code[i] == "{":
                depth += 1
            elif code[i] == "}":
                depth -= 1
            i += 1
        blocks.append(code[start : i - 1])
    return blocks


def top_level_statements(body: str) -> list[str]:
    """The declarations at the block's own scope: up to a `;` at depth 0, or a
    braced body (a function's, a class's) that returns to depth 0."""
    statements = []
    current = []
    depth = 0
    for ch in body:
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                current.append(ch)
                statements.append("".join(current))
                current = []
                continue
        if depth == 0 and ch == ";":
            statements.append("".join(current))
            current = []
            continue
        current.append(ch)
    statements.append("".join(current))
    return [s.strip() for s in statements if s.strip()]


def skip_balanced(text: str, opening: str, closing: str) -> str:
    """The text after the balanced group it starts with, or "" if it does not close."""
    depth = 0
    for j, ch in enumerate(text):
        if ch == opening:
            depth += 1
        elif ch == closing:
            depth -= 1
            if depth == 0:
                return text[j + 1 :].lstrip()
    return ""


def strip_template_header(statement: str) -> str:
    """The statement without a leading `template <...>` and its `requires (...)` clause."""
    text = statement.lstrip()
    while text.startswith("template"):
        text = skip_balanced(text[text.index("<") :], "<", ">")
        while text.startswith("requires"):
            text = text[len("requires") :].lstrip()
            text = skip_balanced(text, "(", ")") if text.startswith("(") else ""
    return text


NOT_A_FUNCTION = re.compile(
    r"^(?:using|struct|class|enum|union|concept|namespace|static_assert|friend|typedef|export|template)\b"
)


def function_name(statement: str) -> str | None:
    """The name of the function a top-level statement declares, or None."""
    text = strip_template_header(statement)
    text = re.sub(r"\[\[.*?\]\]", " ", text)
    text = re.sub(r"\bdecltype\s*\(\s*auto\s*\)", "auto", text).strip()
    if not text or NOT_A_FUNCTION.match(text):
        return None
    paren = text.find("(")
    if paren < 0:
        return None
    head = text[:paren]
    if "=" in head or "operator" in head:
        return None  # a variable, or an operator, which is not named
    match = re.search(r"([A-Za-z_]\w*)\s*$", head)
    if not match or match.group(1).isupper():
        return None  # no name, or a macro
    return match.group(1)


def exported_functions(root: Path) -> dict[str, dict[str, set[str]]]:
    """Per library, every exported function name and the interface files that declare it."""
    names: dict[str, dict[str, set[str]]] = {library: {} for library in LIBRARIES}
    for library in LIBRARIES:
        for path in sorted((root / "src" / library).glob("*.cppm")):
            code = blank_comments_and_strings(path.read_text(encoding="utf-8"))
            for block in export_blocks(code):
                for statement in top_level_statements(block):
                    name = function_name(statement)
                    if name:
                        names[library].setdefault(name, set()).add(path.relative_to(root).as_posix())
    return names


def declared_operations(root: Path) -> dict[str, list[str]]:
    """Per library, the names its capability catalogue declares, in order."""
    declared: dict[str, list[str]] = {}
    for library in LIBRARIES:
        path = root / "src" / library / "Capabilities.cppm"
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        # Comments may quote a Declare call; only code counts.
        code = re.sub(r"//[^\n]*", "", text)
        declared[library] = [match.group("name") for match in DECLARE.finditer(code)]
    return declared


def check(root: Path, not_operations: dict[str, str]) -> list[str]:
    """Every violation, as messages."""
    exported = exported_functions(root)
    declared = declared_operations(root)
    errors = []
    for library in LIBRARIES:
        catalogue = f"src/{library}/Capabilities.cppm"
        names = declared[library]
        for name in sorted({n for n in names if names.count(n) > 1}):
            errors.append(f"{catalogue}: {name} is declared more than once")
        for name in sorted(set(names) - set(exported[library])):
            errors.append(f"{catalogue}: {name} is declared but the {library} library exports no such function")
        for name in sorted(set(names) & set(not_operations)):
            errors.append(f"{catalogue}: {name} is declared and also listed as not an operation")
        for name, files in sorted(exported[library].items()):
            if name not in names and name not in not_operations:
                errors.append(
                    f"{sorted(files)[0]}: {name} is exported but has no capability declaration in "
                    f"{catalogue} and is not listed in NOT_OPERATIONS (ADR-0021, section 5)"
                )
    everywhere = set().union(*(set(exported[library]) for library in LIBRARIES))
    for name in sorted(set(not_operations) - everywhere):
        errors.append(f"tools/check_capabilities.py: NOT_OPERATIONS lists {name}, which no library exports")
    return errors


SELF_TEST_FILES = {
    "src/curve/Line.cppm": """
export module Stereon.Curve:Line;
export namespace Stereon {
struct Line { double Eval(double t) const; };
[[nodiscard]] Result<Line> MakeLine(const Point3& origin, const Vector3& direction);
template <class C>
    requires (Curve<C>)
[[nodiscard]] Point3 Eval(const C& curve, double t);
[[nodiscard]] Line Reversed(const Line& line);
[[nodiscard]] Line Reversed(const Line& line, int unused);
inline constexpr int Count = Size(3);
} // namespace Stereon
namespace Stereon::Detail { void Hidden(); }
""",
    "src/curve/Capabilities.cppm": """
export module Stereon.Curve:Capabilities;
export namespace Stereon {
// Declare("Commented", ...) is not a declaration.
inline constexpr OperationCapability CurveCatalogue[] = {
    Declare("MakeLine", {}, {}, "a"),
    Declare("Reversed", {}, {}, "b"),
};
}
""",
}


def self_test() -> int:
    import tempfile

    failures = []
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for relative, text in SELF_TEST_FILES.items():
            (root / relative).parent.mkdir(parents=True, exist_ok=True)
            (root / relative).write_text(text, encoding="utf-8")
        exported = exported_functions(root)["curve"]
        if set(exported) != {"MakeLine", "Eval", "Reversed"}:
            failures.append(f"exported functions: {sorted(exported)}")
        cases = [
            ("a complete catalogue", {"Eval": "x"}, []),
            ("an undeclared function", {}, ["Eval is exported"]),
            ("a stale exception", {"Eval": "x", "Gone": "x"}, ["NOT_OPERATIONS lists Gone"]),
            ("a contradiction", {"Eval": "x", "MakeLine": "x"}, ["MakeLine is declared and also"]),
        ]
        for label, not_operations, expected in cases:
            errors = check(root, not_operations)
            if len(errors) != len(expected) or not all(e in m for e, m in zip(expected, errors)):
                failures.append(f"{label}: {errors}")
        catalogue = root / "src/curve/Capabilities.cppm"
        catalogue.write_text(
            catalogue.read_text(encoding="utf-8").replace('"b"),', '"b"),\n    Declare("Reversed", {}, {}, "c"),\n    Declare("Missing", {}, {}, "d"),'),
            encoding="utf-8",
        )
        errors = check(root, {"Eval": "x"})
        for expected in ("Reversed is declared more than once", "Missing is declared but"):
            if not any(expected in error for error in errors):
                failures.append(f"expected '{expected}' in {errors}")
    for failure in failures:
        print("self-test: " + failure, file=sys.stderr)
    if failures:
        return 1
    print(f"check_capabilities self-test: {len(cases) + 2} cases pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--list", action="store_true", help="print every exported function and exit")
    parser.add_argument("--self-test", action="store_true", help="run the built-in cases and exit")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if args.list:
        for library, names in exported_functions(args.root).items():
            for name, files in sorted(names.items()):
                print(f"{library} {name}: {', '.join(sorted(files))}")
        return 0
    errors = check(args.root, NOT_OPERATIONS)
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        print(f"check_capabilities: {len(errors)} violation(s)", file=sys.stderr)
        return 1
    declared = declared_operations(args.root)
    print(
        f"check_capabilities: {sum(len(names) for names in declared.values())} declared operations, "
        f"{len(NOT_OPERATIONS)} listed non-operations, no violations"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())