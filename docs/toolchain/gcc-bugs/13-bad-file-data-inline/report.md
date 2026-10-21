# Not reproduced on its own; probably the same bug as 06. Do not file separately.

Rows: the second and third "Bad file data" entries ("failed to load pendings
for 'Stereon::BasicVec3'" from inline functions in surface partitions).

## What was tried
- The failing code (inline `RequireFinite<double>`/`RequireInDomain` calls in
  `Stereon.Surface:Erased`, inline checked factories, inline
  `BoundingBox(const Torus&)` indexing `BasicVec3<double>::operator[]`) was
  never committed: a9636c1 (#32) and 6deec29 (#39) are squash merges with the
  workarounds in place.
- Synthetic variants of the 06 testcase on CE GCC 16.2 (A, B, C modules, B
  instantiating `V<int>`):
  - C's inline function calls a non-friend member `V<int>::operator[]`: builds.
  - C's inline function calls a function template from A returning `V<T>`: builds.
  - `operator-` as a member instead of a hidden friend: builds.
  - hidden friend returning `T` instead of `V<T>`: builds.
  Only the hidden-friend form fails. All of Stereon's `Vec3`/`Point3`
  arithmetic is hidden friends, and the failing inline functions almost
  certainly used them besides `operator[]`.
- With the importer also evaluating `p - p`, the CMI loads but the friend is
  emitted twice (assembler "already defined"), noted in the 06 report.

## To confirm locally
After an upstream fix of 06, move one surface factory (e.g. `MakeCylinder`)
back inline into its partition and rebuild a clean tree; if "failed to load
pendings" still appears, reduce it again.
