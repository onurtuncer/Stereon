# Not reproduced on its own; probably the same bug as 06. Do not file separately.

Row: "failed to read compiled module cluster: Bad file data" (Surface:Cylinder
building Stereon.Curve's BasicCircle).

## What was tried
- The failing pre-workaround code (`BasicCylinder` building `BasicCircle<T, 3>`
  and calling its `Eval`, `Deriv`, `Project`, plus `tests/support/Formatting.hpp`'s
  partial specialisations) was never committed: PR #26 (cf1dbaa) was squashed
  with the workaround in place, and `git log --all -S BasicCircle --
  src/surface/Cylinder.cppm` is empty.
- Building the real core/curve/surface modules of 74fc6cc^ on Compiler
  Explorer (16.2, Linux) produced this row's error text, raised in
  Stereon.Surface reading its Cylinder partition, and reduced to the
  three-module hidden-friend pattern of ../06-bad-file-data-hidden-friend.
- Synthetic variants with a partial specialisation on an imported template in
  the importing TU, and a constexpr function template using the hidden friend
  in a `static_assert` in the third module: neither failed.

## To confirm locally
Temporarily make `BasicCylinder`'s `Eval` build a `BasicCircle<T, 3>` and call
its `Eval`, rebuild a clean tree. If it still fails after 06 is fixed
upstream, reduce it again.
