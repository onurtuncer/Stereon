# Ellipse-plane CI replay

Source: PR #272, macOS Actions run 43273260370, job 131374874702.
RapidCheck seed: 17041817369036985825. The shrunk draws were
`0,-16,-40,-40,-13,15,4,2,0,0,-40,-40,-13`.

The ellipse center is (-40/32,-40/32,-13/32), major radius 1.5,
minor radius 1.5*4/10, X=(1,0,0), Y=(0,cos(phi),sin(phi)), with
phi=-16*(2*pi)/64. The plane has the same center and axes X and (0,1,0).
The failure was `Ellipse x Plane: unresolved root existence or multiplicity
in requested domain` in the linearization property.

Independent oracle: the plane equation reduces to minor*sin(phi)*sin(t)=0.
Both intersections are transverse, at center +/- major*X. The deterministic
test requires successful decided crossings, both spatial roots, parameter-domain
membership and evaluation on both supports. A whole closed turn may retain
both parameter representatives of the seam, so two or three entries are accepted,
but no third spatial root. No unresolved result is discarded or accepted.

The campaign tests powers-of-two scales 1/4, 1, 4 with scaled linear resolution
and tolerance ceiling, x translations 0 and 8, and curve reversal. The default
angular policy is unchanged and execution is serial. It uses the production
default curve domain, including its seam, as did the failing property.

This is stacked on #283's conic arc-cut fix. The previous CI failure is observed;
passing the new deterministic replay on #283 remains to be established by CI.
The ledger therefore stays `fixture captured`, not verified. The broader
linearization property discards IntersectionFailed on #283; this fixed test
specifically prevents that discard from hiding this known transverse fixture.

Replay after building `stereon_intersect_tests`:

```sh
ctest --preset gcc16-debug -R 'intersect: MR macOS ellipse plane shrunk linearization fixture resolves' --output-on-failure
```
