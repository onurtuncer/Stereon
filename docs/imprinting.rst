..
   SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
   SPDX-License-Identifier: MPL-2.0

Imprinting
==========

Imprinting writes the intersection of two solids into each of them. It is
steps 3 and 4 of the Boolean pipeline (:doc:`modeling`): every face of one
body is met with every face of the other, each meeting is restricted to the
two faces (*clipping*), and the pieces that lie in a face are written into
it as new edges, splitting the face (*imprinting*). After it, every curve
along which the boundaries of the two solids meet is an edge of both, and
each body's faces fall into regions that lie wholly inside or wholly outside
the other solid. Classification (:doc:`modeling`, *Classification*) and
assembly (*Assembly*) work on those regions and never look at an
intersection again.

This chapter describes the algorithm as ``src/boolean/Clip.cpp``,
``src/boolean/Imprint.cpp`` and ``src/boolean/Bodies.cpp`` implement it, in
the order the work is done: which face pairs are met, how a branch is
clipped to a face, how the pieces become a network of curves on each face,
how faces that share a plane are handled, how a network is written into the
topology through the Euler operators, how the two bodies' edges are paired
afterwards, and what is decided exactly, what within a tolerance, what is
deterministic and what fails. The governing decisions are ADR-0001
(tolerances), ADR-0002 (exact predicates), ADR-0009 (edges and their
pcurves), ADR-0010 (surface intersection), ADR-0027 (controlled topology
edits), ADR-0029 (intersection outcomes) and ADR-0030 (lineage).

.. contents:: On this page
   :local:
   :depth: 2

Vocabulary
----------

``branch``
   One connected component of the intersection of two surfaces
   (``IntersectionBranch``): a 3D curve, its pcurve on each surface on one
   shared parameter, a parameter range, the contact type (``Crossing``,
   ``Tangential``, ``Overlap``), the side transition across it
   (``OppositeSide``, ``SameSide``, ``Unresolved``) and a certified bound on
   its deviation from both surfaces (ADR-0009, section 4).

``pcurve``
   A curve in a face's own parameter space :math:`(u, v)`. The surface
   :math:`S` maps it to 3D: an edge's points are :math:`S(p(t))`. Every use
   of an edge by a face carries one, and all of them run on the edge's
   parameter (ADR-0009).

``clip``
   A branch restricted to one face (``FaceClip``): where it crosses the
   face's boundary, the stretches of it inside the face, and the stretches
   where it runs along the boundary.

``piece``
   A stretch of a branch inside a face (``FacePiece``): a parameter range,
   and at each end the boundary crossing it ends at, or nothing for a free
   end.

``crossing``
   A point where a branch meets a face's boundary (``LoopCrossing``): the
   coedge whose edge it meets, the parameter on that edge and on the
   branch, the point in the face's parameter space, the contact, whether it
   was decided exactly, and the loop vertex it is at when it lies within that
   vertex's tolerance.

``junction``
   A network vertex not on the face's boundary: where a branch passes from
   one face of the other body to the next, or where it ends inside the face.

``network``
   What one face receives (``FaceNetwork``): vertices (crossings and
   junctions) and edges (pieces) between them.

``imprint``
   The network written into the face (``Imprint``): an edge per piece, the
   vertices made for the network's vertices, and the vertices the splits
   of crossed edges made.

``pair``
   An edge of A and an edge of B that are one curve of the result
   (``SharedEdge``): a piece of the intersection, imprinted into both
   bodies, or two edges that coincide on faces sharing a plane.

The pipeline
------------

``ImprintBodies(a, b, ctx)`` runs the whole of this chapter:

#. **Gather** the networks. Every face of ``a`` is met with every face of
   ``b`` whose bounding box it meets, in input order: the surfaces are
   intersected over the faces' parameter boxes, each branch is clipped to
   both faces, and the stretches inside both are added to both faces'
   networks as a piece of the intersection. Faces whose surfaces coincide
   are handled as described in `Coincident faces`_.
#. **Imprint** each body: the body is opened in a ``ShapeBuilder``, every
   network is written into it inside one transaction, and the result is
   committed (``Commit``), which validates it and derives its History.
#. **Pair** the edges: each piece of the intersection, now an edge in
   each body, becomes a ``SharedEdge`` with the branch's side transition;
   then the edges of the two bodies that coincide on shared planes are
   paired as well.

The result, ``ImprintedBodies``, holds both committed bodies (each an
``OperationResult`` with its History and its topology events) and the
pairs.

Candidate face pairs
--------------------

Two faces are met only when their bounding boxes meet. A face's box is the
hull of its vertices, each widened by its own tolerance, and of the boxes
of its edges' own curves where those are bounded (a circle's is; a line's
is not, and its vertices already bound it). The test is conservative: a
box that misses the other proves the faces apart, and nothing else is
skipped.

The surfaces are then intersected (``Intersect``) over the faces'
*parameter boxes*, not over the whole surfaces. A face's parameter box is
the hull of its pcurves over their ranges (a line pcurve by its ends, any
other by its bounding box), widened by the Context's resolution so that a
branch reaching the face's boundary is not cut short by rounding. The
intersection is complete over those boxes (ADR-0010): every branch,
isolated point and overlap region of the two surfaces inside them is
reported, with its contact type and side transition.

Clipping a branch to a face
---------------------------

``ClipToFace(shape, face, branch, side, ctx)`` restricts one branch to one
face, working entirely in the face's parameter space on the branch's
pcurve for that face (``branch.OnA`` or ``branch.OnB``). The branch's range
must be bounded.

Crossings
~~~~~~~~~

The face's loops are read as a list of *loop edges*: each coedge's pcurve
over its edge's range, with the edge's vertices at the range's ends. The
branch's pcurve is intersected with every loop edge's pcurve
(``Intersect`` on plane curves, :doc:`intersection`):

- every intersection point becomes a crossing, with the contact type the
  plane curve intersection reports and whether it decided the point
  exactly. Lines, circles and ellipses are decided exactly on the input
  doubles (ADR-0002); a B-spline pcurve's roots are bracketed to rounding
  and the clip is then marked ``Decided = false``;
- every overlap, a stretch where the two pcurves coincide, is recorded as
  a stretch of the branch *on the boundary* (``OnBoundary``).

**Snapping to vertices.** A crossing whose 3D point lies within a loop
vertex's tolerance of that vertex is *at* the vertex (ADR-0001): its edge
parameter is snapped to that end of the edge's range and the crossing
records the vertex. Two edges meeting at the vertex both report the
crossing; ordering keeps one.

**End crossings.** A branch end lying within an edge's tolerance of that
edge is a crossing of the edge, whether or not the plane curve intersection
found a root there: a B-spline pcurve's root at the very end of its domain
can be missed by subdivision. The tolerance is the larger of the edge's and
the branch's certified bound. Such a crossing is bracketed, not decided,
and it too snaps to a vertex.

**Ordering.** Crossings are sorted by their parameter on the branch, and a
vertex reported by both of its edges at one parameter is kept once.

Pieces
~~~~~~

The branch's range is cut at its crossings into stretches, and each
stretch's middle point decides whether the stretch lies inside the face:

- only crossings of contact type ``Crossing`` cut the branch; a tangential
  touch of the boundary does not change sides;
- a crossing within :math:`10^{-9}` of the branch's range of the previous
  cut, or of either end of the branch, cuts no stretch of its own: a closed
  branch starting on a seam meets that seam at both of its ends, found to
  rounding. A crossing at the branch's own end *is* that end, recorded as
  the piece's crossing there;
- a stretch whose middle lies on an ``OnBoundary`` stretch is along the
  boundary, neither inside nor outside, and is not a piece;
- otherwise the middle point's side is the parity of a ray's crossings of
  the loops (below), and a stretch inside becomes a piece, its ends the
  crossings that cut it (nothing at a free end of the branch).

The inside test
~~~~~~~~~~~~~~~

A point :math:`q` of the parameter space lies inside the face when a ray
from it crosses the face's loops an odd number of times. The ray is a line
pcurve from :math:`q` at an angle of :math:`k \varphi \pi` for
:math:`k = 0, 1, \dots, 7`, :math:`\varphi` the golden ratio: irrational
fractions of a turn, so that no two directions align with the loops' edges
in the same way. Each ray is intersected exactly with every loop edge's
pcurve, and a direction is abandoned, the next tried, when the ray

- runs along a loop edge (an overlap),
- meets an edge at an end of its range, a loop vertex, where the
  neighbouring edge reports the same point and the count would be
  ambiguous, or
- starts on an edge.

A tangential touch of an edge does not count as a crossing. When all eight
directions are ambiguous the clip fails with ``IntersectionFailed``: the
side cannot be settled, and nothing is guessed.

``FaceContains(shape, face, uv, ctx)`` uses the same test for a single
point, after one more check: a point within an edge's tolerance of that
edge, measured in 3D between :math:`S(uv)` and the nearest point of the
edge's pcurve, is *on the boundary*, and the answer is nothing rather than
a side (ADR-0001).

Pieces of the intersection
--------------------------

A branch of two faces' surfaces is clipped to both faces, ``fa`` of A on
side A and ``fb`` of B on side B. A piece of the intersection is a stretch
of the branch inside *both* faces: every pair of an inside piece of each
clip contributes the intersection of their ranges, unless it is empty or
narrower than :math:`10^{-9}` of the branch's range. A stretch along one
face's boundary is inside neither, so a branch running along an edge of
one face adds nothing; the edge is already there.

Each such stretch is added to both faces' networks as an edge from a start
vertex to an end vertex:

- an end at which the stretch ends at a clip's crossing is that crossing's
  vertex, made once per crossing;
- any other end is a **junction**, where the branch passes from this face
  of the other body to the next one or ends inside the face. Junctions are
  merged by position: a junction within the larger of the two tolerances
  of an existing junction of the same face's network *is* that junction.
  Crossings and junctions never merge with each other;
- a stretch covering the whole branch with no crossing at either end is a
  loop, its start and end one vertex, only when the branch is a closed
  curve: its ends within the tolerance of each other (``IsLoop``). A
  straight stretch lying wholly inside a face also covers its whole branch
  without a crossing, and is no loop.

The network's vertices take the tolerance of the piece, the larger of the
Context's resolution and the branch's certified bound. Each edge records
the sense of the boundary coedge at each crossing end, which matters on a
seam (below). The pair of edges, one in each network, is remembered with
the branch's side transition, to become a ``SharedEdge`` once both bodies
are imprinted.

Coincident faces
----------------

When two faces' surfaces coincide, the surface intersection reports an
*overlap region* (``OverlapRegion``, ADR-0029, section 1) rather than
branches: its domain on each surface, the affine correspondence of their
parameters, whether their normals agree, how the coincidence was
established (exactly in the input doubles, or tolerantly: exactly parallel
and within the pair tolerance, with the distance recorded), and its own
boundary branches, which bound the *parameter boxes* the surfaces were
intersected over, not the faces. Coincident supports do not make
overlapping faces: the faces' trims decide that.

The prototype imprints coincident **planar** faces. Coincident curved
faces, two cylinders sharing an axis and radius, fail with
``NotSupported`` naming both faces.

Boundaries imprinted on each other
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a coincident pair ``fa`` and ``fb``, each face's boundary is imprinted
on the other. Every coedge of every loop of ``fb`` gives a branch on
``fa``'s plane:

- its 3D curve is the edge's own curve (a line, circle or ellipse; an edge
  without one of its own fails with ``NotSupported``);
- its pcurve on ``fa`` is that curve expressed in the plane's frame, exactly
  and on the curve's own parameter. With the plane
  :math:`P(u, v) = O + u X + v Y` and its axes orthonormal,

  .. math::

     \mathrm{uv}(p) = \big((p - O)\cdot X,\; (p - O)\cdot Y\big),
     \qquad \mathrm{dir}(d) = (d\cdot X,\; d\cdot Y),

  a line :math:`o + t d` maps to the line :math:`\mathrm{uv}(o) + t\,
  \mathrm{dir}(d)`, and a circle or ellipse of centre :math:`c` and axes
  :math:`X', Y'` to the plane curve of centre :math:`\mathrm{uv}(c)` and
  axes :math:`\mathrm{dir}(X'), \mathrm{dir}(Y')`, with the same radii;
- its second pcurve is the coedge's own pcurve on ``fb``;
- its range is the edge's range, its contact ``Overlap``, its transition
  ``Unresolved``, and its bound the edge's tolerance plus the distance
  between the coincident planes.

The branch is clipped to ``fa`` exactly as an intersection branch is, and
each stretch inside ``fa`` is added to ``fa``'s network *alone*
(``AddOneSided``): ``fb`` already has the edge. A stretch along ``fa``'s own
boundary, where the two faces share a stretch of edge, is not added; the
pairing below joins those. The same is done with ``fa``'s boundary on
``fb``. The overlap's own boundary branches are skipped: they bound the
parameter boxes, and the faces' boundaries take their place.

A face of one body meeting a coplanar face of the other also meets the
other body's faces *across* that plane, along the edges of the coplanar
face. Those transverse branches lie on the coplanar face's boundary, so on
that side the clip reports them along the boundary and they add no piece:
the boundary imprint has already drawn them, once.

Pairing coincident edges
~~~~~~~~~~~~~~~~~~~~~~~~~

After both bodies are imprinted, every edge of A that lies on one of the
coincident planes is paired with the edge of B it coincides with. An edge
lies on a plane when its ends and its middle (on its own curve, or for an
edge without one, on its primary side :math:`S(p(t))`) are within the
edge's tolerance of the plane. Two such edges are one when their ends
coincide, either way round, and their middles coincide, all within the
larger of their tolerances. Each edge of B pairs at most once, with the
first edge of A that matches, and edges already paired as pieces of the
intersection are left out. The pairs are ``SharedEdge``\ s with transition
``Unresolved``: the boundary of an overlap separates regions without saying
which side either is on, and classification decides them on their own.

Imprinting a network
--------------------

``ImprintNetworks(builder, shape, networks, ctx)`` writes a body's networks
into a builder opened on the body, inside the caller's transaction, and
returns what each network became. Everything is done through the
controlled operators of ADR-0027, so every intermediate state is a valid
topology except where an operator's contract allows otherwise, and a
failure anywhere rolls the whole body back.

Splitting the crossed edges
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, every edge that some network's vertices cross is split at the
crossings:

#. the crossings of all of the body's networks are grouped by the edge they
   cross; a crossing *at a loop vertex* needs no split, its vertex is the
   existing one;
#. each edge's crossings are sorted by their parameter on the edge;
#. the edge is split at the first crossing with SEMV
   (``SplitEdgeMakeVertex``), which keeps the lower stretch in the edge and
   makes the upper one a new edge and the split point a new vertex at the
   crossing's position, its tolerance the larger of the edge's and the
   crossing's;
#. every following crossing of the same edge within that tolerance of the
   split point is the same vertex. An edge two faces share is crossed by
   both faces' networks at the same point; this is how they come to meet
   at one vertex;
#. the next crossing splits the upper stretch, and so on.

The SEMV events (an edge ``Split`` into its upper stretch, a vertex
``Created``) are recorded in the builder's topology events, which the
assembly later reads to trace each piece of an edge back to the edge it
came from.

Writing the pieces
~~~~~~~~~~~~~~~~~~

Then each network's pieces are written in rounds, until every piece is
done:

- **Between two made vertices.** A piece whose start and end vertices both
  exist is written between them. The coedges of the face's loops that end
  at each vertex are found (``UseEndingAt``); when both lie on one loop the
  piece splits the face (MEF, ``MakeEdgeFace``): the loop is cut in two and
  a new face takes the part from the start coedge round to the end coedge.
  When they lie on two loops of one face, the piece bridges them (MEKR,
  ``MakeEdgeKillRing``) and the two loops become one.
- **From a made vertex to a new one.** A piece whose start exists but whose
  end does not is grown as a strut (MEV, ``MakeEdgeVertex``): an edge from
  the existing vertex to a new vertex at the end's position. After each
  strut the round starts again, since the new vertex may complete other
  pieces first.
- **A component touching no boundary.** When no piece can be written, some
  piece has neither end made: a component of the network lies inside the
  face without touching its boundary (a hole's outline, say). It is *hung*
  on a helper strut: a ``Spatial`` straight edge from a vertex of the face's
  first loop to the component's first unmade vertex. The component is then
  written from it, its last piece closing it with MEF. When every piece is
  done, each helper strut is removed (KEMR, ``KillEdgeMakeRing``): the
  component, cut free, is left as a separate loop of the face, a *ring*,
  and the face MEF cut off is its *disc*.

A network whose pieces can never all be written (a piece runs into the
face from a junction no other piece reaches) fails with
``IntersectionFailed``; so does a piece whose two ends are not on the
boundary of one face.

Rings and their orientation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ring left in the original face must be a *hole*. Loops run
counter-clockwise about their face's normal, so a hole runs clockwise about
:math:`S_u \times S_v` in a face oriented with its surface, and
counter-clockwise in a reversed face. MEF has no loop order to decide which
of its two loops stays in the original face: for a single closed piece it
makes an edge from a coedge to itself, and for several pieces it joins two
struts. It may leave the disc's loop in the original face and the hole's
in the new one.

``OrientRing`` checks every ring after its helper is removed. KEMR keeps
the helper's anchor side in the original loop and returns the other, and
which side the component lay on depends on the way the strut was hung, so
either of the two loops it leaves may be the ring: the ring is the one made
of the network's pieces (``RingAmong``). ``OrientRing`` samples the pcurves
of the ring's uses, in the order and direction the loop runs (64 samples
each), and takes the polygon's signed area in the face's parameter space,

.. math::

   A = \tfrac12 \sum_i \left(u_i v_{i+1} - u_{i+1} v_i\right),

positive when the ring turns counter-clockwise. When the ring turns the
wrong way for the face, the loops change faces: the disc's loop moves into
the ring's face (KFMRH, ``KillFaceMakeRingHole``), which keeps the outer
boundary whether it is the original face or the one MEF made, and the ring
moves out into a new face (MFKRH, ``MakeFaceKillRingHole``), which then is
the face of the piece that closed the component. Before the ring was told
from the face's own boundary this check never ran, and a ring imprinted on
a reversed face, a cavity's, kept the turn of an outer loop.

Seams
~~~~~

A periodic surface's seam is one edge used twice by its face, once from
either side of the chart: on a cylinder, at :math:`u = 0` and at
:math:`u = 2\pi`. When a piece ends on a seam, both uses of the edge end at
the vertex the crossing made, and the clip names one of them without
telling which side the piece arrives from. The imprint takes the use the
piece's *own* pcurve reaches: of the face's uses of the crossed edge, the
one whose pcurve, at the crossing's parameter on the edge, lies nearest the
piece's end in the face's parameter space (``SenseAt``). Writing the piece
against the other use would leave the face's loop discontinuous in
parameter space, a jump of a whole period, which the winding number of
classification (``WindingNumberOf``) refuses as ``InvalidResult``.

The edges the pieces become
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A piece becomes an ``Exact`` edge, carrying the branch's 3D curve as its
own, when that curve is a line, a circle or an ellipse. Any other branch,
the quartic of two skew cylinders or a fitted B-spline, becomes a ``Trim``
edge (ADR-0009, section 8): the 3D curve is derived from the face's surface
through the pcurve, and a cached curve is built from the primary side
(``CacheEdgeCurve``) and certified within the edge's tolerance where
interval evaluation covers the types. Both uses of the new edge, one on
each side of it in the face, take the piece's pcurve. The edge's tolerance,
like its vertices', is the larger of the Context's resolution and the
branch's certified bound.

After the imprint
-----------------

Each body's builder is committed on its own (``Commit``): the result is
validated, Structural and, in checked builds or when the Context asks,
Geometric, and its History is derived from the body it was opened on. Its
topology events, the SEMV and MEF splits among them, are published with it.
A result that fails validation is not published: the imprint fails.

The pieces of the intersection are then resolved to the edges the two
imprints made of them, and paired with their branches' side transitions;
the coincident edges are paired as above. A piece of the intersection is
always one curve over one range in both bodies, so its two edges run the
same way; a coincident pair may run either way, and the assembly handles
both.

What is exact, what is tolerant
-------------------------------

Every *decision* the imprint takes is one of three kinds:

- **Exact, on the input doubles** (ADR-0002): whether a line, circle or
  ellipse pcurve meets a loop's line, circle or ellipse pcurve, and where;
  whether a ray crosses, touches or runs along an edge; the order of
  crossings along a branch. These never depend on rounding, evaluation
  order or thread count.
- **Within a tolerance** (ADR-0001): whether a crossing is at a loop vertex
  (the vertex's tolerance), whether a branch end lies on an edge (the
  edge's tolerance or the branch's bound), whether a point is on a face's
  boundary (the edge's tolerance), whether two junctions or two crossings
  of one edge are one vertex (the larger tolerance), whether a branch is a
  closed curve, and whether two edges coincide on a shared plane. No
  tolerance is a literal: each is an entity's own, the branch's certified
  bound or the Context's resolution, and the epsilon lint
  (``tools/check_epsilon.py``) holds ``src/boolean/`` to that.
- **A parameter gap.** Two crossings closer along a branch than
  :math:`10^{-9}` of the branch's range are one cut, and a common stretch of
  two clips narrower than that is no piece. This is not a geometric
  tolerance: it absorbs the rounding of roots found twice, at the two ends
  of a closed branch or from the two edges at a vertex, and it is marked in
  the source as a parameter gap.

A crossing the plane curve intersection could only bracket (a B-spline
pcurve's) clears the clip's ``Decided`` flag. The imprint still uses it;
the flag tells a caller the result rests on a bracketed root.

Determinism
-----------

The imprint gives the same result for the same input on every run and
thread count (ADR-0007, G1). Face pairs are met in input order, A's faces
outermost; crossings are sorted by parameter; rays turn through a fixed
sequence of directions; junctions merge with the first match in creation
order; crossed edges are split in parameter order; network pieces are
written in a fixed round order; and coincident edges pair with the first
match. Nothing depends on a hash, a pointer or completion order.

Failures
--------

Every failure is a ``KernelError``, and the imprint, done inside one
transaction per body, publishes nothing when it fails.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Code
     - When
   * - ``InvalidInput``
     - A face's coedge has no pcurve; a branch's range is unbounded or
       empty; a network's face is not live in the builder; the face has no
       loop to hang a component on.
   * - ``IntersectionFailed``
     - Every one of the eight rays from a point meets a loop vertex, so its
       side cannot be settled; a piece's ends are not on the boundary of
       one face; a piece runs into the face from a junction no other piece
       reaches, so the network never closes.
   * - ``NotSupported``
     - Coincident curved faces; a coincident face bounded by an edge with
       no line, circle or ellipse of its own.
   * - ``InvalidResult``
     - An imprinted body fails Structural or Geometric validation at
       ``Commit`` (a kernel bug, reported rather than published).
   * - (from the operators)
     - Any refusal of an Euler operator's contract (ADR-0027), naming the
       operator and the entities.

Worked examples
---------------

Two boxes overlapping at a corner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`A = [0, 10]^3` and :math:`B = [5, 15]^3`. Three faces of each body
meet three faces of the other, each pair of perpendicular planes in a line.
A's face :math:`x = 10` meets B's face :math:`y = 5` along the line
:math:`x = 10, y = 5`; on A's face it runs inside from :math:`z = 0` to
:math:`z = 10`, on B's face from :math:`z = 5` to :math:`z = 15`, so the
piece of the intersection is :math:`z \in [5, 10]`. Its end at
:math:`z = 10` is a crossing of A's top edge, its end at :math:`z = 5` a
crossing of B's bottom edge; on each face the other end is a junction,
where the intersection turns onto the next face. Each of A's three faces
receives two pieces meeting at a junction, an L from one edge to another:
the crossed edges are split, the first piece is a strut from its crossing
to the junction and the second closes the face with MEF, cutting off a
:math:`5 \times 5` square. Six pieces of the intersection, six pairs; each
body ends with 9 faces in two regions, the three squares and the rest.

A cylinder through a box
~~~~~~~~~~~~~~~~~~~~~~~~

:math:`A = [0, 10]^3`, and B a cylinder of radius 3 about the line
:math:`x = y = 5` from :math:`z = -2` to :math:`z = 13`. A's top and bottom
faces each meet B's lateral face in a full circle, a closed branch that
touches no edge of A's square: a component touching no boundary. It is
hung on a helper strut from a corner of the square, closed with MEF from
its coedge to itself, the helper removed, and the ring oriented as a hole:
A's top face gets a circular hole and a new disc face. On B's lateral face
the same circle runs from the seam all the way round to the seam: its
crossings are the seam's, at :math:`u = 0` and :math:`u = 2\pi`, and both
ends take the use of the seam their own pcurve reaches, so the face's
loops stay continuous. B's lateral face is cut into three bands.

A smaller box against a face
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`A = [0, 10]^3` and :math:`B = [10, 20] \times [2, 8]^2`. A's face
:math:`x = 10` and B's face :math:`x = 10` coincide, with opposite normals.
B's square boundary is imprinted on A's face: four straight branches,
each wholly inside A's face, so four pieces between four junctions, each a
segment and not a loop. They form a component touching no boundary: hung
on a helper, written as three struts and an MEF, the helper removed, and
the four-piece ring oriented as a hole by its polygon's signed area. A's
face becomes a frame and a square. A's boundary lies wholly outside B's
face and imprints nothing on it. B's side faces meet A's face along the
square's sides, which lie on B's boundary, so they add no piece. The
pairing then pairs the four sides of A's square with B's four edges, which
run either way. Classification finds A's square ``OnOpposite`` B's face and
the frame outside.

Limits of the prototype
-----------------------

- Coincident faces are imprinted for planes only.
- A coincident face's boundary edges must carry a line, circle or ellipse
  of their own; a ``Boundary`` or ``Trim`` edge is refused.
- Coincident edges are paired by their ends and middle; two different arcs
  of one circle sharing both ends would be taken for one, a case two
  coincident planar faces of valid solids do not produce.
- The inside test turns its ray through eight directions; a point every
  one of which meets a vertex fails rather than being decided.
