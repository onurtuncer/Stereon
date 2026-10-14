# ADR-0031: Specify Attribute Propagation and Conflict Policies

- **Status:** Proposed
- **Date:** 2026-10-14
- **Owner:** Onur Tuncer
- **Phase:** 3 (v0.3 lineage/attributes); Booleans in v0.4–v0.6 and tagged meshing in v0.7
- **Extends:** ADR-0006, ADR-0011, ADR-0014, ADR-0026, ADR-0030

## Context

A geometric edit can leave a structurally valid model with incorrect engineering meaning. Copying a colour to both halves of a face is reasonable; duplicating a total force or a unique persistent-ID value may not be. Merging two faces tagged as different boundary conditions requires an explicit choice, and preserving a threaded-hole label after changing its supporting surface may be invalid.

ADR-0011 already requires history-driven propagation and conflict warnings. This ADR specifies the per-key semantics, split/merge distinctions, callback behavior and publication guarantees needed to make that policy usable.

Stroud §8.3, particularly §§8.3.2–8.3.4, describes copying information after splits, checking clashes after merges and revalidating features after shape changes. The application policy interface below is Stereon's design rather than a generic semantic interpretation of arbitrary user data.

## Options

1. **Blindly copy every attribute to every descendant.** Simple but can duplicate unique identifiers, conserved quantities and obsolete feature information.
2. **Drop all metadata on modification.** Avoids stale values but destroys useful labels and forces geometric matching in every application.
3. **Propagate through semantic lineage with a per-key transfer/conflict policy and structured diagnostics.** More work, but passive data has a useful default and engineering data can require stronger handling.

## Decision

Adopt option 3. Sparse attribute storage remains as in ADR-0004. ADR-0030 provides the source/target relation; this ADR determines the meaning of transferring a particular key along that relation. Neither geometric proximity nor equality of local indices establishes inheritance.

### 1. Declare policies per attribute key

A registered key has a value schema/version, applicable entity kinds, deterministic equality rule, transfer policy, merge policy and required/optional behavior. Where necessary it declares dependencies on geometry, orientation, parameterization or a support region. Policy configuration belongs to the operation's Context/application configuration; there is no mutable global registry.

Ordinary unregistered passive keys use ADR-0011's same-kind inheritance with deterministic `KeepFirst` conflict handling and warnings. This is a data-preservation convention, not a claim that the kernel understands an arbitrary key. Known semantic keys must not silently fall back to passive handling when their registered policy is missing. Callers must register policies for values such as total loads, unique IDs or geometry-dependent feature labels before relying on their semantics.

| Attribute family | Typical policy obligation |
| --- | --- |
| Passive colour or descriptive label | Copy to same-kind continuations; merge equal values or resolve differences explicitly. |
| Region/boundary-condition tag | Preserve supported regions and orientation assumptions; require revalidation or reject a merge that cannot preserve the application's region distinction. |
| Intensive field, e.g. pressure | Define restriction and orientation/frame handling; a shared numeric value alone may be insufficient. |
| Extensive quantity, e.g. total force | Application-defined distribution/conservation across a split; never generic duplication presented as conservation. |
| Feature/shape modifier | Revalidate dependencies after changes, invalidate/drop with a diagnostic, or fail if required. |
| Application persistent/unique ID | Application-defined regeneration or split/merge handling. The reserved ID slot remains uninterpreted by the kernel; passive copying does not promise uniqueness. |

The kernel implements transfer mechanics and diagnostics, not a universal solver for these semantics.

### 2. Apply transfer rules to events, not raw provenance counts

| Event | Default passive behavior | Semantic-key obligation |
| --- | --- | --- |
| Unchanged continuation | Preserve the value. | Preserve only while registered dependencies remain valid; an orientation-dependent key must check its relevant use/orientation. |
| Modified continuation | Copy to a same-kind continuation. | Revalidate, transform, restrict or explicitly drop/fail according to the key. |
| Split | Copy a passive value to each same-kind retained child. | Apply declared restriction/distribution/ID policy to the complete split event. |
| Merge | Collect candidate values from all same-kind contributors; equal values can collapse to one. | Validate compatibility even for equal values when supports, frames or dependencies differ; apply the configured conflict rule. |
| Generated entity | Inherit eligible same-kind generator values; no automatic cross-kind copying. | Select declared generators/roles or compute a new value. A face made from an edge does not acquire every edge attribute automatically. |
| Deleted entity | No value in the output for that entity. | Preserve provenance/diagnostic access where requested, without attaching data to a nonexistent target. |

The event manifest distinguishes an absent value from an explicit null value. By default, absence contributes no passive candidate and one present candidate may be retained. A registered required key can reject missing contributors or missing output coverage. Creation roots with no input entities receive only explicit constructor/import values or registered defaults; they do not invent inherited attributes.

A key representing one uniform value on an entity cannot preserve two conflicting subregions after a face merge merely by choosing one value. The application must authorize a resolution, use a supported region-valued representation, or fail; future topology-preserving alternatives require their own operation contract. Attribute propagation does not silently rewrite topology to satisfy a tag.

### 3. Conflicts are explicit and ordered

Preserve ADR-0011's `MergePolicy` choices:

- **KeepFirst:** select the lowest deterministic source key `(input ordinal, entity kind, input index, generator role/event key)` among eligible candidates. Document operand-order dependence; thread scheduling never determines the winner.
- **Drop:** remove the key from the target and report why.
- **Custom:** a registered resolver receives the key, target, ordered source values including missingness, qualified lineage and relevant event/dependency information. It returns a validated value, an explicit drop, or a structured failure.

Every conflict reports the key, source/target references, conflicting candidates, chosen policy and disposition in `OperationResult`. Every nontrivial invalidation/drop/revalidation failure is also a diagnostic, even if no two values conflict. Standard passive copies need no warning. Diagnostic data must remain interpretable while the result owns its manifests; logging alone is insufficient.

A strict application can promote conflict/drop/missing-required-key warnings to errors under ADR-0006 before result publication. Required output attributes must satisfy their policy postconditions; a missing required key is a failure even if geometry passes validation. Error reporting distinguishes attribute/policy failures through stable diagnostic reasons; any needed `ErrorCode` extension is added with the implementation, not by reusing `InvalidInput` for valid geometric inputs.

### 4. Callbacks preserve determinism and failure guarantees

Resolvers operate on immutable inputs and a private result candidate. They cannot mutate input shapes, publish the candidate, change Context settings or call back into the same active builder. Results are validated against the key schema and applicability constraints.

Invoke resolvers in a deterministic serialized order in the first implementation. Callbacks must be deterministic functions of supplied data, without externally visible side effects; the kernel cannot roll back arbitrary application effects. They report expected failures through the structured result interface and obey ADR-0006's exception contract (`std::bad_alloc` may propagate). A thrown exception outside that contract is programmer misuse, not a successful fallback to KeepFirst. Python callbacks in GIL-released workers remain excluded under ADR-0020; this ADR does not introduce that facility.

Compute lineage, transfer values, collect diagnostics, apply strict warning/coverage policies and perform required validation before atomically publishing shape/history/metadata. Failure discards the private candidate under ADR-0026/0027. Required metadata policy failure cannot leave a successful geometry-only result in the ordinary operation API.

### 5. Persistence and mesh export preserve the declared scope

Persist key identifiers/schema versions, values and serializable support/ownership references when saving attributes. Persist policy identifiers and effective configuration in replay provenance where needed; executable callbacks are not serialized. Replaying a semantic transfer without its required policy/version fails explicitly rather than defaulting to passive copying.

ADR-0012 still permits unknown optional attribute chunks to be skipped with warnings. An application requiring those keys promotes the resulting loss/missing-coverage diagnostic to failure; a successful geometry load alone is not proof that its engineering metadata survived.

For the ADR-0014 duct workflow, B-rep face tags map to mesh elements through their face associations. Every output triangle has exactly one required inlet/outlet/wall tag and all lost/conflicting coverage is detected before the analysis workflow accepts it. This requirement does not assert that all possible attributes are uniform per face or that arbitrary field interpolation is implemented.

## Consequences

Passive labels remain easy to use. Engineering consumers gain explicit control over semantics and can refuse lossy results. Policy registration, complete lineage, deterministic resolver execution and metadata validation add work; expensive application-specific transfer remains outside the core geometry algorithms.

## Verification

1. Split a coloured/tagged face: passive copies reach all retained children. A declared total-load policy conserves the total across children; a missing required semantic policy fails instead of duplicating the load.
2. Merge equal and unequal values; test absent versus null, same local indices in different inputs, reversed operand order and several generator roles. KeepFirst follows the documented total order and reports every conflict.
3. Changing a feature's supporting geometry triggers its registered validation; an invalid label is dropped with diagnostics or causes failure if required. Test orientation/frame-dependent data as well as scalar labels.
4. Cross-kind generation does not copy unrelated keys. A deleted source can supply legitimate Generated provenance without leaving dangling attribute targets. Persistent-ID tests distinguish passive copying from an application uniqueness policy.
5. Inject resolver failure, allocation failure and strict warning promotion: no partial shape/history/metadata is published and immutable inputs remain unchanged. G1 checks include values, conflicts and dispositions.
6. Native round trips preserve supported key schemas and values; unknown optional metadata produces visible loss, and strict required-key applications reject incomplete loads. Replay without a necessary resolver version fails explicitly.
7. The duct cut → validation → mesh workflow preserves required tags, labels the new hole wall as wall, records resolved conflicts and rejects missing/ambiguous triangle coverage before external analysis.

## References and scope

- Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006, §8.3 and §§8.3.2–8.3.4, printed pp. 262–269 (reviewed in the supplied book).
- [ADR-0006](0006-error-handling.md), [ADR-0011](0011-history-and-naming.md), [ADR-0012](0012-native-file-format.md), [ADR-0014](0014-mesh-generation-boundary.md), [ADR-0020](0020-cpp-c-abi-and-python-bindings.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0030](0030-entity-lineage-in-operation-results.md).

This proposal supplements the existing attribute policy and supersedes none. Its acceptance criteria are staged with implementation; no runtime attribute or callback support is claimed by this document.
