// SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
// SPDX-License-Identifier: MPL-2.0
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

module;

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

export module Stereon.Core:Sensitivity;

import :Error;
import :Handle;
import :Vector;

/// The foundation of shape sensitivities (ADR-0021, stage 1): the types an
/// operation uses to say what its result depends on and what it can
/// differentiate, and an evaluation uses to say what it computed. No
/// operation computes a derivative yet; every one is declared and reports
/// NotSupported, never a silent zero. Names are provisional until the first
/// derivative slice applies ADR-0017.
///
/// - Parameter identities (section 3, 7): a ParameterId with its quantity
///   and unit (lengths in millimetres, angles in radians, ADR-0001), and a
///   ParameterSet that fixes the order of every parameter vector, seed and
///   gradient over it and each parameter's reference scale.
/// - Statuses and results (section 5): DerivativeStatus, and the
///   SensitivityResult a completed evaluation returns as
///   Result<SensitivityResult>, which carries a derivative only for Valid
///   and NearEvent.
/// - Capability declarations (section 5, Decision): OperationCapability,
///   per operation, queryable before any evaluation; each library keeps a
///   catalogue (Capabilities partitions), and tools/check_capabilities.py
///   fails the build's checks when a public operation has none.
/// - Correspondence (section 4): a TrackedSample names its entity by the
///   snapshot it belongs to and its slot there, its orientation, its local
///   coordinates and its tracking convention. It is bound to that snapshot,
///   so Compact() or any renumbering cannot redefine it: a sample of one
///   snapshot is not a sample of another.
/// - The sensitivity context (section 3): caller-owned, one evaluation's
///   seeds, adjoints and workspace, separate from the operation Context and
///   never stored in a Shape; not copyable, so two contexts share nothing.

export namespace Stereon {

// --- Parameter identities -------------------------------------------------------

/// What a design parameter measures (ADR-0021, section 3).
enum class ParameterQuantity : std::uint8_t {
    /// A length, in millimetres (ADR-0001).
    Length,
    /// An angle, in radians (ADR-0001).
    Angle,
    /// A number without a unit: a weight, a ratio, a parameter value.
    Dimensionless,
};

/// The unit of a quantity, as ADR-0001 fixes it: "mm", "rad" or "1".
[[nodiscard]] constexpr std::string_view UnitOf(const ParameterQuantity quantity) noexcept
{
    switch (quantity) {
    case ParameterQuantity::Length:
        return "mm";
    case ParameterQuantity::Angle:
        return "rad";
    case ParameterQuantity::Dimensionless:
        return "1";
    }
    return "?";
}

/// A design parameter: its position in the ParameterSet that made it, and
/// what it measures. Default-constructed it names nothing.
class ParameterId {
public:
    static constexpr std::uint32_t Invalid = 0xFFFFFFFFu;

    constexpr ParameterId() = default;
    constexpr ParameterId(const std::uint32_t index, const ParameterQuantity quantity) noexcept
        : m_Index(index)
        , m_Quantity(quantity)
    {}

    [[nodiscard]] constexpr bool IsValid() const noexcept { return m_Index != Invalid; }
    /// The parameter's position in every vector over its set.
    [[nodiscard]] constexpr std::uint32_t GetIndex() const noexcept { return m_Index; }
    [[nodiscard]] constexpr ParameterQuantity GetQuantity() const noexcept { return m_Quantity; }

    constexpr bool operator==(const ParameterId&) const = default;

private:
    std::uint32_t m_Index = Invalid;
    ParameterQuantity m_Quantity = ParameterQuantity::Dimensionless;
};

/// The parameters of one model, in order. The order is the convention of
/// every vector over them (ADR-0021, section 7): component k of a seed, a
/// gradient or a row of a Jacobian belongs to the parameter with index k.
/// Each parameter has a reference scale in its unit: seeds and derivatives
/// are per reference, so a step of 1 in a seed moves the parameter by its
/// reference, which keeps millimetre and radian parameters comparable.
class ParameterSet {
public:
    /// Adds a parameter at the end. InvalidInput for an empty or repeated
    /// name, or a reference that is not finite and positive.
    [[nodiscard]] Result<ParameterId> Add(std::string_view name, ParameterQuantity quantity,
                                          double reference = 1.0);

    [[nodiscard]] std::size_t GetSize() const noexcept { return m_Names.size(); }
    /// True when the id is one this set made.
    [[nodiscard]] bool Contains(ParameterId id) const noexcept;
    /// The parameter's name, unit and reference scale; the id must be one
    /// this set made.
    [[nodiscard]] std::string_view GetName(ParameterId id) const;
    [[nodiscard]] std::string_view GetUnit(ParameterId id) const;
    [[nodiscard]] double GetReference(ParameterId id) const;
    /// The parameter of that name, invalid if there is none.
    [[nodiscard]] ParameterId Find(std::string_view name) const noexcept;

private:
    std::vector<std::string> m_Names;
    std::vector<ParameterQuantity> m_Quantities;
    std::vector<double> m_References;
};

// --- Statuses ---------------------------------------------------------------------

/// The validity of a completed derivative evaluation (ADR-0021, section 5).
enum class DerivativeStatus : std::uint8_t {
    /// Computed under the declared regularity assumptions and the checks
    /// that apply; not a global mathematical certificate.
    Valid,
    /// A local derivative, but near a branch event or poorly conditioned:
    /// usable, and qualified as such.
    NearEvent,
    /// The requested ordinary derivative is not defined under the selected
    /// semantics.
    Nondifferentiable,
    /// The requested product or operation has no implementation.
    NotSupported,
    /// Correspondence, convergence or regularity could not be established
    /// well enough to return a trustworthy derivative.
    Indeterminate,
};

[[nodiscard]] constexpr std::string_view ToString(const DerivativeStatus status) noexcept
{
    switch (status) {
    case DerivativeStatus::Valid:
        return "Valid";
    case DerivativeStatus::NearEvent:
        return "NearEvent";
    case DerivativeStatus::Nondifferentiable:
        return "Nondifferentiable";
    case DerivativeStatus::NotSupported:
        return "NotSupported";
    case DerivativeStatus::Indeterminate:
        return "Indeterminate";
    }
    return "DerivativeStatus";
}

/// True for the statuses whose result carries a derivative: Valid and the
/// qualified NearEvent. The others carry none (ADR-0021, section 5).
[[nodiscard]] constexpr bool HasDerivative(const DerivativeStatus status) noexcept
{
    return status == DerivativeStatus::Valid || status == DerivativeStatus::NearEvent;
}

// --- Capability declarations --------------------------------------------------------

/// The derivative products an operation can implement (ADR-0021, section 7),
/// discoverable independently.
enum class DerivativeProduct : std::uint8_t {
    /// Jacobian-vector product, J p': forward mode.
    Jvp = 1,
    /// Vector-Jacobian product, J^T x': reverse mode.
    Vjp = 2,
    /// The dense Jacobian, for small problems; never the only interface.
    DenseJacobian = 4,
};

/// A set of products as bits; empty for an operation that implements none.
class ProductSet {
public:
    constexpr ProductSet() = default;
    constexpr ProductSet(const std::initializer_list<DerivativeProduct> products) noexcept
    {
        for (const DerivativeProduct p : products) {
            m_Bits |= static_cast<std::uint8_t>(p);
        }
    }
    [[nodiscard]] constexpr bool Contains(const DerivativeProduct p) const noexcept
    {
        return (m_Bits & static_cast<std::uint8_t>(p)) != 0;
    }
    [[nodiscard]] constexpr bool IsEmpty() const noexcept { return m_Bits == 0; }
    constexpr bool operator==(const ProductSet&) const = default;

private:
    std::uint8_t m_Bits = 0;
};

/// What an operation's result depends on (ADR-0021, section 3): the inputs
/// a derivative would have to follow, as bits.
enum class Dependency : std::uint8_t {
    /// The values that define new geometry: a centre, radii and axes; poles,
    /// weights and knots.
    DefiningValues = 1,
    /// Geometry given as input: the curve converted, the profile swept, the
    /// surface trimmed or intersected.
    InputGeometry = 2,
    /// Points given as input: the points projected or fitted.
    InputPoints = 4,
    /// The operation's own numeric arguments: a range, a direction, an axis,
    /// a parameter value, a length.
    OperationValues = 8,
    /// Tolerances and budgets, which steer discrete choices (ADR-0001,
    /// ADR-0021 section 8): a change between two evaluations is a branch
    /// event.
    Tolerances = 16,
    /// The input shape's connectivity, which a derivative never changes.
    Topology = 32,
};

/// A set of dependencies as bits.
class DependencySet {
public:
    constexpr DependencySet() = default;
    constexpr DependencySet(const std::initializer_list<Dependency> dependencies) noexcept
    {
        for (const Dependency d : dependencies) {
            m_Bits |= static_cast<std::uint8_t>(d);
        }
    }
    [[nodiscard]] constexpr bool Contains(const Dependency d) const noexcept
    {
        return (m_Bits & static_cast<std::uint8_t>(d)) != 0;
    }
    [[nodiscard]] constexpr bool IsEmpty() const noexcept { return m_Bits == 0; }
    constexpr bool operator==(const DependencySet&) const = default;

private:
    std::uint8_t m_Bits = 0;
};

/// How a sample follows the geometry when a parameter moves (ADR-0021,
/// section 4).
enum class TrackingConvention : std::uint8_t {
    /// A point at fixed (u, v) of a surface or face.
    FixedSurfaceCoordinates = 1,
    /// A point at a fixed parameter of an edge or a curve.
    FixedEdgeParameter = 2,
    /// A point on a trim that moves with it: its surface coordinates follow
    /// the moving pcurve or intersection.
    MovingTrimCoordinates = 4,
    /// The point a uniquely defined projection gives.
    UniqueProjection = 8,
};

[[nodiscard]] constexpr std::string_view ToString(const TrackingConvention convention) noexcept
{
    switch (convention) {
    case TrackingConvention::FixedSurfaceCoordinates:
        return "FixedSurfaceCoordinates";
    case TrackingConvention::FixedEdgeParameter:
        return "FixedEdgeParameter";
    case TrackingConvention::MovingTrimCoordinates:
        return "MovingTrimCoordinates";
    case TrackingConvention::UniqueProjection:
        return "UniqueProjection";
    }
    return "TrackingConvention";
}

/// A set of tracking conventions as bits; empty for an operation whose
/// result has no samples to track (a functional's value).
class ConventionSet {
public:
    constexpr ConventionSet() = default;
    constexpr ConventionSet(const std::initializer_list<TrackingConvention> conventions) noexcept
    {
        for (const TrackingConvention c : conventions) {
            m_Bits |= static_cast<std::uint8_t>(c);
        }
    }
    [[nodiscard]] constexpr bool Contains(const TrackingConvention c) const noexcept
    {
        return (m_Bits & static_cast<std::uint8_t>(c)) != 0;
    }
    [[nodiscard]] constexpr bool IsEmpty() const noexcept { return m_Bits == 0; }
    constexpr bool operator==(const ConventionSet&) const = default;

private:
    std::uint8_t m_Bits = 0;
};

/// One operation's derivative capability, declared before any evaluation
/// (ADR-0021, section 5 and Decision): its name (the public function, all
/// its overloads), what its result depends on, the conventions its samples
/// may be tracked by, the assumptions under which a derivative would be
/// valid, and the products it implements. An empty product set means
/// NotSupported for every product, which is the declaration every
/// operation has today.
struct OperationCapability {
    std::string_view Operation;
    DependencySet Dependencies;
    ConventionSet Correspondence;
    std::string_view Assumptions;
    ProductSet Products;

    [[nodiscard]] constexpr bool Supports(const DerivativeProduct product) const noexcept
    {
        return Products.Contains(product);
    }
};

/// A declaration, for the capability catalogues. tools/check_capabilities.py
/// reads the `Declare("Name", ...)` calls of each library's Capabilities
/// partition and checks them against the library's public functions.
[[nodiscard]] constexpr OperationCapability Declare(const std::string_view operation,
                                                    const DependencySet dependencies,
                                                    const ConventionSet correspondence,
                                                    const std::string_view assumptions,
                                                    const ProductSet products = {}) noexcept
{
    OperationCapability capability;
    capability.Operation = operation;
    capability.Dependencies = dependencies;
    capability.Correspondence = correspondence;
    capability.Assumptions = assumptions;
    capability.Products = products;
    return capability;
}

/// The declaration of `operation` in a catalogue, or null when it has none.
[[nodiscard]] constexpr const OperationCapability*
FindCapability(const std::span<const OperationCapability> catalogue,
               const std::string_view operation) noexcept
{
    for (const OperationCapability& capability : catalogue) {
        if (capability.Operation == operation) {
            return &capability;
        }
    }
    return nullptr;
}

// --- Correspondence -----------------------------------------------------------------

/// Whether a sample uses its entity as stored or reversed (ADR-0037).
enum class SampleOrientation : std::uint8_t {
    AsStored,
    Reversed,
};

/// A point an analysis follows as parameters move (ADR-0021, section 4): on
/// a face at (u, v), or on an edge at t, of one frozen snapshot, under a
/// tracking convention. The snapshot's store and the entity's slot identify
/// it together, never the slot alone: a sample of one snapshot is no sample
/// of the compacted or edited one, whose slots mean something else; carrying
/// it across is the History's business (ADR-0011), not a renumbering's.
class TrackedSample {
public:
    /// A sample on a face. InvalidInput for a zero store or slot index of
    /// no entity, coordinates that are not finite, or a convention that does
    /// not track a point of a face (FixedEdgeParameter).
    [[nodiscard]] static Result<TrackedSample>
    OnFace(StoreId snapshot, std::uint32_t face, const Point2& uv, TrackingConvention convention,
           SampleOrientation orientation = SampleOrientation::AsStored);
    /// A sample on an edge. InvalidInput as above, or for a convention that
    /// does not track a point of an edge (FixedSurfaceCoordinates).
    [[nodiscard]] static Result<TrackedSample>
    OnEdge(StoreId snapshot, std::uint32_t edge, double t, TrackingConvention convention,
           SampleOrientation orientation = SampleOrientation::AsStored);

    [[nodiscard]] EntityKind GetKind() const noexcept { return m_Kind; }
    [[nodiscard]] StoreId GetSnapshot() const noexcept { return m_Snapshot; }
    [[nodiscard]] std::uint32_t GetIndex() const noexcept { return m_Index; }
    [[nodiscard]] SampleOrientation GetOrientation() const noexcept { return m_Orientation; }
    [[nodiscard]] TrackingConvention GetConvention() const noexcept { return m_Convention; }
    /// (u, v) on a face; (t, 0) on an edge.
    [[nodiscard]] Point2 GetCoordinates() const noexcept { return m_Coordinates; }
    /// True when the sample belongs to the snapshot with this store.
    [[nodiscard]] bool IsOf(const StoreId snapshot) const noexcept
    {
        return m_Snapshot == snapshot;
    }

private:
    TrackedSample() = default;

    EntityKind m_Kind = EntityKind::Face;
    StoreId m_Snapshot = 0;
    std::uint32_t m_Index = 0;
    SampleOrientation m_Orientation = SampleOrientation::AsStored;
    TrackingConvention m_Convention = TrackingConvention::FixedSurfaceCoordinates;
    Point2 m_Coordinates;
};

// --- The sensitivity context ----------------------------------------------------------

/// What one derivative evaluation owns (ADR-0021, section 3): the parameters
/// it is over, the seed of a forward product, the adjoint of a reverse one,
/// and scratch space for solvers. Owned by the caller, separate from the
/// operation Context (which keeps tolerances and schedulers) and never
/// stored in a Shape, which stays immutable (ADR-0005). It cannot be copied,
/// so two contexts share no mutable state and evaluations in separate
/// contexts can run concurrently; it can be moved.
class SensitivityContext {
public:
    explicit SensitivityContext(ParameterSet parameters);
    SensitivityContext(const SensitivityContext&) = delete;
    SensitivityContext& operator=(const SensitivityContext&) = delete;
    SensitivityContext(SensitivityContext&&) noexcept = default;
    SensitivityContext& operator=(SensitivityContext&&) noexcept = default;
    ~SensitivityContext() = default;

    [[nodiscard]] const ParameterSet& GetParameters() const noexcept { return m_Parameters; }

    /// The seed p' of a forward product, one component per parameter in the
    /// set's order, each per the parameter's reference. InvalidInput for a
    /// length other than the set's size or a component that is not finite
    /// (ADR-0021, section 5: invalid seeds are KernelErrors).
    [[nodiscard]] Result<void> SetSeed(std::span<const double> seed);
    [[nodiscard]] std::span<const double> GetSeed() const noexcept { return m_Seed; }

    /// The adjoint x' of a reverse product: a covector over `size` sample
    /// coordinates under the operation's declared pairing. InvalidInput for
    /// a component that is not finite.
    [[nodiscard]] Result<void> SetAdjoint(std::span<const double> adjoint);
    [[nodiscard]] std::span<const double> GetAdjoint() const noexcept { return m_Adjoint; }

    /// Scratch space of at least `size` doubles for one evaluation's
    /// solvers, reused between calls; its contents are the evaluation's.
    [[nodiscard]] std::span<double> Workspace(std::size_t size);

    /// Drops the seed, the adjoint and the workspace, keeping the parameters.
    void Clear() noexcept;

private:
    ParameterSet m_Parameters;
    std::vector<double> m_Seed;
    std::vector<double> m_Adjoint;
    std::vector<double> m_Workspace;
};

// --- Results --------------------------------------------------------------------------

/// What a completed derivative evaluation returns, inside a
/// Result<SensitivityResult> (ADR-0021, section 5; name provisional): its
/// status; the derivative, only for Valid and NearEvent; and what is known
/// of the regime: the branch and topology signature, the entities it
/// concerns, the tolerance policy it was computed under, the solver residual
/// and the conditioning, each NaN where unknown.
class SensitivityResult {
public:
    /// A computed derivative. `status` must be Valid or NearEvent; any other
    /// status gives InvalidInput, since it would carry a payload it must not.
    [[nodiscard]] static Result<SensitivityResult>
    Computed(DerivativeStatus status, std::vector<double> values, std::uint64_t signature);
    /// No derivative: NotSupported, Nondifferentiable or Indeterminate, with
    /// the reason. InvalidInput for Valid or NearEvent, which need values.
    [[nodiscard]] static Result<SensitivityResult> Without(DerivativeStatus status,
                                                           std::string reason);
    /// The answer for a product an operation does not implement:
    /// NotSupported, naming the operation and the product.
    [[nodiscard]] static SensitivityResult NotSupported(const OperationCapability& capability,
                                                        DerivativeProduct product);

    [[nodiscard]] DerivativeStatus GetStatus() const noexcept { return m_Status; }
    /// The derivative for Valid and NearEvent; nothing otherwise, so that no
    /// caller can mistake an absent derivative for a zero one.
    [[nodiscard]] std::optional<std::span<const double>> GetValues() const noexcept;
    /// True for NearEvent: the derivative is usable but qualified.
    [[nodiscard]] bool IsQualified() const noexcept
    {
        return m_Status == DerivativeStatus::NearEvent;
    }
    [[nodiscard]] std::string_view GetReason() const noexcept { return m_Reason; }

    /// A digest of the branch and topology the derivative belongs to: two
    /// evaluations with different signatures were not in one regime.
    [[nodiscard]] std::uint64_t GetSignature() const noexcept { return m_Signature; }
    /// The entities of the input the derivative concerns, by snapshot and
    /// slot, as tracked samples name theirs.
    [[nodiscard]] std::span<const std::pair<StoreId, std::uint32_t>> GetEntities() const noexcept
    {
        return m_Entities;
    }
    void AddEntity(StoreId snapshot, std::uint32_t index);
    /// The resolution and tolerance ceiling held fixed during the evaluation
    /// (ADR-0021, section 8); NaN until set.
    [[nodiscard]] double GetResolution() const noexcept { return m_Resolution; }
    [[nodiscard]] double GetToleranceCeiling() const noexcept { return m_ToleranceCeiling; }
    void SetTolerancePolicy(double resolution, double toleranceCeiling) noexcept;
    /// The residual of the solved equations and a condition estimate; NaN
    /// where the evaluation has none.
    [[nodiscard]] double GetResidual() const noexcept { return m_Residual; }
    [[nodiscard]] double GetConditioning() const noexcept { return m_Conditioning; }
    void SetDiagnostics(double residual, double conditioning) noexcept;

private:
    SensitivityResult() = default;

    DerivativeStatus m_Status = DerivativeStatus::NotSupported;
    std::vector<double> m_Values;
    std::string m_Reason;
    std::uint64_t m_Signature = 0;
    std::vector<std::pair<StoreId, std::uint32_t>> m_Entities;
    double m_Resolution = std::numeric_limits<double>::quiet_NaN();
    double m_ToleranceCeiling = std::numeric_limits<double>::quiet_NaN();
    double m_Residual = std::numeric_limits<double>::quiet_NaN();
    double m_Conditioning = std::numeric_limits<double>::quiet_NaN();
};

} // namespace Stereon
