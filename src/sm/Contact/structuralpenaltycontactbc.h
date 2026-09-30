/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef structuralpenaltycontactbc_h
#define structuralpenaltycontactbc_h


#include "Contact/contactbc.h"
#include "Contact/contactpoint.h"
#include "ContactSurface/structuralfecontactsurface.h"
#include <map>
#include <string>
#include <iosfwd>

///@name Input fields for _IFT_ContactElement
//@{
#define _IFT_StructuralPenaltyContactBoundaryCondition_Name "structuralpenaltycontactbc"
#define _IFT_StructuralPenaltyContactBoundaryCondition_penaltyNormal "pn"
#define _IFT_StructuralPenaltyContactBoundaryCondition_penaltyTangential "pt"
#define _IFT_StructuralPenaltyContactBoundaryCondition_friction "friction"
#define _IFT_StructuralPenaltyContactBoundaryCondition_frictionTransition "frictiontransition"
#define _IFT_StructuralPenaltyContactBoundaryCondition_frictionHardening "frictionhardening"
////

#define _IFT_StructuralPenaltyContactBoundaryCondition_masterSurfaceNum "mastersurface"
#define _IFT_StructuralPenaltyContactBoundaryCondition_slaveSurfaceNum "slavesurface"
///
#define _IFT_StructuralPenaltyContactBoundaryCondition_nsd "nsd"
///
#define _IFT_StructuralPenaltyContactBoundaryCondition_algo "algo"
#define _IFT_StructuralPenaltyContactBoundaryCondition_searchPadding "searchpadding"
#define _IFT_StructuralPenaltyContactBoundaryCondition_facetOwnershipHysteresis "facethysteresis"
#define _IFT_StructuralPenaltyContactBoundaryCondition_searchTol "searchtol"
#define _IFT_StructuralPenaltyContactBoundaryCondition_generalizedFeatures "generalizedfeatures"
#define _IFT_StructuralPenaltyContactBoundaryCondition_directionalProjection "directionalprojection"
#define _IFT_StructuralPenaltyContactBoundaryCondition_autoPenalty "autopenalty"
#define _IFT_StructuralPenaltyContactBoundaryCondition_tangentMode "tangentmode"
#define _IFT_StructuralPenaltyContactBoundaryCondition_fdCheck "fdcheck"
#define _IFT_StructuralPenaltyContactBoundaryCondition_fdPerturbation "fdperturbation"
#define _IFT_StructuralPenaltyContactBoundaryCondition_fdOutputPrefix "fdoutputprefix"
#define _IFT_StructuralPenaltyContactBoundaryCondition_fdTolerance "fdtolerance"

//@}

namespace oofem {
class Domain;
class SparseMtrx;
class TimeStep;
class DofManager;
class UnknownNumberingScheme;
class FloatMatrix;

/**
 * Geometrically exact surface-to-surface penalty contact with Coulomb friction.
 *
 * The formulation follows A. Konyukhov and K. Schweizerhof, Computational
 * Contact Mechanics: Geometrically Exact Theory for Arbitrary Shaped Bodies,
 * Springer, LNACM 67 (2013): closest-point projection in Chapter 3
 * (pp. 35-62), penalty/friction evolution and return mapping in Section 6.1
 * (pp. 136-143), and consistent tangents in Section 7.1 (pp. 185-194).
 * Equation references in the implementation refer to this book.
 *
 * @warning The frictional branch (\c friction != 0, and the optional
 * \c frictionTransition / \c frictionHardening regularizations) is an
 * experimental, development-phase feature that has not yet been fully
 * verified. Use \c friction 0 (the default assumption for production
 * analyses) unless you are actively validating the friction path. The
 * normal-contact (frictionless) path is unaffected and is exercised by the
 * regular regression/benchmark suite; the frictional cases are intentionally
 * kept out of that automatically-run suite (see
 * tests/regression/benchmark/sm/contact/friction_wip/).
 *
 * Contact pairs and history are owned by the selected search algorithm.
 */


enum class ContactProcess{ None, Sticking, Sliding};

struct ContactTractionState
{
    double normalTraction = 0.0;
    FloatArray tangentialTraction;
    FloatArray trialTangentialTraction;
    ContactProcess mode = ContactProcess::None;
    double tangentialNorm = 0.0;
    double frictionBound = 0.0;
    double accumulatedPlasticSlip = 0.0;
};
  

class OOFEM_EXPORT StructuralPenaltyContactBoundaryCondition : public ContactBoundaryCondition
{
private:
    double penalty_normal;  ///< penalty in the normal direction.
    double penalty_tangential; ///< penalty in the tangent directions.
    double friction; ///< coefficient of friction
    /// Smoothness of the differentiable stick/slip projection. Zero recovers
    /// the sharp Coulomb return map; smaller positive values are sharper.
    double frictionTransition = 0.0;
    /// Dimensionless post-yield tangent ratio in [0,1). Zero is perfect
    /// Coulomb friction; positive values enable isotropic slip hardening.
    double frictionHardening = 0.0;
    /// Absolute broad-phase padding; negative selects the geometry-based default.
    double searchPadding = -1.0;
    /// Relative distance-squared margin required to unseat a pair's currently
    /// owned master facet in favor of a new candidate (see
    /// isBetterContactProjection in contactsearch.C). Zero (default) disables
    /// hysteresis. Nonzero damps Newton chattering caused by a slave point
    /// sitting on a near-degenerate shared edge between adjacent facets,
    /// where the true closest facet flips every iteration.
    double facetOwnershipHysteresis = 0.0;
    /// Parametric-domain margin passed to the contact surfaces' "still inside
    /// this facet" test (FEContactSurface::domainTolerance). Comparable to
    /// other codes' sliding-elastic-interface search tolerance (typical
    /// default 0.01). Default here preserves the previous unconditional
    /// near-zero tolerance.
    double searchTol = 1.e-10;
    /// Extends the closest-point search from plain facet projection to also
    /// consider edge and vertex features of the master surface. Mutually
    /// exclusive with directionalProjection.
    bool generalizedFeatures = false;
    /// Projects each slave point along the slave surface's own normal
    /// direction instead of a general closest-point search. Requires nsd 3.
    /// Mutually exclusive with generalizedFeatures.
    bool directionalProjection = false;
    /// When true, pn/pt are ignored and the penalty stiffnesses are instead
    /// computed per slave contact element from the contacting materials'
    /// initial Young's modulus and element geometry.
    bool autoPenalty = false;
    /// 0: consistent automatic selection, 1: rate-form analytical diagnostic,
    /// 2: branch-frozen finite-difference tangent, 3: exact finite-step
    /// analytical tangent on a smooth surface branch, including facet-history
    /// columns (FD fallback for unsupported projection features).
    int tangentMode = 0;
    /// E_n A/V factor indexed by the slave contact element.
    std::map<int, double> automaticPenaltyFactors;
    int surface_dimension; ///< dimension of the surface, i.e., nsd - 1
    /// 0: plain surface-to-surface search, 1: sweep-and-prune broad-phase
    /// search (3D only).
    int algo;
    /// Development/debugging aid: verify the analytical contact tangent
    /// against a finite-difference tangent on every call. Adds significant
    /// runtime cost.
    bool finiteDifferenceCheck = false;
    /// Relative perturbation size used by finiteDifferenceCheck.
    double finiteDifferencePerturbation = 1.e-7;
    /// Relative tolerance used by finiteDifferenceCheck.
    double finiteDifferenceRelativeTolerance = 0.0;
    /// Filename prefix for finiteDifferenceCheck diagnostic output.
    std::string finiteDifferenceOutputPrefix = "contact_fd";
    bool finiteDifferenceCheckInProgress = false;
    ContactPair *finiteDifferenceContactPair = nullptr;
    double finiteDifferenceContactArea = 0.0;
    
    //@todo: move to contact search algorithm
    int masterSurfaceNumber;
    int slaveSurfaceNumber;
    IntArray masterSurfaceElements;
    IntArray slaveSurfaceElements;
    // should be generalized to structural contact surface, i.e., a surface that is not necessarly discretized into finite elements
    StructuralFEContactSurface *masterContactSurface;
    StructuralFEContactSurface *slaveContactSurface;

public:
    /// Constructor.
    StructuralPenaltyContactBoundaryCondition(int n, Domain *d) : ContactBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~StructuralPenaltyContactBoundaryCondition() {};
    // initialization, i.e., reading input filex
    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    void postInitialize() override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    virtual const char *giveClassName() const override { return "StructuralPenaltyContactBoundaryCondition"; }
    virtual const char *giveInputRecordName() const override { return _IFT_StructuralPenaltyContactBoundaryCondition_Name; }
    //
    void setupContactSearchAlgorithm() override;
private:
    // compute tangent stiffness matrix
    void  computeTangentFromContact(FloatMatrix &answer, ContactPair *cp, TimeStep *tStep) override;
    // compute internal forces
    void computeInternalForcesFromContact(FloatArray &answer, ContactPair *cp, TimeStep *tStep) override;
    ContactTractionState computeTractionState(ContactPair* contactPair, TimeStep *tStep);
    // location arrays
    void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
 private:
    FloatMatrix computeContravariantMetric(const std::vector<FloatArray> &tangent_vectors);
    FloatMatrix computeCovariantMetric(const std::vector<FloatArray> &tangent_vectors);
    /**
     * Reconstructs the committed physical friction traction and, when the
     * closest-point carrier changes facet or feature, parallel-transports it
     * by the minimal rotation between the committed and current normals.
     */
    FloatArray computeTransportedPreviousPhysicalTraction(ContactPair *contactPair);
    double computeTangentialNorm(const FloatArray &tangentialTraction, const FloatMatrix &contravariantMetric) const;
    void computeFrictionReturnMagnitude(double &returnedMagnitude,
                                        double &derivativeWrtTrialMagnitude,
                                        double &derivativeWrtFrictionBound,
                                        double trialMagnitude,
                                        double frictionBound) const;
    double giveContactAreaForLinearization(ContactPair *contactPair, TimeStep *tStep) const;
    void initializeAutomaticPenaltyFactors();
    double giveAutomaticPenaltyFactor(ContactPair *contactPair) const;
    double giveNormalPenalty(ContactPair *contactPair) const;
    double giveTangentialPenalty(ContactPair *contactPair) const;
    bool computeDiscreteTangent(FloatMatrix &answer,
                                ContactPair *contactPair,
                                TimeStep *tStep);
    bool updateFiniteDifferencePairGeometry(ContactPair *contactPair, int masterElementId,
                                            ContactFeatureType featureType, int featureIndex,
                                            TimeStep *tStep);
    void dumpFiniteDifferenceCheck(const FloatMatrix &analyticalTangent, ContactPair *contactPair, TimeStep *tStep);
    void dumpFiniteDifferenceKinematicCheck(std::ostream &output,
                                            const FloatMatrix &suppliedTangent,
                                            const FloatMatrix &numericalTangent,
                                            ContactPair *contactPair,
                                            TimeStep *tStep);
    FloatMatrix computeFiniteDifferenceTangent(std::size_t pairIndex, TimeStep *tStep);
    std::size_t findContactPairIndex(const ContactPair *contactPair) const;

    
};
} // end namespace oofem
#endif // node2segmentpenaltycontact_h
