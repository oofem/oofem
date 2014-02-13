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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "mixedgradientpressureweakperiodic.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "engngm.h"
#include "node.h"
#include "element.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "masterdof.h"
#include "classfactory.h" // For sparse matrix creation.
#include "sparsemtrxtype.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "set.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"

namespace oofem {
REGISTER_BoundaryCondition(MixedGradientPressureWeakPeriodic);

MixedGradientPressureWeakPeriodic :: MixedGradientPressureWeakPeriodic(int n, Domain *d) : MixedGradientPressureBC(n, d)
{
    this->voldman = new Node(0, d); // Node number lacks meaning here.
    this->voldman->appendDof( new MasterDof( 1, voldman, ( DofIDItem ) d->giveNextFreeDofID() ) );
    this->tractionsdman = new Node(0, this->domain); // Node number lacks meaning here.
}


MixedGradientPressureWeakPeriodic :: ~MixedGradientPressureWeakPeriodic()
{
    delete voldman;
    delete tractionsdman;
}


int MixedGradientPressureWeakPeriodic :: giveNumberOfInternalDofManagers()
{
    return 2;
}


DofManager *MixedGradientPressureWeakPeriodic :: giveInternalDofManager(int i)
{
    if ( i == 1 ) {
        return this->tractionsdman;
    } else {
        return this->voldman;
    }
}


IRResultType MixedGradientPressureWeakPeriodic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    MixedGradientPressureBC :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->order, _IFT_MixedGradientPressureWeakPeriodic_order);
    if ( this->order < 0 ) {
        OOFEM_ERROR("MixedGradientPressureWeakPeriodic :: initializeFrom - order must be at least 0");
    }

    int nsd = this->domain->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * ( int ) pow(double ( order + 1 ), nsd - 1);
    this->tractionsdman->setNumberOfDofs(total);
    for ( int i = 1; i <= total; ++i ) {
        // Simply use t_i = S_i . n, where S_1 = [1,0,0;0,0,0;0,0,0], S_2 = [0,1,0;0,0,0;0,0,0], etc.
        // then the linear terms, [x,0,0], [0,x,0] ... and so on
        this->tractionsdman->setDof( i, new MasterDof( i, tractionsdman, ( DofIDItem ) this->domain->giveNextFreeDofID() ) );
    }

    return IRRT_OK;
}


void MixedGradientPressureWeakPeriodic :: setPrescribedDeviatoricGradientFromVoigt(const FloatArray &t)
{
    this->constructFullMatrixForm(this->devGradient, t);
    this->volGradient = 0.;
    for ( int i = 1; i <= this->devGradient.giveNumberOfRows(); i++ ) {
        this->volGradient += this->devGradient.at(i, i);
    }
}


void MixedGradientPressureWeakPeriodic :: constructFullMatrixForm(FloatMatrix &d, const FloatArray &d_voigt) const
{
    int n = d_voigt.giveSize();
    if ( n == 6 ) { // Then 3D
        d.resize(3, 3);
        d.at(1, 1) = d_voigt.at(1);
        d.at(2, 2) = d_voigt.at(2);
        d.at(3, 3) = d_voigt.at(3);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        d.at(1, 2) = d.at(2, 1) = d_voigt.at(6) * 0.5;
        d.at(1, 3) = d.at(3, 1) = d_voigt.at(5) * 0.5;
        d.at(2, 3) = d.at(3, 2) = d_voigt.at(4) * 0.5;
    } else if ( n == 3 ) { // Then 2D
        d.resize(2, 2);
        d.at(1, 1) = d_voigt.at(1);
        d.at(2, 2) = d_voigt.at(2);
        d.at(1, 2) = d.at(2, 1) = d_voigt.at(3);
    } else if ( n == 1 ) { // Then 1D
        d.resize(1, 1);
        d.at(1, 1) = d_voigt.at(1);
    } else {
        OOFEM_ERROR("MixedGradientPressureWeakPeriodic :: setPrescribedDeviatoricGradientFromVoigt: Tensor is in strange voigt format. Should be 3 or 6. Use setPrescribedTensor directly if needed.");
    }
}


void MixedGradientPressureWeakPeriodic :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                                                             const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }

    if ( eid != EID_MomentumBalance || type != TangentStiffnessMatrix ) {
        return;
    }

    IntArray loc_r, loc_c, t_loc_r, t_loc_c, e_loc_r, e_loc_c;

    // Fetch the columns/rows for the tractions;
    this->tractionsdman->giveCompleteLocationArray(t_loc_r, r_s);
    this->tractionsdman->giveCompleteLocationArray(t_loc_c, c_s);

    // Fetch the columns/rows for dvol;
    this->voldman->giveCompleteLocationArray(e_loc_r, r_s);
    this->voldman->giveCompleteLocationArray(e_loc_c, c_s);

    Set *set = this->giveDomain()->giveSet(this->set);
    IntArray bNodes;
    const IntArray &boundaries = set->giveBoundaryList();

    rows.resize(boundaries.giveSize() + 2);
    cols.resize(boundaries.giveSize() + 2);
    int i = 0;
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
        e->giveBoundaryLocationArray(loc_r, bNodes, eid, r_s);
        e->giveBoundaryLocationArray(loc_c, bNodes, eid, c_s);
        // For most uses, *loc_r == *loc_c
        rows [ i ] = loc_r;
        cols [ i ] = t_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = t_loc_r;
        cols [ i ] = loc_c;
        i++;
    }

    // The volumetric part:
    rows [ i ] = t_loc_r;
    cols [ i ] = e_loc_c;
    i++;

    rows [ i ] = e_loc_r;
    cols [ i ] = t_loc_c;
    i++;
}


void MixedGradientPressureWeakPeriodic :: integrateTractionVelocityTangent(FloatMatrix &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dv dA
    FloatArray normal, n, m, t, surfCoords, coords, contrib;
    FloatMatrix nMatrix, mMatrix;

    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    ///@todo Check the order here:
    int maxorder = this->order + interp->giveInterpolationOrder() * 3 - 2;
    IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);
    int nsd = el->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * ( int ) pow(double ( order + 1 ), nsd - 1);

    surfCoords.resize(nsd - 1);
    mMatrix.resize(nsd, total);

    answer.clear();
    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = * gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        interp->boundaryEvalN(n, boundary, lcoords, cellgeo);
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);

        // Construct the basis functions for the tractions:
        mMatrix.zero();
        int pos = 0;
        for ( int kj = 0; kj < nsd; ++kj ) {
            // First we compute the surface coordinates. Its the global coordinates without the kj component.
            for ( int c = 0, q = 0; c < nsd; ++c ) {
                if ( c != kj ) {
                    surfCoords(q) = coords(c);
                    q++;
                }
            }
            // Evaluate all the basis functions for a given column.
            this->evaluateTractionBasisFunctions(t, surfCoords);
            for ( int ti = 0; ti < t.giveSize(); ++ti ) {
                for ( int ki = 0; ki < nsd; ++ki ) {
                    mMatrix(ki, pos) = t(ti) * normal(kj);
                    pos++;
                }
            }
        }
        nMatrix.beNMatrixOf(n, nsd);

        answer.plusProductUnsym( mMatrix, nMatrix, detJ * gp->giveWeight() );
    }
    delete ir;
}


void MixedGradientPressureWeakPeriodic :: integrateTractionXTangent(FloatMatrix &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dx_m dA
    FloatArray normal, t, surfCoords, coords, contrib;

    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder() * 3 - 2;
    IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);
    int nsd = el->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * ( int ) pow(double ( order + 1 ), nsd - 1);

    surfCoords.resize(nsd - 1);

    FloatArray tmpAnswer;
    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = * gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
        contrib.resize(total);
        int pos = 0;
        for ( int kj = 0; kj < nsd; ++kj ) {
            // First we compute the surface coordinates. Its the global coordinates without the kj component.
            for ( int c = 0, q = 0; c < nsd; ++c ) {
                if ( c != kj ) {
                    surfCoords(q) = coords(c);
                    q++;
                }
            }
            // Evaluate all the basis functions for a given column.
            this->evaluateTractionBasisFunctions(t, surfCoords);
            for ( int ti = 0; ti < t.giveSize(); ++ti ) {
                for ( int ki = 0; ki < nsd; ++ki ) {
                    contrib(pos) = t(ti) * normal(kj) * coords(ki) / 3.0;
                    pos++;
                }
            }
        }

        tmpAnswer.add(detJ * gp->giveWeight(), contrib);
    }
    delete ir;
    answer.resize(tmpAnswer.giveSize(), 1);
    answer.setColumn(tmpAnswer, 1);
}


void MixedGradientPressureWeakPeriodic :: integrateTractionDev(FloatArray &answer, Element *el, int boundary, const FloatMatrix &ddev)
{
    // Computes the integral: int dt . dx dA
    FloatArray normal, t, surfCoords, coords, vM_dev, contrib;

    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder() * 3 - 2;
    IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);
    int nsd = el->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * ( int ) pow(double ( order + 1 ), nsd - 1);

    surfCoords.resize(nsd - 1);
    contrib.resize(total);
    answer.clear();

    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = * gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        // Compute v_m = d_dev . x
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
        vM_dev.beProductOf(ddev, coords);

        // "pos" loops over all the traction lagrange multipliers.
        int pos = 0;
        for ( int kj = 0; kj < nsd; ++kj ) {
            // First we compute the surface coordinates. Its the global coordinates without the kj component.
            for ( int c = 0, q = 0; c < nsd; ++c ) {
                if ( c != kj ) {
                    surfCoords(q) = coords(c);
                    q++;
                }
            }
            // Evaluate all the basis functions for a given column.
            this->evaluateTractionBasisFunctions(t, surfCoords);
            for ( int ti = 0; ti < t.giveSize(); ++ti ) {
                for ( int ki = 0; ki < nsd; ++ki ) {
                    contrib(pos) = t(ti) * normal(kj) * vM_dev(ki);
                    pos++;
                }
            }
        }

        answer.add(detJ * gp->giveWeight(), contrib);
    }
    delete ir;
}


void
MixedGradientPressureWeakPeriodic :: evaluateTractionBasisFunctions(FloatArray &answer, const FloatArray &coords)
{
    // This evaluates the function x^a * y^b, for a, b in [0,order]
    int nsd = coords.giveSize() + 1;
    int total = ( int ) pow(double ( order + 1 ), nsd - 1);
    answer.resize(total);
    int pos = 0;
    double tx = 1.;
    for ( int xi = 0; xi < ( this->order + 1 ); ++xi ) {
        double ty = 1.;
        for ( int yi = 0; yi < ( this->order + 1 ); ++yi ) {
            answer(pos) = tx * ty;
            pos++;
            ty *= coords(1);
        }
        tx *= coords(0);
    }
}


void MixedGradientPressureWeakPeriodic :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                                         CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    // Boundary condition only acts on the momentumbalance part.
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }

    if ( eid != EID_MomentumBalance ) {
        return;
    }

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    IntArray v_loc, t_loc, e_loc;  // For the velocities and stress respectively
    IntArray velocityDofIDs, tractionDofIDs, dvolDofID, bNodes;
    this->tractionsdman->giveCompleteLocationArray(t_loc, s);
    this->voldman->giveCompleteLocationArray(e_loc, s);
    this->tractionsdman->giveCompleteMasterDofIDArray(tractionDofIDs);
    this->voldman->giveCompleteMasterDofIDArray(dvolDofID);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. on the traction and on dvol.
        double rve_size = this->domainSize();

        if ( e_loc.at(1) ) {
            answer.at( e_loc.at(1) ) -= rve_size * pressure; // Note the negative sign (pressure as opposed to mean stress)
            if ( eNorms ) {
                eNorms->at( dvolDofID.at(1) ) = rve_size * pressure * rve_size * pressure;
            }
        }

        // The second contribution is on the momentumbalance equation; int t . [[ d_dev . x ]] dA = int t . [[ d_dev . x ]] dA
        FloatArray fe;
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            this->integrateTractionDev(fe, e, boundary, this->devGradient);
            fe.negated();

            answer.assemble(fe, t_loc);
            if ( eNorms ) {
                eNorms->assembleSquared(fe, velocityDofIDs);
            }
        }
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ke_v, Ke_e;
        FloatArray fe_v, fe_t, fe_t2, fe_e(1);
        FloatArray t, e, v;

        // Fetch the current values of internal dofs and their master dof ids;
        this->tractionsdman->giveCompleteUnknownVector(t, mode, tStep);
        this->voldman->giveCompleteUnknownVector(e, mode, tStep);

        // Assemble: -int       t . [[ delta_v ]] dA
        //            int delta_t . [[ e.x - v ]] dA
        //            int       t . [[ x ]]       dA delta_e
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            el->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
            el->giveBoundaryLocationArray(v_loc, bNodes, eid, s, & velocityDofIDs);
            el->computeBoundaryVectorOf(bNodes, eid, mode, tStep, v);

            // Integrate the tangents;
            this->integrateTractionVelocityTangent(Ke_v, el, boundary);
            this->integrateTractionXTangent(Ke_e, el, boundary);

            // We just use the tangent, less duplicated code
            fe_v.beTProductOf(Ke_v, t);
            fe_v.negated();
            fe_t.beProductOf(Ke_v, v);
            fe_t.negated();
            fe_t2.beProductOf(Ke_e, e);
            fe_t.add(fe_t2);
            fe_e.beTProductOf(Ke_e, t);

            answer.assemble(fe_v, v_loc); // Contributions to delta_v equations
            answer.assemble(fe_t, t_loc); // Contribution to delta_t equations
            answer.assemble(fe_e, e_loc); // Contribution to delta_e equations
            if ( eNorms ) {
                eNorms->assembleSquared(fe_v, velocityDofIDs);
                eNorms->assembleSquared(fe_t, tractionDofIDs);
                eNorms->assembleSquared(fe_e, dvolDofID);
            }
        }
    }
}


void MixedGradientPressureWeakPeriodic :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                                                   CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }

    if ( eid != EID_MomentumBalance ) {
        return;
    }

    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == StiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke_v, Ke_vT, Ke_e, Ke_eT;
        IntArray v_loc_r, v_loc_c, t_loc_r, t_loc_c, e_loc_r, e_loc_c;
        IntArray bNodes;
        Set *set = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = set->giveBoundaryList();

        // Fetch the columns/rows for the stress contributions;
        this->tractionsdman->giveCompleteLocationArray(t_loc_r, r_s);
        this->tractionsdman->giveCompleteLocationArray(t_loc_c, c_s);

        this->voldman->giveCompleteLocationArray(e_loc_r, r_s);
        this->voldman->giveCompleteLocationArray(e_loc_c, c_s);

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            el->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
            el->giveBoundaryLocationArray(v_loc_r, bNodes, eid, r_s);
            el->giveBoundaryLocationArray(v_loc_c, bNodes, eid, c_s);

            this->integrateTractionVelocityTangent(Ke_v, el, boundary);
            this->integrateTractionXTangent(Ke_e, el, boundary);

            Ke_v.negated();
            Ke_vT.beTranspositionOf(Ke_v);
            Ke_eT.beTranspositionOf(Ke_e);

            answer->assemble(t_loc_r, v_loc_c, Ke_v);
            answer->assemble(v_loc_r, t_loc_c, Ke_vT);

            answer->assemble(t_loc_r, e_loc_c, Ke_e);
            answer->assemble(e_loc_r, t_loc_c, Ke_eT);
        }
    }
}


void MixedGradientPressureWeakPeriodic :: computeFields(FloatArray &sigmaDev, double &vol, EquationID eid, TimeStep *tStep)
{
    int nsd = domain->giveNumberOfSpatialDimensions();
    double rve_size = this->domainSize();
    FloatArray tractions, t, normal, surfCoords, coords;
    FloatMatrix sigma;
    this->tractionsdman->giveCompleteUnknownVector(tractions, VM_Total, tStep);
    ///@todo This would be a lot simpler if I bothered to create orthogonal basis functions for the tractions. If so, it would just be the constant terms.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    // Reminder: sigma = int t * n dA, where t = sum( C_i * n t_i ).
    // This loop will construct sigma in matrix form.
    sigma.resize(3, 3);
    sigma.zero();
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

        int maxorder = this->order + interp->giveInterpolationOrder() * 3 - 2;
        IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);

        surfCoords.resize(nsd - 1);
        for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
            GaussPoint *gp = ir->getIntegrationPoint(i);
            FloatArray &lcoords = * gp->giveCoordinates();
            FEIElementGeometryWrapper cellgeo(el);

            double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
            // Compute v_m = d_dev . x
            interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
            // "pos" loops over all the traction lagrange multipliers.
            int pos = 0;
            for ( int kj = 0; kj < nsd; ++kj ) { // Columns of the basis functions for the tractions
                // First we compute the surface coordinates. Its the global coordinates without the kj component.
                for ( int c = 0, q = 0; c < nsd; ++c ) {
                    if ( c != kj ) {
                        surfCoords(q) = coords(c);
                        q++;
                    }
                }
                // Evaluate all the (scalar) basis functions for a given column.
                this->evaluateTractionBasisFunctions(t, surfCoords);
                for ( int ti = 0; ti < t.giveSize(); ++ti ) {
                    for ( int ki = 0; ki < nsd; ++ki ) { // Columns of the basis functions for the tractions
                        // Order of tractions is very important here! It must in the right order here (consistent with other loops)
                        for ( int kJ = 0; kJ < nsd; ++kJ ) {
                            sigma(ki, kJ) += tractions(pos) * t(ti) * normal(kj) * coords(kJ) * detJ * gp->giveWeight();
                        }
                        //printf("pos = %d, traction = %f, t = %f, normal_j = %f, normal_i = %f\n", pos, tractions(pos), t(ti), normal(kj), normal(ki));
                        pos++;
                    }
                }
            }
        }
        delete ir;
    }
    sigma.times(1. / rve_size);

    double pressure = 0.;
    for ( int i = 1; i <= nsd; i++ ) {
        pressure += sigma.at(i, i);
    }
    pressure /= 3; // Not 100% sure about this for 2D cases.
    if ( nsd == 3 ) {
        sigmaDev.resize(6);
        sigmaDev.at(1) = sigma.at(1, 1) - pressure;
        sigmaDev.at(2) = sigma.at(2, 2) - pressure;
        sigmaDev.at(3) = sigma.at(3, 3) - pressure;
        sigmaDev.at(4) = 0.5 * ( sigma.at(2, 3) + sigma.at(3, 2) );
        sigmaDev.at(5) = 0.5 * ( sigma.at(1, 3) + sigma.at(3, 1) );
        sigmaDev.at(6) = 0.5 * ( sigma.at(1, 2) + sigma.at(2, 1) );
    } else if ( nsd == 2 ) {
        sigmaDev.resize(3);
        sigmaDev.at(1) = sigma.at(1, 1) - pressure;
        sigmaDev.at(2) = sigma.at(2, 2) - pressure;
        sigmaDev.at(3) = 0.5 * ( sigma.at(1, 2) + sigma.at(2, 1) );
    } else {
        sigmaDev.resize(1);
        sigmaDev.at(1) = sigma.at(1, 1) - pressure;
    }

    // And the volumetric part is much easier.
    vol = this->voldman->giveDof(1)->giveUnknown(VM_Total, tStep);
    vol /= rve_size;
    vol -= volGradient; // This is needed for consistency; We return the volumetric "residual" if a gradient with volumetric contribution is set.
}


void MixedGradientPressureWeakPeriodic :: computeTangents(FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, EquationID eid, TimeStep *tStep)
{
    //double size = this->domainSize();
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    SparseLinearSystemNM *solver = classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ); // = rve->giveLinearSolver();
    SparseMtrx *Kff;
    SparseMtrxType stype = SMT_PetscMtrx; // = rve->giveSparseMatrixType();
    EModelDefaultEquationNumbering fnum;
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    double rve_size = this->domainSize();

    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    Kff = classFactory.createSparseMtrx(stype);
    if ( !Kff ) {
        OOFEM_ERROR2("MixedGradientPressureWeakPeriodic :: computeTangents - Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure(rve, this->domain->giveNumber(), eid, fnum);
    rve->assemble(Kff, tStep, eid, StiffnessMatrix, fnum, fnum, this->domain);

    // Setup up indices and locations
    int neq = Kff->giveNumberOfRows();

    // Indices and such of internal dofs
    int nsd = this->devGradient.giveNumberOfColumns();
    int nd = nsd * ( 1 + nsd ) / 2;

    IntArray t_loc, e_loc;

    // Fetch the positions;
    this->tractionsdman->giveCompleteLocationArray( t_loc, EModelDefaultEquationNumbering() );
    this->voldman->giveCompleteLocationArray( e_loc, EModelDefaultEquationNumbering() );

    // Matrices and arrays for sensitivities
    FloatMatrix ddev_pert(neq, nd);
    FloatArray p_pert(neq);
    FloatMatrix ddev_tmp(nsd, nsd);

    // Solutions for the pertubations:
    FloatMatrix s_d(neq, nd);
    FloatArray s_p(neq);

    // Unit pertubations for d_dev
    FloatArray tmp(neq), fe, ddevVoigt(nd);
    for ( int i = 1; i <= nd; ++i ) {
        ddevVoigt.zero();
        ddevVoigt.at(i) = 1.0 * rve_size;
        this->constructFullMatrixForm(ddev_tmp, ddevVoigt);
        for ( int k = 1; k <= boundaries.giveSize() / 2; ++k ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(k * 2 - 1) );
            int boundary = boundaries.at(k * 2);

            this->integrateTractionDev(fe, e, boundary, ddev_tmp);
            fe.negated();

            tmp.assemble(fe, t_loc);
        }
        ddev_pert.setColumn(tmp, i);
    }

    // Unit pertubation for p
    p_pert.zero();
    p_pert.at( e_loc.at(1) ) = 1.0 * rve_size;

    // Solve all sensitivities
    solver->solve(Kff, ddev_pert, s_d);
    solver->solve(Kff, & p_pert, & s_p);

    // Extract the tractions from the sensitivity solutions s_d and s_p:
    FloatArray tractions_p( t_loc.giveSize() );
    FloatMatrix tractions_d(t_loc.giveSize(), nd);
    tractions_p.beSubArrayOf(s_p, t_loc);
    for ( int i = 1; i <= nd; ++i ) {
        for ( int k = 1; k <= nd; ++k ) {
            tractions_d.at(k, i) = s_d.at(t_loc.at(k), i);
        }
    }

    // The de/dp tangent:
    Cp = s_p.at( e_loc.at(1) ) / rve_size;
    // The de/dd tangent:
    Cd.resize(nd);
    if ( nsd == 3 ) {
        Cd.at(1) = s_d.at(e_loc.at(1), 1);
        Cd.at(2) = s_d.at(e_loc.at(1), 2);
        Cd.at(3) = s_d.at(e_loc.at(1), 3);
        Cd.at(4) = s_d.at(e_loc.at(1), 4);
        Cd.at(5) = s_d.at(e_loc.at(1), 5);
        Cd.at(6) = s_d.at(e_loc.at(1), 6);
    } else if ( nsd == 2 ) {
        Cd.at(1) = s_d.at(e_loc.at(1), 1);
        Cd.at(2) = s_d.at(e_loc.at(1), 2);
        Cd.at(3) = s_d.at(e_loc.at(1), 3);
    }
    // The ds/de tangent:
    Ep = Cd; // This is much simpler, but it's "only" correct if there is a potential (i.e. symmetric problems). This could be generalized if needed.
    // The ds/dd tangent:
    FloatMatrix Ed_tmp(nd, nd);
    FloatArray surfCoords, normal, coords, t;
    FloatMatrix sigma;
    for ( int dpos = 1; dpos <= nd; ++dpos ) {
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

            int maxorder = this->order + interp->giveInterpolationOrder() * 3 - 2;
            IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);

            surfCoords.resize(nsd - 1);
            for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
                GaussPoint *gp = ir->getIntegrationPoint(i);
                FloatArray &lcoords = * gp->giveCoordinates();
                FEIElementGeometryWrapper cellgeo(el);

                double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
                // Compute v_m = d_dev . x
                interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);

                // "pos" loops over all the traction lagrange multipliers.
                int pos = 0;
                for ( int kj = 0; kj < nsd; ++kj ) {
                    // First we compute the surface coordinates. Its the global coordinates without the kj component.
                    for ( int c = 0, q = 0; c < nsd; ++c ) {
                        if ( c != kj ) {
                            surfCoords(q) = coords(c);
                            q++;
                        }
                    }
                    // Evaluate all the basis functions for a given column.
                    this->evaluateTractionBasisFunctions(t, surfCoords);
                    for ( int ti = 0; ti < t.giveSize(); ++ti ) {
                        for ( int ki = 0; ki < nsd; ++ki ) {
                            // Order of tractions is very important here! It must in the right order here (consistent with other loops)
                            for ( int kJ = 0; kJ < nsd; ++kJ ) {
                                sigma(ki, kJ) += tractions_d(pos, dpos) * t(ti) * normal(kj) * coords(kJ) * detJ * gp->giveWeight();
                            }
                            pos++;
                        }
                    }
                }
            }
            delete ir;
        }
        sigma.times(1. / rve_size);

        double pressure = 0.;
        for ( int i = 1; i <= nsd; i++ ) {
            pressure += sigma.at(i, i);
        }
        pressure /= 3; // Not 100% sure about this for 2D cases.
        if ( nsd == 3 ) {
            Ed_tmp.at(1, dpos) = sigma.at(1, 1) - pressure;
            Ed_tmp.at(2, dpos) = 0.5 * ( sigma.at(2, 3) + sigma.at(3, 2) );
            Ed_tmp.at(3, dpos) = 0.5 * ( sigma.at(1, 3) + sigma.at(3, 1) );
            Ed_tmp.at(4, dpos) = sigma.at(2, 2) - pressure;
            Ed_tmp.at(5, dpos) = 0.5 * ( sigma.at(1, 2) + sigma.at(2, 1) );
            Ed_tmp.at(6, dpos) = sigma.at(3, 3) - pressure;
        } else if ( nsd == 2 ) {
            Ed_tmp.at(1, dpos) = sigma.at(1, 1) - pressure;
            Ed_tmp.at(2, dpos) = 0.5 * ( sigma.at(1, 2) + sigma.at(2, 1) );
            Ed_tmp.at(3, dpos) = sigma.at(2, 2) - pressure;
        } else {
            Ed_tmp.at(1, dpos) = sigma.at(1, 1) - pressure;
        }
        dpos++;
    }

    // Reorder to voigt form:
    IntArray indx(nd);
    if ( nsd == 3 ) {
        indx.setValues(6,  1, 4, 6, 5, 3, 2);
    } else if ( nsd == 2 ) {
        indx.setValues(3,  1, 3, 2);
    } else if ( nsd == 1 ) {
        indx.setValues(1,  1);
    }
    Ed.beSubMatrixOf(Ed_tmp, indx, indx);

    delete Kff;
    delete solver;
}

void MixedGradientPressureWeakPeriodic :: giveInputRecord(DynamicInputRecord &input)
{
    MixedGradientPressureBC :: giveInputRecord(input);
    input.setField(this->pressure, _IFT_MixedGradientPressure_pressure);
    OOFEM_ERROR("MixedGradientPressureWeakPeriodic :: giveInputRecord - Not supported yet\n");
    //FloatArray devGradientVoigt;
    //input.setField(devGradientVoigt, _IFT_MixedGradientPressure_devGradient);
}

void MixedGradientPressureWeakPeriodic :: scale(double s)
{
    devGradient.times(s);
    pressure *= s;
}
} // end namespace oofem
