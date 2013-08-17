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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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

REGISTER_BoundaryCondition( MixedGradientPressureWeakPeriodic );

MixedGradientPressureWeakPeriodic :: MixedGradientPressureWeakPeriodic(int n, Domain *d) : MixedGradientPressureBC(n,d)
{
    this->voldman = new Node(0, d); // Node number lacks meaning here.
    this->voldman->appendDof(new MasterDof(1, voldman, (DofIDItem)d->giveNextFreeDofID()));
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
    int total = nsd * nsd * pow(order + 1, nsd - 1);
    this->tractionsdman->setNumberOfDofs( total );
    for ( int i = 1; i <= total; ++i ) {
        // Simply use t_i = S_i . n, where S_1 = [1,0,0;0,0,0;0,0,0], S_2 = [0,1,0;0,0,0;0,0,0], etc.
        // then the linear terms, [x,0,0], [0,x,0] ... and so on
        this->tractionsdman->setDof(i, new MasterDof(i, tractionsdman, (DofIDItem)this->domain->giveNextFreeDofID() ) );
    }

    return IRRT_OK;
}


void MixedGradientPressureWeakPeriodic :: setPrescribedDeviatoricGradientFromVoigt(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 1 ) { // Then 1D
        this->devGradient.resize(1, 1);
        this->devGradient.at(1, 1) = t.at(1);
    } else if ( n == 3 ) { // Then 2D
        this->devGradient.resize(2, 2);
        this->devGradient.at(1, 1) = t.at(1);
        this->devGradient.at(2, 2) = t.at(2);
        this->devGradient.at(1, 2) = this->devGradient.at(2, 1) = t.at(3);
    } else if ( n == 6 ) { // Then 3D
        this->devGradient.resize(3, 3);
        this->devGradient.at(1, 1) = t.at(1);
        this->devGradient.at(2, 2) = t.at(2);
        this->devGradient.at(3, 3) = t.at(3);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->devGradient.at(1, 2) = this->devGradient.at(2, 1) = t.at(6) * 0.5;
        this->devGradient.at(1, 3) = this->devGradient.at(3, 1) = t.at(5) * 0.5;
        this->devGradient.at(2, 3) = this->devGradient.at(3, 2) = t.at(4) * 0.5;
    } else {
        OOFEM_ERROR("MixedGradientPressureWeakPeriodic :: setPrescribedDeviatoricGradientFromVoigt: Tensor is in strange voigt format. Should be 3 or 6. Use setPrescribedTensor directly if needed.");
    }
}


void MixedGradientPressureWeakPeriodic :: giveLocationArrays(std::vector<IntArray> &rows, std::vector<IntArray> &cols, EquationID eid, CharType type,
    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance || type != TangentStiffnessMatrix) {
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

    rows.resize(boundaries.giveSize()+2);
    cols.resize(boundaries.giveSize()+2);
    int i = 0;
    for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
        int boundary = boundaries.at(pos*2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
        e->giveBoundaryLocationArray(loc_r, bNodes, eid, r_s);
        e->giveBoundaryLocationArray(loc_c, bNodes, eid, c_s);
        // For most uses, *loc_r == *loc_c
        rows[i] = loc_r;
        cols[i] = t_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows[i] = t_loc_r;
        cols[i] = loc_c;
        i++;
    }

    // The volumetric part:
    rows[i] = t_loc_r;
    cols[i] = e_loc_c;
    i++;

    rows[i] = e_loc_r;
    cols[i] = t_loc_c;
    i++;
}


void MixedGradientPressureWeakPeriodic :: integrateTractionVelocityTangent(FloatMatrix &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dv dA
    FloatArray normal, n, m, t, surfCoords, coords, v_m, contrib;
    FloatMatrix nMatrix, mMatrix;
    IntArray boundaryNodes;
    
    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    ///@todo Check the order here:
    int maxorder = this->order + interp->giveInterpolationOrder()*3 - 2;
    IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);
    int nsd = el->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * pow(order + 1, nsd - 1);

    surfCoords.resize( nsd - 1 );
    mMatrix.resize( nsd, total );

    answer.resize(0,0);
    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = *gp->giveCoordinates();
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
            // Then pass those coordinates to the function which evaluates all the basis functions at that point.
            this->evaluateTractionBasisFunctions(t, surfCoords);
            for ( int ti = 0; ti < t.giveSize(); ++ti ) {
                for ( int ki = 0; ki < nsd; ++ki ) {
                    mMatrix(ki, pos) = t(ti) * normal(kj);
                    pos++;
                }
            }
        }
        nMatrix.beNMatrixOf(n, nsd);

        answer.plusProductUnsym(mMatrix, nMatrix, detJ * gp->giveWeight());
    }
    delete ir;
}


void MixedGradientPressureWeakPeriodic :: integrateTractionXTangent(FloatMatrix &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dx_m dA
    FloatArray normal, t, surfCoords, coords, v_m, contrib;
    IntArray boundaryNodes;
    
    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder()*3 - 2;
    IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);
    int nsd = el->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * pow(order + 1, nsd - 1);

    surfCoords.resize( nsd - 1 );

    FloatArray tmpAnswer;
    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = *gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
        contrib.resize( total );
        int pos = 0;
        for ( int kj = 0; kj < nsd; ++kj ) {
            // First we compute the surface coordinates. Its the global coordinates without the kj component.
            for ( int c = 0, q = 0; c < nsd; ++c ) {
                if ( c != kj ) {
                    surfCoords(q) = coords(c);
                    q++;
                }
            }
            // Then pass those coordinates to the function which evaluates all the basis functions at that point.
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


void MixedGradientPressureWeakPeriodic :: integrateTractionDev(FloatArray &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dx dA
    FloatArray normal, t, surfCoords, coords, vM_dev, contrib;
    IntArray boundaryNodes;
    
    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder()*3 - 2;
    IntegrationRule *ir = interp->giveBoundaryIntegrationRule(maxorder, boundary);
    int nsd = el->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * pow(order + 1, nsd - 1);

    surfCoords.resize(nsd - 1);
    contrib.resize( total );
    answer.resize(0);

    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i ) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = *gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        // Compute v_m = d_dev . x
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
        vM_dev.beProductOf(this->devGradient, coords);

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
            // Then pass those coordinates to the function which evaluates all the basis functions at that point.
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
    int total = pow(order + 1, nsd - 1);
    answer.resize(total);
    int pos = 0;
    double tx = 1.;
    for ( int xi = 0; xi < (this->order + 1); ++xi ) {
        double ty = 1.;
        for ( int yi = 0; yi < (this->order + 1); ++yi ) {
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
    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance)
        return;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    IntArray v_loc, t_loc, e_loc;  // For the velocities and stress respectively
    IntArray velocityDofIDs, tractionDofIDs, dvolDofID, bNodes;
    this->tractionsdman->giveCompleteLocationArray(t_loc, s);
    this->voldman->giveCompleteLocationArray(e_loc, s);
    this->tractionsdman->giveCompleteMasterDofIDArray(tractionDofIDs);
    this->voldman->giveCompleteMasterDofIDArray(dvolDofID);

    if (type == ExternalForcesVector) {
        // The external forces have two contributions. on the traction and on dvol.
        double rve_size = this->domainSize();

        if ( e_loc.at(1) ) {
            answer.at(e_loc.at(1)) -= rve_size*pressure; // Note the negative sign (pressure as opposed to mean stress)
            if ( eNorms ) eNorms->at(dvolDofID.at(1)) = rve_size*pressure*rve_size*pressure;
        }

        // The second contribution is on the momentumbalance equation; int t . [[ d_dev . x ]] dA = int t . [[ d_dev . x ]] dA
        FloatArray fe;
        for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
            int boundary = boundaries.at(pos*2);

            e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
            e->giveBoundaryLocationArray(v_loc, bNodes, eid, s, &velocityDofIDs);
            this->integrateTractionDev(fe, e, boundary);
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
        for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
            int boundary = boundaries.at(pos*2);
            
            // Fetch the element information;
            el->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
            el->giveBoundaryLocationArray(v_loc, bNodes, eid, s, &velocityDofIDs);
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
    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance)
        return;

    if (type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == StiffnessMatrix || type == ElasticStiffnessMatrix) {
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

        for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
            int boundary = boundaries.at(pos*2);
            
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
#if 0
    int nsd = this->giveDomain()->giveNumberOfSpatialDimensions();
    FloatArray sigmaDevBase;
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    // Fetch the current values of the stress in deviatoric base;
    sigmaDevBase.resize(this->sigmaDev->giveNumberOfDofs());
    for ( int i = 1; i <= this->sigmaDev->giveNumberOfDofs(); i++ ) {
        sigmaDevBase.at(i) = this->sigmaDev->giveDof(i)->giveUnknown(VM_Total, tStep);
    }
    // Convert it back from deviatoric base:
    if (nsd == 3) {
        this->fromDeviatoricBase3D(sigmaDev, sigmaDevBase);
    } else if (nsd == 2) {
        this->fromDeviatoricBase2D(sigmaDev, sigmaDevBase);
    } else {
        sigmaDev.resize(0);
    }

    // Postprocessing; vol = int v . n dA
    FloatArray unknowns, fe;
    vol = 0.;
    for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
        int boundary = boundaries.at(pos*2);
        
        e->computeVectorOf(eid, VM_Total, tStep, unknowns);
        this->integrateVolTangent(fe, e, boundary);
        vol += fe.dotProduct(unknowns);
    }
    double rve_size = this->domainSize();
    vol /= rve_size;
    vol -= volGradient; // This is needed for consistency; We return the volumetric "residual" if a gradient with volumetric contribution is set.
#endif
}


void MixedGradientPressureWeakPeriodic :: computeTangents(
    FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, EquationID eid, TimeStep *tStep)
{
#if 0
    //double size = this->domainSize();
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    SparseLinearSystemNM *solver = classFactory.createSparseLinSolver(ST_Petsc, this->domain, this->domain->giveEngngModel());// = rve->giveLinearSolver();
    SparseMtrx *Kff;
    SparseMtrxType stype = SMT_PetscMtrx;// = rve->giveSparseMatrixType();
    EModelDefaultEquationNumbering fnum;
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    double rve_size = this->domainSize();
    
    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    Kff = classFactory.createSparseMtrx(stype);
    if ( !Kff ) {
        OOFEM_ERROR2("MixedGradientPressureWeakPeriodic :: computeTangents - Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure( rve, this->domain->giveNumber(), eid, fnum );
    rve->assemble(Kff, tStep, eid, StiffnessMatrix, fnum, fnum, this->domain );

    // Setup up indices and locations
    int neq = Kff->giveNumberOfRows();

    // Indices and such of internal dofs
    int ndev = this->sigmaDev->giveNumberOfDofs();

    // Matrices and arrays for sensitivities
    FloatMatrix ddev_pert(neq,ndev); // In fact, npeq should most likely equal ndev
    FloatArray p_pert(neq); // RHS for d_dev [d_dev11, d_dev22, d_dev12] in 2D

    FloatMatrix s_d(neq,ndev); // Sensitivity fields for d_dev
    FloatArray s_p(neq); // Sensitivity fields for p

    // Unit pertubations for d_dev
    ddev_pert.zero();
    for (int i = 1; i <= ndev; ++i) {
        int eqn = this->sigmaDev->giveDof(i)->giveEquationNumber(fnum);
        ddev_pert.at(eqn, i) = -1.0*rve_size;
    }
    
    // Unit pertubation for d_p
    p_pert.zero();
    FloatArray fe;
    IntArray loc;
    for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
        int boundary = boundaries.at(pos*2);

        e->giveBoundaryLocationArray(loc, bNodes, eid, fnum);
        this->integrateVolTangent(fe, e, boundary);
        fe.times(-1.0); // here d_p = 1.0
        p_pert.assemble(fe, loc);
    }

    // Solve all sensitivities
    solver->solve(Kff,ddev_pert,s_d);
    solver->solve(Kff,&p_pert,&s_p);

    // Extract the stress response from the solutions
    FloatArray sigma_p(ndev);
    FloatMatrix sigma_d(ndev,ndev);
    for (int i = 1; i <= ndev; ++i) {
        int eqn = this->sigmaDev->giveDof(i)->giveEquationNumber(fnum);
        sigma_p.at(i) = s_p.at(eqn);
        for (int j = 1; j <= ndev; ++j) {
            sigma_d.at(i,j) = s_d.at(eqn,j);
        }
    }

    // Post-process the volumetric rate of deformations in the sensitivity fields;
    FloatArray e_d(ndev); e_d.zero();
    double e_p = 0.0;
    for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
        int boundary = boundaries.at(pos*2);
    
        this->integrateVolTangent(fe, e, boundary);
        e->giveBoundaryLocationArray(loc, bNodes, eid, fnum);
        
        // Using "loc" to pick out the relevant contributions. This won't work at all if there are local coordinate systems in these nodes
        // or slave nodes etc. The goal is to compute the velocity from the sensitivity field, but we need to avoid going through the actual
        // engineering model. If this ever becomes an issue it needs to perform the same steps as Element::giveUnknownVector does.
        for (int i = 1; i <= fe.giveSize(); ++i) {
            if (loc.at(i) > 0) {
                e_p += fe.at(i) * s_p.at(loc.at(i));
                for (int j = 1; j <= ndev; ++j) {
                    e_d.at(j) += fe.at(i) * s_d.at(loc.at(i), j);
                }
            }
        }
    }
    e_p /= rve_size;
    e_d.times(1./rve_size);

    // Now we need to express the tangents in the normal cartesian coordinate system (as opposed to the deviatoric base we use during computations
    Cp = e_p; // Scalar components are of course the same
    double nsd = this->giveDomain()->giveNumberOfSpatialDimensions();
    if (nsd == 3) {
        this->fromDeviatoricBase3D(Cd, e_d);
        this->fromDeviatoricBase3D(Ep, sigma_p);
        this->fromDeviatoricBase3D(Ed, sigma_d);
        
    } else if (nsd == 2) {
        this->fromDeviatoricBase2D(Cd, e_d);
        this->fromDeviatoricBase2D(Ep, sigma_p);
        this->fromDeviatoricBase2D(Ed, sigma_d);

    } else { // For 1D case, there simply are no deviatoric components!
        Cd.resize(0);
        Ep.resize(0);
        Ed.beEmptyMtrx();
    }

    delete Kff;
    delete solver; ///@todo Remove this when solver is taken from engngmodel
#endif
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

