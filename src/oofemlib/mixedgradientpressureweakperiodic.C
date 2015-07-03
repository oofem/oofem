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
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_BoundaryCondition(MixedGradientPressureWeakPeriodic);

MixedGradientPressureWeakPeriodic :: MixedGradientPressureWeakPeriodic(int n, Domain *d) : MixedGradientPressureBC(n, d),
    voldman( new Node(-1, d) ),
    tractionsdman( new Node(-1, d) )
{
    int dofid = this->domain->giveNextFreeDofID();
    v_id.followedBy(dofid);
    this->voldman->appendDof( new MasterDof( voldman.get(), ( DofIDItem ) dofid ) );
}


MixedGradientPressureWeakPeriodic :: ~MixedGradientPressureWeakPeriodic()
{
}


int MixedGradientPressureWeakPeriodic :: giveNumberOfInternalDofManagers()
{
    return 2;
}


DofManager *MixedGradientPressureWeakPeriodic :: giveInternalDofManager(int i)
{
    if ( i == 1 ) {
        return this->tractionsdman.get();
    } else {
        return this->voldman.get();
    }
}


IRResultType MixedGradientPressureWeakPeriodic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;


    IR_GIVE_FIELD(ir, this->order, _IFT_MixedGradientPressureWeakPeriodic_order);
    if ( this->order < 0 ) {
        OOFEM_ERROR("order must be at least 0");
    }

    int nsd = this->domain->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * ( int ) pow(double ( order + 1 ), nsd - 1);
    this->tractionsdman->setNumberOfDofs(0);
    t_id.clear();
    for ( int i = 1; i <= total; ++i ) {
        int dofid = this->domain->giveNextFreeDofID();
        t_id.followedBy(dofid);
        // Simply use t_i = S_i . n, where S_1 = [1,0,0;0,0,0;0,0,0], S_2 = [0,1,0;0,0,0;0,0,0], etc.
        // then the linear terms, [x,0,0], [0,x,0] ... and so on
        this->tractionsdman->appendDof( new MasterDof( tractionsdman.get(), ( DofIDItem ) dofid ) );
    }

    return MixedGradientPressureBC :: initializeFrom(ir);
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
        d.at(2, 3) = d.at(3, 2) = d_voigt.at(4) * 0.5;
        d.at(1, 3) = d.at(3, 1) = d_voigt.at(5) * 0.5;
        d.at(1, 2) = d.at(2, 1) = d_voigt.at(6) * 0.5;
    } else if ( n == 3 ) { // Then 2D
        d.resize(2, 2);
        d.at(1, 1) = d_voigt.at(1);
        d.at(2, 2) = d_voigt.at(2);
        d.at(1, 2) = d.at(2, 1) = d_voigt.at(3) * 0.5;
    } else if ( n == 1 ) { // Then 1D
        d.resize(1, 1);
        d.at(1, 1) = d_voigt.at(1);
    } else {
        OOFEM_ERROR("Tensor is in strange voigt format.");
    }
}


void MixedGradientPressureWeakPeriodic :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                             const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( type != TangentStiffnessMatrix ) {
        return;
    }

    IntArray loc_r, loc_c, t_loc_r, t_loc_c, e_loc_r, e_loc_c;

    // Fetch the columns/rows for the tractions;
    this->tractionsdman->giveLocationArray(t_id, t_loc_r, r_s);
    this->tractionsdman->giveLocationArray(t_id, t_loc_c, c_s);

    // Fetch the columns/rows for dvol;
    this->voldman->giveLocationArray(v_id, e_loc_r, r_s);
    this->voldman->giveLocationArray(v_id, e_loc_c, c_s);

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
        e->giveBoundaryLocationArray(loc_r, bNodes, this->dofs, r_s);
        e->giveBoundaryLocationArray(loc_c, bNodes, this->dofs, c_s);
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


void MixedGradientPressureWeakPeriodic :: evaluateTractionBasisFunctions(FloatArray &answer, const FloatArray &surfCoords)
{
    // This evaluates the function x^a * y^b, for a, b in [0,order]
    int total = ( int ) pow(double ( order + 1 ), surfCoords.giveSize());
    answer.resize(total);
    int pos = 0;
    double tx = 1.;
    for ( int xi = 0; xi < ( this->order + 1 ); ++xi ) {
        double ty = 1.;
        for ( int yi = 0; yi < ( this->order + 1 ); ++yi ) {
            answer(pos) = tx * ty;
            pos++;
            ty *= surfCoords(1);
        }
        tx *= surfCoords(0);
    }
}


void MixedGradientPressureWeakPeriodic :: constructMMatrix(FloatMatrix &mMatrix, FloatArray &coords, FloatArray &normal)
{
    FloatArray t, surfCoords;
    int nsd = this->giveDomain()->giveNumberOfSpatialDimensions();
    int total = nsd * nsd * ( int ) pow(double ( this->order + 1 ), nsd - 1);

    surfCoords.resize(nsd - 1);
    mMatrix.resize(nsd, total);
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
}


void MixedGradientPressureWeakPeriodic :: integrateTractionVelocityTangent(FloatMatrix &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dv dA
    FloatArray normal, n, coords;
    FloatMatrix nMatrix, mMatrix;

    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder() * 3;
    std :: unique_ptr< IntegrationRule >ir( interp->giveBoundaryIntegrationRule(maxorder, boundary) );
    int nsd = this->giveDomain()->giveNumberOfSpatialDimensions();

    answer.clear();
    for ( GaussPoint *gp: *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        interp->boundaryEvalN(n, boundary, lcoords, cellgeo);
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);

        // Construct the basis functions for the tractions:
        this->constructMMatrix(mMatrix, coords, normal);
        nMatrix.beNMatrixOf(n, nsd);

        answer.plusProductUnsym( mMatrix, nMatrix, detJ * gp->giveWeight() );
    }
}


void MixedGradientPressureWeakPeriodic :: integrateTractionXTangent(FloatMatrix &answer, Element *el, int boundary)
{
    // Computes the integral: int dt . dx_m dA
    FloatMatrix mMatrix;
    FloatArray normal, coords, vM_vol;

    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder() * 3;
    std :: unique_ptr< IntegrationRule >ir( interp->giveBoundaryIntegrationRule(maxorder, boundary) );

    FloatArray tmpAnswer;
    for ( GaussPoint *gp: *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);

        vM_vol.beScaled(1.0/3.0, coords);
        this->constructMMatrix(mMatrix, coords, normal);

        tmpAnswer.plusProduct(mMatrix, vM_vol, detJ * gp->giveWeight());
    }
    answer.resize(tmpAnswer.giveSize(), 1);
    answer.setColumn(tmpAnswer, 1);
}


void MixedGradientPressureWeakPeriodic :: integrateTractionDev(FloatArray &answer, Element *el, int boundary, const FloatMatrix &ddev)
{
    // Computes the integral: int dt . dx dA
    FloatMatrix mMatrix;
    FloatArray normal, coords, vM_dev;

    FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

    int maxorder = this->order + interp->giveInterpolationOrder() * 3;
    std :: unique_ptr< IntegrationRule >ir( interp->giveBoundaryIntegrationRule(maxorder, boundary) );
    answer.clear();

    for ( GaussPoint *gp: *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(el);

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        // Compute v_m = d_dev . x
        interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
        vM_dev.beProductOf(ddev, coords);

        this->constructMMatrix(mMatrix, coords, normal);

        answer.plusProduct(mMatrix, vM_dev, detJ * gp->giveWeight());
    }
}


void MixedGradientPressureWeakPeriodic :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                         CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    IntArray v_loc, t_loc, e_loc;  // For the velocities and stress respectively
    IntArray velocityDofIDs, bNodes;
    this->tractionsdman->giveLocationArray(t_id, t_loc, s);
    this->voldman->giveLocationArray(v_id, e_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. on the traction and on dvol.
        double rve_size = this->domainSize();

        if ( e_loc.at(1) ) {
            answer.at( e_loc.at(1) ) -= rve_size * pressure; // Note the negative sign (pressure as opposed to mean stress)
            if ( eNorms ) {
                eNorms->at( v_id.at(1) ) += rve_size * pressure * rve_size * pressure;
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
        this->tractionsdman->giveUnknownVector(t, t_id, mode, tStep);
        this->voldman->giveUnknownVector(e, v_id, mode, tStep);

        // Assemble: -int       t . [[ delta_v ]] dA
        //            int delta_t . [[ e.x - v ]] dA
        //            int       t . [[ x ]]       dA delta_e
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            el->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
            el->giveBoundaryLocationArray(v_loc, bNodes, this->dofs, s, & velocityDofIDs);
            el->computeBoundaryVectorOf(bNodes, this->dofs, mode, tStep, v);

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
                eNorms->assembleSquared(fe_t, t_id);
                eNorms->assembleSquared(fe_e, v_id);
            }
        }
    }
}


void MixedGradientPressureWeakPeriodic :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                                   CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke_v, Ke_vT, Ke_e, Ke_eT;
        IntArray v_loc_r, v_loc_c, t_loc_r, t_loc_c, e_loc_r, e_loc_c;
        IntArray bNodes;
        Set *set = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = set->giveBoundaryList();

        // Fetch the columns/rows for the stress contributions;
        this->tractionsdman->giveLocationArray(t_id, t_loc_r, r_s);
        this->tractionsdman->giveLocationArray(t_id, t_loc_c, c_s);

        this->voldman->giveLocationArray(v_id, e_loc_r, r_s);
        this->voldman->giveLocationArray(v_id, e_loc_c, c_s);

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            el->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
            el->giveBoundaryLocationArray(v_loc_r, bNodes, this->dofs, r_s);
            el->giveBoundaryLocationArray(v_loc_c, bNodes, this->dofs, c_s);

            this->integrateTractionVelocityTangent(Ke_v, el, boundary);
            this->integrateTractionXTangent(Ke_e, el, boundary);

            Ke_v.negated();
            Ke_vT.beTranspositionOf(Ke_v);
            Ke_eT.beTranspositionOf(Ke_e);

            answer.assemble(t_loc_r, v_loc_c, Ke_v);
            answer.assemble(v_loc_r, t_loc_c, Ke_vT);

            answer.assemble(t_loc_r, e_loc_c, Ke_e);
            answer.assemble(e_loc_r, t_loc_c, Ke_eT);
        }
    }
}


void MixedGradientPressureWeakPeriodic :: computeFields(FloatArray &sigmaDev, double &vol, TimeStep *tStep)
{
    double rve_size = this->domainSize();
    FloatArray tractions;
    this->tractionsdman->giveUnknownVector(tractions, t_id, VM_Total, tStep);
    this->computeStress(sigmaDev, tractions, rve_size);

    // And the volumetric part is much easier.
    vol = (*this->voldman->begin())->giveUnknown(VM_Total, tStep);
    vol /= rve_size;
    vol -= volGradient; // This is needed for consistency; We return the volumetric "residual" if a gradient with volumetric contribution is set.
}


void MixedGradientPressureWeakPeriodic :: computeStress(FloatArray &sigmaDev, FloatArray &tractions, double rve_size)
{
    FloatMatrix mMatrix;
    FloatArray normal, coords, t;

    int nsd = domain->giveNumberOfSpatialDimensions();
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    // Reminder: sigma = int t * n dA, where t = sum( C_i * n t_i ).
    // This loop will construct sigma in matrix form.

    FloatMatrix sigma;

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *el = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        FEInterpolation *interp = el->giveInterpolation(); // Geometry interpolation. The displacements or velocities must have the same interpolation scheme (on the boundary at least).

        int maxorder = this->order + interp->giveInterpolationOrder() * 3;
        std :: unique_ptr< IntegrationRule >ir( interp->giveBoundaryIntegrationRule(maxorder, boundary) );

        for ( GaussPoint *gp: *ir ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            FEIElementGeometryWrapper cellgeo(el);

            double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
            // Compute v_m = d_dev . x
            interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);

            this->constructMMatrix(mMatrix, coords, normal);
            t.beProductOf(mMatrix, tractions);
            sigma.plusDyadUnsym(t, coords, detJ * gp->giveWeight());
        }
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
}

void MixedGradientPressureWeakPeriodic :: computeTangents(FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, TimeStep *tStep)
{
    //double size = this->domainSize();
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM > solver( 
        classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    EModelDefaultEquationNumbering fnum;
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    double rve_size = this->domainSize();

    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx(stype) );
    if ( !Kff ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure(rve, this->domain->giveNumber(), fnum);
    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, fnum, this->domain);

    // Setup up indices and locations
    int neq = Kff->giveNumberOfRows();

    // Indices and such of internal dofs
    int nsd = this->devGradient.giveNumberOfColumns();
    int nd = nsd * ( 1 + nsd ) / 2;

    IntArray t_loc, e_loc;

    // Fetch the positions;
    this->tractionsdman->giveLocationArray( t_id, t_loc, fnum );
    this->voldman->giveLocationArray( v_id, e_loc, fnum );

    // Matrices and arrays for sensitivities
    FloatMatrix ddev_pert(neq, nd);
    FloatArray p_pert(neq);

    // Solutions for the pertubations:
    FloatMatrix s_d(neq, nd);
    FloatArray s_p(neq);

    // Unit pertubations for d_dev
    for ( int i = 1; i <= nd; ++i ) {
        FloatMatrix ddev_tmp;
        FloatArray tmp(neq), fe, fetot, ddevVoigt(nd);
        ddevVoigt.at(i) = 1.0;
        this->constructFullMatrixForm(ddev_tmp, ddevVoigt);
        for ( int k = 1; k <= boundaries.giveSize() / 2; ++k ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(k * 2 - 1) );
            int boundary = boundaries.at(k * 2);

            this->integrateTractionDev(fe, e, boundary, ddev_tmp);
            fetot.subtract(fe);
        }
        tmp.assemble(fetot, t_loc);
        ddev_pert.setColumn(tmp, i);
    }

    // Unit pertubation for p
    p_pert.zero();
    p_pert.at( e_loc.at(1) ) = - 1.0 * rve_size;

    // Solve all sensitivities
    solver->solve(*Kff, ddev_pert, s_d);
    solver->solve(*Kff, p_pert, s_p);

    // Extract the tractions from the sensitivity solutions s_d and s_p:
    FloatArray tractions_p( t_loc.giveSize() );
    FloatMatrix tractions_d(t_loc.giveSize(), nd);
    tractions_p.beSubArrayOf(s_p, t_loc);
    for ( int i = 1; i <= nd; ++i ) {
        for ( int k = 1; k <= t_loc.giveSize(); ++k ) {
            tractions_d.at(k, i) = s_d.at(t_loc.at(k), i);
        }
    }

    // The de/dp tangent:
    Cp = - s_p.at( e_loc.at(1) );
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
    Ed.resize(nd, nd);
    for ( int dpos = 1; dpos <= nd; ++dpos ) {
        FloatArray tractions, sigmaDev;
        tractions.beColumnOf(tractions_d, dpos);
        this->computeStress(sigmaDev, tractions, rve_size);
        Ed.setColumn(sigmaDev, dpos);
    }
}

void MixedGradientPressureWeakPeriodic :: giveInputRecord(DynamicInputRecord &input)
{
    MixedGradientPressureBC :: giveInputRecord(input);
    input.setField(this->pressure, _IFT_MixedGradientPressure_pressure);
    OOFEM_ERROR("Not supported yet");
    //FloatArray devGradientVoigt;
    //input.setField(devGradientVoigt, _IFT_MixedGradientPressure_devGradient);
}

void MixedGradientPressureWeakPeriodic :: scale(double s)
{
    devGradient.times(s);
    pressure *= s;
}
} // end namespace oofem
