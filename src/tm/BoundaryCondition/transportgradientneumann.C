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

#include "tm/BoundaryCondition/transportgradientneumann.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "tm/Elements/transportelement.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition(TransportGradientNeumann);

TransportGradientNeumann :: TransportGradientNeumann(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    //PrescribedGradientHomogenization(),
    mpFluxHom( std::make_unique<Node>(0, d) )
{
    int nsd = d->giveNumberOfSpatialDimensions();
    for ( int i = 0; i < nsd; i++ ) {
        // Just putting in X_i id-items since they don't matter.
        int dofId = d->giveNextFreeDofID();
        mFluxIds.followedBy(dofId);
        mpFluxHom->appendDof( new MasterDof( mpFluxHom.get(), ( DofIDItem ) ( dofId ) ) );
    }
}


void TransportGradientNeumann :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, mGradient, _IFT_TransportGradientNeumann_gradient);

    IR_GIVE_FIELD(ir, surfSets, _IFT_TransportGradientNeumann_surfSets);
    this->mCenterCoord.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, mCenterCoord, _IFT_TransportGradientNeumann_centerCoords)
    
    this->dispControl = ir.hasField(_IFT_TransportGradientNeumann_dispControl);    
}


void TransportGradientNeumann :: giveInputRecord(DynamicInputRecord &input)
{
    //ActiveBoundaryCondition :: giveInputRecord(input);
    //PrescribedGradientHomogenization :: giveInputRecord(input);
    input.setField(mGradient, _IFT_TransportGradientNeumann_gradient);
}


void TransportGradientNeumann :: postInitialize()
{
    ActiveBoundaryCondition :: postInitialize();
    
    if ( this->dispControl ) this->computeEta();
}


DofManager *TransportGradientNeumann :: giveInternalDofManager(int i)
{
    return mpFluxHom.get();
}

void TransportGradientNeumann :: scale(double s)
{
    this->mGradient.times(s);
}

double TransportGradientNeumann :: domainSize()
{
    Domain *domain = this->giveDomain();
    int nsd = domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    for ( auto &setNum : this->surfSets ) {
        Set *set = domain->giveSet(setNum);
        const IntArray &boundaries = set->giveBoundaryList();

        for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
            Element *e = domain->giveElement( boundaries[pos * 2] );
            int boundary = boundaries[pos * 2 + 1];
            FEInterpolation *fei = e->giveInterpolation();
            domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
        }
    }
    return fabs(domain_size / nsd);
}

void TransportGradientNeumann :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                   CharType type, ValueModeType mode,
                                                   const UnknownNumberingScheme &s, 
                                                   FloatArray *eNorm,
                                                   void* lock)
{
    IntArray flux_loc;  // For the displacements and flux respectively
    mpFluxHom->giveLocationArray(mFluxIds, flux_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for flux, the load is simply the prescribed gradient.
        double rve_size = this->domainSize();
        FloatArray fluxLoad;

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        fluxLoad.beScaled(-rve_size*loadLevel, mGradient);
#ifdef _OPENMP
        if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
        answer.assemble(fluxLoad, flux_loc);
#ifdef _OPENMP
        if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ke;
        FloatArray fe_v, fe_s;
        FloatArray fluxHom, e_u;
        IntArray loc, masterDofIDs, fluxMasterDofIDs;

        // Fetch the current values of the flux;
        mpFluxHom->giveUnknownVector(fluxHom, mFluxIds, mode, tStep);
        // and the master dof ids for flux used for the internal norms
        mpFluxHom->giveMasterDofIDArray(mFluxIds, fluxMasterDofIDs);

        // Assemble
        for ( int i = 0; i < this->surfSets.giveSize(); ++i ) {
            int surfSet = this->surfSets[i];
            Set *setPointer = this->giveDomain()->giveSet(surfSet);
            const IntArray &boundaries = setPointer->giveBoundaryList();
            for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
                Element *e = this->giveDomain()->giveElement( boundaries[pos * 2] );
                int boundary = boundaries[pos * 2 + 1];

                // Fetch the element information;
                const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary);
                e->giveBoundaryLocationArray(loc, bNodes, this->dofs, s, & masterDofIDs);
                e->computeBoundaryVectorOf(bNodes, this->dofs, mode, tStep, e_u);
                this->integrateTangent(Ke, e, boundary, i, pos);

                // We just use the tangent, less duplicated code (the addition of flux is linear).
                fe_v.beProductOf(Ke, e_u);
                fe_s.beTProductOf(Ke, fluxHom);

                // Note: The terms appear negative in the equations:
                fe_v.negated();
                fe_s.negated();
#ifdef _OPENMP
                if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
                answer.assemble(fe_s, loc); // Contributions to delta_v equations
                answer.assemble(fe_v, flux_loc); // Contribution to delta_s_i equations
                if ( eNorm != NULL ) {
                    eNorm->assembleSquared(fe_s, masterDofIDs);
                    eNorm->assembleSquared(fe_v, fluxMasterDofIDs);
                }
#ifdef _OPENMP
                if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
            }
        }
    }
}

void TransportGradientNeumann :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                          CharType type, const UnknownNumberingScheme &r_s, 
                                          const UnknownNumberingScheme &c_s, double scale,
                                          void* lock)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        IntArray loc_r, loc_c, flux_loc_r, flux_loc_c;

        // Fetch the columns/rows for the flux contributions;
        mpFluxHom->giveLocationArray(mFluxIds, flux_loc_r, r_s);
        mpFluxHom->giveLocationArray(mFluxIds, flux_loc_c, c_s);

        for ( int i = 0; i < this->surfSets.giveSize(); ++i ) {
            int surfSet = this->surfSets[i];
            Set *set = this->giveDomain()->giveSet(surfSet);
            const IntArray &boundaries = set->giveBoundaryList();
            for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
                Element *e = this->giveDomain()->giveElement( boundaries[pos * 2] );
                int boundary = boundaries[pos * 2 + 1];

                const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary);
                e->giveBoundaryLocationArray(loc_r, bNodes, this->dofs, r_s);
                e->giveBoundaryLocationArray(loc_c, bNodes, this->dofs, c_s);

                this->integrateTangent(Ke, e, boundary, i, pos);
                Ke.times(-scale);
                KeT.beTranspositionOf(Ke);
#ifdef _OPENMP
                if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
                answer.assemble(flux_loc_r, loc_c, Ke); // Contribution to delta_s_i equations
                answer.assemble(loc_r, flux_loc_c, KeT); // Contributions to delta_v equations
#ifdef _OPENMP
                if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
            }
        }
    } else   {
        printf("Skipping assembly in TransportGradientNeumann::assemble().\n");
    }
}

void TransportGradientNeumann :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                       const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc_r, loc_c, flux_loc_r, flux_loc_c;

    // Fetch the columns/rows for the flux contributions;
    mpFluxHom->giveLocationArray(mFluxIds, flux_loc_r, r_s);
    mpFluxHom->giveLocationArray(mFluxIds, flux_loc_c, c_s);

    rows.clear();
    cols.clear();
    
    int i = 0;    
    for ( auto &surfSet : this->surfSets ) {
        Set *set = this->giveDomain()->giveSet(surfSet);
        const IntArray &boundaries = set->giveBoundaryList();

        rows.resize( rows.size() + boundaries.giveSize() );
        cols.resize( cols.size() + boundaries.giveSize() );
        for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries[pos * 2] );
            int boundary = boundaries[pos * 2 + 1];

            const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary);
            e->giveBoundaryLocationArray(loc_r, bNodes, this->dofs, r_s);
            e->giveBoundaryLocationArray(loc_c, bNodes, this->dofs, c_s);

            // For most uses, loc_r == loc_c, and flux_loc_r == flux_loc_c.
            rows [ i ] = loc_r;
            cols [ i ] = flux_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows [ i ] = flux_loc_r;
            cols [ i ] = loc_c;
            i++;
        }
    }
}

void TransportGradientNeumann :: computeField(FloatArray &flux, TimeStep *tStep)
{
    mpFluxHom->giveUnknownVector(flux, mFluxIds, VM_Total, tStep);
}


void TransportGradientNeumann :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM > solver( 
        classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    double rve_size = this->domainSize();

    // Solving the system:
    // [Kuu Kus] [u] = [0]  => Kuu*u + Kus*s = 0
    // [Ksu 0  ] [s] = [I]  => Ksu*u         = I
    // 1. u = -Kuu^(-1)*Kus*s. Denote us = -Kuu^(-1)*Kus
    // 2. -[Ksu*us]*s = I. Denote Ks = Ksu*us
    // 3. s = Ks^(-1)*I
    // Here, s = tangent as rhs is identity.
    
    // 1.
    // This is not very good. We have to keep Kuu and Kff in memory at the same time. Not optimal
    // Consider changing this approach.
    EModelDefaultEquationNumbering fnum;
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx(stype) );
    if ( !Kff ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type  %d\n", stype);
    }
    Kff->buildInternalStructure(rve, this->domain->giveNumber(), fnum);

    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    
    IntArray loc_u, loc_s;
    this->mpFluxHom->giveLocationArray(this->mFluxIds, loc_s, fnum);
    int neq = Kff->giveNumberOfRows();
    loc_u.resize(neq - loc_s.giveSize());
    int k = 0;
    for ( int i = 1; i <= neq; ++i ) {
        if ( !loc_s.contains(i) ) {
            loc_u.at(++k) = i;
        }
    }

    std :: unique_ptr< SparseMtrx > Kuu = Kff->giveSubMatrix(loc_u, loc_u);
    FloatMatrix KusD;
#if 0
    // NOTE: Kus is actually a dense matrix, but we have to make it a dense matrix first
    std :: unique_ptr< SparseMtrx > Kus = Kff->giveSubMatrix(loc_u, loc_s);
    Kus->toFloatMatrix(KusD);
    Kus.reset();
#else
    // Ugly code-duplication, but worth it for performance.
    // PETSc doesn't allow you to take non-symmetric submatrices from a SBAIJ-matrix.
    // Which means extracting "Kus" puts limitations on Kuu.
    FloatMatrix Ke, KeT;
    IntArray loc_r;

    KusD.resize(neq - loc_s.giveSize(), loc_s.giveSize());
    for ( int i = 0; i < this->surfSets.giveSize(); ++i ) {
        int surfSet = this->surfSets[i];
        Set *set = this->giveDomain()->giveSet(surfSet);
        const IntArray &boundaries = set->giveBoundaryList();
        for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries[pos * 2] );
            int boundary = boundaries[pos * 2 + 1];

            const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary);
            e->giveBoundaryLocationArray(loc_r, bNodes, this->dofs, fnum);

            this->integrateTangent(Ke, e, boundary, i, pos);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            KusD.assemble(KeT, loc_r, {1, 2, 3}); // Contributions to delta_v equations
        }
    }
#endif
    // Release a large chunk of redundant memory early.
    Kff.reset();

    // 1.
    FloatMatrix us;
#if 1
    // Initial guess can help significantly for KSP solver,
    // However, it is difficult to construct a cheap, reliable estimate for Neumann b.c.
    // We need the tangent to relate q-bar to g-bar. Since this would be meaningless 
    // we shall instead assume isotropic response, and 
    int nsd = domain->giveNumberOfSpatialDimensions();
#if 0
    FloatMatrix Di(nsd, nsd), tmpD, tmpDi;
    double vol = 0;
    for ( auto &e : domain->giveElements() ) {
        static_cast< TransportElement* >(e.get())->computeConstitutiveMatrixAt(tmpD, Capacity, 
            e->giveDefaultIntegrationRulePtr()->getIntegrationPoint(1), tStep);
        double tmpVol = e->computeVolumeAreaOrLength();
        vol += tmpVol;
        tmpDi.beInverseOf(tmpD);
        Di.add(tmpVol, tmpDi);
    }
    Di.times(1/vol);
#endif
    us.resize(KusD.giveNumberOfRows(), KusD.giveNumberOfColumns());
    for ( auto &n : domain->giveDofManagers() ) {
        int k1 = n->giveDofWithID( this->dofs[0] )->__giveEquationNumber();
        if ( k1 ) {
            const auto &coords = n->giveCoordinates();
            for ( int i = 1; i <= nsd; ++i ) {
                us.at(k1, i) += -(coords.at(i) - mCenterCoord.at(i));
            }
        }
    }
    // Now "us" will have roughly the correct shape, but the wrong magnitude, so we rescale it:
    FloatMatrix KusD0;
    Kuu->times(us, KusD0);
    us.times( KusD.computeFrobeniusNorm() / KusD0.computeFrobeniusNorm() );
#endif

    if ( solver->solve(*Kuu, KusD, us) & NM_NoSuccess ) {
        OOFEM_ERROR("Failed to solve Kuu");
    }
    us.negated();

    // 2.
    FloatMatrix Ks;
    Ks.beTProductOf(KusD, us);

    // 3.
    tangent.beInverseOf(Ks);
    tangent.times(rve_size);
    tangent.negated();
}


void TransportGradientNeumann :: giveFluxLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s)
{
    mpFluxHom->giveLocationArray(mFluxIds, oCols, r_s);
}


void TransportGradientNeumann :: computeEta()
{
    OOFEM_LOG_INFO("Computing eta for Neumann b.c.\n");
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    eta.resize(this->surfSets.giveSize());
    for ( int i = 0; i < this->surfSets.giveSize(); ++i ) {
        // Compute the coordinate indices based on what surface we're on.
        int i_r;
        int i_t;
        if ( i % 3 == 0 ) { // Plane facing +- x
            i_r = 1;
            i_t = 2;
        } else if ( i % 3 == 1 ) { // Plane facing +- y
            i_r = 0;
            i_t = 2;
        } else { // Plane facing +- z
            i_r = 0;
            i_t = 1;
        }
        FloatArray coords, normal, tmp;
        FloatMatrix d;
        Set *setPointer = this->giveDomain()->giveSet(surfSets[i]);
        const IntArray &boundaries = setPointer->giveBoundaryList();

        eta[i].resize(boundaries.giveSize() / 2);

        // Set up and solve: c * eta_c = f
        // which represents the equations:
        // int (eta - 1) dA = 0
        // int (eta - 1) r dA = 0
        // int (eta - 1) t dA = 0
        //
        // int eta_0 + eta_r * r + eta_t * t dA = A
        // int (eta_0 + eta_r * r + eta_t * t) r dA = int r dA
        // int (eta_0 + eta_r * r + eta_t * t) t dA = int t dA
        FloatArray f(3), eta_c;
        FloatMatrix c(3, 3);
        for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries[pos * 2] );
            int boundary = boundaries[pos * 2 + 1];

            FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryIntegrationRule(order, boundary) );
            static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(d, Capacity, e->giveDefaultIntegrationRulePtr()->getIntegrationPoint(1), tStep);

            for ( auto &gp: *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);

                double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
                double dA = detJ * gp->giveWeight();
                interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
                coords.subtract(this->mCenterCoord);
                double r = coords[i_r], t = coords[i_t];
                
                // Compute material property
                ///@todo Can't do this yet, problem with material model interface (doesn't know the GP material mode). This should be changed.
                //static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(d, Capacity, gp, tStep);
                tmp.beProductOf(d, normal);
                double k = tmp.dotProduct(normal);

                f[0] += dA;
                f[1] += dA * r;
                f[2] += dA * t;

                c(0, 0) += dA * k;
                c(0, 1) += dA * k * r;
                c(0, 2) += dA * k * t;
                c(1, 0) += dA * k * r;
                c(1, 1) += dA * k * r * r;
                c(1, 2) += dA * k * t * r;
                c(2, 0) += dA * k * t;
                c(2, 1) += dA * k * r * t;
                c(2, 2) += dA * k * t * t;
            }
        }
        c.solveForRhs(f, eta_c);

        for ( int pos = 0; pos < boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries[pos * 2] );
            int boundary = boundaries[pos * 2 + 1];

            FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryIntegrationRule(order, boundary) );
            static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(d, Capacity, e->giveDefaultIntegrationRulePtr()->getIntegrationPoint(1), tStep);

            eta[i][pos].resize(ir->giveNumberOfIntegrationPoints());
            for ( auto &gp: *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);

                interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
                interp->boundaryLocal2Global(coords, boundary, lcoords, cellgeo);
                coords.subtract(this->mCenterCoord);
                double r = coords[i_r], t = coords[i_t];
                
                // Compute material property
                //static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(d, Capacity, gp, tStep);
                tmp.beProductOf(d, normal);
                double k = tmp.dotProduct(normal);

                eta[i][pos].at(gp->giveNumber()) = (eta_c[0] + eta_c[1] * r + eta_c[2] * t) * k;
            }
        }
    }
}


void TransportGradientNeumann :: integrateTangent(FloatMatrix &oTangent, Element *e, int boundary, int surfSet, int pos)
{
    FloatArray normal, n;

    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation
    //FEInterpolation *interpUnknown = e->giveInterpolation(this->dofs[0]);

    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryIntegrationRule(order, boundary) );

    oTangent.clear();

    for ( auto &gp: *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);
        // If eta isn't set, assume that classical Neumann b.c. is used (eta = 1)
        double scale = 1.0;
        if ( !eta.empty() ) {
            //printf(" surfset %d / %d, pos %d, number %d\n", surfSet, pos, gp->giveNumber(), eta.size());
            scale = eta[surfSet][pos].at(gp->giveNumber());
        }

        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        interp->boundaryEvalN(n, boundary, lcoords, cellgeo);
        oTangent.plusDyadUnsym(normal, n, scale * detJ * gp->giveWeight());
    }
}
} /* namespace oofem */
