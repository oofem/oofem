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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "prescribeddispslipbcneumannrc.h"
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
#include "crosssection.h"
#include "fei2dquadlin.h"


namespace oofem {
REGISTER_BoundaryCondition(PrescribedDispSlipBCNeumannRC);

PrescribedDispSlipBCNeumannRC :: PrescribedDispSlipBCNeumannRC(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedDispSlipHomogenization(),
    mpSigmaHom( new Node(0, d) )
{
    int nsd = d->giveNumberOfSpatialDimensions();
    for ( int i = 0; i < nsd * nsd; i++ ) {
        // Just putting in X_i id-items since they don't matter.
        int dofId = d->giveNextFreeDofID();
        mSigmaIds.followedBy(dofId);
        mpSigmaHom->appendDof( new MasterDof( mpSigmaHom.get(), ( DofIDItem ) ( dofId ) ) );
    }
}

PrescribedDispSlipBCNeumannRC :: ~PrescribedDispSlipBCNeumannRC()
{
}


void PrescribedDispSlipBCNeumannRC :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    PrescribedDispSlipHomogenization :: initializeFrom(ir);
}


void PrescribedDispSlipBCNeumannRC :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedDispSlipHomogenization :: giveInputRecord(input);
}


DofManager *PrescribedDispSlipBCNeumannRC :: giveInternalDofManager(int i)
{
    return mpSigmaHom.get();
}

void PrescribedDispSlipBCNeumannRC :: scale(double s)
{
    this->dispGradient.times(s);
}

void PrescribedDispSlipBCNeumannRC :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                   CharType type, ValueModeType mode,
                                                   const UnknownNumberingScheme &s, 
                                                   FloatArray *eNorm, 
                                                   void* lock)
{
    IntArray sigma_loc;  // For the displacements and stress respectively
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for sigma, the load is simply the prescribed gradient.
        double rve_size = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());
        FloatArray stressLoad;
        FloatArray gradVoigt;
        giveDispGradient(gradVoigt);

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        stressLoad.beScaled(-rve_size*loadLevel, gradVoigt);

        answer.assemble(stressLoad, sigma_loc);

    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ke;
        FloatArray fe_v, fe_s;
        FloatArray sigmaHom, e_u;
        IntArray loc, masterDofIDs, sigmaMasterDofIDs;

        // Fetch the current values of the stress;
        mpSigmaHom->giveUnknownVector(sigmaHom, mSigmaIds, mode, tStep);
        // and the master dof ids for sigmadev used for the internal norms
        mpSigmaHom->giveMasterDofIDArray(mSigmaIds, sigmaMasterDofIDs);

        // Assemble
        Set *setPointer = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = setPointer->giveBoundaryList();
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            e->giveLocationArray(loc, s, & masterDofIDs);
            e->computeVectorOf(mode, tStep, e_u);
            this->integrateTangent(Ke, e, boundary);

            // We just use the tangent, less duplicated code (the addition of sigmaDev is linear).
            fe_v.beProductOf(Ke, e_u);
            fe_s.beTProductOf(Ke, sigmaHom);

            // temporary hack for inconsistent normal vectors
            auto interp = dynamic_cast< FEI2dQuadLin * >( e->giveInterpolation());
            if ( interp != nullptr ) {
            } else {
                fe_v.negated();
                fe_s.negated();
            }

            answer.assemble(fe_s, loc); // Contributions to delta_v equations
            answer.assemble(fe_v, sigma_loc); // Contribution to delta_s_i equations
            if ( eNorm != NULL ) {
                eNorm->assembleSquared(fe_s, masterDofIDs);
                eNorm->assembleSquared(fe_v, sigmaMasterDofIDs);
            }
        }
    }
}

void PrescribedDispSlipBCNeumannRC :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                             CharType type, 
                                             const UnknownNumberingScheme &r_s, 
                                             const UnknownNumberingScheme &c_s, 
                                             double scale,
                                             void* lock)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;
        Set *set = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = set->giveBoundaryList();

        // Fetch the columns/rows for the stress contributions;
        mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_r, r_s);
        mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_c, c_s);

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            e->giveLocationArray(loc_r, r_s);
            e->giveLocationArray(loc_c, c_s);

            this->integrateTangent(Ke, e, boundary);

            //temporary hack for inconsistent normal vectors
            auto interp = dynamic_cast< FEI2dQuadLin * >( e->giveInterpolation() );
            if ( interp != nullptr ) {
            } else {
                Ke.negated();
            }

            Ke.times(scale);
            KeT.beTranspositionOf(Ke);

            answer.assemble(sigma_loc_r, loc_c, Ke); // Contribution to delta_s_i equations
            answer.assemble(loc_r, sigma_loc_c, KeT); // Contributions to delta_v equations
        }
    } else   {
        OOFEM_LOG_DEBUG("Skipping assembly in PrescribedDispSlipBCNeumann::assemble().");
    }
}

void PrescribedDispSlipBCNeumannRC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                       const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;

    // Fetch the columns/rows for the stress contributions;
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_r, r_s);
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_c, c_s);

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    rows.resize( boundaries.giveSize() );
    cols.resize( boundaries.giveSize() );
    int i = 0;
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );

        // Here, we could use only the nodes actually located on the boundary, but we don't.
        // Instead, we use all nodes belonging to the element, which is allowed because the
        // basis functions related to the interior nodes will be zero on the boundary.
        // Obviously, this is less efficient, so why do we want to do it this way?
        // Because it is easier when XFEM enrichments are present. /ES
        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        // For most uses, loc_r == loc_c, and sigma_loc_r == sigma_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = sigma_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = sigma_loc_r;
        cols [ i ] = loc_c;
        i++;
    }
}

void PrescribedDispSlipBCNeumannRC :: computeStress(FloatArray &sigma, TimeStep *tStep)
{
    double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
    double areaRVE = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());
    double thick = volRVE/areaRVE;
    mpSigmaHom->giveUnknownVector(sigma, mSigmaIds, VM_Total, tStep);
    sigma.times(1/thick);
}


void PrescribedDispSlipBCNeumannRC :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    std :: unique_ptr< SparseLinearSystemNM > solver( 
        classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    double rve_size = this->domainSize(this->giveDomain(), this->giveSetNumber());

    // 1. Kuu*us = -Kus*s   =>  us = -Kuu\Ku  where u = us*s
    // 2. Ks = Kus'*us
    // 3. Ks*lambda = I
    
    // 1.
    // This is not very good. We have to keep Kuu and Kff in memory at the same time. Not optimal
    // Consider changing this approach.
    EModelDefaultEquationNumbering fnum;
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx(stype) );
    if ( !Kff ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure(rve, this->domain->giveNumber(), fnum);

    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    
    IntArray loc_u, loc_s;
    this->mpSigmaHom->giveLocationArray(this->mSigmaIds, loc_s, fnum);
    int neq = Kff->giveNumberOfRows();
    loc_u.resize(neq - loc_s.giveSize());
    int k = 0;
    for ( int i = 1; i <= neq; ++i ) {
        if ( !loc_s.contains(i) ) {
            loc_u.at(k++) = i;
        }
    }

    std :: unique_ptr< SparseMtrx > Kuu = Kff->giveSubMatrix(loc_u, loc_u);
    // NOTE: Kus is actually a dense matrix, but we have to make it a dense matrix first
    std :: unique_ptr< SparseMtrx > Kus = Kff->giveSubMatrix(loc_u, loc_s);
    FloatMatrix eye(Kus->giveNumberOfColumns(), Kus->giveNumberOfColumns());
    eye.beUnitMatrix();
    FloatMatrix KusD;
    Kus->times(KusD, eye);

    // Release a large chunk of redundant memory early.
    Kus.reset();
    Kff.reset();

    // 1.
    FloatMatrix us;
    solver->solve(*Kuu, KusD, us);
    us.negated();

    // 2.
    FloatMatrix Ks;
    Ks.beTProductOf(KusD, us);
    Kus->times(Ks, us);

    // 3.
    tangent.beInverseOf(Ks);
    tangent.times(rve_size);
}


void PrescribedDispSlipBCNeumannRC :: giveStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s)
{
    mpSigmaHom->giveLocationArray(mSigmaIds, oCols, r_s);
}

void PrescribedDispSlipBCNeumannRC :: integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex)
{
    FloatArray normal, n;
    FloatMatrix nMatrix, E_n;
    FloatMatrix contrib;

    Domain *domain = e->giveDomain();

    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();

    // Interpolation order
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;

    XfemElementInterface *xfemElInt = dynamic_cast< XfemElementInterface * >( e );
    if ( xfemElInt && domain->hasXfemManager() ) {
        FEInterpolation2d *interp2d = dynamic_cast< FEInterpolation2d * >( interp );
        if ( !interp2d ) {
            OOFEM_ERROR("failed to cast to FEInterpolation2d.")
        }
        const auto &edgeNodes = interp2d->computeLocalEdgeMapping(iBndIndex);

//        const auto &xS = * ( e->giveDofManager( edgeNodes.at(1) )->giveCoordinates() );
//        const auto &xE = * ( e->giveDofManager( edgeNodes.at( edgeNodes.giveSize() ) )->giveCoordinates() );

        std :: vector< Line >segments;
        std :: vector< FloatArray >intersecPoints;
        xfemElInt->partitionEdgeSegment(iBndIndex, segments, intersecPoints);
        MaterialMode matMode = e->giveMaterialMode();
        ir = std::make_unique<DiscontinuousSegmentIntegrationRule>(1, e, segments);
        int numPointsPerSeg = 1;
        ir->SetUpPointsOnLine(numPointsPerSeg, matMode);
    } else {
        ir = interp->giveBoundaryIntegrationRule(order, iBndIndex);
    }

    oTangent.clear();

    for ( auto &gp: *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        // Evaluate the normal;
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        // Compute global coordinates of Gauss point
        FloatArray globalCoord;

        interp->boundaryLocal2Global(globalCoord, iBndIndex, lcoords, cellgeo);

        // Compute local coordinates on the element
        FloatArray bulkElLocCoords;
        e->computeLocalCoordinates(bulkElLocCoords, globalCoord);

        // If cracks cross the edge, special treatment is necessary.
        // Exploit the XfemElementInterface to minimize duplication of code.
        if ( xfemElInt != NULL && domain->hasXfemManager() ) {

            xfemElInt->XfemElementInterface_createEnrNmatrixAt(nMatrix, bulkElLocCoords, * e, false);
        } else {
            // Evaluate the velocity/displacement coefficients
            interp->evalN(n, bulkElLocCoords, cellgeo);
            nMatrix.beNMatrixOf(n, nsd);
        }

        E_n.resize(nsd*nsd, nsd);
        E_n.zero();

        if ( nsd == 3 ) {
            E_n.at(1, 1) = normal.at(1);
            E_n.at(2, 2) = normal.at(2);
            E_n.at(3, 3) = normal.at(3);
            E_n.at(4, 1) = normal.at(2);
            E_n.at(5, 1) = normal.at(3);
            E_n.at(6, 2) = normal.at(3);
            E_n.at(7, 2) = normal.at(1);
            E_n.at(8, 3) = normal.at(1);
            E_n.at(9, 3) = normal.at(2);
        } else if ( nsd == 2 ) {
            E_n.at(1, 1) = normal.at(1);
            E_n.at(2, 2) = normal.at(2);
            E_n.at(3, 1) = normal.at(2);
            E_n.at(4, 2) = normal.at(1);
        } else {
            E_n.at(1, 1) = normal.at(1);
        }

        contrib.beProductOf(E_n, nMatrix);

        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}

double PrescribedDispSlipBCNeumannRC::domainSize( Domain *d, int set )
{
    double omegaBox = PrescribedDispSlipHomogenization::domainSize(d, this->giveSetNumber());

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant in 2D
        Element *e = this->giveDomain()->giveElement( this->giveDomain()->giveSet( this->giveSetNumber() )->giveBoundaryList().at(1) );
        std::unique_ptr<IntegrationRule> ir = e->giveInterpolation()->giveIntegrationRule( e->giveInterpolation()->giveInterpolationOrder() );
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(0);
        double thickness = cs->give(CS_Thickness, gp);
        return omegaBox * thickness;
    } else {
        return omegaBox;
    }
}

} /* namespace oofem */
