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

#include "transportgradientdirichlet.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "function.h"
#include "engngm.h"
#include "set.h"
#include "node.h"
#include "element.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "assemblercallback.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "feinterpol3d.h"
#include "transportelement.h"
#include "gausspoint.h"
#include "sparselinsystemnm.h"

namespace oofem {
REGISTER_BoundaryCondition(TransportGradientDirichlet);


IRResultType TransportGradientDirichlet :: initializeFrom(InputRecord *ir)
{
    IRResultType result;           // Required by IR_GIVE_FIELD macro
    IR_GIVE_FIELD(ir, mGradient, _IFT_TransportGradientDirichlet_gradient);

    mCenterCoord.resize(3);
    IR_GIVE_OPTIONAL_FIELD(ir, mCenterCoord, _IFT_TransportGradientDirichlet_centerCoords);

    if ( ir->hasField(_IFT_TransportGradientDirichlet_usePsi) ) {
        IR_GIVE_FIELD(ir, surfSets, _IFT_TransportGradientDirichlet_surfSets);
        IR_GIVE_FIELD(ir, edgeSets, _IFT_TransportGradientDirichlet_edgeSets);
    }

    return GeneralBoundaryCondition :: initializeFrom(ir);
}


void TransportGradientDirichlet :: giveInputRecord(DynamicInputRecord &input)
{
    input.setField(mGradient, _IFT_TransportGradientDirichlet_gradient);
    input.setField(mCenterCoord, _IFT_TransportGradientDirichlet_centerCoords);
    if ( this->surfSets.giveSize() > 0 ) {
        input.setField(_IFT_TransportGradientDirichlet_usePsi);
        input.setField(surfSets, _IFT_TransportGradientDirichlet_surfSets);
        input.setField(edgeSets, _IFT_TransportGradientDirichlet_edgeSets);
    }

    return GeneralBoundaryCondition :: giveInputRecord(input);
}


double TransportGradientDirichlet :: give(Dof *dof, ValueModeType mode, double time)
{
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    if ( coords->giveSize() != this->mCenterCoord.giveSize() ) {
        OOFEM_ERROR("Size of coordinate system different from center coordinate in b.c.");
    }

    double factor = 0;
    if ( mode == VM_Total ) {
        factor = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        factor = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        factor = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }

    // Reminder: u_i = g_i . (x_i - xc_i + psi_i)
    FloatArray dx;
    dx.beDifferenceOf(* coords, this->mCenterCoord);
    // Add "psi" if it is defined. Classical Dirichlet b.c. is retained if this isn't defined (or set to zero).
    if ( !psis.empty() ) {
        dx.add(psis[dof->giveDofManager()->giveNumber()]);
    }

    return mGradient.dotProduct(dx) * factor;
}


double TransportGradientDirichlet :: domainSize()
{
    Domain *domain = this->giveDomain();
    int nsd = domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = domain->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = domain->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
}


void TransportGradientDirichlet :: updateCoefficientMatrix(FloatMatrix &C)
// v_prescribed = C.g = (x-xbar + psi).g;
// C = [x-psi_x y-psi_y]
//     [ .. ] in 2D, voigt form [g_1, g_2]
// C = [x-psi_x y-psi_y z-psi_z]
//     [ ....] in 3D, voigt form [g_1, g_2, g_3]
{
    Domain *domain = this->giveDomain();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    C.resize(npeq, nsd);
    C.zero();

    for ( auto &n : domain->giveDofManagers() ) {
        // Add "psi" if it is defined. Classical Dirichlet b.c. is retained if this isn't defined (or set to zero).
        FloatArray psi(nsd);
        if ( !psis.empty() ) {
            psi = psis[n->giveNumber()];
        }
        FloatArray *coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID( this->dofs(0) );
        int k1 = d1->__givePrescribedEquationNumber();
        if ( k1 ) {
            for ( int i = 1; i <= nsd; ++i ) {
                C.at(k1, i) = coords->at(i) - mCenterCoord.at(i) + psi.at(i);
            }
        }
    }
}


void TransportGradientDirichlet :: computeField(FloatArray &sigma, TimeStep *tStep)
{
    EngngModel *emodel = this->domain->giveEngngModel();
    int npeq = emodel->giveNumberOfDomainEquations( this->giveDomain()->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    FloatArray R_c(npeq), R_ext(npeq);

    R_c.zero();
    R_ext.zero();
    emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    R_c.subtract(R_ext);

    // Condense it;
    FloatMatrix C;
    this->updateCoefficientMatrix(C);
    sigma.beTProductOf(C, R_c);
    sigma.times( 1. / this->domainSize() );
}


void TransportGradientDirichlet :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
// a = [a_c; a_f];
// K.a = [R_c,0];
// [K_cc, K_cf; K_fc, K_ff].[a_c;a_f] = [R_c; 0];
// a_c = d.[x-x_b] = [x-x_b].d = C.d
// E = C'.(K_cc - K_cf.K_ff^(-1).K_fc).C
//   = C'.(K_cc.C - K_cf.(K_ff^(-1).(K_fc.C)))
//   = C'.(K_cc.C - K_cf.a)
//   = C'.X
{
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM >solver( classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    EModelDefaultEquationNumbering fnum;
    EModelDefaultPrescribedEquationNumbering pnum;

    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    std :: unique_ptr< SparseMtrx >Kff( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx >Kfp( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx >Kpf( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx >Kpp( classFactory.createSparseMtrx(stype) );
    if ( !Kff ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure(rve, 1, fnum);
    Kfp->buildInternalStructure(rve, 1, fnum, pnum);
    Kpf->buildInternalStructure(rve, 1, pnum, fnum);
    Kpp->buildInternalStructure(rve, 1, pnum);
    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    rve->assemble(*Kfp, tStep, TangentAssembler(TangentStiffness), fnum, pnum, this->domain);
    rve->assemble(*Kpf, tStep, TangentAssembler(TangentStiffness), pnum, fnum, this->domain);
    rve->assemble(*Kpp, tStep, TangentAssembler(TangentStiffness), pnum, this->domain);

    FloatMatrix C, X, Kpfa, KfpC, a;

    this->updateCoefficientMatrix(C);
    Kpf->timesT(C, KfpC);
    solver->solve(*Kff, KfpC, a);
    Kpp->times(C, X);
    Kpf->times(a, Kpfa);
    X.subtract(Kpfa);
    tangent.beTProductOf(C, X);
    tangent.times( 1. / this->domainSize() );
}

void TransportGradientDirichlet :: computePsi()
{
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    FloatMatrix D, Dred, B, Bred, DB, Ke;
    FloatArray f, q, fx, qx, x;
    FloatArray n, b, fe, Ne, cvec;
    IntArray loc;
    // As defined for the hex element:
    IntArray surfaceOrder{3, 3, 1, 2, 1, 2};
    IntArray edgeOrder{2, 1, 2, 1, 3, 3, 3, 3, 2, 1, 2, 1};

    std :: unique_ptr< SparseLinearSystemNM > solver( classFactory.createSparseLinSolver(ST_Petsc, this->giveDomain(), this->giveDomain()->giveEngngModel()) );

    // One "psi" per node on the boundary. Initially all to zero.
    this->psis.clear();
    const IntArray &totalNodes = domain->giveSet(this->giveSetNumber())->giveNodeList();
    for ( int node : totalNodes ) {
        this->psis.emplace(node, FloatArray(3));
    }

    // Identify corner nodes and all the total edge nodes (needed for equation numbering)
    IntArray totalCornerNodes;
    IntArray totalEdgeNodes;
    for ( int i = 0; i < this->edgeSets.giveSize(); ++i ) {
        Set *setPointer = this->giveDomain()->giveSet(edgeSets[i]);
        const IntArray &nodes = setPointer->giveNodeList();
        for ( int n : nodes ) {
            if ( !totalEdgeNodes.insertSortedOnce(n, 10) ) {
                totalCornerNodes.insertSortedOnce(n);
            }
        }
    }

    // First we must determine the values along the edges, which become boundary conditions for the surfaces
    for ( int i = 0; i < this->edgeSets.giveSize(); ++i ) {
        Set *setPointer = this->giveDomain()->giveSet(edgeSets[i]);
        const IntArray &edges = setPointer->giveEdgeList();
        int t_index = edgeOrder[i];
        
        // Number the equations along this edge set
        const IntArray &edgeNodes = setPointer->giveNodeList();
        IntArray eqs(edgeNodes.giveSize());
        int eq_count = 0;
        for ( int n : edgeNodes ) {
            if ( totalCornerNodes.containsSorted(n) ) {
                eqs[n] = eq_count++;
            }
        }

        std :: unique_ptr< SparseMtrx > K( classFactory.createSparseMtrx(solver->giveRecommendedMatrix(true)) );
        ///@todo Preallocation
        for ( int pos = 0; pos < edges.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( edges[pos * 2] );
            int edge = edges[pos * 2 + 1];
            
            FEInterpolation3d *interp = static_cast< FEInterpolation3d* >( e->giveInterpolation() );
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryEdgeIntegrationRule(order, edge) );

            // Need the coordinates arranged as [x1, y1, z1, x2, y2, z2, ...] for the RHS.
            IntArray bNodes;
            interp->boundaryEdgeGiveNodes(bNodes, edge);
            cvec.resize(bNodes.giveSize());
            int count = 0;
            for ( auto &node : bNodes ) {
                Node *n = e->giveNode(node);
                cvec[count++] = n->giveCoordinate(t_index);
                loc.followedBy(eqs[n->giveNumber()]);
            }

            for ( auto &gp: *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);

                interp->boundaryEdgeEvalN(n, edge, lcoords, cellgeo);

                double detJ = interp->boundaryEdgeGiveTransformationJacobian(edge, lcoords, cellgeo);
                interp->edgeEvaldNdxi(b, edge, lcoords, cellgeo);
                b.times(1. / detJ);
                double dL = detJ * gp->giveWeight();
                
                // Compute material property
                static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(D, Capacity, gp, tStep);
                double k = D.at(t_index, t_index);

                Ke.plusDyadSymmUpper(b, - k * dL);
                Ne.add(dL, n);
            }
            Ke.symmetrized();
            fe.beProductOf(Ke, cvec);

            K->assemble(loc, Ke);
            f.assemble(fe, loc);
            q.assemble(Ne, loc);
        }

        solver->solve(*K, q, qx);
        solver->solve(*K, f, fx);
        double lambda = (q.dotProduct(fx)) / (q.dotProduct(qx));
        x = fx;
        x.add(-lambda, q);

        for ( int n : edgeNodes ) {
            if ( eqs[n] > 0 ) {
                this->psis[n].at(t_index) = x.at(eqs[n]);
            }
        }
    }
    
    // Surfaces use the edge solutions are boundary conditions:
    for ( int i = 0; i < this->surfSets.giveSize(); ++i ) {
        Set *setPointer = this->giveDomain()->giveSet(surfSets[i]);
        const IntArray &surfs = setPointer->giveBoundaryList();
        int t_index = surfaceOrder[i];
        
        // Number the equations along this surface set
        const IntArray &surfNodes = setPointer->giveNodeList();
        IntArray eqs(surfNodes.giveSize());
        int eq_count = 0;
        for ( int n : surfNodes ) {
            if ( totalEdgeNodes.containsSorted(n) ) {
                eqs[n] = eq_count++;
            }
        }

        std :: unique_ptr< SparseMtrx > K( classFactory.createSparseMtrx(solver->giveRecommendedMatrix(true)) );
        ///@todo Preallocation
        for ( int pos = 0; pos < surfs.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( surfs[pos * 2] );
            int surf = surfs[pos * 2 + 1];

            FEInterpolation3d *interp = static_cast< FEInterpolation3d* >( e->giveInterpolation() );
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryIntegrationRule(order, surf) );

            // Need [x1 + psi1, x2 + psi2, x3 + psi3, ...] for the RHS.
            IntArray bNodes;
            interp->boundaryGiveNodes(bNodes, surf);
            cvec.resize(bNodes.giveSize());
            int count = 0;
            for ( auto &node : bNodes ) {
                Node  *n = e->giveNode(node);
                cvec[count++] = n->giveCoordinate(t_index) + this->psis[node].at(t_index);
                loc.followedBy(eqs[n->giveNumber()]);
            }

            for ( auto &gp: *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);

                interp->boundaryEvalN(n, surf, lcoords, cellgeo);

                double detJ = interp->boundaryGiveTransformationJacobian(surf, lcoords, cellgeo);
                interp->surfaceEvaldNdx(B, surf, lcoords, cellgeo);
                double dA = detJ * gp->giveWeight();
                
                // Compute material property
                static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(D, Capacity, gp, tStep);

                IntArray subindx, allindx;
                if ( t_index == 1 ) {
                    subindx = {2, 3};
                } else if ( t_index == 2 ) {
                    subindx = {1, 3};
                } else {
                    subindx = {1, 2};
                }
                allindx.enumerate(Bred.giveNumberOfColumns());
                Bred.beSubMatrixOf(B, subindx, allindx);
                Dred.beSubMatrixOf(D, subindx, subindx);


                DB.beProductOf(Dred, Bred);
                Ke.plusProductSymmUpper(Bred, DB, - dA);
                Ne.add(dA, n);
            }
            Ke.symmetrized();
            fe.beProductOf(Ke, cvec);

            K->assemble(loc, Ke);
            f.assemble(fe, loc);
            q.assemble(Ne, loc);
        }
        solver->solve(*K, q, qx);
        solver->solve(*K, f, fx);
        double lambda = (q.dotProduct(fx)) / (q.dotProduct(qx));
        x = fx;
        x.add(-lambda, q);
        
        for ( int n : surfNodes ) {
            if ( eqs[n] > 0 ) {
                this->psis[n].at(t_index) = x.at(eqs[n]);
            }
        }
    }
}

} // end namespace oofem
