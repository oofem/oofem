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

#include "tm/BoundaryCondition/transportgradientdirichlet.h"
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
#include "tm/Elements/transportelement.h"
#include "gausspoint.h"
#include "sparselinsystemnm.h"

#include "timer.h" // Benchmarking

#include <tuple>
#include <vector>
#include <algorithm>
#include <iterator>

namespace oofem {
REGISTER_BoundaryCondition(TransportGradientDirichlet);


void TransportGradientDirichlet :: initializeFrom(InputRecord &ir)
{
    GeneralBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, mGradient, _IFT_TransportGradientDirichlet_gradient);

    mCenterCoord.resize(3);
    IR_GIVE_OPTIONAL_FIELD(ir, mCenterCoord, _IFT_TransportGradientDirichlet_centerCoords);

    this->tractionControl = ir.hasField(_IFT_TransportGradientDirichlet_tractionControl);
    if ( this->tractionControl ) {
        IR_GIVE_FIELD(ir, surfSets, _IFT_TransportGradientDirichlet_surfSets);
        //IR_GIVE_FIELD(ir, edgeSets, _IFT_TransportGradientDirichlet_edgeSets);
    }
}


void TransportGradientDirichlet :: giveInputRecord(DynamicInputRecord &input)
{
    input.setField(mGradient, _IFT_TransportGradientDirichlet_gradient);
    input.setField(mCenterCoord, _IFT_TransportGradientDirichlet_centerCoords);
    input.setField(surfSets, _IFT_TransportGradientDirichlet_surfSets);
    //input.setField(edgeSets, _IFT_TransportGradientDirichlet_edgeSets);
    if ( this->tractionControl ) {
        input.setField(_IFT_TransportGradientDirichlet_tractionControl);
    }

    return GeneralBoundaryCondition :: giveInputRecord(input);
}


void TransportGradientDirichlet :: postInitialize()
{
    BoundaryCondition :: postInitialize();
    
    if ( this->tractionControl ) this->computeXi();
}


double TransportGradientDirichlet :: give(Dof *dof, ValueModeType mode, double time)
{
    const auto &coords = dof->giveDofManager()->giveCoordinates();

    if ( coords.giveSize() != this->mCenterCoord.giveSize() ) {
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

    // Reminder: u_i = g_i . (x_i - xc_i + xi_i)
    auto dx = coords - this->mCenterCoord;
    // Add "xi" if it is defined. Classical Dirichlet b.c. is retained if this isn't defined (or set to zero).
    if ( !xis.empty() ) {
        dx.add(xis[dof->giveDofManager()->giveNumber()]);
    }

    return mGradient.dotProduct(dx) * factor;
}


double TransportGradientDirichlet :: domainSize()
{
    Domain *domain = this->giveDomain();
    int nsd = domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    if ( this->tractionControl ) {
        for ( auto &surf : this->surfSets ) {
            const IntArray &boundaries = domain->giveSet(surf)->giveBoundaryList();

            for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
                Element *e = domain->giveElement( boundaries.at(pos * 2 - 1) );
                int boundary = boundaries.at(pos * 2);
                FEInterpolation *fei = e->giveInterpolation();
                domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
            }
        }
    } else {
        const IntArray &boundaries = domain->giveSet(this->set)->giveBoundaryList();

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = domain->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);
            FEInterpolation *fei = e->giveInterpolation();
            domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
        }
    }
    return fabs(domain_size / nsd);
}


void TransportGradientDirichlet :: computeCoefficientMatrix(FloatMatrix &C)
// v_prescribed = C.g = (x-xbar + xi).g;
// C = [x-xi_x y-xi_y]
// C = [x-xi_x y-xi_y z-xi_z]
{
    Domain *domain = this->giveDomain();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    C.resize(npeq, nsd);
    C.zero();

    for ( auto &n : domain->giveDofManagers() ) {
        const auto &coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID( this->dofs[0] );
        int k1 = d1->__givePrescribedEquationNumber();
        if ( k1 ) {
            // Add "xi" if it is defined. Classical Dirichlet b.c. is retained if this isn't defined (or set to zero).
            FloatArray xi(nsd);
            if ( this->tractionControl ) {
                xi = xis[n->giveNumber()];
            }
            for ( int i = 1; i <= nsd; ++i ) {
                C.at(k1, i) = coords.at(i) - mCenterCoord.at(i) + xi.at(i);
            }
            //printf("C.at(%d, :) = %e, %e, %e\n", k1, C.at(k1, 1), C.at(k1, 2), C.at(k1, 3));
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
    this->computeCoefficientMatrix(C);
    sigma.beTProductOf(C, R_c);
    sigma.times( -1. / this->domainSize() );
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
    //std :: unique_ptr< SparseMtrx >Kfp( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx >Kpf( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx >Kpp( classFactory.createSparseMtrx(stype) );
    if ( !Kff ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure(rve, 1, fnum);
    //Kfp->buildInternalStructure(rve, 1, fnum, pnum);
    Kpf->buildInternalStructure(rve, 1, pnum, fnum);
    Kpp->buildInternalStructure(rve, 1, pnum);

    Timer t;
    t.startTimer();
#if 0
    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    //rve->assemble(*Kfp, tStep, TangentAssembler(TangentStiffness), fnum, pnum, this->domain);
    rve->assemble(*Kpf, tStep, TangentAssembler(TangentStiffness), pnum, fnum, this->domain);
    rve->assemble(*Kpp, tStep, TangentAssembler(TangentStiffness), pnum, this->domain);
#else
    auto ma = TangentAssembler(TangentStiffness);
    IntArray floc, ploc;
    FloatMatrix mat, R;

    int nelem = domain->giveNumberOfElements();
#ifdef _OPENMP
 #pragma omp parallel for shared(Kff, Kpf, Kpp) private(mat, R, floc, ploc)
#endif
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement(ielem);
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote || !element->isActivated(tStep) ) {
            continue;
        }

        ma.matrixFromElement(mat, *element, tStep);

        if ( mat.isNotEmpty() ) {
            ma.locationFromElement(floc, *element, fnum);
            ma.locationFromElement(ploc, *element, pnum);
            ///@todo This rotation matrix is not flexible enough.. it can only work with full size matrices and doesn't allow for flexibility in the matrixassembler.
            if ( element->giveRotationMatrix(R) ) {
                mat.rotatedWith(R);
            }

#ifdef _OPENMP
 #pragma omp critical
#endif
            {
                Kff->assemble(floc, mat);
                Kpf->assemble(ploc, floc, mat);
                Kpp->assemble(ploc, mat);
            }
        }
    }
    Kff->assembleBegin();
    Kpf->assembleBegin();
    Kpp->assembleBegin();

    Kff->assembleEnd();
    Kpf->assembleEnd();
    Kpp->assembleEnd();
#endif
    t.stopTimer();
    OOFEM_LOG_INFO("Assembly time %.3f s\n", t.getUtime());

    FloatMatrix C, X, Kpfa, KfpC, sol;

    this->computeCoefficientMatrix(C);
    Kpf->timesT(C, KfpC);
    
    // Initial guess (Taylor assumption) helps KSP-iterations
    sol.resize(KfpC.giveNumberOfRows(), KfpC.giveNumberOfColumns());
    int nsd = domain->giveNumberOfSpatialDimensions();
    for ( auto &n : domain->giveDofManagers() ) {
        int k1 = n->giveDofWithID( this->dofs[0] )->__giveEquationNumber();
        if ( k1 ) {
            const auto &coords = n->giveCoordinates();
            for ( int i = 1; i <= nsd; ++i ) {
                sol.at(k1, i) = -(coords.at(i) - mCenterCoord.at(i));
            }
        }
    }

    t.startTimer();
    if ( solver->solve(*Kff, KfpC, sol) != CR_CONVERGED ) {
        OOFEM_ERROR("Failed to solve Kff");
    }
    t.stopTimer();
    Kpp->times(C, X);
    Kpf->times(sol, Kpfa);
    X.subtract(Kpfa);
    tangent.beTProductOf(C, X);
    tangent.times( 1. / this->domainSize() );
    OOFEM_LOG_INFO("Tangent problem solve time %.3f s\n", t.getUtime());
}

void TransportGradientDirichlet :: computeXi()
{
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();

    std :: unique_ptr< SparseLinearSystemNM > solver( classFactory.createSparseLinSolver(ST_Petsc, this->giveDomain(), this->giveDomain()->giveEngngModel()) );

    // One "xi" per node on the boundary. Initially all to zero.
    this->xis.clear();
    for ( int node : domain->giveSet(this->giveSetNumber())->giveNodeList() ) {
        this->xis.emplace(node, FloatArray(3));
    }
    for ( auto &surf : this->surfSets ) {
        for ( int node : domain->giveSet(surf)->giveNodeList() ) {
            this->xis.emplace(node, FloatArray(3));
        }
    }

    OOFEM_LOG_INFO("Computing edge sets from surfaces\n");
    // Instead of requiring the input file to specify the edges, this code will automatically detect them
    std :: vector< std :: vector< int > > surf2edges = {{1, 2, 3, 4},
        {9, 10, 11, 12},
        {1, 5, 6, 9},
        {2, 6, 7, 10},
        {3, 7, 8, 11},
        {4, 5, 8, 12}};

    // Edge sets generated from surface sets:
    // 1-2, 1-3, 1-4, .. 2-3, 2-4, ..., 5-6 (skipping the empty sets)
    // In case we ever want to use more arbitrary RVE "cutouts", this could be used.
    std :: vector< Set > edgeSets;
    std :: vector< std :: vector< std :: tuple< int, int > > > surfedges( this->surfSets.giveSize() );
    for ( int i = 0; i < this->surfSets.giveSize(); ++i ) {
        const IntArray &surfs = this->giveDomain()->giveSet(surfSets[i])->giveBoundaryList();
        surfedges[i].reserve( 4 * surfs.giveSize() / 2 );
        for ( int pos = 0; pos < surfs.giveSize() / 2; ++pos ) {
            for ( int edgenum : surf2edges[surfs[pos * 2 + 1]-1] ) {
                surfedges[i].emplace_back( std :: make_tuple(surfs[pos * 2], edgenum) );
            }
        }
    }

    for ( int i = 0; i < this->surfSets.giveSize() - 1; ++i ) {
        for ( int j = i+1; j < this->surfSets.giveSize(); ++j ) {
            std :: vector< std :: tuple< int, int > > ijEdgeSet;
            std :: set_intersection(surfedges[i].begin(), surfedges[i].end(), 
                                    surfedges[j].begin(), surfedges[j].end(), std::back_inserter(ijEdgeSet));
            IntArray edgelist;
            edgelist.preallocate(ijEdgeSet.size() * 2);
            for ( auto &edge : ijEdgeSet ) {
                edgelist.followedBy( std :: get<0>(edge) );
                edgelist.followedBy( std :: get<1>(edge) );
            }

            if ( edgelist.giveSize() > 0 ) {
                Set s(0, this->giveDomain());
                s.setEdgeList(edgelist);
                edgeSets.emplace_back(std :: move(s));
            }
        }
    }
    // END OF EDGE-SET GENERATION

    // Identify corner nodes and all the total edge nodes (needed for equation numbering)
    IntArray totalCornerNodes;
    IntArray totalEdgeNodes;
    for ( auto &setPointer : edgeSets ) {
        const IntArray &nodes = setPointer.giveNodeList();
        for ( int n : nodes ) {
            if ( !totalEdgeNodes.insertSortedOnce(n, 10) ) {
                totalCornerNodes.insertSortedOnce(n);
            }
        }
    }

#if 1
    //IntArray edgeOrder{2, 1, 2, 1, 3, 3, 3, 3, 2, 1, 2, 1};
    // First we must determine the values along the edges, which become boundary conditions for the surfaces
    OOFEM_LOG_INFO("Computing xi on edges\n");
    for ( auto &setPointer : edgeSets ) {
        const IntArray &edges = setPointer.giveEdgeList();

        // Number the equations along this edge set
        std :: map< int, int > eqs;
        int eq_count = 0;
        for ( int n : setPointer.giveNodeList() ) {
            if ( totalCornerNodes.containsSorted(n) ) {
                eqs[n] = 0;
            } else {
                eqs[n] = ++eq_count;
            }
        }
        
        FloatMatrix f;
        FloatArray q;
        ///@todo Preallocation(?)
        std :: unique_ptr< SparseMtrx > K( classFactory.createSparseMtrx(solver->giveRecommendedMatrix(true)) );
        K->buildInternalStructure(domain->giveEngngModel(), eq_count, eq_count, {}, {});
        f.resize(eq_count, 3);

        K->assembleBegin();
        for ( int pos = 0; pos < edges.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( edges[pos * 2] );
            int edge = edges[pos * 2 + 1];

            FloatMatrix D, Ke, fe;
            FloatArray b;
            
            FEInterpolation3d *interp = static_cast< FEInterpolation3d* >( e->giveInterpolation() );
            const auto &bNodes = interp->boundaryEdgeGiveNodes(edge, e->giveGeometryType());
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryEdgeIntegrationRule(order, edge, e->giveGeometryType()) );
            static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(D, Capacity, e->giveDefaultIntegrationRulePtr()->getIntegrationPoint(1), tStep);

            // Compute integral of B'*D*B and N:
            for ( auto &gp: *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);

                double detJ = interp->boundaryEdgeGiveTransformationJacobian(edge, lcoords, cellgeo);
                interp->edgeEvaldNdxi(b, edge, lcoords, cellgeo);
                b.times(1. / detJ);
                double dL = detJ * gp->giveWeight();
                
                // Compute material property
                ///@todo Can't do this yet, problem with material model interface (doesn't know the GP material mode). This should be changed.
                //static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(D, Capacity, gp, tStep);
#if 1
                FloatArray t;
                for ( int i = 1; i <= b.giveSize(); ++i ) {
                    t.add(b.at(i), e->giveNode(bNodes.at(i))->giveCoordinates());
                }
                t.normalize();
                FloatArray tmp;
                tmp.beProductOf(D, t);
                double k = tmp.dotProduct(t);
#else
                double k = D.at(1,1);
#endif

                Ke.plusDyadSymmUpper(b, k * dL);
            }
            Ke.symmetrized();
            
            // Need the element-nodal coordinates for the RHS, as the associated location array:
            IntArray loc(bNodes.giveSize());
            FloatMatrix cvec(bNodes.giveSize(), 3);
            for ( int i = 1; i <= bNodes.giveSize(); ++i ) {
                int enode = bNodes.at(i);
                Node *n = e->giveNode(enode);
                const auto &x = n->giveCoordinates();
                cvec.at(i, 1) = x.at(1);
                cvec.at(i, 2) = x.at(2);
                cvec.at(i, 3) = x.at(3);
                loc.at(i) = eqs[n->giveNumber()];
            }
            
            fe.beProductOf(Ke, cvec);
            fe.negated();
            f.assemble(fe, loc, {1, 2, 3});
            K->assemble(loc, Ke);
        }
        K->assembleEnd();

        FloatMatrix x;
        solver->solve(*K, f, x);

        for ( int n : setPointer.giveNodeList() ) {
            int eq = eqs[n];
            if ( eq > 0 ) {
                this->xis[n] = {x.at(eq, 1), x.at(eq, 2), x.at(eq, 3)};
            }
        }
    }
#endif

    OOFEM_LOG_INFO("Computing xi on surface sets\n");
#if 1
    // Surfaces use the edge solutions are boundary conditions:
    for ( auto &setNum : surfSets ) {
        Set *setPointer = this->giveDomain()->giveSet(setNum);
        const IntArray &surfs = setPointer->giveBoundaryList();

        // Number the equations along this surface set
        std :: map< int, int > eqs;
        int eq_count = 0;
        for ( int n : setPointer->giveNodeList() ) {
            if ( totalEdgeNodes.containsSorted(n) ) {
                eqs[n] = 0;
            } else {
                eqs[n] = ++eq_count;
            }
        }

        FloatMatrix f;
        ///@note We can use the single constraint, but this assumes flat surfaces
        FloatArray q;
        std :: unique_ptr< SparseMtrx > K( classFactory.createSparseMtrx(solver->giveRecommendedMatrix(true)) );
        K->buildInternalStructure(domain->giveEngngModel(), eq_count, eq_count, {}, {});
        f.resize(eq_count, 3);
        q.resize(eq_count);

        K->assembleBegin();
        ///@todo Preallocation
        for ( int pos = 0; pos < surfs.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( surfs[pos * 2] );
            int surf = surfs[pos * 2 + 1];

            FloatMatrix D, B, dNdx, DB, Ke, fe;
            FloatArray n, qe, normal;

            FEInterpolation3d *interp = static_cast< FEInterpolation3d* >( e->giveInterpolation() );
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir( interp->giveBoundaryIntegrationRule(order, surf, e->giveGeometryType()) );
            static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(D, Capacity, e->giveDefaultIntegrationRulePtr()->getIntegrationPoint(1), tStep);

            for ( auto &gp: *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);

                interp->boundaryEvalN(n, surf, lcoords, cellgeo);

                double detJ = interp->boundaryEvalNormal(normal, surf, lcoords, cellgeo);
                interp->surfaceEvaldNdx(dNdx, surf, lcoords, cellgeo);
                B.beTranspositionOf(dNdx);
                double dA = detJ * gp->giveWeight();
                
                for (int i = 1; i <= B.giveNumberOfRows(); ++i ) {
                    for (int j = 1; j <= B.giveNumberOfColumns(); ++j ) {
                        double tmp = 0;
                        for (int k = 1; k <= 3; ++k ) {
                            tmp += normal.at(k) * B.at(k,j);
                        }
                        B.at(i,j) -= normal.at(i) * tmp;
                    }
                }

                // Compute material property
                ///@todo Can't do this yet, problem with material model interface (doesn't know the GP material mode). This should be changed.
                //static_cast< TransportElement* >(e)->computeConstitutiveMatrixAt(D, Capacity, gp, tStep);

                // Vector:
                DB.beProductOf(D, B);
                Ke.plusProductSymmUpper(B, DB, dA);
                qe.add(dA, n);
            }
            Ke.symmetrized();

            const auto &bNodes = interp->boundaryGiveNodes(surf, e->giveGeometryType());
            IntArray loc(bNodes.giveSize());
            FloatMatrix cvec(bNodes.giveSize(), 3);
            for ( int i = 1; i <= bNodes.giveSize(); ++i ) {
                int enode = bNodes.at(i);
                Node *n = e->giveNode(enode);
                auto x = n->giveCoordinates() + this->xis[n->giveNumber()];
                cvec.at(i, 1) = x.at(1);
                cvec.at(i, 2) = x.at(2);
                cvec.at(i, 3) = x.at(3);
                loc.at(i) = eqs[n->giveNumber()];
            }
            fe.beProductOf(Ke, cvec);
            fe.negated();

            K->assemble(loc, Ke);
            f.assemble(fe, loc, {1, 2, 3});
            q.assemble(qe, loc);
        }
        K->assembleEnd();

        // Solve with constraints:
        // [K   Q] [xi    ] = [f]
        // [Q^T 0] [lambda ]   [0]
        // ==>
        // [Q^T . K^(-1) . Q] . lamda = Q^T . K^(-1) . f
        // xi = K^(-1) . f - K^(-1) . Q . lambda
        // alternatively:
        // [Q^T . Qx] . lamda = Q^T . fx
        // xi = fx - Qx . lambda

        double qTKiq;
        FloatMatrix x;
        FloatArray qx, lambda, tmp, qTKif;
        solver->solve(*K, q, qx);
        solver->solve(*K, f, x);
        
        qTKif.beTProductOf(x, q);
        qTKiq = q.dotProduct(qx);
        lambda.beScaled(1./qTKiq, qTKif);
        x.plusDyadUnsym(qx, lambda, -1.0);

        for ( int n : setPointer->giveNodeList() ) {
            int eq = eqs[n];
            if ( eq > 0 ) {
                this->xis[n] = {x.at(eq, 1), x.at(eq, 2), x.at(eq, 3)};
            }
        }

    }
#endif
}

} // end namespace oofem
