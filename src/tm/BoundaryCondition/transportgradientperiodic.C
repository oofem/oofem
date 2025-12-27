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

#include "tm/BoundaryCondition/transportgradientperiodic.h"
#include "dofiditem.h"
#include "node.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "engngm.h"
#include "node.h"
#include "activedof.h"
#include "masterdof.h"
#include "classfactory.h"
#include "sparsemtrxtype.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "dynamicinputrecord.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "feinterpol.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "timestep.h"
#include "mathfem.h"

namespace oofem {
REGISTER_BoundaryCondition(TransportGradientPeriodic);

TransportGradientPeriodic :: TransportGradientPeriodic(int n, Domain *d) : ActiveBoundaryCondition(n, d), //PrescribedGradientHomogenization(),
    grad( new Node(1, d) )
{
    int nsd = d->giveNumberOfSpatialDimensions();
    // The prescribed strains.
    for ( int i = 0; i < nsd; i++ ) {
        int dofid = d->giveNextFreeDofID();
        grad_ids.followedBy(dofid);
        // Just putting in X_i id-items since they don't matter.
        // These don't actually need to be active, they are masterdofs with prescribed values, its
        // easier to just have them here rather than trying to make another Dirichlet boundary condition.
        grad->appendDof( new ActiveDof( grad.get(), (DofIDItem)dofid, this->giveNumber() ) );
        //grad->appendDof( new MasterDof(grad.get(), this->giveNumber(), 0, (DofIDItem)dofid ) );
    }
}


int TransportGradientPeriodic :: giveNumberOfInternalDofManagers()
{
    return 1;
}


DofManager *TransportGradientPeriodic :: giveInternalDofManager(int i)
{
    return this->grad.get();
}


void TransportGradientPeriodic :: findSlaveToMasterMap()
{
    FloatArray coord;
    SpatialLocalizer *sl = this->domain->giveSpatialLocalizer();
    //Set *masterSet = this->domain->giveSet(2);
    const IntArray &nodes = this->domain->giveSet(this->set)->giveNodeList(); // Split into slave set and master set?
    int nsd = jump.giveSize();
    std :: vector< FloatArray > jumps;
    // Construct all the possible jumps;
    jumps.reserve((2 << (nsd-1)) - 1);
    if ( nsd != 3 ) {
        OOFEM_ERROR("Only 3d implemented yet!");
    }
    jumps.emplace_back(jump);
    jumps.emplace_back(FloatArray{jump.at(1), jump.at(2), 0.});
    jumps.emplace_back(FloatArray{jump.at(1), 0., jump.at(3)});
    jumps.emplace_back(FloatArray{0., jump.at(2), jump.at(3)});
    jumps.emplace_back(FloatArray{jump.at(1), 0., 0.});
    jumps.emplace_back(FloatArray{0., jump.at(2), 0.});
    jumps.emplace_back(FloatArray{0., 0., jump.at(3)});

    double maxdist = jump.computeNorm()*1e-5;

    this->slavemap.clear();
    for ( int inode : nodes ) {
        Node *masterNode = nullptr;
        Node *node = this->domain->giveNode(inode);
        const auto &masterCoord = node->giveCoordinates();
        //printf("node %d\n", node->giveLabel()); masterCoord.printYourself();
        // The difficult part, what offset to subtract to find the master side;
        for ( FloatArray &testJump : jumps ) {
            coord.beDifferenceOf(masterCoord, testJump);
            masterNode = sl->giveNodeClosestToPoint(coord, maxdist);
            if ( masterNode ) {
                //printf("Found master (%d) to node %d (distance = %e)\n",  masterNode->giveNumber(), node->giveNumber(),
                //       distance(*masterNode->giveCoordinates(), coord));
                break;
            }
        }
        if ( masterNode ) {
            this->slavemap.insert({node->giveNumber(), masterNode->giveNumber()});
        } else {
            OOFEM_ERROR("Couldn't find master node!");
        }
    }
}


double TransportGradientPeriodic :: domainSize(Domain *d, int setNum)
{
    int nsd = d->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = d->giveSet(setNum);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = d->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
}


int TransportGradientPeriodic :: giveNumberOfMasterDofs(ActiveDof *dof)
{
    if ( this->isGradDof(dof) ) {
        return 1;
    }
    return this->giveDomain()->giveNumberOfSpatialDimensions() + 1;
}


Dof *TransportGradientPeriodic :: giveMasterDof(ActiveDof *dof, int mdof)
{
    if ( this->isGradDof(dof) ) {
        return nullptr;
    }
    if ( mdof == 1 ) {
        int node = this->slavemap[dof->giveDofManager()->giveNumber()];
        return this->domain->giveDofManager(node)->giveDofWithID(dof->giveDofID());
    } else {
        return this->grad->giveDofWithID(this->grad_ids[mdof-2]);
    }
}


void TransportGradientPeriodic :: computeField(FloatArray &flux, TimeStep *tStep)
{
    DofIDEquationNumbering pnum(true, grad_ids);
    EngngModel *emodel = this->giveDomain()->giveEngngModel();
    FloatArray tmp;
    int npeq = grad_ids.giveSize();
    // sigma = residual (since we use the slave dofs) = f_ext - f_int
    flux.resize(npeq);
    flux.zero();
    emodel->assembleVector(flux, tStep, InternalForceAssembler(), VM_Total, pnum, this->domain);
    tmp.resize(npeq);
    tmp.zero();
    emodel->assembleVector(tmp, tStep, ExternalForceAssembler(), VM_Total, pnum, this->domain);
    flux.subtract(tmp);
    // Divide by the RVE-volume
    flux.times(1.0 / ( this->domainSize(this->giveDomain(), this->set) + this->domainSize(this->giveDomain(), this->masterSet) ));
}


void TransportGradientPeriodic :: computeTangent(FloatMatrix &k, TimeStep *tStep)
{
    EModelDefaultEquationNumbering fnum;
    //DofIDEquationNumbering pnum(true, this->grad_ids);
    EModelDefaultPrescribedEquationNumbering pnum;
    int nsd = this->domain->giveNumberOfSpatialDimensions();

    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM > solver( classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx( stype ) );
    std :: unique_ptr< SparseMtrx > Kfp( classFactory.createSparseMtrx( stype ) );
    std :: unique_ptr< SparseMtrx > Kpp( classFactory.createSparseMtrx( stype ) );

    Kff->buildInternalStructure(rve, this->domain->giveNumber(), fnum);
    int neq = Kff->giveNumberOfRows();
    Kfp->buildInternalStructure(rve, this->domain->giveNumber(), fnum, pnum);
    Kpp->buildInternalStructure(rve, this->domain->giveNumber(), pnum);
    //Kfp->buildInternalStructure(rve, neq, nsd, {}, {});
    //Kpp->buildInternalStructure(rve, nsd, nsd, {}, {});
#if 1
    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    rve->assemble(*Kfp, tStep, TangentAssembler(TangentStiffness), fnum, pnum, this->domain);
    rve->assemble(*Kpp, tStep, TangentAssembler(TangentStiffness), pnum, this->domain);
#else
    auto ma = TangentAssembler(TangentStiffness);
    IntArray floc, ploc;
    FloatMatrix mat, R;

    int nelem = domain->giveNumberOfElements();
#ifdef _OPENMP
 #pragma omp parallel for shared(Kff, Kfp, Kpp) private(mat, R, floc, ploc)
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
                Kfp->assemble(floc, ploc, mat);
                Kpp->assemble(ploc, mat);
            }
        }
    }
    Kff->assembleBegin();
    Kfp->assembleBegin();
    Kpp->assembleBegin();

    Kff->assembleEnd();
    Kfp->assembleEnd();
    Kpp->assembleEnd();
#endif

    FloatMatrix grad_pert(nsd, nsd), rhs, sol(neq, nsd);
    grad_pert.resize(nsd, nsd);
    grad_pert.beUnitMatrix();
    // Workaround since the matrix size is inflexible with custom dof numbering (so far, planned to be fixed).
    IntArray grad_loc;
    this->grad->giveLocationArray(this->grad_ids, grad_loc, pnum);
    FloatMatrix pert(Kpp->giveNumberOfRows(), nsd);
    pert.assemble(grad_pert, grad_loc, {1,2,3});
    //pert.printYourself("pert");

    //printf("Kfp = %d x %d\n", Kfp->giveNumberOfRows(), Kfp->giveNumberOfColumns());
    //printf("Kff = %d x %d\n", Kff->giveNumberOfRows(), Kff->giveNumberOfColumns());
    //printf("Kpp = %d x %d\n", Kpp->giveNumberOfRows(), Kpp->giveNumberOfColumns());

    // Compute the solution to each of the pertubation of eps
    Kfp->times(pert, rhs);
    //rhs.printYourself("rhs");

    // Initial guess (Taylor assumption) helps KSP-iterations
    for ( auto &n : domain->giveDofManagers() ) {
        int k1 = n->giveDofWithID( this->dofs[0] )->__giveEquationNumber();
        if ( k1 ) {
            const auto &coords = n->giveCoordinates();
            for ( int i = 1; i <= nsd; ++i ) {
                sol.at(k1, i) = -(coords.at(i) - mCenterCoord.at(i));
            }
        }
    }

    if ( solver->solve(*Kff, rhs, sol) != CR_CONVERGED ) {
        OOFEM_ERROR("Failed to solve Kff");
    }
    // Compute the solution to each of the pertubation of eps
    Kfp->timesT(sol, k); // Assuming symmetry of stiffness matrix
    // This is probably always zero, but for generality
    FloatMatrix tmpMat;
    Kpp->times(pert, tmpMat);
    k.subtract(tmpMat);
    k.times( - 1.0 / ( this->domainSize(this->giveDomain(), this->set) + this->domainSize(this->giveDomain(), this->masterSet) ));
    
    // Temp workaround on sizing issue mentioned above:
    FloatMatrix k2 = k;
    k.beSubMatrixOf(k2, grad_loc, {1,2,3});
}


void TransportGradientPeriodic :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    const auto &coords = dof->giveDofManager()->giveCoordinates();
    const auto &masterCoords = master->giveCoordinates();

    FloatArray dx = coords - masterCoords;

    masterContribs.resize(dx.giveSize() + 1);

    masterContribs.at(1) = 1.; // Master dof is always weight 1.0
    for ( int i = 1; i <= dx.giveSize(); ++i ) {
        masterContribs.at(i+1) = dx.at(i);
    }
}


double TransportGradientPeriodic :: giveUnknown(double val, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    FloatArray dx, g;
    dx.beDifferenceOf(dof->giveDofManager()->giveCoordinates(), master->giveCoordinates());
    this->grad->giveUnknownVector(g, this->grad_ids, mode, tStep);
    return val + g.dotProduct(dx);
}


double TransportGradientPeriodic :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if ( this->isGradDof(dof) ) {
        int ind = grad_ids.findFirstIndexOf(dof->giveDofID()) - 1;
        return this->mGradient(ind) * this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
    }

    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(field, mode, tStep);
    return this->giveUnknown(val, mode, tStep, dof);
}


double TransportGradientPeriodic :: giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if ( this->isGradDof(dof) ) {
        int ind = grad_ids.findFirstIndexOf(dof->giveDofID()) - 1;
        return this->mGradient(ind) * this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
    }

    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);

    if ( mode == VM_Incremental ) {
        double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(mode, tStep);
        return this->giveUnknown(val, mode, tStep, dof);
    }
    double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(mode, tStep);
    return this->giveUnknown(val, mode, tStep, dof);
}


bool TransportGradientPeriodic :: isPrimaryDof(ActiveDof *dof)
{
    return this->isGradDof(dof);
}


double TransportGradientPeriodic :: giveBcValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    int index = grad_ids.findFirstIndexOf(dof->giveDofID()) - 1;
    return this->mGradient(index) * this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
}


bool TransportGradientPeriodic :: hasBc(Dof *dof, TimeStep *tStep)
{
    return this->isGradDof(dof);
}


bool TransportGradientPeriodic :: isGradDof(Dof *dof)
{
    return this->grad.get() == dof->giveDofManager();
}


void TransportGradientPeriodic :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    //PrescribedGradientHomogenization::initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->mGradient, _IFT_TransportGradientPeriodic_gradient)
    this->mCenterCoord = {0., 0., 0.};
    IR_GIVE_OPTIONAL_FIELD(ir, this->mCenterCoord, _IFT_TransportGradientPeriodic_centerCoords)

    IR_GIVE_FIELD(ir, this->masterSet, _IFT_TransportGradientPeriodic_masterSet)
    IR_GIVE_FIELD(ir, this->jump, _IFT_TransportGradientPeriodic_jump)
}


void TransportGradientPeriodic :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    //PrescribedGradientHomogenization :: giveInputRecord(input);
    input.setField(this->mGradient, _IFT_TransportGradientPeriodic_gradient);
    input.setField(this->mCenterCoord, _IFT_TransportGradientPeriodic_centerCoords);
    
    input.setField(this->masterSet, _IFT_TransportGradientPeriodic_masterSet);
    input.setField(this->jump, _IFT_TransportGradientPeriodic_jump);
}


void TransportGradientPeriodic :: postInitialize()
{
    this->findSlaveToMasterMap();
}
} // end namespace oofem
