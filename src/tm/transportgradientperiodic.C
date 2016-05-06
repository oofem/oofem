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

#include "transportgradientperiodic.h"
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

namespace oofem {
REGISTER_BoundaryCondition(TransportGradientPeriodic);

TransportGradientPeriodic :: TransportGradientPeriodic(int n, Domain *d) : ActiveBoundaryCondition(n, d)//, PrescribedGradientHomogenization(),
    mGradient( new Node(1, d) )
{
    int nsd = d->giveNumberOfSpatialDimensions();
    // The prescribed strains.
    for ( int i = 0; i < nsd; i++ ) {
        int dofid = d->giveNextFreeDofID();
        grad_ids.followedBy(dofid);
        // Just putting in X_i id-items since they don't matter.
        // These don't actually need to be active, they are masterdofs with prescribed values, its
        // easier to just have them here rather than trying to make another Dirichlet boundary condition.
        //strain->appendDof( new ActiveDof( mGradient.get(), (DofIDItem)dofid, this->giveNumber() ) );
        mGradient->appendDof( new MasterDof(mGradient.get(), this->giveNumber(), 0, (DofIDItem)dofid ) );
    }
}


TransportGradientPeriodic :: ~TransportGradientPeriodic()
{
}


int TransportGradientPeriodic :: giveNumberOfInternalDofManagers()
{
    return 1;
}


DofManager *TransportGradientPeriodic :: giveInternalDofManager(int i)
{
    return this->mGradient.get();
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

    this->slavemap.clear();
    for ( int inode : nodes ) {
        Node *masterNode = NULL;
        Node *node = this->domain->giveNode(inode);
        const FloatArray &masterCoord = *node->giveCoordinates();
        //printf("node %d\n", node->giveLabel()); masterCoord.printYourself();
        // The difficult part, what offset to subtract to find the master side;
        for ( FloatArray &testJump : jumps ) {
            coord.beDifferenceOf(masterCoord, testJump);
            masterNode = sl->giveNodeClosestToPoint(coord, fabs(jump.at(1))*1e-5);
            if ( masterNode != NULL ) {
                //printf("Found master (%d) to node %d (distance = %e)\n",  masterNode->giveNumber(), node->giveNumber(),
                //       masterNode->giveCoordinates()->distance(coord));
                break;
            }
        }
        if ( masterNode != NULL ) {
            this->slavemap.insert({node->giveNumber(), masterNode->giveNumber()});
        } else {
            OOFEM_ERROR("Couldn't find master node!");
        }
    }
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
        return NULL;
    }
    if ( mdof == 1 ) {
        int node = this->slavemap[dof->giveDofManager()->giveNumber()];
        return this->domain->giveDofManager(node)->giveDofWithID(dof->giveDofID());
    } else {
        return this->mGradient->giveDofWithID(this->grad_ids[mdof-2]);
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
    DofIDEquationNumbering pnum(true, this->grad_ids);
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM > solver( classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx( stype ) );
    std :: unique_ptr< SparseMtrx > Kfp( classFactory.createSparseMtrx( stype ) );
    std :: unique_ptr< SparseMtrx > Kpp( classFactory.createSparseMtrx( stype ) );

    Kff->buildInternalStructure(rve, this->domain->giveNumber(), fnum);
    Kfp->buildInternalStructure(rve, this->domain->giveNumber(), fnum, pnum);
    Kpp->buildInternalStructure(rve, this->domain->giveNumber(), pnum);
    rve->assemble(*Kff, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    rve->assemble(*Kfp, tStep, TangentAssembler(TangentStiffness), fnum, pnum, this->domain);
    rve->assemble(*Kpp, tStep, TangentAssembler(TangentStiffness), pnum, this->domain);

    int neq = Kfp->giveNumberOfRows();
    int nsd = this->domain->giveNumberOfSpatialDimensions();

    FloatMatrix grad_pert(nsd, nsd), rhs, sol(neq, nsd);
    grad_pert.resize(nsd, nsd);
    grad_pert.beUnitMatrix();

    // Compute the solution to each of the pertubation of eps
    Kfp->times(grad_pert, rhs);
    solver->solve(*Kff, rhs, sol);

    // Compute the solution to each of the pertubation of eps
    Kfp->timesT(sol, k); // Assuming symmetry of stiffness matrix
    // This is probably always zero, but for generality
    FloatMatrix tmpMat;
    Kpp->times(grad_pert, tmpMat);
    k.subtract(tmpMat);
    k.times( - 1.0 / ( this->domainSize(this->giveDomain(), this->set) + this->domainSize(this->giveDomain(), this->masterSet) ));
}


void TransportGradientPeriodic :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();
    FloatArray *masterCoords = master->giveCoordinates();

    FloatArray dx;
    dx.beDifferenceOf(* coords, * masterCoords );

    masterContribs.resize(dx.giveSize() + 1);

    masterContribs.at(1) = 1.; // Master dof is always weight 1.0
    for ( int i = 1; i <= dx.giveSize(); ++i ) {
        masterContribs.at(i+1) = dx.at(i);
    }
}


double TransportGradientPeriodic :: giveUnknown(double val, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    FloatArray dx, grad;
    dx.beDifferenceOf(* dof->giveDofManager()->giveCoordinates(), * master->giveCoordinates());
    this->mGradient->giveUnknownVector(grad, mode, this->grad_ids, tStep);
    return val + grad.dotProduct(dx);
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
    return this->mGradient.get() == dof->giveDofManager();
}


IRResultType TransportGradientPeriodic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, this->mGradient, _IFT_TransportGradientPeriodic_gradient)
    IR_GIVE_FIELD(ir, this->mCenterCoords, _IFT_TransportGradientPeriodic_centerCoords)

    IR_GIVE_FIELD(ir, this->masterSet, _IFT_TransportGradientPeriodic_masterSet)
    IR_GIVE_FIELD(ir, this->jump, _IFT_TransportGradientPeriodic_jump)

    //ActiveBoundaryCondition :: initializeFrom(ir);
    return PrescribedGradientHomogenization::initializeFrom(ir);
}


void TransportGradientPeriodic :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    //PrescribedGradientHomogenization :: giveInputRecord(input);
    input.setField(this->mGradient, _IFT_TransportGradientPeriodic_gradient);
    input.setField(this->mCenterCoords, _IFT_TransportGradientPeriodic_centerCoords);
    
    input.setField(this->masterSet, _IFT_TransportGradientPeriodic_masterSet);
    input.setField(this->jump, _IFT_TransportGradientPeriodic_jump);
}


void TransportGradientPeriodic :: postInitialize()
{
    this->findSlaveToMasterMap();
}
} // end namespace oofem
