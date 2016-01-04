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

#include "prescribedgradientbcperiodic.h"
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
REGISTER_BoundaryCondition(PrescribedGradientBCPeriodic);

PrescribedGradientBCPeriodic :: PrescribedGradientBCPeriodic(int n, Domain *d) : ActiveBoundaryCondition(n, d), PrescribedGradientHomogenization(),
    strain( new Node(1, d) )
{
    // The unknown volumetric strain
    int nsd = d->giveNumberOfSpatialDimensions();
    int components = nsd * nsd;
    // The prescribed strains.
    for ( int i = 0; i < components; i++ ) {
        int dofid = d->giveNextFreeDofID();
        strain_id.followedBy(dofid);
        // Just putting in X_i id-items since they don't matter.
        // These don't actually need to be active, they are masterdofs with prescribed values, its
        // easier to just have them here rather than trying to make another Dirichlet boundary condition.
        strain->appendDof( new ActiveDof( strain.get(), (DofIDItem)dofid, this->giveNumber() ) );
        //strain->appendDof( new MasterDof(strain.get(), this->giveNumber(), 0, (DofIDItem)dofid ) );
    }
}


PrescribedGradientBCPeriodic :: ~PrescribedGradientBCPeriodic()
{
}


int PrescribedGradientBCPeriodic :: giveNumberOfInternalDofManagers()
{
    return 1;
}


DofManager *PrescribedGradientBCPeriodic :: giveInternalDofManager(int i)
{
    return this->strain.get();
}


void PrescribedGradientBCPeriodic :: findSlaveToMasterMap()
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


int PrescribedGradientBCPeriodic :: giveNumberOfMasterDofs(ActiveDof *dof)
{
    if ( this->isStrainDof(dof) ) {
        return 1;
    }
    return this->giveDomain()->giveNumberOfSpatialDimensions() + 1;
}


Dof *PrescribedGradientBCPeriodic :: giveMasterDof(ActiveDof *dof, int mdof)
{
    if ( this->isStrainDof(dof) ) {
        return NULL;
    }
    if ( mdof == 1 ) {
        int node = this->slavemap[dof->giveDofManager()->giveNumber()];
        //printf("dofid = %d, slave node = %d, master node = %d\n", dof->giveDofID(),dof->giveDofManager()->giveNumber(), node );
        //this->domain->giveDofManager(node)->printYourself();
        return this->domain->giveDofManager(node)->giveDofWithID(dof->giveDofID());
    } else {
        DofIDItem dofid = dof->giveDofID();
        FloatArray *coords = dof->giveDofManager()->giveCoordinates();
        int nsd = coords->giveSize();
        if ( dofid == D_u || dofid == V_u ) {
            return this->strain->giveDofWithID(strain_id[nsd*(mdof-2)]);
        } else if ( dofid == D_v || dofid == V_v ) {
            return this->strain->giveDofWithID(strain_id[nsd*(mdof-2)+1]);
        } else /* if ( dofid == D_u || dofid == V_u ) */ {
            return this->strain->giveDofWithID(strain_id[nsd*(mdof-2)+2]);
        }
    }
}


void PrescribedGradientBCPeriodic :: computeField(FloatArray &sigma, TimeStep *tStep)
{
    DofIDEquationNumbering pnum(true, strain_id);
    EngngModel *emodel = this->giveDomain()->giveEngngModel();
    FloatArray tmp, sig_tmp;
    int npeq = strain_id.giveSize();
    // sigma = residual (since we use the slave dofs) = f_ext - f_int
    sig_tmp.resize(npeq);
    sig_tmp.zero();
    emodel->assembleVector(sig_tmp, tStep, InternalForceAssembler(), VM_Total, pnum, this->domain);
    tmp.resize(npeq);
    tmp.zero();
    emodel->assembleVector(tmp, tStep, ExternalForceAssembler(), VM_Total, pnum, this->domain);
    sig_tmp.subtract(tmp);
    // Divide by the RVE-volume
    sig_tmp.times(1.0 / ( this->domainSize(this->giveDomain(), this->set) + this->domainSize(this->giveDomain(), this->masterSet) ));

    sigma.resize(sig_tmp.giveSize());
    if ( sig_tmp.giveSize() == 9 ) {
        sigma.assemble(sig_tmp, {1, 9, 8, 6, 2, 7, 5, 4, 3});
    } else if ( sig_tmp.giveSize() == 4 ) {
        sigma.assemble(sig_tmp, {1, 4, 3, 2});
    } else {
        sigma = sig_tmp;
    }
}


void PrescribedGradientBCPeriodic :: computeTangent(FloatMatrix &E, TimeStep *tStep)
{
    EModelDefaultEquationNumbering fnum;
    DofIDEquationNumbering pnum(true, strain_id);
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
    int ncomp = nsd * nsd;

    FloatMatrix grad_pert(ncomp, ncomp), rhs, sol(neq, ncomp);
    grad_pert.resize(ncomp, ncomp); // In fact, npeq should most likely equal ndev
    grad_pert.beUnitMatrix();

    // Compute the solution to each of the pertubation of eps
    Kfp->times(grad_pert, rhs);
    solver->solve(*Kff, rhs, sol);

    // Compute the solution to each of the pertubation of eps
    FloatMatrix E_tmp;
    Kfp->timesT(sol, E_tmp); // Assuming symmetry of stiffness matrix
    // This is probably always zero, but for generality
    FloatMatrix tmpMat;
    Kpp->times(grad_pert, tmpMat);
    E_tmp.subtract(tmpMat);
    E_tmp.times( - 1.0 / ( this->domainSize(this->giveDomain(), this->set) + this->domainSize(this->giveDomain(), this->masterSet) ));
    
    E.resize(E_tmp.giveNumberOfRows(), E_tmp.giveNumberOfColumns());
    if ( nsd == 3 ) {
        E.assemble(E_tmp, {1, 9, 8, 6, 2, 7, 5, 4, 3});
    } else if ( nsd == 2 ) {
        E.assemble(E_tmp, {1, 4, 3, 2});
    } else {
        E = E_tmp;
    }
}


void PrescribedGradientBCPeriodic :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();
    FloatArray *masterCoords = master->giveCoordinates();

    FloatArray dx;
    dx.beDifferenceOf(* coords, * masterCoords );

    int nsd = dx.giveSize(); // Number of spatial dimensions

    masterContribs.resize(nsd + 1);

    masterContribs.at(1) = 1.; // Master dof is always weight 1.0
    for ( int i = 1; i <= dx.giveSize(); ++i ) {
        masterContribs.at(i+1) = dx.at(i);
    }
}


double PrescribedGradientBCPeriodic :: giveUnknown(double val, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();
    FloatArray *masterCoords = master->giveCoordinates();
    FloatArray dx, uM;
    dx.beDifferenceOf(* coords, * masterCoords );

    int ind;
    if ( id == D_u || id == V_u || id == P_f || id == T_f ) {
        ind = 1;
    } else if ( id == D_v || id == V_v ) {
        ind = 2;
    } else { /*if ( id == D_w || id == V_w )*/   // 3D only:
        ind = 3;
    }

    FloatMatrix grad(3, 3);
    for ( int i = 0; i < this->strain_id.giveSize(); ++i ) {
        Dof *dof = this->strain->giveDofWithID(strain_id[i]);
        grad(i % 3, i / 3) = dof->giveUnknown(mode, tStep);
    }
    uM.beProductOf(grad, dx); // The "jump" part of the unknown ( u^+ = [[u^M]] + u^- )

    return val + uM.at(ind);
}


double PrescribedGradientBCPeriodic :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if ( this->isStrainDof(dof) ) {
        int ind = strain_id.findFirstIndexOf(dof->giveDofID()) - 1;
        return this->mGradient(ind % 3, ind / 3) * this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
    }

    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(field, mode, tStep);
    return this->giveUnknown(val, mode, tStep, dof);
}


double PrescribedGradientBCPeriodic :: giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if ( this->isStrainDof(dof) ) {
        int ind = strain_id.findFirstIndexOf(dof->giveDofID()) - 1;
        return this->mGradient(ind % 3, ind / 3) * this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
    }

    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);

    if ( mode == VM_Incremental ) {
        double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(mode, tStep);
        return this->giveUnknown(val, mode, tStep, dof);
    }
    double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(mode, tStep);
    return this->giveUnknown(val, mode, tStep, dof);
}


bool PrescribedGradientBCPeriodic :: isPrimaryDof(ActiveDof *dof)
{
    return this->isStrainDof(dof);
}


double PrescribedGradientBCPeriodic :: giveBcValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    if ( this->isStrainDof(dof) ) {
        int index = strain_id.findFirstIndexOf(dof->giveDofID()) - 1;
        return this->mGradient( index % 3, index / 3 ) * this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());;
    }
    OOFEM_ERROR("Has no prescribed value from bc.");
    return 0.0;
}


bool PrescribedGradientBCPeriodic :: hasBc(Dof *dof, TimeStep *tStep)
{
    return this->isStrainDof(dof);
}


bool PrescribedGradientBCPeriodic :: isStrainDof(Dof *dof)
{
    return this->strain.get() == dof->giveDofManager();
}


IRResultType PrescribedGradientBCPeriodic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, this->masterSet, _IFT_PrescribedGradientBCPeriodic_masterSet)
    IR_GIVE_FIELD(ir, this->jump, _IFT_PrescribedGradientBCPeriodic_jump)

    ActiveBoundaryCondition :: initializeFrom(ir);
    return PrescribedGradientHomogenization::initializeFrom(ir);
}


void PrescribedGradientBCPeriodic :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedGradientHomogenization :: giveInputRecord(input);
    input.setField(this->masterSet, _IFT_PrescribedGradientBCPeriodic_masterSet);
    input.setField(this->jump, _IFT_PrescribedGradientBCPeriodic_jump);
}


void PrescribedGradientBCPeriodic :: postInitialize()
{
    this->findSlaveToMasterMap();
}
} // end namespace oofem
