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

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCPeriodic);

PrescribedGradientBCPeriodic :: PrescribedGradientBCPeriodic(int n, Domain *d) : PrescribedGradientBC(n, d)
{
    // The unknown volumetric strain
    int nsd = d->giveNumberOfSpatialDimensions();
    int components = nsd * nsd;
    // The prescribed strains.
    strain.reset(new Node(1, d));
    strain_id.clear();
    for ( int i = 0; i < components; i++ ) {
        int dofid = d->giveNextFreeDofID();
        strain_id.followedBy(dofid);
        // Just putting in X_i id-items since they don't matter.
        // These don't actually need to be active, they are masterdofs with prescribed values, its
        // easier to just have them here rather than trying to make another Dirichlet boundary condition.
        strain->appendDof( new ActiveDof( strain.get(), (DofIDItem)dofid, this->giveNumber() ) );
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


double PrescribedGradientBCPeriodic :: domainSize()
{
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundaries.at(pos * 2), FEIElementGeometryWrapper(e) );
    }

    Set *masterSet = this->giveDomain()->giveSet(this->masterSet);
    const IntArray &masterBoundaries = masterSet->giveBoundaryList();
    for ( int pos = 1; pos <= masterBoundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( masterBoundaries.at(pos * 2 - 1) );
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( masterBoundaries.at(pos * 2), FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
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
        // The difficult part, what offset to subtract to find the master side;
        for ( FloatArray &testJump : jumps ) {
            coord.beDifferenceOf(masterCoord, testJump);
            masterNode = sl->giveNodeClosestToPoint(coord);
            if ( masterNode->giveCoordinates()->distance(coord) <= fabs(jump.at(1))*1e-6) {
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


void PrescribedGradientBCPeriodic :: computeFields(FloatArray &sigma, TimeStep *tStep)
{
    EngngModel *emodel = this->giveDomain()->giveEngngModel();
    FloatArray tmp;
    int npeq = emodel->giveNumberOfDomainEquations( this->giveDomain()->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    // sigma = residual (since we use the slave dofs) = f_ext - f_int
    sigma.resize(npeq);
    sigma.zero();
    emodel->assembleVector(sigma, tStep, InternalForcesVector, VM_Total, EModelDefaultPrescribedEquationNumbering(), this->domain);
    tmp.resize(npeq);
    tmp.zero();
    emodel->assembleVector(tmp, tStep, ExternalForcesVector, VM_Total, EModelDefaultPrescribedEquationNumbering(), this->domain);
    sigma.subtract(tmp);
    // Divide by the RVE-volume
    sigma.times( 1.0 / this->domainSize() );
}


void PrescribedGradientBCPeriodic :: computeTangents(FloatMatrix &E, TimeStep *tStep)
{
    ///@todo Implement tangent computations
    EModelDefaultEquationNumbering fnum;
    EModelDefaultPrescribedEquationNumbering pnum;
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    SparseLinearSystemNM *solver = classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx( stype ) );

    Kff->buildInternalStructure(rve, 1, fnum);
    rve->assemble(Kff.get(), tStep, StiffnessMatrix, fnum, this->domain);

    OOFEM_ERROR("Not implemented yet");

    delete solver; ///@todo Remove this when solver is taken from engngmodel
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
    uM.beProductOf(this->mGradient, dx); // The "jump" part of the unknown ( u^+ = [[u^M]] + u^- )

    int ind;
    if ( id == D_u || id == V_u || id == P_f || id == T_f ) {
        ind = 1;
    } else if ( id == D_v || id == V_v ) {
        ind = 2;
    } else { /*if ( id == D_w || id == V_w )*/   // 3D only:
        ind = 3;
    }
    return val + uM.at(ind);
}


double PrescribedGradientBCPeriodic :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if ( this->isStrainDof(dof) ) {
        int ind = strain_id.findFirstIndexOf(dof->giveDofID());
        return this->mGradient(ind % 3, ind / 3);
    }

    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(field, mode, tStep);
    return this->giveUnknown(val, mode, tStep, dof);
}


double PrescribedGradientBCPeriodic :: giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
    if ( this->isStrainDof(dof) ) {
        int ind = strain_id.findFirstIndexOf(dof->giveDofID());
        return this->mGradient(ind % 3, ind / 3);
    }

    DofManager *master = this->domain->giveDofManager(this->slavemap[dof->giveDofManager()->giveNumber()]);
    double val = master->giveDofWithID(dof->giveDofID())->giveUnknown(mode, tStep);
    return this->giveUnknown(val, mode, tStep, dof);
}


bool PrescribedGradientBCPeriodic :: isPrimaryDof(ActiveDof *dof)
{
    return this->isStrainDof(dof);
}


double PrescribedGradientBCPeriodic :: giveBcValue(ActiveDof *dof, ValueModeType mode, TimeStep *tStep)
{
    if ( this->isStrainDof(dof) ) {
        int index = strain_id.findFirstIndexOf(dof->giveDofID());
        return this->mGradient( index % 3, index / 3 );
    }
    OOFEM_ERROR("Has no prescribed value from bc.");
    return 0.0;
}


bool PrescribedGradientBCPeriodic :: hasBc(ActiveDof *dof, TimeStep *tStep)
{
    return this->isStrainDof(dof);
}


bool PrescribedGradientBCPeriodic :: isStrainDof(ActiveDof *dof)
{
    return this->strain.get() == dof->giveDofManager();
}


IRResultType PrescribedGradientBCPeriodic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, this->masterSet, _IFT_PrescribedGradientBCPeriodic_masterSet)
    IR_GIVE_FIELD(ir, this->jump, _IFT_PrescribedGradientBCPeriodic_jump)

    return PrescribedGradientBC :: initializeFrom(ir);
}


void PrescribedGradientBCPeriodic :: giveInputRecord(DynamicInputRecord &input)
{
    PrescribedGradientBC :: giveInputRecord(input);
    input.setField(this->masterSet, _IFT_PrescribedGradientBCPeriodic_masterSet);
    input.setField(this->jump, _IFT_PrescribedGradientBCPeriodic_jump);
}


void PrescribedGradientBCPeriodic :: postInitialize()
{
    this->findSlaveToMasterMap();
}
} // end namespace oofem