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

#include "Contact/contactelement.h"
#include "dofmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "valuemodetype.h"
#include "dofiditem.h"
#include "timestep.h"
#include "dof.h"
#include "masterdof.h"
#include "unknownnumberingscheme.h"
#include "domain.h"


namespace oofem {

Node2NodeContact :: Node2NodeContact(DofManager *master, DofManager *slave) : ContactElement()
{   
    this->masterNode = master;
    this->slaveNode = slave;
    this->area = 1.0;
    this->eps = 1.0e11; // penalty
};   
  

void
Node2NodeContact :: computeGap(FloatArray &answer, TimeStep *tStep)
{

    FloatArray xs, xm, uS, uM;
    xs = *this->slaveNode->giveCoordinates();
    xm = *this->masterNode->giveCoordinates();
    this->slaveNode->giveUnknownVector(uS, {D_u, D_v, D_w}, VM_Total, tStep, true);
    this->masterNode->giveUnknownVector(uM, {D_u, D_v, D_w}, VM_Total, tStep, true);
    uS.printYourself("uS");
    uM.printYourself("uM");
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs-xm;
    
    // debug - set normal to x-axis - constant
    // need info from surrounding elements in order to compute an average normal at the point
    FloatArray normal = {1.0, 0.0, 0.0};
    answer = {dx.dotProduct(normal), 0.0, 0.0};
    
    printf("normal gap = %e \n", answer.at(1));
    if ( answer.at(1) < 0.0 ) {
        //printf("normal gap = %e \n", answer.at(1));
        this->inContact = true;
    } else {
        //this->inContact = false;
    }
    
}


void
Node2NodeContact :: computeCmatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *TimeStep)
{
  
    //giveNormal -not updated for node2node
    //TODO: debug
    FloatArray normal = {1.0, 0.0, 0.0};
    // C = {n -n}
    answer = {  normal.at(1),  normal.at(2),  normal.at(3),
               -normal.at(1), -normal.at(2), -normal.at(3) };
    
  
}



void
Node2NodeContact :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{
    // call to constitutive model
    // gap should be in a local system
    if ( gap.at(1) < 0.0 ) {
        t = this->eps * gap;
    } else {
        t = {0.0, 0.0, 0.0};
    }
  
}





void
Node2NodeContact :: computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    //Loop through all the master objects and let them do their thing
    FloatArray gap, C, Fc;
      
    // we have a node 2 node contact so we know that both parties should be nodes - use this?
    this->computeGap(gap, tStep);
    
    GaussPoint *gp = NULL; // Need to create integration points later
    // better name?
    FloatArray t;
    this->computeContactTractionAt(gp, t ,gap, tStep);
    // rotate to to global system
    this->computeCmatrixAt(gp, C, tStep);
    
    // compute load vector
    // for penalty: fc = C^T * traction * A, Area - optional par
    answer = t.at(1) * this->area * C;
      
    
  
}

  
void
Node2NodeContact :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
   // Need to set up an integration rule 
    GaussPoint *gp = NULL;
    
    FloatArray C;
    this->computeCmatrixAt(gp, C, tStep);
    answer.beDyadicProductOf(C,C);
    // this is the interface stiffness and should be obtained from that model
    answer.times( this->eps * this->area );
}
  
  
  
  


void
Node2NodeContact :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
    // should return location array for the master and the slave node
    // TODO this whole thing should be rewritten  
  
  
    answer.resize(6);
    answer.zero();
    //TODO should use a proper unknownnumberingscheme
    IntArray dofIdArray = {D_u, D_v, D_w};
    
    // master node
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->masterNode->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->masterNode->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(i) = s.giveDofEquationNumber( dof );
        } 
    }

    // slave node
    for ( int i = 1; i <= dofIdArray.giveSize(); i++ ) {
        if ( this->slaveNode->hasDofID( (DofIDItem)dofIdArray.at(i) ) ) { // add corresponding number
            Dof *dof= this->slaveNode->giveDofWithID( (DofIDItem)dofIdArray.at(i) );
            answer.at(3 + i) = s.giveDofEquationNumber( dof );
        } 
    }    
    
}    


// node 2 node lagrange


Node2NodeContactL :: Node2NodeContactL(DofManager *master, DofManager *slave) : Node2NodeContact(master, slave)
{   
    this->masterNode = master;
    this->slaveNode = slave;
    this->area = 1.0;
};   


int
Node2NodeContactL :: instanciateYourself(DataReader *dr)
{

    // Create dofs for coefficients
   
    //gammaDman = new Node(0, this->domain);
    //gamma_ids.clear();
    //for ( int i = 0; i < ndof; i++ ) {
        //int dofid = this->domain->giveNextFreeDofID();
        //gamma_ids.followedBy(dofid);
        //gammaDman->appendDof( new MasterDof( gammaDman, ( DofIDItem )dofid ) );
    //}
  
    // Add new dof to master dof manager -> new lagrange multiplier
    // TODO fix in a nicer way
    int dofid = this->masterNode->giveDomain()->giveNextFreeDofID();
    this->masterNode->appendDof( new MasterDof( this->masterNode, ( DofIDItem )dofid ) );
    this->lagrangeId = dofid;
    
    return 1;
}



void
Node2NodeContactL :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
  
    Node2NodeContact :: giveLocationArray(answer, s);
    
    //add lagrange dof
    if ( this->masterNode->hasDofID( (DofIDItem)this->lagrangeId ) ) { 
        Dof *dof= this->masterNode->giveDofWithID( (DofIDItem)this->lagrangeId );
        answer.followedBy( s.giveDofEquationNumber(dof) );
    }
    
}    



void
Node2NodeContactL :: computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
  
    //Loop through all the master objects and let them do their thing
    FloatArray gap, C, Fc;
      
    // we have a node 2 node contact so we know that both parties should be nodes - use this?
    this->computeGap(gap, tStep);
    
    GaussPoint *gp = NULL; // Need to create integration points later
    // better name?
    FloatArray t;
    //this->computeContactTractionAt(gp, t ,gap, tStep);
    bool flag =false;
    flag = masterNode->hasDofID( (DofIDItem)this->lagrangeId );
    
    IntArray temp2;
    masterNode->giveCompleteMasterDofIDArray(temp2);
    temp2.printYourself();
    Dof *dof = masterNode->giveDofWithID( this->lagrangeId );
    
    double lambda = dof->giveUnknown(mode, tStep);
    // rotate to to global system
    this->computeCmatrixAt(gp, C, tStep);
    
    // compute load vector
    // for penalty: fc = {C^T * lambda * A, gap_n}
    FloatArray temp = lambda *this->area * C;
    answer.resize( C.giveSize() + 1);
    answer.addSubVector(temp,1);
    answer.at( C.giveSize() + 1 ) = gap.at(1);
        
  
}
    


void
Node2NodeContactL :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
   // Need to set up an integration rule 
    GaussPoint *gp = NULL;
    
    FloatArray C;
    this->computeCmatrixAt(gp, C, tStep);
    int sz = C.giveSize();
    answer.resize(sz+1,sz+1);
    answer.zero();
    
    C.times(this->area);
    
    answer.addSubVectorCol(C, 1, sz + 1);
    answer.addSubVectorCol(C, sz + 1, 1);
}
  
    
    
    
    
    
}












