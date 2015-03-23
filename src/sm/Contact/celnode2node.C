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

#include "Contact/celnode2node.h"
#include "floatmatrix.h"
#include "masterdof.h"
#include "unknownnumberingscheme.h"
#include "gaussintegrationrule.h"

namespace oofem {
  
  
  
Node2NodeContact :: Node2NodeContact(DofManager *master, DofManager *slave) : ContactElement()
{   
    this->masterNode = master;
    this->slaveNode = slave;
    this->area = 1.0;   // should be optional parameter
    this->epsN = 1.0e6; // penalty - should be given by 'contactmaterial'
};   
  
int
Node2NodeContact :: instanciateYourself(DataReader *dr)
{
    // compute normal as direction vector from master node to slave node
    FloatArray xs, xm, _normal;
    xs = *this->slaveNode->giveCoordinates();
    xm = *this->masterNode->giveCoordinates();
    
	_normal = xs-xm;
    double norm = _normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", 
          masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() )
    } else {
        normal = _normal*(1.0/norm);
    }
    return 1;
}


void
Node2NodeContact :: computeGap(FloatArray &answer, TimeStep *tStep)
{

    FloatArray xs, xm, uS, uM;
    xs = *this->slaveNode->giveCoordinates();
    xm = *this->masterNode->giveCoordinates();
    this->slaveNode->giveUnknownVector(uS, {D_u, D_v, D_w}, VM_Total, tStep, true);
    this->masterNode->giveUnknownVector(uM, {D_u, D_v, D_w}, VM_Total, tStep, true);
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs-xm;
    
    FloatArray normal = this->giveNormal();
    answer = {dx.dotProduct(normal), 0.0, 0.0};
    
    //printf("normal gap = %e \n", answer.at(1));
    if ( answer.at(1) < 0.0 ) {
        //printf("normal gap = %e \n", answer.at(1));
        this->inContact = true; // store in gp?
    } else {
        this->inContact = false;
    }
    
}


void
Node2NodeContact :: computeCmatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *TimeStep)
{
  
    // The normal is not updated for node2node which is for small deformations only
    // C = {n -n}
    FloatArray normal = this->giveNormal();
    answer = {  normal.at(1),  normal.at(2),  normal.at(3),
               -normal.at(1), -normal.at(2), -normal.at(3) };
    
}



void
Node2NodeContact :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{
    // should be replaced with a call to constitutive model
    // gap should be in a local system
    if ( gap.at(1) < 0.0 ) {
        t = this->epsN * gap;
    } else {
        t = {0.0, 0.0, 0.0};
    }
  
}





void
Node2NodeContact :: computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    answer.clear();
    FloatArray gap, C;
      
    this->computeGap(gap, tStep);
    if ( gap.at(1) < 0.0 ) {
        GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
        FloatArray t;
        this->computeContactTractionAt(gp, t ,gap, tStep);
        
        this->computeCmatrixAt(gp, C, tStep);
        
        // compute load vector
        // fc = C^T * traction * A, Area - should be optional par
        answer = t.at(1) * this->area * C;
          
    }
  
}

  
void
Node2NodeContact :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
   // Need to set up an integration rule 
    GaussPoint *gp = NULL;
    FloatArray gap;
      
    this->computeGap(gap, tStep);

    FloatArray C;
    this->computeCmatrixAt(gp, C, tStep);
    answer.beDyadicProductOf(C,C);
    // this is the interface stiffness and should be obtained from that model
    answer.times( this->epsN * this->area );
    answer.negated();

    if( gap.at(1) > 0.0 ) {
        answer.zero();
    }
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


void
Node2NodeContact :: setupIntegrationPoints()
{
    // Sets up the integration rule array which contains all the necessary integration points
    if ( this->integrationRule == NULL ) {
        //TODO sets a null pointer for the element in the iRule 
        this->integrationRule = new GaussIntegrationRule(1, NULL) ;
        this->integrationRule->SetUpPoint(_Unknown);
    }
    
  
}













// node 2 node Lagrange


Node2NodeContactL :: Node2NodeContactL(DofManager *master, DofManager *slave) : Node2NodeContact(master, slave)
{   
    this->masterNode = master;
    this->slaveNode = slave;
    this->area = 1.0;
};   



void
Node2NodeContactL :: giveLocationArray(IntArray &answer, const UnknownNumberingScheme &s)
{
  
    Node2NodeContact :: giveLocationArray(answer, s);
    
    // Add one lagrange dof
    if ( this->masterNode->hasDofID( (DofIDItem)this->giveDofIdArray().at(1) ) ) { 
        Dof *dof= this->masterNode->giveDofWithID( (DofIDItem)this->giveDofIdArray().at(1) );
        answer.followedBy( s.giveDofEquationNumber(dof) );
    }
    
}    



void
Node2NodeContactL :: computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
  
    //Loop through all the master objects and let them do their thing
    FloatArray gap, C;
    this->computeGap(gap, tStep);
    
    GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
    FloatArray t;
    this->computeContactTractionAt(gp, t ,gap, tStep);
    this->computeCmatrixAt(gp, C, tStep);
    
    // compute load vector
    // for Lagrange: fc = traction * C^T * A (traction = lambda)
    FloatArray temp = t.at(1) *this->area * C;
    
    answer.resize( C.giveSize() + 1);
    answer.zero();
    if( gap.at(1) < 0.0 ) {
        answer.addSubVector(temp,1);
        answer.at( C.giveSize() + 1 ) = -gap.at(1);
    }
}
    


void
Node2NodeContactL :: computeContactTangent(FloatMatrix &answer, CharType type, TimeStep *tStep)
{
    answer.resize(7,7);
    answer.zero();
    
    FloatArray gap;
    this->computeGap(gap, tStep);
    
    if( gap.at(1) < 0.0 ) {
      
        GaussPoint *gp = this->integrationRule->getIntegrationPoint(0);
        
        FloatArray C;
        this->computeCmatrixAt(gp, C, tStep);
        int sz = C.giveSize();
        C.times(this->area);
        
        answer.addSubVectorCol(C, 1, sz + 1);
        answer.addSubVectorRow(C, sz + 1, 1);
    }
    
    //TODO need to add a small number for the solver
    for ( int i = 1; i <= 7; i++ ) {
        answer.at(i,i) += 1.0e-8;
    }

}
  

  
void
Node2NodeContactL :: computeContactTractionAt(GaussPoint *gp, FloatArray &t, FloatArray &gap, TimeStep *tStep)
{
    // should be replaced with a call to constitutive model
    // gap should be in a local system
    if ( gap.at(1) < 0.0 ) {
        
        Dof *dof = masterNode->giveDofWithID( this->giveDofIdArray().at(1) );
        double lambda = dof->giveUnknown(VM_Total, tStep);
        t = {lambda, 0.0, 0.0};
        //printf("lambda %e \n\n", lambda);
    } else {
        t = {0.0, 0.0, 0.0};
    }
  
} 
  
    
void
Node2NodeContactL :: giveDofManagersToAppendTo(IntArray &answer)
{
    answer = {this->masterNode->giveNumber()};
}
    
    
    
}












