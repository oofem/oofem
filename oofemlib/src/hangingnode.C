/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/
/*
  Contribution by Ladislav Svoboda
*/

#include "hangingnode.h"
#include "hangingdof.h"
#include "nodload.h"
#include "timestep.h"
#include "masterdof.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "debug.h"
#include "verbose.h"
#include "mathfem.h"

#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

#include "fei1dlin.h"
#include "fei2dtrlin.h"
#include "fei2dquadlin.h"
#include "fei3dhexalin.h"


//IMPLEMENT_DYNAMIC_TYPE(HangingNode, Node)


HangingNode :: HangingNode (int n, Domain* aDomain) : Node (n,aDomain), masterDofMngr()
  // Constructor. Creates a node with number n, belonging to aDomain.
{
  locoords.resize (3);
}

HangingNode :: ~HangingNode()
  // Destructor.
{
}


IRResultType
HangingNode :: initializeFrom (InputRecord* ir)
  // Gets from the source line from the data file all the data of the receiver.
{
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro
  
  int i;
  IntArray dofIDArry;
  
  Node::initializeFrom (ir);
  
  IR_GIVE_FIELD (ir, type, IFT_HangingNode_type, "type"); // Macro
  
  IR_GIVE_FIELD (ir, masterDofMngr, IFT_HangingNode_masters, "masters"); // Macro
  countOfMasterDofMngr = masterDofMngr.giveSize();

  IR_GIVE_FIELD          (ir, locoords.at(1), IFT_HangingNode_ksi, "ksi"); // Macro
  IR_GIVE_OPTIONAL_FIELD (ir, locoords.at(2), IFT_HangingNode_eta, "eta"); // Macro
  IR_GIVE_OPTIONAL_FIELD (ir, locoords.at(3), IFT_HangingNode_dzeta, "dzeta"); // Macro
  
  //  compute_naturalcoord ();
  
  if (this->resolveDofIDArray (ir, dofIDArry) != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_Unknown, "", ir, result);
  
  dofArray = new Dof* [this->giveNumberOfDofs()] ;
  for (i=0 ; i<numberOfDofs ; i++)
    dofArray[i] = new HangingDof(i+1,this,(DofID) dofIDArry.at(i+1)) ;
  
  return IRRT_OK;
}

void
HangingNode :: compute_naturalcoord ()
{
  int i,first,x,y,cMDM;
  double locksi,loceta,det;
  FloatArray **mnc;
  
  first = 1;
  cMDM = countOfMasterDofMngr;
  mnc = new FloatArray* [cMDM];

  for (i=1;i<=cMDM;i++)
    mnc[i-1] = domain->giveNode (masterDofMngr.at(i))->giveCoordinates ();
  
  switch (type) {
  case 211:
    for (i=1;i<=3;i++) {
      det = mnc[1]->at(i) - mnc[0]->at(i);
      if (det>0.000000001){
        locksi = (2*coordinates.at(i) - mnc[1]->at(i) - mnc[0]->at(i)) / det;
        if ( !first && (locksi-locoords.at(1))>0.00000001 ) {
          _error ("compute_naturalcoord: internal err: locksi != ksi");
        } else {
          locoords.at(1) = locksi;
          first = 0;
        }
      }
    }
    break;
    
  case 312:
    for (i=1;i<=3;i++) {
      x = i;
      y = i%3 + 1;
      det = (mnc[1]->at(x) - mnc[0]->at(x)) * (mnc[2]->at(y) - mnc[0]->at(y)) - (mnc[2]->at(x) - mnc[0]->at(x)) * (mnc[1]->at(y) - mnc[0]->at(y));
      if (det>0.000000001){
        locksi = ( (mnc[1]->at(x) * mnc[2]->at(y) - mnc[2]->at(x) * mnc[1]->at(y) )
                   +(                mnc[1]->at(y) -                 mnc[2]->at(y) ) * coordinates.at(x)
                   +(mnc[2]->at(x)                 - mnc[1]->at(x)                 ) * coordinates.at(y) ) / det;
        loceta = ( (mnc[2]->at(x) * mnc[0]->at(y) - mnc[0]->at(x) * mnc[2]->at(y) )
                   +(                mnc[2]->at(y) -                 mnc[0]->at(y) ) * coordinates.at(x)
                   +(mnc[0]->at(x)                 - mnc[2]->at(x)                 ) * coordinates.at(y) ) / det;
        if ( !first && (locksi-locoords.at(1))>0.00000001 ) {
          _error ("compute_naturalcoord: internal err: locksi != ksi");
        } else {
          locoords.at(1) = locksi;
          locoords.at(2) = loceta;
          first = 0;
        }
      }
    }
    break;
 
  case 412:
    for (i=1;i<=3;i++) {
      // toto uz je nekde v oofemu udelano asi Danem
      x = i;
      y = i%3 + 1;
      det = ( (mnc[1]->at(x)-mnc[0]->at(x))*(mnc[2]->at(y)-mnc[0]->at(y)) - (mnc[2]->at(x)-mnc[0]->at(x))*(mnc[1]->at(y)-mnc[0]->at(y))
              +(mnc[2]->at(x)-mnc[0]->at(x))*(mnc[3]->at(y)-mnc[0]->at(y)) - (mnc[3]->at(x)-mnc[0]->at(x))*(mnc[2]->at(y)-mnc[0]->at(y)) ) / 2.0;
      if (det>0.000000001){
        int num;
        double xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4;
        double a,b,c,x1,x2,x3;
        
        xx1 = mnc[0]->at(x)+mnc[1]->at(x)+mnc[2]->at(x)+mnc[3]->at(x);
        yy1 = mnc[0]->at(y)+mnc[1]->at(y)+mnc[2]->at(y)+mnc[3]->at(y);
        xx2 = mnc[0]->at(x)-mnc[1]->at(x)-mnc[2]->at(x)+mnc[3]->at(x);
        yy2 = mnc[0]->at(y)-mnc[1]->at(y)-mnc[2]->at(y)+mnc[3]->at(y);
        xx3 = mnc[0]->at(x)+mnc[1]->at(x)-mnc[2]->at(x)-mnc[3]->at(x);
        yy3 = mnc[0]->at(y)+mnc[1]->at(y)-mnc[2]->at(y)-mnc[3]->at(y);
        xx4 = mnc[0]->at(x)-mnc[1]->at(x)+mnc[2]->at(x)-mnc[3]->at(x);
        yy4 = mnc[0]->at(y)-mnc[1]->at(y)+mnc[2]->at(y)-mnc[3]->at(y);
        
        a = yy3*xx4 - yy4*xx3;
        b = yy1*xx4 - yy2*xx3 + yy3*xx2 - yy4*xx1 + yy4*4.*coordinates.at(x) - xx4*4.*coordinates.at(y);
        c = yy1*xx2 - yy2*xx1                     + yy2*4.*coordinates.at(x) - xx2*4.*coordinates.at(y);
        
        cubic (0.,a,b,c,&x1,&x2,&x3,&num);
        if ( (num == 0) || (num == 1 && (fabs (x1) > 1.00001)) ) {
          _error2 ("compute_naturalcoord: internal err on node %d", this->giveNumber());
        }
        if (num == 2 && (fabs (x1) > 1.00001))
          if (fabs(x2) > 1.00001) {
            _error2 ("compute_naturalcoord: internal err on node %d", this->giveNumber());
          } else {x1 = x2;}
        
        loceta = x1;
        locksi = (4.0*coordinates.at(x)-xx1-x1*xx3) / (xx2 + x1*xx4);
        
        if ( !first && (locksi-locoords.at(1))>0.00000001 ) {
          _error ("compute_naturalcoord: internal err: locksi != ksi");
        } else {
          locoords.at(1) = locksi;
          locoords.at(2) = loceta;
          first = 0;
        }
      }
    }
    break;
    
  case 322:
  case 622:
  case 822:
    _error ("compute_naturalcoord: unsupported element type");
    break;
    
  default:
    _error ("compute_naturalcoord: unsupported element type");
    break;
    
  }
}


int
HangingNode :: checkConsistency ()
{
  /*
    Checks internal data consistency in node. 
    Current implementation checks (when receiver has slave dofs) if receiver has the same 
    coordinate system as master dofManager of slave dof.
  */
  
  int result = 1;
  int i,j;
  
  result = result && Node::checkConsistency();
  
  // check if master is RigidArmNode or HangingNode - not supported
  for (i=1;i<=countOfMasterDofMngr;i++)
    if (giveMasterDofMngr(i)->giveClassID() == HangingNodeClass ||
        giveMasterDofMngr(i)->giveClassID() == RigidArmNodeClass) {
      _warning ("checkConsistency: chaining of HangingNodes is not allowed");
      result = 0;
    }

  for (i=1; i<= numberOfDofs; i++)
    for (j=1;j<=countOfMasterDofMngr;j++) {
      if (this->giveDof (i)->giveDofID() != giveMasterDofMngr(j)->giveDof (i)->giveDofID()) {
        _warning ("checkConsistency: dofID mismatch on HangingNode and master");
        return 0;
      }
    }
  
  return result;
}

/*
void 
HangingNode :: giveUnknownVector (FloatArray& answer, IntArray& dofMask, UnknownType type, UnknownTypeMode mode, TimeStep* stepN)
{
  int i,j,size;
  IntArray dofArray;
  FloatArray masterUnknowns(countOfMasterDofMngr);
  
  size = dofMask.giveSize();
  answer.resize (size * countOfMasterDofMngr);
  this-> giveDofArray(dofMask, dofArray);
  
  for (i=1;i<=size;i++) {
    ((HangingDof*)(this->giveDof(dofArray.at(i))))->giveMasterUnknowns (masterUnknowns,type,mode,stepN) ;
    for (j=1;j<=countOfMasterDofMngr;j++)
      answer.at(i+size*(j-1)) = masterUnknowns.at(j);
  }
}
*/

void 
HangingNode :: giveUnknownVector (FloatArray& answer, const IntArray& dofMask, EquationID type, ValueModeType mode, TimeStep* stepN)
{
  int i,j,k,size;
  FloatArray masterUnknowns;
  
  size = dofMask.giveSize();
  masterUnknowns.resize (size);
  answer.resize (size * countOfMasterDofMngr);
  
  for (i=1,k=1;i<=countOfMasterDofMngr;i++) {
    (giveMasterDofMngr(i))->giveUnknownVector(masterUnknowns,dofMask,type,mode,stepN);
    for (j=1;j<=size;j++)
      answer.at(k++) = masterUnknowns.at(j);
  }
}


void 
HangingNode :: givePrescribedUnknownVector (FloatArray& answer, const IntArray& dofMask, ValueModeType mode, TimeStep* stepN)
{
  int i,j,k,size;
  FloatArray masterPrescribedUnknowns;
  
  size = dofMask.giveSize();
  masterPrescribedUnknowns.resize (size);
  answer.resize (size * countOfMasterDofMngr);
  
  for (i=1,k=1;i<=countOfMasterDofMngr;i++) {
    (giveMasterDofMngr(i))->givePrescribedUnknownVector(masterPrescribedUnknowns,dofMask,mode,stepN);
    for (j=1;j<=size;j++)
      answer.at(k++) = masterPrescribedUnknowns.at(j);
  }
}


void
HangingNode :: computeDofTransformation (FloatMatrix& answer, const IntArray* dofIDArry, DofManTrasfType mode)
{
  // computes trasformation matrix of receiver.
  // transformation should include trasformation from global cs to nodal cs,
  // as well as further necessary transformations (for example in case 
  // rigid arms this must include transformation to master dofs).
  int i,j,size;
  FloatArray N;
  
  if (dofIDArry) size = dofIDArry->giveSize(); else size = numberOfDofs;
  
  
  if (mode != _toGlobalCS) { _error ("computeDofTransformation: unknown mode");  return; }
  
  if (dofIDArry)
    for (i=1;i<=size;i++)
      if (this->findDofWithDofId ((DofID) dofIDArry->at(i)) == 0) {
        _error2 ("computeDofTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
        return;
      }
  
  N.resize (countOfMasterDofMngr);
  
  switch (type){
  case 211:{ // linear truss
    N.at(1) = 0.5 * (1. - locoords.at(1));
    N.at(2) = 0.5 * (1. + locoords.at(1));
    break;
  }
  case 312:{ // linear triangle
    N.at(1) = locoords.at(1);
    N.at(2) = locoords.at(2);
    N.at(3) = 1.0 - locoords.at(1) - locoords.at(2);
    break;
  }
  case 412:{ // linear rectangle
    FEI2dQuadLin(0,0).evalN (N,locoords,0.0);
    break;
  }
  case 813:{ // linear hexahedron
    FEI3dHexaLin().evalN (N,locoords,0.0);
    break;
  }
  case 321:{ // quadratic truss
    N.at(1) = 0.5 * locoords.at(1) * (locoords.at(1) - 1.);
    N.at(2) = 0.5 * locoords.at(1) * (locoords.at(1) + 1.);
    N.at(3) = 1. - locoords.at(1) * locoords.at(1);
    break;
  }
  case 622:{ // quadratic triangle
    locoords.at(3) = 1. - locoords.at(1) - locoords.at(2);
    N.at(1) = (2. * locoords.at(1) - 1.) * locoords.at(1);   
    N.at(2) = (2. * locoords.at(2) - 1.) * locoords.at(2);
    N.at(3) = (2. * locoords.at(3) - 1.) * locoords.at(3);
    N.at(4) =  4. * locoords.at(1) * locoords.at(2);
    N.at(5) =  4. *                  locoords.at(2) * locoords.at(3);
    N.at(6) =  4. *                                   locoords.at(3) * locoords.at(1);
    break;
  }
  case 822:{ // quadratic rectangle
    N.at(1) = (1. + locoords.at(1)) * (1. + locoords.at(2)) * 0.25 * ( locoords.at(1) + locoords.at(2) - 1.) ;
    N.at(2) = (1. - locoords.at(1)) * (1. + locoords.at(2)) * 0.25 * (-locoords.at(1) + locoords.at(2) - 1.) ;
    N.at(3) = (1. - locoords.at(1)) * (1. - locoords.at(2)) * 0.25 * (-locoords.at(1) - locoords.at(2) - 1.) ;
    N.at(4) = (1. + locoords.at(1)) * (1. - locoords.at(2)) * 0.25 * ( locoords.at(1) - locoords.at(2) - 1.) ;
    N.at(5) = 0.5 * (1.-locoords.at(1)*locoords.at(1)) * (1.+locoords.at(2));
    N.at(6) = 0.5 * (1.-locoords.at(1)) * (1.-locoords.at(2)*locoords.at(2));
    N.at(7) = 0.5 * (1.-locoords.at(1)*locoords.at(1)) * (1.-locoords.at(2));
    N.at(8) = 0.5 * (1.+locoords.at(1)) * (1.-locoords.at(2)*locoords.at(2));
    break;
  }
  default:{
    _error ("computeDofTransformation: unknown type");
    break;
  }
  }
  
  answer.resize (size,size*countOfMasterDofMngr);
  answer.zero ();
  
  for (i=1;i<=size;i++)
    for (j=1;j<=countOfMasterDofMngr;j++)
      answer.at(i,i+size*(j-1)) = N.at(j);
  
  //answer.printYourself();   //termit
  
}


void
HangingNode :: computeLoadTransformation (FloatMatrix& answer, const IntArray* dofIDArry, DofManTrasfType mode)
{
  // computes trasformation matrix of receiver.
  // transformation should include trasformation from global cs to nodal cs,
  // as well as further necessary transformations (for example in case 
  // rigid arms this must include transformation to master dofs).
  FloatMatrix t;
  
  if (mode != _toNodalCS) { _error ("computeLoadTransformation: unsupported mode");  return; }
  
  this->computeDofTransformation (t, dofIDArry, _toGlobalCS);
  answer.beTranspositionOf (t);
}


void
HangingNode :: computeLoadVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)
  // Computes the vector of the nodal loads of the receiver.
{
  FloatMatrix masterTransf;
  //FloatArray localAnswer;
  
  // assemble answer of receiver for "local dofs"
  Node :: computeLoadVectorAt (answer,stepN, mode);

  // transform "local dofs" to master dofs
  if (answer.isNotEmpty ()) {
    // assemble transformation contributions from local dofs
    this->computeLoadTransformation (masterTransf, NULL, _toGlobalCS);
    answer.rotatedWith (masterTransf,'n');
  }
}


void 
HangingNode :: giveLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const
  // Returns the location array of the receiver. Creates this array if it
  // does not exist yet. The location array contains the equation number of
  // every  requested degree of freedom of the receiver.
  // In dofIDArry are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
  //
  // Rigid arm return always full location array
{
  int i,j,k,size;
  Node* master;
  // prevents some size problem when connecting different elements with 
  // different number of dofs
  size = dofIDArry.giveSize();
  locationArray.resize(size*countOfMasterDofMngr) ;
  
  for (i=1,k=1;i<=countOfMasterDofMngr;i++) {
    master = giveMasterDofMngr(i);
    for (j=1;j<=size;j++)
      locationArray.at(k++) = master->giveDof(dofIDArry.at(j))->giveEquationNumber();         // predpoklad - hangingnode je stale konzistentni !!!

      //if ((index=master->findDofWithDofId ((DofID) dofIDArry.at(j))) == 0) {            
      //char buff [40];                                                                   
      //sprintf(buff,"DofManager %d : uncompatible dof requested",master->giveNumber());  
      //_error (buff);
      //return;}
      //locationArray.at(k++) = master->giveDof(index)->giveEquationNumber();
  }
  return  ;
}


void 
HangingNode :: givePrescribedLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const
  // Returns the location array of prescribed equations of the receiver. Creates this array if it
  // does not exist yet. The location array contains the equation number of
  // every  requested degree of freedom of the receiver.
  // In dofIDArry are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
  //
  // Rigid arm return always full location array
{
  int i,j,k,size,indx ;
  Node* master;
  // prevents some size problem when connecting different elements with 
  // different number of dofs
  size = dofIDArry.giveSize();
  locationArray.resize(size*countOfMasterDofMngr) ;
  
  for (i=1,k=1;i<=countOfMasterDofMngr;i++) {
    master = giveMasterDofMngr(i);
    for (j=1;j<=size;j++) {
      if ((indx=master->findDofWithDofId ((DofID) dofIDArry.at(j))) == 0) {
	char buff [40];
	sprintf(buff,"DofManager %d : uncompatible dof requested",master->giveNumber());
	_error (buff);
	return;
      }
      locationArray.at(k++) = master->giveDof(indx)->givePrescribedEquationNumber();
    }
  }
  return  ;
}


void
HangingNode :: giveCompleteLocationArray (IntArray& locationArray) const
  // Returns the complete location array of the receiver.
  // including all available dofs
{
  int i,j,k;
  Node* master;
  // prevents some size problem when connecting different elements with 
  // different number of dofs
  locationArray.resize (numberOfDofs*countOfMasterDofMngr) ;
  
  for (i=1,k=1;i<=countOfMasterDofMngr;i++) {
    master = giveMasterDofMngr(i);
    for (j=1;j<=numberOfDofs;j++) {
      locationArray.at(k++) = master->giveDof(j)->giveEquationNumber();
    }
  }
  return ;
}


void
HangingNode :: giveCompletePrescribedLocationArray (IntArray& locationArray) const
  // Returns the complete location array of the receiver.
  // including all available dofs
{
  int i,j,k;
  Node* master;
  // prevents some size problem when connecting different elements with 
  // different number of dofs
  locationArray.resize (numberOfDofs*countOfMasterDofMngr) ;
  
  for (i=1,k=1;i<=countOfMasterDofMngr;i++) {
    master = giveMasterDofMngr(i);
    for (j=1;j<=numberOfDofs;j++) {
      locationArray.at(k++) = master->giveDof(j)->givePrescribedEquationNumber();
    }
  }
  return ;
}

