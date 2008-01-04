
//
//  termitovo - zaloha
//

#include "slavedof.h"
#include "node.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <ctype.h>
#endif



/**
 */
SlaveDof :: SlaveDof (int n, DofManager* aNode, DofID id) : Dof (n, aNode, id), masterContribution()
{
  countOfPrimaryMasterDofs = -1;
}


/**
 */
void
SlaveDof :: initialize (int cntOfMstrDfMngr, Node** mstrNode, const IntArray* mstrDofID, const FloatArray* mstrContribution)
{
  int i, id;
  bool idSame = false;
  
  
  if (mstrDofID==NULL)  idSame = true;
  else
    if (mstrDofID->giveSize() < cntOfMstrDfMngr)
      _error3 ("initialize: mstrDofID.giveSize %d != cntOfMstrDfMngr %d", mstrDofID->giveSize(), cntOfMstrDfMngr);
  
  if (mstrContribution->giveSize() < cntOfMstrDfMngr)
    _error3 ("initialize: mstrContribution.giveSize %d != cntOfMstrDfMngr %d", mstrContribution->giveSize(), cntOfMstrDfMngr);
  
  
  countOfMasterDofs  = cntOfMstrDfMngr;
  masterContribution = *mstrContribution;
  
  masterDof = new Dof* [countOfMasterDofs];
  
  for (i=1; i<=countOfMasterDofs; i++){
    if (idSame) id = this->dofID;
    else        id = mstrDofID->at(i);
    
    masterDof[i-1] = mstrNode[i-1]->giveDofWithID (id);
  }
}

int
SlaveDof :: giveNumberOfPrimaryMasterDofs (void)
{
  if (countOfPrimaryMasterDofs > 0)
    return countOfPrimaryMasterDofs;
  else
    if (countOfPrimaryMasterDofs==0)
      _error2 ("giveNumberOfPrimaryDofs: slaveDof number %ld is own master", this->giveNumber());
  
  countOfPrimaryMasterDofs = 0;
  
  long i, c=0;
  for (i=1; i<=countOfMasterDofs; i++)
    if (! masterDof[i-1]->isPrimaryDof())
      c += masterDof[i-1]->giveNumberOfPrimaryMasterDofs();
    else
      c += 1;
  
  return countOfPrimaryMasterDofs = c;
}

void
SlaveDof :: giveUnknowns (FloatArray& masterUnknowns, EquationID type, ValueModeType mode, TimeStep* stepN)
{
  int i, k;
  FloatArray mstrUnknwns;
  
  masterUnknowns.resize (this->giveNumberOfPrimaryMasterDofs());
  
  for (k=1, i=1; i<=countOfMasterDofs; i++) {
    if (! masterDof[i-1]->isPrimaryDof()) {
      masterDof[i-1]->giveUnknowns (mstrUnknwns, type, mode, stepN);
      masterUnknowns.copySubVector (mstrUnknwns, k);
      k += mstrUnknwns.giveSize();
    }
    else
      masterUnknowns.at(k++) = masterDof[i-1]->giveUnknown (type, mode, stepN);
  }
}

void
SlaveDof :: giveUnknowns (FloatArray& masterUnknowns, PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{
  int i, k;
  FloatArray mstrUnknwns;
  
  masterUnknowns.resize (this->giveNumberOfPrimaryMasterDofs());
  
  for (k=1, i=1; i<=countOfMasterDofs; i++) {
    if (! masterDof[i-1]->isPrimaryDof()) {
      masterDof[i-1]->giveUnknowns (mstrUnknwns, field, mode, stepN);
      masterUnknowns.copySubVector (mstrUnknwns, k);
      k += mstrUnknwns.giveSize();
    }
    else
      masterUnknowns.at(k++) = masterDof[i-1]->giveUnknown (field, mode, stepN);
  }
}

void
SlaveDof :: giveBcValues (FloatArray& masterBcValues, ValueModeType mode, TimeStep* stepN)
{
  int i, k;
  FloatArray mstrBcVlus;
  Dof *dofJ;
  
  masterBcValues.resize (this->giveNumberOfPrimaryMasterDofs());
  
  for (k=1, i=1; i<=countOfMasterDofs; i++) {
    if (! masterDof[i-1]->isPrimaryDof()) {
      masterDof[i-1]->giveBcValues (mstrBcVlus, mode, stepN);
      masterBcValues.copySubVector (mstrBcVlus, k);
      k += mstrBcVlus.giveSize();
    }
    else {
      dofJ = masterDof[i-1];
      if (dofJ-> hasBc(stepN))
	masterBcValues.at(k++) = dofJ-> giveBcValue(mode,stepN);
      else
	masterBcValues.at(k++) = 0.0;
    }
  }
}

void
SlaveDof :: computeDofTransformation (FloatArray& masterContribs)
{
  int i, k;
  FloatArray mstrContrs;
  
  masterContribs.resize (this->giveNumberOfPrimaryMasterDofs());
  
  for (k=1, i=1; i<=countOfMasterDofs; i++) {
    if (! masterDof[i-1]->isPrimaryDof()) {
      masterDof[i-1]->computeDofTransformation (mstrContrs);
      mstrContrs.times (masterContribution.at(i));
      masterContribs.copySubVector (mstrContrs, k);
      k += mstrContrs.giveSize();
    }
    else {
      masterContribs.at(k++) = masterContribution.at(i);
    }
  }
}

void
SlaveDof :: giveEquationNumbers (IntArray& masterEqNumbers)
{
  int i, k;
  IntArray mstrEqNmbrs;
  
  masterEqNumbers.resize (this->giveNumberOfPrimaryMasterDofs());
  
  for (k=1, i=1; i<=countOfMasterDofs; i++) {
    if (! masterDof[i-1]->isPrimaryDof()) {
      masterDof[i-1]->giveEquationNumbers (mstrEqNmbrs);
      masterEqNumbers.copySubVector (mstrEqNmbrs, k);
      k += mstrEqNmbrs.giveSize();
    }
    else {
      masterEqNumbers.at(k++) = masterDof[i-1]->giveEquationNumber();
    }
  }
}

void
SlaveDof :: givePrescribedEquationNumbers (IntArray& masterEqNumbers)
{
  int i, k;
  IntArray mstrEqNmbrs;
  
  masterEqNumbers.resize (this->giveNumberOfPrimaryMasterDofs());
  
  for (k=1, i=1; i<=countOfMasterDofs; i++) {
    if (! masterDof[i-1]->isPrimaryDof()) {
      masterDof[i-1]->givePrescribedEquationNumbers (mstrEqNmbrs);
      masterEqNumbers.copySubVector (mstrEqNmbrs, k);
      k += mstrEqNmbrs.giveSize();
    }
    else {
      masterEqNumbers.at(k++) = masterDof[i-1]->givePrescribedEquationNumber();
    }
  }
}


/**
   Returns the value of the unknown associated with the receiver at given time step.
   Slave simply asks vector of corresponding master dofs and own transformation
   vector and returns result as dot product of these vectors. Standart element
   services have to transform global unknown vector transform into their local c.s
   before using it (when computing strain vector by \eps=Br, for example, 
   where B is element geometrical matrix). This transformation should contain also
   nodal to global coordinate system transformation. So, this specialized 
   standard method for unknown query returns the corresponding master DOF value.
   @see MasterDof::giveUnknown function
*/
double SlaveDof :: giveUnknown (EquationID type, ValueModeType mode, TimeStep* stepN)
{
  FloatArray masterUnknowns, t;
  
  giveUnknowns (masterUnknowns, type, mode, stepN);
  computeDofTransformation (t);
  
  return dotProduct (masterUnknowns, t, t.giveSize());
}
double SlaveDof :: giveUnknown (PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{
  FloatArray masterUnknowns, t;
  
  giveUnknowns (masterUnknowns, field, mode, stepN);
  computeDofTransformation (t);
  
  return dotProduct (masterUnknowns, t, t.giveSize());
}
