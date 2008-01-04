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


#include "hangingdof.h"
#include "hangingnode.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"

#include "flotarry.h"
#include "dictionr.h"

//#include "string.h"
#include "debug.h"
#include "cltypes.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <ctype.h>
#endif


/**
   Constructor. Creates slave dof vith number i, belonging to aNode dof manager.
   Slave will be linked to master dof with id type belonging to dofManager with 
   number given in master variable.
   @param i dof number
   @param aNode receiver will belong to aNode dof manager
   @param master number of dofManager which contain master dof
     @param id DofID of master dof (and slave too).
*/
HangingDof::HangingDof (int i, DofManager* aNode, DofID id) : Dof (i, aNode, id)
{
}


/**
   Returns the value of the unknown associated with the receiver
   at given time step. Slave simply asks corresponding master dof and 
   returns master result. Standart element services have to transform 
   global unknown vector transform into their local c.s before using it 
   (when computing strain vector by \eps=Br, for example, 
   where B is element geometrical matrix). This transformation should contain also
   nodal to global coordinate system transformation. So, this specialized 
   standard method for unknown query returns the corresponding master DOF value.
   @see MasterDof::giveUnknown function
*/
double HangingDof::giveUnknown (EquationID type, ValueModeType mode, TimeStep* stepN)
  // The key method of class Dof. Returns the value of the unknown 'u'
  // (e.g., the displacement) of the receiver, at stepN. This value may,
  // or may not, be already available. It may depend on a boundary (if it
  // is not a predicted unknown) or initial condition. stepN is not the
  // current time step n, it is assumed to be the previous one (n-1).
{
  int i,index,cMDM;
  FloatArray masterUnknowns;
  Node* master;
  
  cMDM = ((HangingNode*)dofManager)->giveCountOfMasterDofMngr();
  masterUnknowns.resize (cMDM);
  
  for (i=1;i<=cMDM;i++) {
    master = ((HangingNode*)dofManager)->giveMasterDofMngr (i);
    if ( (index=master->findDofWithDofId(this->giveDofID())) == 0 ) {
      _error2("hangingdof on node %d : uncompatible dof requested",master->giveNumber());
    }
    masterUnknowns.at(i) = master->giveDof(index)->giveUnknown (type, mode, stepN);
  }
  
  return giveMasterUnknowns (masterUnknowns);
}


/**
 */
double HangingDof::giveUnknown (PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{
  int i,index,cMDM;
  FloatArray masterUnknowns;
  Node* master;
  
  cMDM = ((HangingNode*)dofManager)->giveCountOfMasterDofMngr();
  masterUnknowns.resize (cMDM);
  
  for (i=1;i<=cMDM;i++) {
    master = ((HangingNode*)dofManager)->giveMasterDofMngr (i);
    if ( (index=master->findDofWithDofId(this->giveDofID())) == 0 ) {
      _error2 ("hangingdof on node %d : uncompatible dof requested",master->giveNumber());
    }
    masterUnknowns.at(i) = master->giveDof(index)->giveUnknown (field, mode, stepN);
  }
  
  return giveMasterUnknowns (masterUnknowns);
}


/**
 */
double HangingDof::giveMasterUnknowns (FloatArray &masterUnknowns)
{
  long i,cMDM;
  double answer;
  FloatMatrix t;
  IntArray dofIDArry(1);
  
  dofIDArry.at(1) = this->giveDofID();
  ((HangingNode*)dofManager)->computeDofTransformation (t,&dofIDArry,_toGlobalCS);
  
  cMDM = masterUnknowns.giveSize();
  answer = 0.0;
  for (i=1;i<=cMDM;i++)
    answer += t.at(1,i) * masterUnknowns.at(i);
  
  return answer;
}
