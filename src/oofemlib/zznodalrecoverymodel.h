/* $Header: /home/cvs/bp/oofem/oofemlib/src/zznodalrecoverymodel.h,v 1.2 2003/04/06 14:08:26 bp Exp $ */
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

//   **************************************************
//   *** CLASS Zienkiewicz-Zhu NODAL RECOVERY MODEL ***
//   **************************************************

#ifndef zznodalrecoverymodel_h 

#include "compiler.h"
#include "cltypes.h"
#include "nodalrecoverymodel.h"
#include "interface.h"

class GaussPoint;
class ZZNodalRecoveryModelInterface;

/**
 The nodal recovery model based on paper of Zienkiewicz and Zhu "A Simple Estimator and Adaptive
 Procedure for Practical Engineering Analysis". The recovery is based 
 on nodal averaging or projection process in which it is assumed that the stress sigma_star is 
 interpolated by the same function as the displacement.
*/
class ZZNodalRecoveryModel : public NodalRecoveryModel {

protected:

public:
 /// Constructor
 ZZNodalRecoveryModel (Domain *d);
 /// Destructor
 ~ZZNodalRecoveryModel();
 /** Recovers the nodal values for all regions of given Domain.
  @param d domain of interest
  @param type determines the type of internal variable to be recovered
  @param tStep time step
  */
 int recoverValues (InternalStateType type, TimeStep* tStep);
private:
 /**
  Initializes the region table indicating regions to skip.
  @param regionMap region tabl, the nonzero entry for region indicates region to skip due to 
  unsupported elements or incompatible value size
  @param regionValSize contains the record size for each region
  @param type determines the type of internal variable to be recovered
  */
 void initRegionMap (IntArray& regionMap, IntArray& regionValSize, InternalStateType type);
 
};

/**
 The element interface required by ZZNodalRecoveryModel.
*/
class ZZNodalRecoveryModelInterface : public Interface {

public:
 /// Constructor
 ZZNodalRecoveryModelInterface () {}
 
 /**
  Computes the element contribution to \f$\int_\Omega N^T\alpha\;d\Omega\f$,
  where \f$\alpha\f$ is quantity to be recovered (for example stress or strain vector).
  The size of answer should be recordSize*numberofDofManagers.
  @param answer contains the result
  @param type determines the type of internal variable to be recovered
  @param tStep time step
  */
 // virtual void ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, TimeStep* tStep);
 virtual void ZZNodalRecoveryMI_computeNValProduct (FloatMatrix& answer, InternalStateType type, TimeStep* tStep);
 /**
  Computes the element contribution to \f$\int_\Omega N^TN\;d\Omega\f$ term.
  The size of answer should be [recordSize*numberofDofManagers].
  @param answer contain diagonalized result
  @param type determines the type of internal variable to be recovered
   */
 virtual void ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type);
 /**
  Returns the size of DofManger record required to hold recovered values for given mode.
  @param type determines the type of internal variable to be recovered
  @return size of DofManger record required to hold recovered values
  */
 virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type) = 0;
 /**
  Returns the corresponding element to interface 
  */
 virtual Element* ZZNodalRecoveryMI_giveElement () = 0;
  /**
  Evaluates N matrix (interpolation estimated value matrix).
  */
 virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (FloatMatrix& answer, GaussPoint* aGaussPoint, 
                                  InternalStateType type) = 0;
};


#define zznodalrecoverymodel_h
#endif






