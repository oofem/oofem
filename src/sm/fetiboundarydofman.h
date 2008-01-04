/* $Header: /home/cvs/bp/oofem/sm/src/fetiboundarydofman.h,v 1.2 2003/04/06 14:08:30 bp Exp $ */
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
#ifdef __PARALLEL_MODE
#ifndef fetiboundarydofmanager_h

#include "intarray.h"
#ifndef __MAKEDEPEND
#include <map>
#endif
using namespace std;

/**
 Represent the abstraction for dof manager. This dof manager is on partition boudary and 
 influences to master equations. It is typically shared by two or more partitions.
 It keeps its associated global number, number of partitions sharing it, and corresponding
 number of DOFs (only those with associate equation are taken into account).
 If n partitions share the dof manager, then on master level n-1 different code numbers for each DOF
 must be maintained. The compatibility is then enforced using lagrange multipliers.
 It is necessary to select one partition (reference one), to which other partitions DOFS are 
 "linked" through lagrange multipliers. Such partition contributes to all sharing partitions DOFs.
 Such partition contribute to all allocated cod numbers. The partition with lowest rank is 
 selected as refence one. 
 */
class FETIBoundaryDofManager {
protected:
 /// Associated global number of dofManager
 int globalNumber;
 /// Total number of partitions sharing receiver
 int numberOfPartitions;
 // number of nonprescribed dofs, i.e, those, for which equation is necessary
 int ndofs;
 // reference partition is partition to which other partitions sharing the dof manager
 // are linked using lagrange multipliers. We use the partition, which has its number 
 // the lowest from all sharing partitions.
 int referencePartition;
 /// List of partitions sharing dof manager
 IntArray partitions;
 /**
  Contains code numbers for each linked partition for each DOF (ndofs*(numberOfPartitions-1) DOFs).
  The reference partition contributes to all code numbers.
  The code numbers for n-th linked partition (its number is determined by index 
  of corresponding entry in partitions array) are stored at ((n-1)*ndofs+1, ..., n*ndofs) positions
  in codeNumbers array.
  */
 IntArray codeNumbers;
public:
 FETIBoundaryDofManager() ;
 FETIBoundaryDofManager(int, int, int) ;
 FETIBoundaryDofManager(const FETIBoundaryDofManager&);
 
 /// Returns number of partitions sharing receiver
 int giveNumberOfSharedPartitions() {return numberOfPartitions;}
 /// Returns number of DOFs (with associated equation) of receiver
 int giveNumberOfDofs() {return ndofs;}
 /// Returns correcponding global number of receiver
 int giveGlobalNumber() {return globalNumber;}
 /// Returns reference partition number of receiver
 int giveReferencePratition() {return referencePartition;}
 /// Returns number of i-th shared partition of receiver
 int giveSharedPartition (int i) {return partitions.at(i);}
 /**
  Returns code number corresponding to partition number partition_num and to dof_num-th DOF 
  @param partition_num partition number for which code number is required
  @param dof_num the specifies the particular DOF
  @return value of correspong code number, zero if such partition does not share 
  the receiver or if code number for reference partition is requested.
  */
 int giveCodeNumber (int partition_num, int dof_num);
 /**
  Returns code numbers for all DOFs associated with shared partition.
  @param rank partition number
  @param locArray the location array of size ndof
  @return nonzero if o.k, zero if no such partition shared or if 
  code numbers for reference partition required.
  */
 int giveCompleteLocationArray (int rank, IntArray& locationArray);
 /**
  Adds partition to list of partitions, sharing this dof manager.
  The referencePartition is updated if necessary.
  @param partitionNumber new partition number (0..size-1)
  */
 void addPartition (int partitionNumber) ;
 /**
  Associates the eqution numbers to particular DOFs. 
  @param equationCounter current equation counter, updated
  @return new value of equationCounter
  */
 int  setCodeNumbers (int &equationCounter) ;
} ;

#define fetiboundarydofmanager_h
#endif
#endif
