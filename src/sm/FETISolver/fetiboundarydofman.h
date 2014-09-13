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

#ifndef fetiboundarydofman_h
#define fetiboundarydofman_h

#include "intarray.h"

namespace oofem {
/**
 * Represent the abstraction for DOF manager. This DOF manager is on partition boundary and
 * influences to master equations. It is typically shared by two or more partitions.
 * It keeps its associated global number, number of partitions sharing it, and corresponding
 * number of DOFs (only those with associate equation are taken into account).
 * If n partitions share the DOF manager, then on master level n-1 different code numbers for each DOF
 * must be maintained. The compatibility is then enforced using Lagrange multipliers.
 * It is necessary to select one partition (reference one), to which other partitions DOFs are
 * "linked" through Lagrange multipliers. Such partition contributes to all sharing partitions DOFs.
 * Such partition contribute to all allocated cod numbers. The partition with lowest rank is
 * selected as reference one.
 */
class FETIBoundaryDofManager
{
protected:
    /// Associated global number of dofManager.
    int globalNumber;
    /// Total number of partitions sharing receiver.
    int numberOfPartitions;
    /// Number of nonprescribed dofs, i.e, those, for which equation is necessary
    int ndofs;
    /**
     * Reference partition is partition to which other partitions sharing the dof manager
     * are linked using lagrange multipliers. We use the partition, which has its number
     * the lowest from all sharing partitions.
     */
    int referencePartition;
    /// List of partitions sharing dof manager.
    IntArray partitions;
    /**
     * Contains code numbers for each linked partition for each DOF (ndofs*(numberOfPartitions-1) DOFs).
     * The reference partition contributes to all code numbers.
     * The code numbers for n-th linked partition (its number is determined by index
     * of corresponding entry in partitions array) are stored at ((n-1)*ndofs+1, ..., n*ndofs) positions
     * in codeNumbers array.
     */
    IntArray codeNumbers;

public:
    FETIBoundaryDofManager();
    FETIBoundaryDofManager(int num, int part, int ndof);
    FETIBoundaryDofManager(const FETIBoundaryDofManager &);

    /// Returns number of partitions sharing receiver.
    int giveNumberOfSharedPartitions() { return numberOfPartitions; }
    /// Returns number of DOFs (with associated equation) of receiver.
    int giveNumberOfDofs() { return ndofs; }
    /// Returns corresponding global number of receiver.
    int giveGlobalNumber() { return globalNumber; }
    /// Returns reference partition number of receiver.
    int giveReferencePratition() { return referencePartition; }
    /// Returns number of i-th shared partition of receiver.
    int giveSharedPartition(int i) { return partitions.at(i); }
    /**
     * Returns code number corresponding to partition number partition_num and to dof_num-th DOF
     * @param partition_num Partition number for which code number is required.
     * @param dof_num Specifies the particular DOF.
     * @return Value of corresponding code number, zero if such partition does not share.
     * the receiver or if code number for reference partition is requested.
     */
    int giveCodeNumber(int partition_num, int dof_num);
    /**
     * Returns code numbers for all DOFs associated with shared partition.
     * @param rank Partition number.
     * @param locationArray The location array of size ndof.
     * @return Nonzero if o.k, zero if no such partition shared or if
     * code numbers for reference partition required.
     */
    int giveCompleteLocationArray(int rank, IntArray &locationArray);
    /**
     * Adds partition to list of partitions, sharing this dof manager.
     * The referencePartition is updated if necessary.
     * @param partitionNumber New partition number (0..size-1)
     */
    void addPartition(int partitionNumber);
    /**
     * Associates the equation numbers to particular DOFs.
     * @param equationCounter Current equation counter, updated.
     * @return New value of equationCounter.
     */
    int setCodeNumbers(int &equationCounter);
};
} // end namespace oofem

#endif // fetiboundarydofman_h
