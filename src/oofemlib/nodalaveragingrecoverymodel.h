/* $Header: /home/cvs/bp/oofem/oofemlib/src/nodalaveragingrecoverymodel.h,v 1.2 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef nodalaveragingrecoverymodel_h
#define nodalaveragingrecoverymodel_h

#include "compiler.h"

#include "nodalrecoverymodel.h"
#include "interface.h"

class GaussPoint;

/**
 * The nodal recovery model based on nodal averaging". The recovery is based
 * on nodal averaging or projection process in which element contributions are averaged
 * at node with the same weight.
 */
class NodalAveragingRecoveryModel : public NodalRecoveryModel
{
protected:

public:
    /// Constructor
    NodalAveragingRecoveryModel(Domain *d);
    /// Destructor
    ~NodalAveragingRecoveryModel();
    /** Recovers the nodal values for all regions of given Domain.
     * @param d domain of interest
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    int recoverValues(InternalStateType type, TimeStep *tStep);
private:
    /**
     * Initializes the region table indicating regions to skip.
     * @param regionMap region tabl, the nonzero entry for region indicates region to skip due to
     * unsupported elements or incompatible value size
     * @param regionValSize contains the record size for each region
     * @param type determines the type of internal variable to be recovered
     */
    void initRegionMap(IntArray &regionMap, IntArray &regionValSize, InternalStateType type);
};

/**
 * The element interface required by ZZNodalRecoveryModel.
 */
class NodalAveragingRecoveryModelInterface : public Interface
{
public:
    /// Constructor
    NodalAveragingRecoveryModelInterface() { }

    /**
     * Computes the element value in given node.
     * @param answer contains the result
     * @param node element node number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep) = 0;
    /**
     * Computes the element value in given side.
     * @param answer contains the result
     * @param node element side number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                           InternalStateType type, TimeStep *tStep) = 0;
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type) = 0;
};


#endif // nodalaveragingrecoverymodel_h






