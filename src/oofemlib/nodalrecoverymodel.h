/* $Header: /home/cvs/bp/oofem/oofemlib/src/nodalrecoverymodel.h,v 1.10 2003/04/06 14:08:25 bp Exp $ */
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


//   **********************************
//   *** CLASS NODAL RECOVERY MODEL ***
//   **********************************

#ifndef nodalrecoverymodel_h
#define nodalrecoverymodel_h

#include "compiler.h"

#include "tdictionary.h"
#include "intarray.h"
#include "flotarry.h"
#include "alist.h"
#include "interface.h"
#include "internalstatetype.h"
#include "statecountertype.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
#endif

namespace oofem {
class Domain;
class Element;
class CrossSection;
class TimeStep;
/**
 * The base class for all recovery models, which perform nodal averaging or projection
 * processes for internal variables typically stored in integration points.
 * Since typically many problems results in discontinuos  approximation of stresses
 * some recovery designed to yield a realistic field avoiding unwarranted discontinuities.
 *
 * The NodalRecoveryModel class provides common array of nodal dictionaries, where
 * the recovered nodal values are stored for each region.
 */
class NodalRecoveryModel
{
protected:
    typedef TDictionary< int, FloatArray >vectorDictType;
    /**
     * Array of nodal dictionaries, containing nodal values for each region.
     * The region id is dictionary key to corresponding values.
     */
    AList< vectorDictType >nodalValList;
    /** Determines the type of recovered values */
    InternalStateType valType;
    /** Time stamp of recovered values */
    StateCounterType stateCounter;
    Domain *domain;
    /**
     * Array of value record sizes per region.
     * It is typically determined during recovery phase.
     * The purpose of ths attribute is to cache the result for
     * future requests.
     */
    // IntArray regionValSize;
    /**
     * Array of region values masks, containing mask of reduced
     * indexes of Internal Variable components.
     * It is typically determined during recovery phase.
     * The purpose of ths attribute is to cache the result for
     * future requests.
     *
     * @see Element::giveIPValue and Element::giveIntVarCompFullIndx methods
     */
    // AList<IntArray> regionValueMaps;

#ifdef __PARALLEL_MODE
    /// Common Communicator buffer
    CommunicatorBuff *commBuff;
    /// Communicator.
    ProblemCommunicator *communicator;
    /// communicato init flag
    bool initCommMap;
#endif

public:
    /// Constructor
    NodalRecoveryModel(Domain *d);
    /// Destructor
    virtual ~NodalRecoveryModel();
    /** Recovers the nodal values for all regions of given Domain.
     * @param d domain of interest
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    virtual int recoverValues(InternalStateType type, TimeStep *tStep) = 0;
    /**
     * Clears the receiver's nodal table.
     * @return nonzero if o.k.
     */
    virtual int clear();
    /**
     * Returns vector of recovered values for given node and region.
     * @param ptr pointer to recovered values at node, NULL if not present
     * @param node node number
     * @param region region number;
     * @return nonzero if values are defined, zero otherwise
     */
    int giveNodalVector(const FloatArray * &ptr, int node, int region);
    /**
     * Test if recovered values for given node and region exist.
     * @param node node number
     * @param region region number
     * @return nonzero if entry is in table, zero otherwise
     */
    int includes(int node, int region);
    /**
     * Determines the number of material regions of domain.
     * In the current iimplementation the region is associated with cross section model.
     */
    //int giveNumberOfRegions ();
    /**
     * Returns the region id of given element
     * @param element pointer to element which region id is requsted
     * @return region id (number) for this element
     */
    //int giveElementRegion (Element* element);
    /**
     * Returns the region record size. The default implementation scans the elements and once one belonging
     * to given region is found it is requested for the information. Slow.
     * The overloaded instances can chache these results, since they are easily
     * obtainable during recovery.
     * @param reg region id
     * @param type determines the type of variable, for which size is requested. Should be same as used
     * for recovering values.
     */
    virtual int giveRegionRecordSize(int reg, InternalStateType type);
    /**
     * Returns the region values masks, containing mask of reduced
     * indexes of Internal Variable component for given region.
     * The default implementation scans the elements and once one belonging
     * to given region is found it is requested for the information. Slow.
     * The overloaded instances can chache these results, since they are easily
     * obtainable during recovery.
     * @param answer contains result
     * @param reg region id
     * @param type determines the type of variable, for which size is requested. Should be same as used
     * for recovering values.
     */
    virtual void giveRegionRecordMap(IntArray &answer, int reg, InternalStateType type);
protected:
    int init();
    /**
     * Same as public giveNodalVector,but returns non-const pointer.
     */
    FloatArray *giveNodalVectorPtr(int node, int region);
    /**
     * Determine local region node numbering and determine and check nodal values size
     * @param regionNodalNumbers on return array containing for each dofManager its local region number
     * @param regionDofMans on output total number of region dofMans
     * @param ireg region number
     * @returns nonzero if ok, zero if region has to be skipped
     */
    int initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans, int reg);

    /**
     * Update the nodal table acoording to recovered solution for given region
     * @param ireg region number
     * @param regionNodalNumbers array containing for each dofManager its local region number
     * @param regionValSizevalue size of dofMan record
     * @param rhs array with recovered values
     */
    int updateRegionRecoveredValues(const int ireg, const IntArray &regionNodalNumbers,
                                    int regionValSize, const FloatArray &rhs);
};
} // end namespace oofem
#endif // nodalrecoverymodel_h
