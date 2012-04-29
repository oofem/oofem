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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
 *
 * The recovery can be performed independently on regions of the domain to account for potential discontinuity
 * of recovered variable over region boundaries. To make this concept more general, the
 * virtual regions are introduced. The idea is that several real regions can map to a single virtual region.
 * The mapping is defined by so called virtualRegionMap. The whole domain recovery can be
 * obtained by a single virtual region to which all real regions map.
 *
 * The NodalRecoveryModel class provides common array of nodal dictionaries, where
 * the recovered nodal values are stored for each region.
 */
class NodalRecoveryModel
{
public:
    enum NodalRecoveryModelType { NRM_NodalAveraging = 0, NRM_ZienkiewiczZhu = 1,  NRM_SPR = 2 };

protected:
    typedef TDictionary< int, FloatArray > vectorDictType;
    /**
     * Array of nodal dictionaries, containing nodal values for each region.
     * The region id is dictionary key to corresponding values.
     */
    AList< vectorDictType >nodalValList;
    /// Determines the type of recovered values.
    InternalStateType valType;
    /// Time stamp of recovered values.
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
    /**
     * Number of virtual regions, if positive.
     * When equals to zero a single virtual region, to which all real regions map, is assumed (whole domain recovery)
     * If negative, real regions are used instead.
     */
    int numberOfVirtualRegions;
    /**
     * Mapping from real regions to virtual regions.
     * The size of this array should be equal to number of real regions.
     */
    IntArray virtualRegionMap;


#ifdef __PARALLEL_MODE
    /// Common Communicator buffer.
    CommunicatorBuff *commBuff;
    /// Communicator.
    ProblemCommunicator *communicator;
    /// Communication init flag.
    bool initCommMap;
#endif

public:
    /// Constructor
    NodalRecoveryModel(Domain *d);
    /// Destructor
    virtual ~NodalRecoveryModel();
    /**
     * Recovers the nodal values for all regions.
     * @param type Determines the type of internal variable to be recovered.
     * @param tStep Time step.
     */
    virtual int recoverValues(InternalStateType type, TimeStep *tStep) = 0;
    /**
     * Clears the receiver's nodal table.
     * @return nonzero if o.k.
     */
    virtual int clear();
    /**
     * Initializes the receiver. Called form constructor, but when domain changes,
     * init call necessary to update data structucture.
     * @return nonzero if o.k.
     */
    int init();

    /**
     * Returns vector of recovered values for given node and region.
     * @param ptr Pointer to recovered values at node, NULL if not present.
     * @param node Node number.
     * @param region Region number.
     * @return Nonzero if values are defined, zero otherwise.
     */
    int giveNodalVector(const FloatArray * &ptr, int node, int region);
    /**
     * Test if recovered values for given node and region exist.
     * @param node Node number.
     * @param region Region number.
     * @return Nonzero if entry is in table, zero otherwise.
     */
    int includes(int node, int region);
    /**
     * Returns the region record size. The default implementation scans the elements and once one belonging
     * to given region is found it is requested for the information. Slow.
     * The overloaded instances can cache these results, since they are easily
     * obtainable during recovery.
     * @param reg Virtual region id.
     * @param type Determines the type of variable, for which size is requested. Should be same as used
     * for recovering values.
     */
    virtual int giveRegionRecordSize(int reg, InternalStateType type);
    /**
     * Returns the region values masks, containing mask of reduced
     * indexes of Internal Variable component for given region.
     * The default implementation scans the elements and once one belonging
     * to given region is found it is requested for the information. Slow.
     * The overloaded instances can cache these results, since they are easily
     * obtainable during recovery.
     * @param answer Contains result.
     * @param reg Virtual region id.
     * @param type Determines the type of variable, for which size is requested. Should be same as used
     * for recovering values.
     */
    virtual void giveRegionRecordMap(IntArray &answer, int reg, InternalStateType type);
    /**
     * Sets recovery mode (region by region or whole domain mode) and allows to specify virtual region mapping.
     * @param nvr Specifies number of regions. If set to zero, the recovery is performed over whole domain (single virtual region),
     * when positive, it should be equal to number of virtual regions and vrmap should be provided. If negative, then recovery over
     * real regions is used.
     * @param vrmap When nvr positive, it defines the mapping from real (true) regions to virtual regions. The size of this array
     * should equal to number of true regions and values should be in range <1, nvr>. When nvr is zero or negative, this parameter
     * is ignored. The i-th value defines mapping of true region number to virtual region number. Zero or negative value causes
     * region to be ignored.
     */
    void setRecoveryMode(int nvr, const IntArray &vrmap);
    /**
     * Returns element region number. If virtual region mapping is active, the element virtual region
     * is returned (using map) instead of real element region number.
     */
    int giveElementVirtualRegionNumber(int ielem);
    /**
     * Returns number of recovery regions.
     */
    int giveNumberOfVirtualRegions() { return this->numberOfVirtualRegions; }

protected:
    /**
     * Same as public giveNodalVector,but returns non-const pointer.
     */
    FloatArray *giveNodalVectorPtr(int node, int region);
    /**
     * Determine local region node numbering and determine and check nodal values size.
     * @param regionNodalNumbers on Return array containing for each dofManager its local region number.
     * @param regionDofMans On output total number of region dofMans.
     * @param reg Virtual region number.
     * @returns Nonzero if ok, zero if region has to be skipped.
     */
    int initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans, int reg);

    /**
     * Update the nodal table according to recovered solution for given region.
     * @param ireg Virtual region number.
     * @param regionNodalNumbers Array containing for each dofManager its local region number.
     * @param regionValSize Size of dofMan record.
     * @param rhs Array with recovered values.
     */
    int updateRegionRecoveredValues(const int ireg, const IntArray &regionNodalNumbers,
                                    int regionValSize, const FloatArray &rhs);
};
} // end namespace oofem
#endif // nodalrecoverymodel_h
