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

#ifndef nodalrecoverymodel_h
#define nodalrecoverymodel_h

#include "oofemenv.h"
#include "intarray.h"
#include "floatarray.h"
#include "interface.h"
#include "internalstatetype.h"
#include "statecountertype.h"
#include "set.h"

#include <map>
#include <vector>

namespace oofem {
class Domain;
class Element;
class CrossSection;
class TimeStep;
class CommunicatorBuff;
class ProblemCommunicator;
/**
 * The base class for all recovery models, which perform nodal averaging or projection
 * processes for internal variables typically stored in integration points.
 *
 * The recovery can be performed on one region. Region is defined by given set containing elements defing the region volume.
 * The element set can be set up in the input file, or created dynamically.
 * If Averaging over multiple regions is needed, multiple instances of nodal recovery model should be created,
 * o single instance reused.
 */
class OOFEM_EXPORT NodalRecoveryModel
{
public:
    enum NodalRecoveryModelType { NRM_NodalAveraging = 0, NRM_ZienkiewiczZhu = 1,  NRM_SPR = 2 };

protected:
    /**
     * Map of nodal values. Only nodes for active region are determined and stored.
     * The node number is dictionary key to corresponding values.
     */
    std :: map< int, FloatArray >nodalValList;
    /// Determines the type of recovered values.
    InternalStateType valType;
    /// Time stamp of recovered values.
    StateCounterType stateCounter;
    Domain *domain;

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
    NodalRecoveryModel(Domain * d);
    /// Destructor
    virtual ~NodalRecoveryModel();

    void setDomain(Domain *ipDomain) { domain = ipDomain; }

    /**
     * Recovers the nodal values for all regions.
     * @param type Determines the type of internal variable to be recovered.
     * @param tStep Time step.
     */
    virtual int recoverValues(Set elementSet, InternalStateType type, TimeStep *tStep) = 0;
    /**
     * Clears the receiver's nodal table.
     * @return nonzero if o.k.
     */
    virtual int clear();
    /**
     * Returns vector of recovered values for given node and region.
     * @param ptr Pointer to recovered values at node, NULL if not present.
     * @param node Node number.
     * @return Nonzero if values are defined, zero otherwise.
     */
    int giveNodalVector(const FloatArray * &ptr, int node);
    /**
     * Returns the region record size. Available after recovery.
     * @param reg Virtual region id.
     * @param type Determines the type of variable, for which size is requested. Should be same as used
     * for recovering values.
     */
    virtual int giveRegionRecordSize();

    virtual const char *giveClassName() const = 0;
    std :: string errorInfo(const char *func) { return std :: string(this->giveClassName()) + func; }

protected:
    /**
     * Determine local region node numbering and determine and check nodal values size.
     * @param regionNodalNumbers on Return array containing for each dofManager its local region number.
     * @param regionDofMans On output total number of region dofMans.
     * @param reg Virtual region number.
     * @returns Nonzero if ok, zero if region has to be skipped.
     */
    int initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans, Set &region);

    /**
     * Update the nodal table according to recovered solution for given region.
     * @param ireg Virtual region number.
     * @param regionNodalNumbers Array containing for each dofManager its local region number.
     * @param regionValSize Size of dofMan record.
     * @param rhs Array with recovered values.
     */
    int updateRegionRecoveredValues(const IntArray &regionNodalNumbers,
                                    int regionValSize, const FloatArray &rhs);
};
} // end namespace oofem
#endif // nodalrecoverymodel_h
