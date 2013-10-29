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

#ifndef mmashapefunctprojection_h
#define mmashapefunctprojection_h

#include "alist.h"
#include "materialmappingalgorithm.h"
#include "nodalrecoverymodel.h"
#include "interface.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * Element interface related to the MMAShapeFunctProjection.
 * Declares the required services at element level, which
 * are required by MMAShapeFunctProjection algorithm.
 */
class OOFEM_EXPORT MMAShapeFunctProjectionInterface : public Interface
{
public:
    /**
     * Typedefs to introduce the container type for nodal numbers
     */
    typedef AList< FloatArray >nodalValContainerType;

    enum coordType { coordType_local, coordType_global };

    /// Constructor
    MMAShapeFunctProjectionInterface() { }

    /**
     * @name The element interface required by MMAShapeFunctProjectionInterface
     */
    //@{
    /**
     * Interpolates the internal variables, using given nodal values
     * to given point (given by global coordinates) using element shape functions.
     * @param answer Computed internal variable.
     * @param coords Global/local coordinates of point of interest, see param ct.
     * @param ct Determines type of coordinates system used.
     * @param list Container of nodal values.
     * @param type Internal variable type.
     * @param tStep Solution step.
     */
    virtual void MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep) = 0;
    //@}
};



/**
 * The class implements the transfer of state variables based on
 * projection using shape functions.
 *
 * The algorithm of projecting internal vars (q) can be summarized as follows:
 * -# Gauss point components @f$ q^h @f$ are projected to nodes @f$ q^h_N @f$ of old mesh.
 * -# Nodal components are transferred from the old to new mesh (h+1):
 *   -# For each node A of new mesh (h+1) the so-called background element (in the old mesh (h)) is found.
 *   -# The local coordinates within the background element that correspond to the position of
 *      the node A are computed.
 *   -# By using the shape functions of background element, the state variables q are mapped
 *      from the nodes B of the old mesh (h) to the nodes A of the new mes (+1).
 * -# The state variables at the integration points of the new mesh (h+1) are  obtained by
 *    employing the shape functions of elements of the new mesh (h+1).
 *
 * It is obvious, that for efficiency, it is necessary to compute the nodal values A only once and
 * use them for all values. Therefore, the mmapper is typically declared as static material member,
 * and is used by all IPs. Also it should be used only for mapping of one internal variable.
 * For each internal variable there should be corresponding mapper instance.
 */
class OOFEM_EXPORT MMAShapeFunctProjection : public MaterialMappingAlgorithm
{
protected:
    /// Smother.
    AList< NodalRecoveryModel >smootherList;
    /// Solution state counter.
    StateCounterType stateCounter;
    /// Internal variables in list.
    IntArray intVarTypes;
    /// Source domain.
    Domain *domain;

public:
    /// Constructor
    MMAShapeFunctProjection();
    /// Destructor
    virtual ~MMAShapeFunctProjection();

    virtual void __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep);

    virtual void finish(TimeStep *tStep);

    virtual int mapVariable(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);

    virtual const char *giveClassName() const { return "MMAShapeFunctProjectionInterface"; }
};
} // end namespace oofem
#endif // mmashapefunctprojection_h
