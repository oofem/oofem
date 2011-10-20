/* $Header: /home/cvs/bp/oofem/sm/src/mmashapefunctprojection.h,v 1.7 2003/05/19 13:04:00 bp Exp $ */
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

//   *********************************************************
//   *** CLASS SHAPE FUNCTION PROJECTION MAPPING ALGORITHM ***
//   *********************************************************

#ifndef mmashapefunctprojection_h
#define mmashapefunctprojection_h

#include "alist.h"
#include "materialmappingalgorithm.h"
#include "nodalrecoverymodel.h"
#include "interface.h"
#include "compiler.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * Element interface related to the MMAShapeFunctProjection.
 * Declares the required services at element level, which
 * are required by MMAShapeFunctProjection algorithm.
 */
class MMAShapeFunctProjectionInterface : public Interface
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
     * @param answer computed internal variable
     * @param coords global/local coordinates of point of interest, see param ct
     * @param ct determines type of coordinates system used (0=
     * @param list container of nodal values
     * @param type internal variable type
     * @param tStep solution step
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
 * 1) Gauss point components q^h are projected to nodes q^h_N of old mesh.
 * 2) Nodal components are transfered from the old to new mesh (h+1):
 *  2a) For each node A of new mesh (h+1) the so-called background element (in the old mesh (h))
 *     is found.
 * 2b) The local coordinates within the background element that correspond to the position of
 *     the node A are computed.
 * 2c) By using the shape functions of background element, the state variables q are mapped
 *     from the nodes B of the old mesh (h) to the nodes A of the new mes (+1).
 * 3) The state variables at the integration points of the new mesh (h+1) are  obtained by
 *  employing the shape functions of elements of the new mesh (h+1).
 *
 * It is obvious, that for efficiency, it is necessary to compute the nodal values A only once and
 * use them for all values. Therefore, the mmapper is typically declared as static material member,
 * and is used by all IPs. Also it should be used only for mapping of one internal varialble.
 * For each internal variable there should be corresponding mapper instance.
 *
 *
 */
class MMAShapeFunctProjection : public MaterialMappingAlgorithm
{
protected:
    /// smother
    AList< NodalRecoveryModel >smootherList;
    /// Solution state counter.
    StateCounterType stateCounter;
    /// Internal variables in list
    IntArray intVarTypes;
    /// Source domain
    Domain *domain;
public:
    /** Constructor
     * @param oldd old mesh reference
     */
    MMAShapeFunctProjection();
    /// Destructor
    ~MMAShapeFunctProjection();
    /**
     * Initializes the receiver state before mapping. The idea is to place some
     * comon global oprerations before mapping particular IP's if necessary.
     * Stores Times stamp of last initialization, so multiple calls for same step
     * do not initialize receiver again.
     * @param dold old domain
     * @param varTypes array of InternalStateType values, identifiing all vars to be mapped
     * @param coords coordinates of the receiver point
     * @param region if > 0 region id of receiver point,, if < 0 ignore regions.
     * @param tStep time step
     */
    virtual void __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep);
    /**
     * Finishes the mapping for given time step. Used to perform cleanup.
     * Typically some mappers reguire to compute some global mesh data related to
     * current step, which are valid for all IPs - so they are computed only once for
     * all IPs.
     */
    virtual void finish(TimeStep *tStep);
    /** Maps and update the unknown of given type from
     * old mesh oldd to new mesh to which gp belongs to. The result is stored in answer array.
     * The closest IP point to specified one is found and its state variable is
     * returned as a result.
     * @param answer contains result
     * @param type determines the type of internal variable
     * @param gp Integration point belonging to new domain to which mapping occur
     * @param oldd old mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int mapVariable(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    /** Maps and update the unknown of given type from
     * old mesh oldd to new mesh to which gp belongs to. The result is stored in answer array.
     * @param answer contains result
     * @param type determines the type of internal variable
     * @param coords coordinates of receiver point to which mapping occur
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);

    /** Returns class name of the receiver */
    const char *giveClassName() const { return "MMAShapeFunctProjectionInterface"; }
protected:
};
} // end namespace oofem
#endif // mmashapefunctprojection_h
