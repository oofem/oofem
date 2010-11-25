/* $Header: /home/cvs/bp/oofem/sm/src/mmaclosestiptransfer.h,v 1.6 2003/05/19 13:04:00 bp Exp $ */
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

//   **************************************************************
//   *** CLASS CONTAINING ELEMENT IP TRANSFER MAPPING ALGORITHM ***
//   **************************************************************

#ifndef mmacontainingelementprojection_h
#define mmacontainingelementprojection_h

#include "compiler.h"

#include "materialmappingalgorithm.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * The class implements the transfer of internal variabless based on containing element
 * The element containing the receiving gauss point is localized and its closest gauss point
 * is used as a source gauss point that determines the values of receiver.
 */
class MMAContainingElementProjection : public MaterialMappingAlgorithm
{
protected:

    GaussPoint *source;
    // GaussPoint* recv;
    // FloatArray rc;
public:
    /** Constructor
     * @param gp Integration point belonging to new domain to which mapping occur
     * @param oldd old mesh reference
     */
    MMAContainingElementProjection();
    /**
     * Initializes the receiver state before mapping. The idea is to place some
     * comon global oprerations before mapping particular IP's if necessary.
     * Stores Times stamp of last initialization, so multiple calls for same step
     * do not initialize receiver again.
     * @note new domain can be obtained from given ip.
     * @param dold old domain
     * @param varTypes array of InternalStateType values, identifiing all vars to be mapped
     * @param coords coordinates of the receiver point
     * @param region if > 0 region id of receiver point,, if < 0 ignore regions.
     * @param tStep time step
     */
    void __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep);
    /**
     * Finishes the mapping for given time step. Used to perform cleanup.
     * Typically some mappers reguire to compute some global mesh data related to
     * current step, which are valid for all IPs - so they are computed only once for
     * all IPs.
     */
    void finish(TimeStep *tStep) { };
    /** Maps and update the unknown of given type from
     * old mesh oldd to new mesh to which gp belongs to. The result is stored in answer array.
     * The closest IP point to specified one is found and its state variable is
     * returned as a result.
     * @param answer contains result
     * @param type determines the type of internal variable
     * @param coords coordinates of receiver to which mapping occur
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);
    /** Returns class name of the receiver */
    const char *giveClassName() const { return "MMAContainingElementProjection"; }

protected:
};
} // end namespace oofem
#endif // mmacontainingelementprojection_h
