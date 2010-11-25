/* $Header: /home/cvs/bp/oofem/sm/src/Attic/symmetrybarrier.h,v 1.1.2.1 2004/04/05 15:19:47 bp Exp $ */
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


//
// class PolylineNonlocalBarrier
//
#ifndef symmetrybarrier_h
#define symmetrybarrier_h

#include "nonlocalbarrier.h"
#include "domain.h"

#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
/**
 * Implementation of symmetry nonlocal barier.
 * It allows to specify up to three planes (orthogonal ones) of symmetry
 * It then modifies the integration weights of source points to take into account
 * symmetry of the averaged field.
 * @see NonlocalBarrier class.
 */
class SymmetryBarrier : public NonlocalBarrier
{
protected:
    FloatArray origin;
    FloatArray normals;
    IntArray mask;
    FloatMatrix lcs;

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    SymmetryBarrier(int n, Domain *aDomain);  // constructors
    /// Virtual destructor.
    virtual ~SymmetryBarrier();                // destructor

    /**
     * Abstract method modifying the integration weight between master (c1) and source (c2) point.
     * @param c1 coordinates of master point
     * @param c2 coordinates of source point
     * @param weight original integration weight; on output modified weight
     * @param shieldFlag (output param; set to true if shielding is activated)
     */
    virtual void applyConstraint(const FloatArray &c1, const FloatArray &c2, double &weight,
                                 bool &shieldFlag, NonlocalMaterialExtensionInterface *nei);
    /** Initializes receiver acording to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. Input record can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "SymmetryBarrier"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return NonlocalBarrierClass; }
};
} // end namespace oofem
#endif // symmetrybarrier_h
