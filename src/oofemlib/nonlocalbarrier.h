/* $Header: /home/cvs/bp/oofem/oofemlib/src/nonlocalbarrier.h,v 1.3.4.1 2004/04/05 15:19:43 bp Exp $ */
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
// class NonlocalBarrier
//
#ifndef nonlocalbarrier_h
#define nonlocalbarrier_h

#include "femcmpnn.h"
#include "domain.h"

#include "flotarry.h"
#include "nonlocalmaterialext.h"

namespace oofem {

/**
 * Abstract base class for all nonlocal barriers. The purpose of this class is to
 * model barrier for nonlocal averaging process (visibility criterion, symmetry condition).
 * Usually, the given remote integration point influences to the source point
 * nonlocal average if the averaging function at source point and evaluated for
 * remote point has nonzero value. The barrier allows to exclude additional points,
 * which may be close enough, but due to several reasons there is no influence
 * between these points (for example, they can be  separated by a notch).
 *
 * @see NonlocalMaterialStatusExtensionInterface class.
 */
class NonlocalBarrier : public FEMComponent
{
protected:
public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    NonlocalBarrier(int n, Domain *aDomain);   // constructors
    /// Virtual destructor.
    virtual ~NonlocalBarrier() { };               // destructor

    /**
     * Abstract method returning true if the barrier is activated
     * by interaction of two given points. In this case the nonlocal influence
     * is not considered. Otherwise returns false.
     * @param c1 coordinates of first point
     * @param c2 coordinates of second point
     * @return true if barrier is activated, false otherwise
     */
    // virtual bool isActivated (const FloatArray& c1, const FloatArray& c2) = 0;

    /**
     * Abstract method modifying the integration weight between master (c1) and source (c2) point.
     * @param c1 coordinates of master point
     * @param c2 coordinates of source point
     * @param weight original integration weight; on output modified weight
     * @param shieldFlag (output param; set to true if shielding is activated)
     */
    virtual void applyConstraint(const FloatArray &c1, const FloatArray &c2, double &weight,
                                 bool &shieldFlag, NonlocalMaterialExtensionInterface *nei) = 0;


    /// Returns class name of the receiver.
    const char *giveClassName() const { return "NonlocalBarrier"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return NonlocalBarrierClass; }
};

} // end namespace oofem
#endif // nonlocalbarrier_h

