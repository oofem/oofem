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

#ifndef polylinenonlocalbarrier_h
#define polylinenonlocalbarrier_h

#include "nonlocalbarrier.h"
#include "domain.h"

#include "flotarry.h"

namespace oofem {
/**
 * Implemantation of polyline nonlocal barrier.
 * It is a composite one-dimensional cell consisting of one or more connected lines.
 * The polyline is defined by an ordered list of n+1 vertices (nodes),
 * where n is the number of lines in the polyline.
 * Each pair of points (i,i+1) defines a line.
 *
 * The purpose of this class is to
 * model barrier for nonlocal averaging process (visibility criterion).
 * Usually, the given remote integration point influences to the source point
 * nonlocal average if the averaging function at source point and evaluated for
 * remote point has nonzero value. The barrier allows to exclude additional points,
 * which may be close enough, but due to several reasons there is no influence
 * between these points (for example, they can be  separated by a notch).
 */
class PolylineNonlocalBarrier : public NonlocalBarrier
{
protected:
    /// Local x-coordinate index.
    int localXCoordIndx;
    /// Local y-coordinate index.
    int localYCoordIndx;
    /// List of polyline vertices.
    IntArray vertexNodes;

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    PolylineNonlocalBarrier(int n, Domain *aDomain);
    /// Virtual destructor.
    virtual ~PolylineNonlocalBarrier();

    virtual bool isActivated(const FloatArray &c1, const FloatArray &c2);

    virtual void applyConstraint(const FloatArray &c1, const FloatArray &c2, double &weight,
                                 bool &shieldFlag, NonlocalMaterialExtensionInterface *nei);

    virtual IRResultType initializeFrom(InputRecord *ir);

    const char *giveClassName() const { return "PolylineNonlocalBarrier"; }
    classType giveClassID() const { return PolylineNonlocalBarrierClass; }
};
} // end namespace oofem
#endif // polylinenonlocalbarrier_h
