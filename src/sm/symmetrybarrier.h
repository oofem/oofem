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

#ifndef symmetrybarrier_h
#define symmetrybarrier_h

#include "nonlocalbarrier.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for SymmetryBarrier
//@{
#define _IFT_SymmetryBarrier_Name "symmetrybarrier"
#define _IFT_SymmetryBarrier_origin "origin"
#define _IFT_SymmetryBarrier_normals "normals"
#define _IFT_SymmetryBarrier_activemask "activemask"
//@}

namespace oofem {
/**
 * Implementation of symmetry nonlocal barrier.
 * It allows to specify up to three planes (orthogonal ones) of symmetry
 * It then modifies the integration weights of source points to take into account
 * symmetry of the averaged field.
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
     * @param d Pointer to the domain to which element belongs.
     */
    SymmetryBarrier(int n, Domain * d);
    /// Destructor.
    virtual ~SymmetryBarrier();

    virtual void applyConstraint(const FloatArray &c1, const FloatArray &c2, double &weight,
                                 bool &shieldFlag, NonlocalMaterialExtensionInterface *nei);
    virtual double calculateMinimumDistanceFromBoundary(const FloatArray &coords) { return 1.e10; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return _IFT_SymmetryBarrier_Name; }
    virtual const char *giveClassName() const { return "SymmetryBarrier"; }
};
} // end namespace oofem
#endif // symmetrybarrier_h
