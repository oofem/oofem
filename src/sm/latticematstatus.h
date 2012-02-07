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

#ifndef latticematstatus_h
#define latticematstatus_h

#ifndef __MAKEDEPEND

#endif
#include "structuralms.h"
#include "classtype.h"


class GaussPoint;
class Dictionary;
class Domain;
class NonlocalMaterialStatusExtension;

namespace oofem {
/**
 * This class implements a lattice material status.
 */
class LatticeMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Associated intagration point.
    GaussPoint *gp;
public:

    LatticeMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    virtual ~LatticeMaterialStatus() { }
    /// Print receiver's output to given stream.
    void printOutputAt(FILE *, TimeStep *) { }

    virtual void initTempStatus() { }

    virtual void updateYourself(TimeStep *) { } // update after new equilibrium state reached

    // definition
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "LatticeMaterialStatus"; }
    /// Returns classType id of receiver.
    classType giveClassID() const { return LatticeMaterialStatusClass; }

    ///Sets the temp_crack_flag
    virtual void setTempCrackFlag(int val) = 0;

    /// Returns the crack_flag
    virtual int giveCrackFlag() { return 0; }

    virtual double giveCrackWidth() { return 0; }

    virtual double giveDissipation() { return 0; }

    virtual double giveDeltaDissipation() { return 0; }

    IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};
} // end namespace oofem
#endif // matstatus_h
