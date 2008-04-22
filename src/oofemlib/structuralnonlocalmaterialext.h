/* $Header: /home/cvs/bp/oofem/oofemlib/src/structuralnonlocalmaterialext.h,v 1.8 2003/04/06 14:08:26 bp Exp $ */
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
// class StructuralNonlocalMaterialExtension
//

#ifndef structuralnonlocalmaterialext_h
#define structuralnonlocalmaterialext_h

#include "nonlocalmaterialext.h"
#include "cltypes.h"
#include "matstatus.h"


/**
 * Base class for all nonlocal structural material statuses.
 */
class StructuralNonlocalMaterialStatusExtensionInterface : public NonlocalMaterialStatusExtensionInterface
{
protected:

public:
    // StructuralNonlocalMaterialStatus(int n, Domain* d, GaussPoint* g) : NonlocalMaterialStatus (n,d,g) {}
    StructuralNonlocalMaterialStatusExtensionInterface() : NonlocalMaterialStatusExtensionInterface() { }
    ~StructuralNonlocalMaterialStatusExtensionInterface() { }
};




/**
 * Abstract base class for all nonlocal structural materials. Nonlocal in sence, that response in particular
 * point depends not only on state in that point, but also takes into account state of surrounding
 * points. Response typically depends on some nonlocal quantity obtained as nonlocal average over
 * some characteristic volume.
 * This class declares the necessary interface for all nonlocal structural constitutive  models.
 */
class StructuralNonlocalMaterialExtensionInterface : public NonlocalMaterialExtensionInterface
{
protected:

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    //StructuralNonlocalMaterial (int n,Domain* d) : NonlocalMaterial(n,d)
    StructuralNonlocalMaterialExtensionInterface(Domain *d) : NonlocalMaterialExtensionInterface(d)
    { }
    /// Destructor.
    ~StructuralNonlocalMaterialExtensionInterface()                { }



    /**
     * Declares the service updating local variables in given integration points,
     * which take part in nonlocal average process.
     * Because value of single integration point influences nonlocal variables in several near
     * integration points, it is suitable to compute these variables only once. These should be stored
     * in integration point associated statuses.
     * The implementation is left on derived classes.
     * Provide material local strain increment - as is provided to computeRealStresVector.
     * This allows to update internal vars to be averaged to new state
     * @param strainVector total strain vector in given integration point.
     * @param gp integration point to update.
     * @param atTime solution step indicating time of update.
     */
    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime) = 0;
};

#endif // structuralnonlocalmaterialext_h


