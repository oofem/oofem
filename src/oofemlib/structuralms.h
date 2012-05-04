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

#ifndef structuralms_h
#define structuralms_h

#include "matstatus.h"
#include "flotarry.h"

namespace oofem {

class GaussPoint;
class Dictionary;
class Domain;

/**
 * This class implements a structural material status information. It is attribute of
 * gaussPoint. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 *
 * This is a base class for all material statuses corresponding to materials derived from
 * structural material class.
 * It defines stress and strain vectors and their increments.
 * Functions for accessing these components are defined.
 *
 * Tasks:
 * This is abstract class - only basic functionality is supported like:
 * - maintaining and providing access to stress and strain vectors
 *   (including their increments)
 * - storing and restoring status on tape
 * - printingYourself()
 * - updating Yourself after a new equilibrium state has been reached.
 */
class StructuralMaterialStatus : public MaterialStatus
{
protected:
    /// Equilibrated strain vector in reduced form
    FloatArray strainVector;
    /// Equilibrated stress vector in reduced form
    FloatArray stressVector;
    /// Temporary stress vector in reduced form (increments are used mainly in nonlinear analysis)
    FloatArray tempStressVector;
    /// Temporary strain vector in reduced form (to find balanced state)
    FloatArray tempStrainVector;

public:
    /// Constructor. Creates new StructuralMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    StructuralMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~StructuralMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns the const pointer to receiver's strain vector.
    const FloatArray &giveStrainVector() { return strainVector; }
    /// Returns the const pointer to receiver's stress vector.
    const FloatArray &giveStressVector() { return stressVector; }
    /// Returns the const pointer to receiver's temporary strain vector.
    const FloatArray &giveTempStrainVector() { return tempStrainVector; }
    /// Returns the const pointer to receiver's temporary stress vector.
    const FloatArray &giveTempStressVector() { return tempStressVector; }
    /// Assigns strain vector to given vector v.
    void letStrainVectorBe(const FloatArray &v) { strainVector = v; }
    /// Assigns stressVector to given vector v.
    void letStressVectorBe(const FloatArray &v) { stressVector = v; }
    /// Assigns tempStressVector to given vector v.
    void letTempStressVectorBe(const FloatArray &v) { tempStressVector = v; }
    /// Assigns tempStrainVector to given vector v
    void letTempStrainVectorBe(const FloatArray &v) { tempStrainVector = v; }

    virtual const char *giveClassName() const { return "StructuralMaterialStatus"; }
    virtual classType giveClassID() const { return StructuralMaterialStatusClass; }
};
} // end namespace oofem
#endif // structuralms_h
