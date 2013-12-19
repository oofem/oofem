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

#ifndef pfemparticle_h
#define pfemparticle_h

#include "node.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

///@name Input fields for Pfemparticle
//@{
#define _IFT_PFEMParticle_Name "pfemparticle"
#define _IFT_Node_coords "coords"
#define _IFT_Node_lcs "lcs"
//@}

namespace oofem {
class FloatArray;
class IntArray;

/**
 * Particle class being used in PFEM computations
 */
class PFEMParticle : public Node
{
private:
    /// the particle is not building any element
    bool free;

	/// the particle is a part of alpha-shape
	bool alphaShapeFlag;
	
	
	FloatArray coordinatesAtTimeStepBegin;

public:
    /**
     * Constructor. Creates a particle  with number n, belonging to aDomain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    PFEMParticle(int n, Domain * aDomain);
    /**
     * Destructor.
     */
    ~PFEMParticle(void) { }

    /**
     * Initializes receiver acording to object description stored in input record.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /**
     * Checks internal data consistency in node.
     * @return nonzero if receiver check is o.k.
     */
    int checkConsistency();

    virtual void updateYourself(TimeStep *tStep);

    /// Returns the free-propery flag
    bool isFree() { return free; }
    /// Sets the free-property flag
    virtual void setFree(bool newFlag = true) { free = newFlag; }

	/// Returns true if the particle is on alpha shape
	bool isOnAlphaShape() { return alphaShapeFlag; }
	/// Sets the alphaShapeFlag
	virtual void setOnAlphaShape(bool newFlag = true) { alphaShapeFlag = newFlag; }

	void storeCoordinatesTimeStepBegin();
	
	void updateNodalCoordinates(TimeStep* tStep);
	
	void resetNodalCoordinates();


    virtual void printOutputAt(FILE *stream, TimeStep *stepN);

    /**
     * Returns class name of the receiver.
     */
    const char *giveClassName() const { return "PFEMParticle"; }
    /**
     * Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType giveClassID() const { return PFEMParticleClass; }

#ifdef __OOFEG
    virtual void drawScalar(oofegGraphicContext &gc);
#endif
};
} // end namespace oofem
#endif // pfemparticle_h
