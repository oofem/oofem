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

#ifndef interactionpfemparticle_h
#define interactionpfemparticle_h

#include "pfemparticle.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

///@name Input fields for Pfemparticle
//@{
#define _IFT_InteractionPFEMParticle_Name "interactionpfemparticle"
#define _IFT_Node_coords "coords"
#define _IFT_Node_lcs "lcs"
#define _IFT_InteractionPFEMParticle_CoupledNode "couplednode"
//@}

namespace oofem {
class FloatArray;
class IntArray;

/**
 * Particle class being used in PFEM computations
 */
class InteractionPFEMParticle : public PFEMParticle
{
protected:
	int coupledNode;

public:
    /**
     * Constructor. Creates a particle  with number n, belonging to aDomain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    InteractionPFEMParticle(int n, Domain * aDomain);
    /**
     * Destructor.
     */
    ~InteractionPFEMParticle(void) { }

    /**
     * Initializes receiver acording to object description stored in input record.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Checks internal data consistency in node.
     * @return nonzero if receiver check is o.k.
     */
    virtual int checkConsistency();

    virtual void updateYourself(TimeStep *tStep);

    virtual void printOutputAt(FILE *stream, TimeStep *stepN);

    /**
     * Returns class name of the receiver.
     */
    virtual const char *giveClassName() const { return "InteractionPFEMParticle"; }

	virtual const char *giveInputRecordName() const { return _IFT_InteractionPFEMParticle_Name; }
    /**
     * Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    virtual classType giveClassID() const { return InteractionPFEMParticleClass; }

#ifdef __OOFEG
    virtual void drawScalar(oofegGraphicContext &gc);
#endif
};
} // end namespace oofem
#endif // interactionpfemparticle_h
