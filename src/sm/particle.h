/*
 * class Particle added by Milan Jirasek on 1 Feb 2010
 *
 *****    *****   ******  ******  ***   ***
 **   **  **   **  **      **      ** *** **
 **   **  **   **  ****    ****    **  *  **
 **   **  **   **  **      **      **     **
 **   **  **   **  **      **      **     **
 *****    *****   **      ******  **     **
 *****
 *****
 *****         OOFEM : Object Oriented Finite Element Code
 *****
 *****           Copyright (C) 1993 - 2010   Borek Patzak
 *****
 *****
 *****
 *****   Czech Technical University, Faculty of Civil Engineering,
 *****Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *****
 *****This program is free software; you can redistribute it and/or modify
 *****it under the terms of the GNU General Public License as published by
 *****the Free Software Foundation; either version 2 of the License, or
 *****(at your option) any later version.
 *****
 *****This program is distributed in the hope that it will be useful,
 *****but WITHOUT ANY WARRANTY; without even the implied warranty of
 *****MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *****GNU General Public License for more details.
 *****
 *****You should have received a copy of the GNU General Public License
 *****along with this program; if not, write to the Free Software
 *****Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef particle_h
#define particle_h

#include "node.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {
class FloatArray;
class IntArray;

/**
 * Class implementing spherical particles as special nodes having a certain radius.
 * Such particles are used by the cohesive particle model.
 */
class Particle : public Node
{
protected:
    /// particle radius (only spherical particles considered)
    double radius;

public:
    /**
     * Constructor. Creates a particle with number n, belonging to aDomain.
     * @param n particle number in domain aDomain
     * @param aDomain domain to which the particle belongs
     */
    Particle(int n, Domain *aDomain);
    /**
     * Destructor.
     */
    ~Particle(void) {}

    /**
     * Initializes receiver acording to object description stored in input record.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /**
     * Returns class name of the receiver.
     */
    // note: if the following method returned "Particle",
    // the extractor could not read displacements from the output file
    const char *giveClassName() const { return "Node"; }
    /**
     * Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType giveClassID() const { return ParticleClass; }
    /// Returns the radius of the particle
    double giveRadius() const { return radius; }
};
} // end namespace oofem
#endif
