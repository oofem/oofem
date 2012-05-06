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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef particle_h
#define particle_h

#include "node.h"

namespace oofem {

/**
 * Class implementing spherical particles as special nodes having a certain radius.
 * Such particles are used by the cohesive particle model.
 * @author Milan Jirasek
 */
class Particle : public Node
{
protected:
    /// Particle radius (only spherical particles considered).
    double radius;

public:
    /**
     * Constructor. Creates a particle with number n, belonging to aDomain.
     * @param n Particle number in domain aDomain.
     * @param aDomain Domain to which the particle belongs.
     */
    Particle(int n, Domain *aDomain);
    /**
     * Destructor.
     */
    virtual ~Particle(void) {}

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "Particle"; }
    virtual classType giveClassID() const { return ParticleClass; }

    /// Returns the radius of the particle.
    double giveRadius() const { return radius; }
};
} // end namespace oofem
#endif
