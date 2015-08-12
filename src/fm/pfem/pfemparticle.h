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
 *
 * @author David Krybus
 */
class PFEMParticle : public Node
{
protected:
    /// The particle does not compose any element, but still a part of solution domain and moves obeying Newton's laws of motion
    bool freeFlag;
    /// the particle is a part of alpha-shape
    bool alphaShapeFlag;
    /// Too close particles can be deactivated, e.g. removed from meshing and computation.
    bool activeFlag;

public:
    /**
     * Constructor. Creates a particle  with number n, belonging to aDomain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    PFEMParticle(int n, Domain *aDomain);
    /**
     * Destructor.
     */
    ~PFEMParticle(void) { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    virtual void updateYourself(TimeStep *tStep);

    /// Returns the free-propery flag
    virtual bool isFree() { return freeFlag; }
    /// Sets the free-property flag
    virtual void setFree(bool newFlag = true) { freeFlag = newFlag; }

    /// Returns true if the particle is on alpha shape
    virtual bool isOnAlphaShape() { return alphaShapeFlag; }
    /// Sets the alphaShapeFlag
    virtual void setOnAlphaShape(bool newFlag = true) { alphaShapeFlag = newFlag; }

    /// Returns the activeFlag
    virtual bool isActive() { return activeFlag; }
    /// Sets the activeFlag to false
    virtual void deactivate() { activeFlag = false; }

    virtual void printOutputAt(FILE *stream, TimeStep *stepN);

    virtual const char *giveClassName() const { return "PFEMParticle"; }
    virtual const char *giveInputRecordName() const { return _IFT_PFEMParticle_Name; }

#ifdef __OOFEG
    virtual void drawScalar(oofegGraphicContext &gc);
#endif
};
} // end namespace oofem
#endif // pfemparticle_h
