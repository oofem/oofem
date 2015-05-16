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

#ifndef microplanematerial_bazant_h
#define microplanematerial_bazant_h

#include "../sm/Materials/structuralms.h"
#include "microplanematerial.h"


namespace oofem {
/**
 * Abstract base class for all microplane models according to Bazant's approach.
 * Micro strains on microplane are described using magnitude of normal strain,
 * volumetric normal component and by two orthogonal shear components
 * (m and l direction) in microplane.
 */
class MicroplaneMaterial_Bazant : public MicroplaneMaterial
{
public:
    /**
     * Constructor. Creates Abstract Bazant's Microplane Material belonging
     * to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    MicroplaneMaterial_Bazant(int n, Domain * d);
    /// Destructor.
    virtual ~MicroplaneMaterial_Bazant() { }

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);


    /**
     * Updates the volumetric stress component after computing real stress microplane vectors.
     */
    virtual void updateVolumetricStressTo(Microplane *mPlane, double sigv) = 0;

    virtual const char *giveClassName() const { return "MicroplaneMaterial_Bazant"; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new StructuralMaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // microplanematerial_bazant_h
