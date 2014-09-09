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

#ifndef latticestructuralelement_h
#define latticestructuralelement_h

#include "../sm/Elements/structuralelement.h"

namespace oofem {
/**
 * This class implements the base of a special lattice element following
 * the concepts orginally developed by John Bolander. In this lattice
 * framework, elements consists of rigid bodies connected by a set of axial
 * tranversal and rotational springs at the position of the GP.
 * In the case of an irregular arrangement of lattice elements,
 * the position of the GP is not
 * placed at the midpoint of the element, but the midpoint of the element's midcross-section.
 * There is no way to relate the element geometry to the position of the GP. Instead, this information
 * is part of the input.
 * In this base class common interfaces of derived elements are defined.
 */
class LatticeStructuralElement : public StructuralElement
{
public:
    LatticeStructuralElement(int n, Domain * d);
    virtual ~LatticeStructuralElement();

    virtual IRResultType initializeFrom(InputRecord *ir);

    /**
     * Returns the cross-sectional area of the lattice element.
     * @return Cross-section area.
     */
    virtual double giveArea() { return 0; }

    /**
     * Returns the element length
     * @return Element length.
     */
    virtual double giveLength() { return 0; }

    /**
     * Returns the crack flag
     * @return Crack flag.
     */
    virtual int giveCrackFlag() { return 0; }

    /**
     * Returns the number of crossSection nodes
     * @return Number of crosssection nodes.
     */
    virtual int giveNumberOfCrossSectionNodes() { return 0; }

    /**
     * @return Crack width
     */
    virtual double giveCrackWidth() { return 0; }

    /**
     * @return Crack width
     */
    virtual double giveOldCrackWidth() { return 0; }


    /**
     * Returns the energy dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return dissipation
     */
    virtual double giveDissipation() { return 0; }


    /**
     * Usually computes interpolation function, which is not needed for the lattice elements.
     * However, structural element requires implementation.
     */
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) { return; }


    /**
     * Returns the increment of dissipation computed at the GaussPoint of the element.
     * This function is used for the lattice specific vtk export.
     * @return increment of dissipation
     */
    virtual double giveDeltaDissipation() { return 0; }

    /**
     * Returns the coupling flag.
     * @return couplingFlag
     */

    virtual int giveCouplingFlag() { return 0; }

    /**
     * Returns the coupling numbers
     * @return couplingNumbers.
     */
    virtual void giveCouplingNumbers(IntArray &numbers){return;}

    /**
     * Returns the normal stress.
     * @return normalStress
     */
    virtual double giveNormalStress() { return 0; }

    /**
     * Returns the old normal stress.
     * @return oldNormalStress
     */
    virtual double giveOldNormalStress() { return 0; }


    /**
     * Returns flag indicating if status has been updated
     */
    virtual int hasBeenUpdated() { return 0; }


    /**
     * Gives the GP coordinates
     */
    virtual void  giveGpCoordinates(FloatArray &coords) { return; }
};
} // end namespace oofem
#endif
