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

#ifndef transportcrosssection_h
#define transportcrosssection_h

#include "crosssection.h"


namespace oofem {
class TransportMaterial;

/**
 * Transort cross-section. It's functionality is to be a wrapper around the material behavior.
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT TransportCrossSection : public CrossSection
{
public:
    /**
     * Constructor. Creates cross section with number n belonging to domain d.
     * @param n Cross section number.
     * @param d Domain for cross section.
     */
    TransportCrossSection(int n, Domain * d) : CrossSection(n, d) { }

    /// Destructor.
    virtual ~TransportCrossSection() { }

    /// @todo Temporary function that hands out the material. Must be removed for future layered support, but input files will still look the same.
    virtual TransportMaterial *giveMaterial() = 0;
};
} // end namespace oofem
#endif // transportcrosssection_h
