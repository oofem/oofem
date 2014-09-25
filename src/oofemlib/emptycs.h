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

#ifndef emptycs_h
#define emptycs_h

#include "crosssection.h"

#define _IFT_EmptyCS_Name "emptycs"

namespace oofem {
/**
 * Empty cross section model which doesn't have any material models.
 */
class OOFEM_EXPORT EmptyCS : public CrossSection
{
public:
    /**
     * Constructor. Creates cross section with number n belonging to domain d.
     * @param n Cross section number.
     * @param d Domain for cross section.
     */
    EmptyCS(int n, Domain * d);
    /// Destructor.
    virtual ~EmptyCS();

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) { return 1; }
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) { return 1; }
    virtual int estimatePackSize(DataStream &buff, GaussPoint *ip) { return 0; }

    virtual const char *giveClassName() const { return "EmptyCS"; }
    virtual const char *giveInputRecordName() const { return _IFT_EmptyCS_Name; }
};
} // end namespace oofem
#endif // emptycs_h
