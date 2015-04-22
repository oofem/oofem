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

#ifndef warpingcrosssection_h
#define warpingcrosssection_h

#include "simplecrosssection.h"
// #include "structuralmaterial.h"
// #include "floatarray.h"
//  #include "floatmatrix.h"

///@name Input fields for SimpleCrossSection
//@{
#define _IFT_WarpingCrossSection_Name "warpingcs"
#define _IFT_WarpingCrossSection_WarpingNodeNumber "warpingnode"
//@}

namespace oofem {
/**
 * description of warping cross section...
 */
//class OOFEM_EXPORT WarpingCrossSection : public StructuralCrossSection {
class OOFEM_EXPORT WarpingCrossSection : public SimpleCrossSection
{
public:

    WarpingCrossSection(int n, Domain *d) : SimpleCrossSection(n, d) {
        WarpingNodeNumber = 0;
    }

    virtual IRResultType initializeFrom(InputRecord *ir);


    // identification and auxiliary functions

    virtual const char *giveClassName() const {
        return "WarpingCrossSection";
    }

    virtual const char *giveInputRecordName() const {
        return _IFT_WarpingCrossSection_Name;
    }

    int giveWarpingNodeNumber() const {
        return this->WarpingNodeNumber;
    };

protected:
    int WarpingNodeNumber;     // number of the 4rd node
};
} // end namespace oofem
#endif // warpingcrosssection_h
