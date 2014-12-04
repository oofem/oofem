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

#ifndef lspacebb_h
#define lspacebb_h

#include "../sm/Elements/3D/lspace.h"

#define _IFT_LSpaceBB_Name "lspacebb"

namespace oofem {
/**
 * Three dimensional brick with linear approximation, suitable for incompressible settings
 * This is achieved by selective integration of deviatoric (full integration) and
 * volumetric (one point) strain contributions. Implemented using bbar technique.
 */
class LSpaceBB  : public LSpace
{
public:
    LSpaceBB(int n, Domain * d);
    virtual ~LSpaceBB() { }

    virtual const char *giveInputRecordName() const { return _IFT_LSpaceBB_Name; }
    virtual const char *giveClassName() const { return "LSpaceBB"; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
};
} // end namespace oofem
#endif // lspacebb_h
