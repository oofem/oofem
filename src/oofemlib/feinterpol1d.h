/* $Header: /home/cvs/bp/oofem/oofemlib/src/feinterpol1d.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *****************************
//   *** CLASS FEInterpolation1d ***
//   *****************************


#ifndef feinterpol1d_h
#define feinterpol1d_h

#include "feinterpol.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"

namespace oofem {
/**
 * Class representing a general abstraction for finite element interpolation class.
 */
class FEInterpolation1d : public FEInterpolation
{
protected:

public:
    FEInterpolation1d(int o) : FEInterpolation(o) { }
    /**
     * Returns number of spatial dimensions
     */
    int const giveNsd() { return 1; }
};
} // end namespace oofem
#endif // feinterpol1d_h






