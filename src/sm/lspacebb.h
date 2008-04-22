/* $Header: /home/cvs/bp/oofem/sm/src/lspace.h,v 1.6 2003/05/19 13:04:00 bp Exp $ */
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

//   *******************************************************************************
//   ***      CLASS 3D Brick linear Element (with bBar for incompressibility   ) ***
//   *******************************************************************************

#ifndef lspacebb_h
#define lspacebb_h

#include "lspace.h"

/**
 * Three dimensional brick with linear approximation, suitable for incompressible settings
 * This is achieved by selective integration of deviatoric (full integration) and
 * volumetric (one point) strain contributions. Implemented using bbar technique.
 */
class LSpaceBB  : public LSpace
{
protected:
public:
    LSpaceBB(int, Domain *);                        // constructor
    ~LSpaceBB()  { }                                // destructor

    const char *giveClassName() const { return "LSpaceBB"; }
    classType             giveClassID()  const { return LSpaceBBClass; }

protected:
    void               computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
};

#endif // lspacebb_h
