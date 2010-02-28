/* $Header: /home/cvs/bp/oofem/tm/src/quadaxisym1_ht.h,v 1.1 2003/04/14 16:01:40 bp Exp $ */
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

//   ************************************************************************************************
//   *** CLASS QuadAxisym1_ht: Axisymmetric Quad(2d), linear approximation, Heat Transfer element ***
//   ************************************************************************************************

#ifndef quadaxisym1_ht_h
#define quadaxisym1_ht_h

#include "quad1_ht.h"

namespace oofem {
class QuadAxisym1_ht : public Quad1_ht
{
protected:
public:

    // constructor
    QuadAxisym1_ht(int, Domain *, ElementMode em = HeatTransferEM);
    ~QuadAxisym1_ht();                       // destructor

    double                computeVolumeAround(GaussPoint *);
    // definition
    const char *giveClassName() const { return "QuadAxisym1_htElement"; }
    classType                giveClassID() const { return QuadAxisym1_htClass; }

protected:
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    double computeRadiusAt(GaussPoint *gp);
};
} // end namespace oofem
#endif // quadaxisym1_ht_h
