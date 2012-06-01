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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef traxisym1_ht_h
#define traxisym1_ht_h

#include "tr1_ht.h"

namespace oofem {
/**
 * Triangular axisymmetric element with linear approximation for moisture/heat transfer.
 */
class TrAxisym1_ht : public Tr1_ht
{
public:
    TrAxisym1_ht(int n, Domain *d);
    virtual ~TrAxisym1_ht();

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual const char *giveClassName() const { return "TrAxisym1_htElement"; }
    virtual classType giveClassID() const { return TrAxisym1_htClass; }

protected:
    double computeRadiusAt(GaussPoint *gp);
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual int giveApproxOrder(int unknownIndx) { return 2; }
};
} // end namespace oofem
#endif // traxisym1_ht_h
