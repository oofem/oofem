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

#ifndef quadaxisym1_ht_h
#define quadaxisym1_ht_h

#include "quad1_ht.h"

namespace oofem {
/**
 * Quadratic axisymmetric element with linear approximation for heat transfer.
 * @todo Use the interpolation classes.
 */
class QuadAxisym1_ht : public Quad1_ht
{
public:
    QuadAxisym1_ht(int n, Domain *d);
    virtual ~QuadAxisym1_ht();

    virtual double computeVolumeAround(GaussPoint *gp);

    virtual const char *giveClassName() const { return "QuadAxisym1_ht"; }
    virtual classType giveClassID() const { return QuadAxisym1_htClass; }

protected:
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual double computeRadiusAt(GaussPoint *gp);
    virtual int giveApproxOrder(int unknownIndx) { return 2; }
};

/**
 * Same as QuadAxisym1_ht but for heat+mass transfer.
 */
class QuadAxisym1_hmt : public QuadAxisym1_ht
{
public:
    QuadAxisym1_hmt(int n, Domain *d);

    virtual const char *giveClassName() const { return "QuadAxisym1_hmt"; }
    virtual classType giveClassID() const { return QuadAxisym1_hmtClass; }
};


/**
 * Class for mass transfer.
 */
class QuadAxisym1_mt : public QuadAxisym1_ht
{
public:
    QuadAxisym1_mt(int n, Domain *d);
    virtual const char *giveClassName() const { return "QuadAxisym1_mt"; }
    virtual classType giveClassID() const { return QuadAxisym1_mtClass; }
};
} // end namespace oofem
#endif // quadaxisym1_ht_h
