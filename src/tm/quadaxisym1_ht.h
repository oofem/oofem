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

#ifndef quadaxisym1_ht_h
#define quadaxisym1_ht_h

#include "quad1_ht.h"

#define _IFT_QuadAxisym1_ht_Name "quadaxisym1ht"
#define _IFT_QuadAxisym1_hmt_Name "quadaxisym1hmt"
#define _IFT_QuadAxisym1_mt_Name "quadaxisym1mt"

namespace oofem {
/**
 * Quadratic axisymmetric element with linear approximation for heat transfer.
 * @todo Use the interpolation classes.
 */
class QuadAxisym1_ht : public Quad1_ht
{
public:
    QuadAxisym1_ht(int n, Domain * d);
    virtual ~QuadAxisym1_ht();

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual double giveThicknessAt(const FloatArray &gcoords);

    virtual const char *giveClassName() const { return "QuadAxisym1_ht"; }

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
    QuadAxisym1_hmt(int n, Domain * d);

    virtual const char *giveClassName() const { return "QuadAxisym1_hmt"; }
};


/**
 * Class for mass transfer.
 */
class QuadAxisym1_mt : public QuadAxisym1_ht
{
public:
    QuadAxisym1_mt(int n, Domain * d);
    virtual const char *giveClassName() const { return "QuadAxisym1_mt"; }
};
} // end namespace oofem
#endif // quadaxisym1_ht_h
