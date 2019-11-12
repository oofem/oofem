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

#ifndef traxisym1_ht_h
#define traxisym1_ht_h

#include "tm/Elements/tr1_ht.h"

#define _IFT_TrAxisym1_ht_Name "traxisym1ht"
#define _IFT_TrAxisym1_hmt_Name "traxisym1hmt"
#define _IFT_TrAxisym1_mt_Name "traxisym1mt"

namespace oofem {
/**
 * Triangular axisymmetric element with linear approximation for moisture/heat transfer.
 */
class TrAxisym1_ht : public Tr1_ht
{
public:
    TrAxisym1_ht(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;
    const char *giveInputRecordName() const override { return _IFT_TrAxisym1_ht_Name; }
    const char *giveClassName() const override { return "TrAxisym1_htElement"; }

    double giveThicknessAt(const FloatArray &gcoords) override;

protected:
    double computeRadiusAt(GaussPoint *gp);
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    int giveApproxOrder(int unknownIndx) override { return 2; }
};
/**
 * Same as TrAxisym1_ht but for heat+mass transfer.
 */
class TrAxisym1_hmt : public TrAxisym1_ht
{
public:
    TrAxisym1_hmt(int n, Domain * d);

    const char *giveClassName() const override { return "TrAxisym1_hmt"; }
};


/**
 * Class for mass transfer.
 */
class TrAxisym1_mt : public TrAxisym1_ht
{
public:
    TrAxisym1_mt(int n, Domain * d);
    const char *giveClassName() const override { return "TrAxisym1_mt"; }
};
} // end namespace oofem
#endif // traxisym1_ht_h
