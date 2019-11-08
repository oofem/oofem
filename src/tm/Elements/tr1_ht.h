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

#ifndef tr1_ht_h
#define tr1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Tr1_hmt_Name "tr1hmt"
#define _IFT_Tr1_ht_Name "tr1ht"
#define _IFT_Tr1_mt_Name "tr1mt"

namespace oofem {
class FEI2dTrLin;

/**
 * Triangle (2d) element with linear approximation for heat transfer.
 * @todo Use the interpolation classes.
 */
class Tr1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dTrLin interp;

public:
    Tr1_ht(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;

    // definition
    const char *giveInputRecordName() const override { return _IFT_Tr1_ht_Name; }
    const char *giveClassName() const override { return "Tr1_htElement"; }

    int computeNumberOfDofs() override { return ( emode == HeatMass1TransferEM ) ? 6 : 3; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _2dHeat; }
    double giveThicknessAt(const FloatArray &gcoords) override;

    Interface *giveInterface(InterfaceType t) override;

    FEInterpolation *giveInterpolation() const override;

#ifdef __OOFEG
    // Graphics output
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override {}
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override {}
#endif

protected:
    void computeGaussPoints() override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
};

/**
 * Class for mass transfer.
 */
class Tr1_mt : public Tr1_ht
{
public:
    Tr1_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Tr1_mt_Name; }
    const char *giveClassName() const override { return "Tr1_mt"; }
    MaterialMode giveMaterialMode() override { return _2dHeat; }
};


/**
 * Class for heat and mass transfer.
 */
class Tr1_hmt : public Tr1_ht
{
public:
    Tr1_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Tr1_hmt_Name; }
    const char *giveClassName() const override { return "Tr1_hmt"; }
    MaterialMode giveMaterialMode() override { return _2dHeMo; }
};
} // end namespace oofem
#endif // tr1_ht_h
