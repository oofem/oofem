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

#ifndef quad1_ht_h
#define quad1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Quad1_ht_Name "quad1ht"
#define _IFT_Quad1_hmt_Name "quad1hmt"
#define _IFT_Quad1_mt_Name "quad1mt"

namespace oofem {
class FEI2dQuadLin;

/**
 * Quadratic (2d) element with linear approximation for heat transfer.
 */
class Quad1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadLin interpolation;

public:
    Quad1_ht(int n, Domain * d);

    FEInterpolation *giveInterpolation() const override;
    double computeVolumeAround(GaussPoint *gp) override;

    const char *giveInputRecordName() const override { return _IFT_Quad1_ht_Name; }
    const char *giveClassName() const override { return "Quad1_ht"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_quad_1;}


    //int computeNumberOfDofs() override { return ( emode == HeatTransferEM ) ? 4 : 8; }
    int computeNumberOfDofs() override { return 4; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _2dHeat; }
    double giveThicknessAt(const FloatArray &gcoords) override;

    Interface *giveInterface(InterfaceType t) override;

#ifdef __OOFEG
    // Graphics output
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override {}
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override {}
#endif

protected:
    void computeGaussPoints() override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
};

/**
 * Class for heat and mass transfer.
 */
class Quad1_hmt : public Quad1_ht
{
public:
    Quad1_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Quad1_hmt_Name; }
    const char *giveClassName() const override { return "Quad1_hmt"; }
    int computeNumberOfDofs() override { return 8; }
    MaterialMode giveMaterialMode() override { return _2dHeMo; }
};

/**
 * Class for mass transfer.
 */
class Quad1_mt : public Quad1_ht
{
public:
    Quad1_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Quad1_mt_Name; }
    const char *giveClassName() const override { return "Quad1_mt"; }
    int computeNumberOfDofs() override { return 4; }
    MaterialMode giveMaterialMode() override { return _2dHeat; }
};
} // end namespace oofem
#endif // quad1_ht_h
