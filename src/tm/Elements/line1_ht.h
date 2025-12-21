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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef line1_ht_h
#define line1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Line1_hmt_Name "line1hmt"
#define _IFT_Line1_ht_Name "line1ht"
#define _IFT_Line1_mt_Name "line1mt"

namespace oofem {
class FEI3dLineLin;

/**
 * Two node element for heat or moisture transport with linear interpolation.
 */
class Line1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI3dLineLin interp;

public:
    Line1_ht(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;

    // definition
    const char *giveInputRecordName() const override { return _IFT_Line1_ht_Name; }
    const char *giveClassName() const override { return "Line1_htElement"; }

    int computeNumberOfDofs() override { return ( emode == HeatMass1TransferEM ) ? 4 : 2; }
    void initializeFrom(InputRecord &ir, int priority) override;
    MaterialMode giveMaterialMode() override { return _3dHeat; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}


    Interface *giveInterface(InterfaceType t) override;

    FEInterpolation *giveInterpolation() const override;

#ifdef __OOFEG
    // Graphics output
    //void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override {}
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override {}
#endif

protected:
    void computeGaussPoints() override;
};

/**
 * Class for mass transfer.
 */
class Line1_mt : public Line1_ht
{
public:
    Line1_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Line1_mt_Name; }
    const char *giveClassName() const override { return "Line1_mt"; }
    MaterialMode giveMaterialMode() override { return _3dHeat; }
};


/**
 * Class for heat and mass transfer.
 */
class Line1_hmt : public Line1_ht
{
public:
    Line1_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Line1_hmt_Name; }
    const char *giveClassName() const override { return "Line1_hmt"; }
    MaterialMode giveMaterialMode() override { return _3dHeMo; }
};
} // end namespace oofem
#endif // line1_ht_h
