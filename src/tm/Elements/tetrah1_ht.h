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

#ifndef tetrah1_ht_h
#define tetrah1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Tetrah1_ht_Name "tetrah1ht"
#define _IFT_Tetrah1_hmt_Name "tetrah1hmt"

namespace oofem {
class FEI3dTetLin;

/**
 * Tetrahedral (3d) element with linear approximation for heat and mass transfer.
 */
class Tetrah1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI3dTetLin interpolation;

public:
    Tetrah1_ht(int n, Domain * d);

    FEInterpolation *giveInterpolation() const override;

    double computeVolumeAround(GaussPoint *gp) override;

    // definition
    const char *giveInputRecordName() const override { return _IFT_Tetrah1_ht_Name; }
    const char *giveClassName() const override { return "Tetrah1_ht"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_tetra_1;}


    int computeNumberOfDofs() override { return ( emode == HeatTransferEM ) ? 4 : 8; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _3dHeat; }

    Interface *giveInterface(InterfaceType t) override;
    int testElementExtension(ElementExtension ext) override
    { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }

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
    double computeSurfaceVolumeAround(GaussPoint *gp, int iEdge) override;
};

class Tetrah1_hmt : public Tetrah1_ht
{
public:
    Tetrah1_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Tetrah1_hmt_Name; }
    const char *giveClassName() const override { return "Tetrah1_hmt"; }
    MaterialMode giveMaterialMode() override { return _3dHeMo; }
};
} // end namespace oofem
#endif // tetrah1_ht_h
