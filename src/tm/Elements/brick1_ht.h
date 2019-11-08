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

#ifndef brick1_ht_h
#define brick1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Brick1_ht_Name "brick1ht"
#define _IFT_Brick1_hmt_Name "brick1hmt"
#define _IFT_Brick1_mt_Name "brick1mt"

namespace oofem {
class FEI3dHexaLin;

/**
 * Brick (3d) elements with linear approximation for heat and mass transfer.
 */
class Brick1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface
{
protected:
    static FEI3dHexaLin interpolation;

public:
    Brick1_ht(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;
    FEInterpolation *giveInterpolation() const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Brick1_ht_Name; }
    const char *giveClassName() const override { return "Brick1_ht"; }

    int computeNumberOfDofs() override { return 8; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _3dHeat; }

    Interface *giveInterface(InterfaceType t) override;
    int testElementExtension(ElementExtension ext) override
    { return ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ); }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

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

/**
 * Class for heat and mass transfer.
 */
class Brick1_hmt : public Brick1_ht
{
public:
    Brick1_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Brick1_hmt_Name; }
    const char *giveClassName() const override { return "Brick1_hmt"; }
    int computeNumberOfDofs() override { return 16; }
    MaterialMode giveMaterialMode() override { return _3dHeMo; }
};

/**
 * Class for mass transfer.
 */
class Brick1_mt : public Brick1_ht
{
public:
    Brick1_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_Brick1_mt_Name; }
    const char *giveClassName() const override { return "Brick1_mt"; }
    int computeNumberOfDofs() override { return 8; }
    MaterialMode giveMaterialMode() override { return _3dHeat; }
};
} // end namespace oofem
#endif // brick1_ht_h
