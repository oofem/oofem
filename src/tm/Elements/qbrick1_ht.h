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

#ifndef qbrick1_ht_h
#define qbrick1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_QBrick1_ht_Name "qbrick1ht"
#define _IFT_QBrick1_hmt_Name "qbrick1hmt"
#define _IFT_QBrick1_mt_Name "qbrick1mt"


namespace oofem {
class FEI3dHexaQuad;

/**
 * Brick (3d) elements with quadratic approximation for heat and mass transfer. Each node has 1 (heat) or 2 (heat+moisture) degrees of freedom.
 * @author Vit Smilauer
 */
class QBrick1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface
{
protected:
    static FEI3dHexaQuad interpolation;

public:
    QBrick1_ht(int n, Domain * d);

    double computeVolumeAround(GaussPoint *gp) override;
    FEInterpolation *giveInterpolation() const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QBrick1_ht_Name; }
    const char *giveClassName() const override { return "QBrick1_ht"; }

    int computeNumberOfDofs() override { return ( emode == HeatTransferEM ) ? 20 : 40; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _3dHeat; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_hexa_2;}


    Interface *giveInterface(InterfaceType t) override;
    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

protected:
    void computeGaussPoints() override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    double computeSurfaceVolumeAround(GaussPoint *gp, int iEdge) override;
};

class QBrick1_hmt : public QBrick1_ht
{
public:
    QBrick1_hmt(int n, Domain * d);

    MaterialMode giveMaterialMode() override { return _3dHeMo; }
    const char *giveInputRecordName() const override { return _IFT_QBrick1_hmt_Name; }
    const char *giveClassName() const override { return "QBrick1_hmt"; }
};

/**
 * Class for mass transfer.
 */
class QBrick1_mt : public QBrick1_ht
{
public:
    QBrick1_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_QBrick1_mt_Name; }
    const char *giveClassName() const override { return "QBrick1_mt"; }
    MaterialMode giveMaterialMode() override { return _3dHeat; }
};
} // end namespace oofem
#endif // qbrick1_ht_h
