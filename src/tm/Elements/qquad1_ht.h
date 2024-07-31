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

#ifndef qquad1_ht_h
#define qquad1_ht_h

#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_QQuad1_ht_Name "qquad1ht"
#define _IFT_QQuad1_hmt_Name "qquad1hmt"
#define _IFT_QQuad1_mt_Name "qquad1mt"

namespace oofem {
class FEI2dQuadQuad;

/**
 * Quadratic (2d) element with quadratic approximation for heat transfer.
 */
class QQuad1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dQuadQuad interpolation;

public:
    QQuad1_ht(int n, Domain * d);

    FEInterpolation *giveInterpolation() const override;
    double computeVolumeAround(GaussPoint *gp) override;

    const char *giveInputRecordName() const override { return _IFT_QQuad1_ht_Name; }
    const char *giveClassName() const override { return "QQuad1_ht"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_quad_2;}


    int computeNumberOfDofs() override { return 4; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _2dHeat; }
    double giveThicknessAt(const FloatArray &gcoords) override;

    Interface *giveInterface(InterfaceType t) override;
    
    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep) override;

protected:
    void computeGaussPoints() override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
};

/**
 * Class for heat and mass transfer.
 */
class QQuad1_hmt : public QQuad1_ht
{
public:
    QQuad1_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_QQuad1_hmt_Name; }
    const char *giveClassName() const override { return "QQuad1_hmt"; }
    int computeNumberOfDofs() override { return 8; }
    MaterialMode giveMaterialMode() override { return _2dHeMo; }
};

/**
 * Class for mass transfer.
 */
class QQuad1_mt : public QQuad1_ht
{
public:
    QQuad1_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_QQuad1_mt_Name; }
    const char *giveClassName() const override { return "QQuad1_mt"; }
    int computeNumberOfDofs() override { return 4; }
    MaterialMode giveMaterialMode() override { return _2dHeat; }
};
} // end namespace oofem
#endif // quad1_ht_h
