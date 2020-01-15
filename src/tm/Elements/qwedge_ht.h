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
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of.
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef qwedge_ht_h
#define qwedge_ht_h


#include "tm/Elements/transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"


#define _IFT_QWedge_ht_Name "qwedgeht"
#define _IFT_QWedge_hmt_Name "qwedgehmt"
#define _IFT_QWedge_mt_Name "qwedgemt"

namespace oofem {
class FEI3dWedgeQuad;

/**
 * This class implements a Linear 3d  6 - node thermal finite element.
 */
 class QWedge_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dWedgeQuad interpolation;


public:
    QWedge_ht(int, Domain *);

    double computeVolumeAround(GaussPoint *gp) override;
    FEInterpolation *giveInterpolation() const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QWedge_ht_Name; }
    const char *giveClassName() const override { return "QWedge_ht"; }

    int computeNumberOfDofs() override { return 15; }
    void initializeFrom(InputRecord &ir) override;
    MaterialMode giveMaterialMode() override { return _3dHeat; }

    Interface *giveInterface(InterfaceType t) override;
    int testElementExtension(ElementExtension ext) override
    { return ( ext == Element_EdgeLoadSupport );}

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
class QWedge_hmt : public QWedge_ht
{
public:
    QWedge_hmt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_QWedge_hmt_Name; }
    const char *giveClassName() const override { return "QWedge_hmt"; }
    int computeNumberOfDofs() override { return 24; }
    MaterialMode giveMaterialMode() override { return _3dHeMo; }
};

/**
 * Class for mass transfer.
 */
class QWedge_mt : public QWedge_ht
{
public:
    QWedge_mt(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_QWedge_mt_Name; }
    const char *giveClassName() const override { return "QWedge_mt"; }
    int computeNumberOfDofs() override { return 15; }
    MaterialMode giveMaterialMode() override { return _3dHeat; }
};

 
} // end namespace oofem
#endif
