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

#ifndef cct_h
#define cct_h

#include "sm/Elements/nlstructuralelement.h"
#include "sm/CrossSections/layeredcrosssection.h"
#include "sm/ErrorEstimators/zzerrorestimator.h"
#include "mathfem.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_CCTPlate_Name "cctplate"

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements an triangular three-node plate CCT finite element.
 * Each node has 3 degrees of freedom.
 *
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class CCTPlate : public StructuralElement,
    public LayeredCrossSectionInterface, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public ZZErrorEstimatorInterface
{
protected:
    static FEI2dTrLin interp_lin;
    //static FEI2dTrRot interp_rot;

    double area;

public:
    CCTPlate(int n, Domain *d);
    virtual ~CCTPlate() { }

    FEInterpolation *giveInterpolation() const override;
    FEInterpolation *giveInterpolation(DofIDItem id) const override;

    MaterialMode giveMaterialMode() override { return _2dPlate; }
    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    // overloaded to take into account possible element local cs (in derived cct3d)
    double computeArea() override;

protected:
    void computeGaussPoints() override;
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;


    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3,
                                     double &y1, double &y2, double &y3,
                                     double &z1, double &z2, double &z3);


    //virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp) { answer.clear(); }
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    //virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const { answer.clear(); }
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    //virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) { return 0.; }
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp) override;

public:
    // definition & identification
    const char *giveClassName() const override { return "CCTPlate"; }
    const char *giveInputRecordName() const override { return _IFT_CCTPlate_Name; }
    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 9; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp) override;

    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override;
    double computeVolumeAround(GaussPoint *gp) override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { computeLumpedMassMatrix(answer, tStep); }

    Interface *giveInterface(InterfaceType it) override;

    bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override { return 1; }
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    // layered cross section support functions
    void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                    GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) override;
    // boundary load support
    void computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif
};
} // end namespace oofem
#endif // cct_h
