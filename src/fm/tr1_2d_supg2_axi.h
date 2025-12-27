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

#ifndef tr1_2d_supg2_axi_h
#define tr1_2d_supg2_axi_h

#include "tr1_2d_supg.h"

#define _IFT_TR1_2D_SUPG2_AXI_Name "tr1supg2axi"

namespace oofem {
/**
 * Class representing 2d linear axisymmetric triangular element
 * for solving incompressible fluid with SUPG solver
 *
 * This class is similar to TR1_2D_SUPG_AXI, but difference is in handling
 * multiple fluids. This class uses the interface position within an element to
 * perform an integration for each fluid separately when evaluating contributing terms.
 * It does not rely on rule of mixture which interpolates the properties using VOF value,
 * but uses separate integration on each fluid volume.
 */
class TR1_2D_SUPG2_AXI : public TR1_2D_SUPG
{
protected:
    /**
     * myPoly[0] Occupied by reference fluid.
     * myPoly[1] Occupied by second fluid (air).
     */
    Polygon myPoly [ 2 ];
    std::vector< FloatArray > vcoords [ 2 ];

    integrationDomain id [ 2 ];
    /**
     * mat[0] reference fluid
     * mat[1] second fluid
     */
    int mat [ 2 ];

public:
    TR1_2D_SUPG2_AXI(int n, Domain * d);

    void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override;
    void computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep) override;
    void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep) override;
    void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep) override;
    void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep) override;
    void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep) override;

    void updateStabilizationCoeffs(TimeStep *tStep) override;
    void updateElementForNewInterfacePosition(TimeStep *tStep) override { this->updateIntegrationRules(); }
    double computeCriticalTimeStep(TimeStep *tStep) override;

    // definition
    const char *giveClassName() const override { return "TR1_2D_SUPG2_AXI"; }
    const char *giveInputRecordName() const override { return _IFT_TR1_2D_SUPG2_AXI_Name; }
    MaterialMode giveMaterialMode() override { return _2dAxiFlow; }
    void initializeFrom(InputRecord &ir, int priority) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep) override;
    // Graphics output
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) {}
#endif

protected:
    void computeGaussPoints() override;
    void computeDeviatoricStress(FloatArray &answer, const FloatArray &eps, GaussPoint *gp, TimeStep *tStep) override;
    void computeTangent(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void updateVolumePolygons(Polygon &referenceFluidPoly, Polygon &secondFluidPoly, int &rfPoints, int &sfPoints,
                              const FloatArray &normal, const double p, bool updFlag);
    double computeVolumeAroundID(GaussPoint *gp, integrationDomain id, const std::vector< FloatArray > &idpoly);
    double computeRadiusAt(GaussPoint *gp);
    void computeBMtrx(FloatMatrix &answer, GaussPoint *gp);
    void computeNVector(FloatArray &answer, GaussPoint *gp);
    void updateIntegrationRules();
    Material *_giveMaterial(int indx) { return domain->giveMaterial(mat [ indx ]); }

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag) override;
    void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                const FloatArray &normal, const double p, bool updFlag) override;
    void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                 const FloatArray &normal, const double p, bool updFlag) override;
    double truncateMatVolume(const Polygon &matvolpoly, double &volume) override;
    void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool updFlag) override;
    void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag) override;
    Element *giveElement() override { return this; }
    double computeMyVolume(LEPlic *matInterface, bool updFlag) override;
};
} // end namespace oofem
#endif // tr1_2d_supg2_axi_h
