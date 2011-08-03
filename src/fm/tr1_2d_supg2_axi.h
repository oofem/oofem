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
 *               Copyright (C) 1993 - 2011   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef tr1_2d_supg2_axi_h
#define tr1_2d_supg2_axi_h

#include "tr1_2d_supg.h"

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
    const FloatArray **vcoords [ 2 ];

    integrationDomain id [ 2 ];
    /**
     * mat[0] reference fluid
     * mat[1] second fluid
     */
    int mat [ 2 ];

public:
    TR1_2D_SUPG2_AXI(int n, Domain *d);
    ~TR1_2D_SUPG2_AXI();

    void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime);
    void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime);
    void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime);
    void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime);
    void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime);
    void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime);
    void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime);

    void updateStabilizationCoeffs(TimeStep *tStep);
    void updateElementForNewInterfacePosition(TimeStep *tStep) { this->updateIntegrationRules(); }
    double computeCriticalTimeStep(TimeStep *tStep);

    // definition
    const char *giveClassName() const { return "TR1_2D_SUPG_AXI"; }
    classType giveClassID() const { return SUPGElementClass; }
    MaterialMode giveMaterialMode() { return _2dAxiFlow; }
    IRResultType initializeFrom(InputRecord *ir);

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself (oofegGraphicContext&);
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    void computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *);
    void updateVolumePolygons(Polygon &referenceFluidPoly, Polygon &secondFluidPoly, int &rfPoints, int &sfPoints,
                               const FloatArray &normal, const double p, bool updFlag);
    double computeVolumeAround(GaussPoint *gp, integrationDomain id, const FloatArray **idpoly);
    double computeRadiusAt(GaussPoint *gp);
    void computeBMtrx(FloatMatrix &answer, GaussPoint *gp);
    void computeNVector(FloatArray &answer, GaussPoint *gp);
    void updateIntegrationRules();
    Material *_giveMaterial(int indx) { return domain->giveMaterial(mat [ indx ]); }

    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag);
    virtual void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                        const FloatArray &normal, const double p, bool updFlag);
    virtual void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                         const FloatArray &normal, const double p, bool updFlag);
    virtual double truncateMatVolume(const Polygon &matvolpoly, double &volume);
    virtual void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool updFlag);
    virtual void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag);
    virtual Element *giveElement() { return this; }
    virtual double computeMyVolume(LEPlic *matInterface, bool updFlag);
};
} // end namespace oofem
#endif // tr1_2d_supg2_axi_h
