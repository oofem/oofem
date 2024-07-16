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

#ifndef tr1_2d_cbs_h
#define tr1_2d_cbs_h

#include "cbselement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "primaryfield.h"
//<RESTRICTED_SECTION>
#include "leplic.h"
//</RESTRICTED_SECTION>

///@name Input fields for TR12DCBS
//@{
#define _IFT_TR1_2D_CBS_Name "tr1cbs"
#define _IFT_Tr1CBS_vof "vof"
#define _IFT_Tr1CBS_pvof "pvof"
//@}


namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;
class FEI2dTrLin;

/**
 * This class is the implementation of triangular CFD element with linear (and equal order) interpolation of velocity and pressure fields.
 * Should be used with CBS solution algorithm.
 */
class TR1_2D_CBS : public CBSElement, public SpatialLocalizerInterface, public EIPrimaryFieldInterface,
public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface
//<RESTRICTED_SECTION>
, public LEPlicElementInterface
//</RESTRICTED_SECTION>
{
protected:
    static FEI2dTrLin interp;
    //double a[3];
    double b [ 3 ];
    double c [ 3 ];
    double area = 0.;

public:
    TR1_2D_CBS(int n, Domain * aDomain);

    FEInterpolation *giveInterpolation() const override;

    void computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *tStep) override;
    void computeDiagonalMassMtrx(FloatArray &answer, TimeStep *tStep) override;
    void computeConvectionTermsI(FloatArray &answer, TimeStep *tStep) override;
    void computeDiffusionTermsI(FloatArray &answer, TimeStep *tStep) override;
    void computeDensityRhsVelocityTerms(FloatArray &answer, TimeStep *tStep) override;
    void computeDensityRhsPressureTerms(FloatArray &answer, TimeStep *tStep) override;
    void computePrescribedTractionPressure(FloatArray &answer, TimeStep *tStep) override;
    void computeNumberOfNodalPrescribedTractionPressureContributions(FloatArray &answer, TimeStep *tStep) override;
    void computePressureLhs(FloatMatrix &answer, TimeStep *tStep) override;
    void computeCorrectionRhs(FloatArray &answer, TimeStep *tStep) override;
    double computeCriticalTimeStep(TimeStep *tStep) override;

    // definition
    const char *giveClassName() const override { return "TR1_2D_CBS"; }
    const char *giveInputRecordName() const override { return _IFT_TR1_2D_CBS_Name; }
    MaterialMode giveMaterialMode() override { return _2dFlow; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_triangle_1;}


    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int computeNumberOfDofs() override;
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void updateYourself(TimeStep *tStep) override;
    /// Used to check consistency and initialize some element geometry data (area,b,c)
    int checkConsistency() override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    Interface *giveInterface(InterfaceType) override;

    int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                              const FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                              TimeStep *tStep) override;

    //<RESTRICTED_SECTION>
    double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag) override;
    void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                const FloatArray &normal, const double p, bool updFlag) override;
    void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                 const FloatArray &normal, const double p, bool updFlag) override;
    double truncateMatVolume(const Polygon &matvolpoly, double &volume) override;
    void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd) override;
    void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag) override;
    Element *giveElement() override { return this; }
    double computeMyVolume(LEPlic *matInterface, bool updFlag) override;
    double computeCriticalLEPlicTimeStep(TimeStep *tStep) override { return 1.e6; }
    //</RESTRICTED_SECTION>

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;

    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) override;
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) override;
    int SPRNodalRecoveryMI_giveNumberOfIP() override;
    SPRPatchType SPRNodalRecoveryMI_givePatchType() override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep) override;
    // Graphics output
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    //void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override {}
#endif

protected:
    void computeGaussPoints() override;
    void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // tr1_2d_cbs_h
