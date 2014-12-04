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
    double area;

public:
    TR1_2D_CBS(int n, Domain * aDomain);
    virtual ~TR1_2D_CBS();

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *);
    virtual void computeDiagonalMassMtrx(FloatArray &answer, TimeStep *);
    virtual void computeConvectionTermsI(FloatArray &answer, TimeStep *);
    virtual void computeDiffusionTermsI(FloatArray &answer, TimeStep *);
    virtual void computeDensityRhsVelocityTerms(FloatArray &answer, TimeStep *tStep);
    virtual void computeDensityRhsPressureTerms(FloatArray &answer, TimeStep *tStep);
    virtual void computePrescribedTractionPressure(FloatArray &answer, TimeStep *tStep);
    virtual void computeNumberOfNodalPrescribedTractionPressureContributions(FloatArray &answer, TimeStep *tStep);
    virtual void computePressureLhs(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeCorrectionRhs(FloatArray &answer, TimeStep *tStep);
    virtual double computeCriticalTimeStep(TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "TR1_2D_CBS"; }
    virtual const char *giveInputRecordName() const { return _IFT_TR1_2D_CBS_Name; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual int computeNumberOfDofs();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual void updateYourself(TimeStep *tStep);
    /// Used to check consistency and initialize some element geometry data (area,b,c)
    virtual int checkConsistency();

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);

    virtual int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                      FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                      TimeStep *tStep);

    //<RESTRICTED_SECTION>
    virtual double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag);
    virtual void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                        const FloatArray &normal, const double p, bool updFlag);
    virtual void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                         const FloatArray &normal, const double p, bool updFlag);
    virtual double truncateMatVolume(const Polygon &matvolpoly, double &volume);
    virtual void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd);
    virtual void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag);
    virtual Element *giveElement() { return this; }
    virtual double computeMyVolume(LEPlic *matInterface, bool updFlag);
    virtual double computeCriticalLEPlicTimeStep(TimeStep *tStep) { return 1.e6; }
    //</RESTRICTED_SECTION>

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep);
    // Graphics output
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    //virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) {}
#endif

protected:
    virtual void computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // tr1_2d_cbs_h
