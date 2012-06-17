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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef tr1_2d_supg_h
#define tr1_2d_supg_h

#include "supgelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "fei2dtrlin.h"
#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "leplic.h"
#include "levelsetpcs.h"

namespace oofem {
/**
 * Class representing 2d linear triangular element
 * for solving incompressible fluid with SUPG solver
 *
 * This class is similar to TR1_2D_SUPG2, but difference is in handling
 * multiple fluids. This class uses rule of mixture which interpolates the
 * properties using VOF value, requiring the use of twofluidmaterial class
 * as material model for this situation.
 */
class TR1_2D_SUPG : public SUPGElement,
    public SpatialLocalizerInterface, public EIPrimaryFieldInterface,
    public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public LEPlicElementInterface, public LevelSetPCSElementInterface
{
protected:
    static FEI2dTrLin interp;

    //double a[3];
    double b [ 3 ];
    double c [ 3 ];
    double area;

public:
    TR1_2D_SUPG(int n, Domain *d);
    virtual ~TR1_2D_SUPG();

    virtual FEInterpolation *giveInterpolation() { return &interp; }

    virtual void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime);
    virtual void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime);
    virtual void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime) {
        answer.resize(3, 6);
        answer.zero();
    }
    virtual void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime) {
        answer.resize(3);
        answer.zero();
    }
    virtual void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime);

    virtual void computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime);
    virtual void computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime);
    virtual void computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *atTime);

    virtual void updateStabilizationCoeffs(TimeStep *tStep);
    virtual double computeCriticalTimeStep(TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "TR1_2D_SUPG"; }
    virtual classType giveClassID() const { return TR1_2D_SUPGClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int computeNumberOfDofs(EquationID ut);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateYourself(TimeStep *tStep);
    /// Used to check consistency and initialize some element geometry data (area,b,c).
    virtual int checkConsistency();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual int  EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                       FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                       TimeStep *atTime);

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
    virtual double  computeVolumeAround(GaussPoint *aGaussPoint);
    virtual double computeCriticalLEPlicTimeStep(TimeStep *tStep);

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    //virtual void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    virtual void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual double LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep);
    virtual void LS_PCS_computedN(FloatMatrix &answer);
    virtual double LS_PCS_computeVolume() { return area; }
    virtual double LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep);
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself (oofegGraphicContext&);
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

protected:
    virtual void giveLocalVelocityDofMap(IntArray &map);
    virtual void giveLocalPressureDofMap(IntArray &map);

    void computeNMtrx(FloatArray &answer, GaussPoint *gp);
    virtual void computeGaussPoints();

    virtual void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    virtual void initGeometry();
};
} // end namespace oofem
#endif // tr1_2d_supg_h
