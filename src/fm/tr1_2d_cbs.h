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

#ifndef tr1_2d_cbs_h
#define tr1_2d_cbs_h

#include "cbselement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"

#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
//<RESTRICTED_SECTION>
#include "leplic.h"
//</RESTRICTED_SECTION>

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

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
    //double a[3];
    double b [ 3 ];
    double c [ 3 ];
    double area;

public:
    TR1_2D_CBS(int n, Domain *aDomain);
    virtual ~TR1_2D_CBS();

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

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    // definition
    virtual const char *giveClassName() const { return "TR1_2D_CBS"; }
    virtual classType giveClassID() const { return TR1_2D_CBSClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual MaterialMode giveMaterialMode() { return _2dFlow; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int computeNumberOfDofs(EquationID ut);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateYourself(TimeStep *tStep);
    /// Used to check consistency and initialize some element geometry data (area,b,c)
    virtual int checkConsistency();

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                TimeStep *atTime);

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

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
    virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);

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

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);
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
    virtual void computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // tr1_2d_cbs_h
