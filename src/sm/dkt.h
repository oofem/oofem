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

#ifndef dkt_h
#define dkt_h

#include "mathfem.h"
#include "nlstructuralelement.h"
#include "layeredcrosssection.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "fei2dtrlin.h"
#include "zzerrorestimator.h"

#define _IFT_DKTPlate_Name "dktplate"

namespace oofem {

/**
 * This class implements an triangular Discrete Kirchhoff Theory (DKT) element.
 * This element is a plate element suitable for thin plates, as the traswerse shear strain energy is neglected.
 * The element has only 9 DOFs (w-displacement and rotations along coordinate axes) in each node
 * The derivation starts by assuming quadratic variations of rotation field (fi_x, fi_y)
 * The Kirchhoff hypothesis is imposed at vertices and in side mid nodes
 * The cubic variation of transwerse displacement is assumed along the edges, there is no need to define interpolation for w on the element.
 * As w varies cubically along the edges, its derivative along the edge varies quadratically as the normal rotation along the edge. This allows to satisfy the Kirchhoff hypothesis along the edge. The rotation along the edge is assumed to vary linearly. 
 * This allows to express midside rotations (along and normal to the each edge) as a linear combination of nodal rotations and displacements.
 *
 * Tasks:
 * - calculating its B,D,N matrices and dV.
 */
class DKTPlate : public NLStructuralElement,
    public LayeredCrossSectionInterface, public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface,
    public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface
{
protected:
  /// Element geometry approximation
    static FEI2dTrLin interp_lin;
    double area;

public:
    DKTPlate(int n, Domain *d);
    virtual ~DKTPlate() { }

    virtual FEInterpolation *giveInterpolation() const { return &interp_lin; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    virtual MaterialMode giveMaterialMode()  { return _2dPlate; }
    virtual int giveApproxOrder() { return 1; }
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

protected:
    virtual void computeGaussPoints();
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3,
                                     double &y1, double &y2, double &y3,
                                     double *z = NULL);


    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    //virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp) { answer.resize(0, 0); }
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    //virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const { answer.resize(0); }
    //virtual IntegrationRule *GetSurfaceIntegrationRule(int i) { return NULL; }
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    //virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) { return 0.; }
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    //virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf) { answer.resize(0); }
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);

public:
    // definition & identification
    virtual const char *giveClassName() const { return "DKTPlate"; }
    virtual const char *giveInputRecordName() const { return _IFT_DKTPlate_Name; }
    virtual classType giveClassID() const { return DKTPlateClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    virtual Interface *giveInterface(InterfaceType it);

    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP() { return this->numberOfGaussPoints; }
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();
    // ZZErrorEstimatorInterface
    virtual Element *ZZErrorEstimatorI_giveElement() { return this; }

    // ZZRemeshingCriteriaInterface
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize();
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 1; };

    // layered cross section support functions
    virtual void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                    GaussPoint *slaveGp, TimeStep *tStep);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
#endif
};
} // end namespace oofem
#endif // dkt_h
