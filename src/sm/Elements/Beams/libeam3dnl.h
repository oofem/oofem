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

#ifndef libeam3dnl_h
#define libeam3dnl_h

#include "../sm/Elements/nlstructuralelement.h"

///@name Input fields for LIBeam3dNL
//@{
#define _IFT_LIBeam3dNL_Name "libeam3dnl"
#define _IFT_LIBeam3dNL_refnode "refnode"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 * Geometric nonlinearities are taken into account.
 * Based on Element due to Simo and Vu-Quoc, description taken from
 * Crisfield monograph.
 */
class LIBeam3dNL : public NLStructuralElement
{
private:
    /// Initial length.
    double l0;
    /// Last equilibrium triad at the centre.
    FloatMatrix tc;
    // curvature at the centre
    // FloatArray  kappa;
    /// Temporary triad at the centre.
    FloatMatrix tempTc;
    /// Time stamp of temporary centre triad.
    StateCounterType tempTcCounter;
    /// Reference node.
    int referenceNode;

public:
    LIBeam3dNL(int n, Domain * d);
    virtual ~LIBeam3dNL() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    //int computeGtoLRotationMatrix(FloatMatrix &answer);
    //void computeInitialStressMatrix(FloatMatrix& answer, TimeStep* tStep);

    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LIBeam3dNL_Name; }
    virtual const char *giveClassName() const { return "LIBeam3dNL"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
#endif

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual integrationDomain giveIntegrationDomain() const { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _3dBeam; }


protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    virtual void updateYourself(TimeStep *tStep);
    virtual void initForNewStep();
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)
    { OOFEM_ERROR("not implemented"); }
    //int computeGtoLRotationMatrix(FloatMatrix& answer);

    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeGaussPoints();
    virtual double computeLength();
    //double givePitch();
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);
    /**
     * Updates the temporary triad at the centre to the state identified by given solution step.
     * The attribute tempTc is changed to reflect new state and tempTcCounter is set to
     * solution step counter to avoid multiple updates.
     * @param tStep Solution step identifying reached state.
     */
    void updateTempTriad(TimeStep *tStep);
    /**
     * Compute the temporary curvature at the centre to the state identified by given solution step.
     * @param answer
     * @param tStep Solution step identifying reached state.
     */
    void computeTempCurv(FloatArray &answer, TimeStep *tStep);

    /**
     * Evaluates the S matrix from given vector vec.
     * @param answer Assembled result.
     * @param vec Source vector.
     */
    void computeSMtrx(FloatMatrix &answer, FloatArray &vec);
    /**
     * Evaluates the rotation matrix for large rotations according to Rodrigues formula for given
     * pseudovector psi.
     * @param answer Result.
     * @param psi Pseudovector.
     */
    void computeRotMtrx(FloatMatrix &answer, FloatArray &psi);
    /**
     * Computes X mtrx at given solution state.
     * @param answer Returned x matrix.
     * @param tStep Determines solution state.
     */
    void computeXMtrx(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes x_21' vector for given solution state.
     * @param answer Returned x_21'.
     * @param tStep Determines solution state.
     */
    void computeXdVector(FloatArray &answer, TimeStep *tStep);
};
} // end namespace oofem
#endif // libeam3dnl_h
