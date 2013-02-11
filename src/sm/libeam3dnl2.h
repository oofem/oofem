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

#ifndef libeam3dnl2_h
#define libeam3dnl2_h

#include "nlstructuralelement.h"
#include "gausspnt.h"

namespace oofem {
/**
 * This class implements a 3-dimensional Linear Isoparametric
 * Mindlin theory beam element, with reduced integration.
 * Geometric nonlinearities are taken into account.
 * Based on Element due to Simo and Vu-Quoc, description taken from
 * Crisfield monograph.
 * Similar to Libeam3dNL, but rotational update is done using quaternions.
 */
class LIBeam3dNL2 : public NLStructuralElement
{
private:
    /// Initial length.
    double l0;
    /// Quaternion at the center (last equilibrated).
    FloatArray q;
    /// Temporary quaternion at the center.
    FloatArray tempQ;
    // curvature at the centre
    // FloatArray  kappa;
    /// Time stamp of temporary centre quaternion.
    StateCounterType tempQCounter;
    /// Reference node.
    int referenceNode;

public:
    LIBeam3dNL2(int n, Domain *d);
    virtual ~LIBeam3dNL2() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual int computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    // definition & identification
    virtual const char *giveClassName() const { return "LIBeam3dNL2"; }
    virtual classType giveClassID() const { return LIBeam3dNL2Class; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual integrationDomain giveIntegrationDomain() { return _Line; }
    virtual MaterialMode giveMaterialMode() { return _3dBeam; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);

protected:
    // edge load support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    virtual void updateYourself(TimeStep *tStep);
    virtual void initForNewStep();
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int, int)
    { _error("computeBmatrixAt: not implemented"); }
    //int computeGtoLRotationMatrix(FloatMatrix& answer);

    // nonlinear part of geometrical eqs. for i-th component of strain vector.
    // void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int );
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();
    double giveLength();
    //double givePitch ();
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);

    /**
     * Updates the temporary triad at the centre to the state identified by given solution step.
     * The attribute tempQ is changed to reflect new state and tempQCounter is set to
     * solution step counter to avoid multiple updates.
     * @param tStep Solution step identifying reached state.
     */
    void updateTempQuaternion(TimeStep *tStep);
    /**
     * Compute the temporary curvature at the centre to the state identified by given solution step.
     * @param answer Computed curvature.
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
     * Computes X matrix at given solution state.
     * @param answer Returned X matrix.
     * @param tStep Determines solution state.
     */
    void computeXMtrx(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes rotation matrix from given quaternion.
     * @param answer Returned rotation matrix.
     * @param q Input quaternion
     */
    void computeRotMtrxFromQuaternion(FloatMatrix &answer, FloatArray &q);
    /**
     * Computes the normalized quaternion from the given rotation matrix.
     * @param answer Computed quaternion.
     * @param R Source rotation matrix.
     */
    void computeQuaternionFromRotMtrx(FloatArray &answer, FloatMatrix &R);
    /**
     * Computes x_21' vector for given solution state.
     * @param answer Returned x_21'.
     * @param tStep Determines solution state.
     */
    void computeXdVector(FloatArray &answer, TimeStep *tStep);
};
} // end namespace oofem
#endif // libeam3dnl2_h
