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
    ~LIBeam3dNL2() { }

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
    { computeLumpedMassMatrix(answer, tStep); }
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual int computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    double computeVolumeAround(GaussPoint *gp);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    // definition & identification
    const char *giveClassName() const { return "LIBeam3dNL"; }
    classType giveClassID() const { return LIBeam3dNLClass; }
    IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    integrationDomain giveIntegrationDomain() { return _Line; }
    MaterialMode giveMaterialMode() { return _3dBeam; }

protected:
    // edge load support
    void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer);
    void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);

    void updateYourself(TimeStep *tStep);
    void initForNewStep();
    //void computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int, int)
    { _error("computeBmatrixAt: not implemented"); }
    //int computeGtoLRotationMatrix(FloatMatrix& answer);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj);

    // nonlinear part of geometrical eqs. for i-th component of strain vector.
    // void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int );
    void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeGaussPoints();
    double giveLength();
    //double givePitch ();
    int giveLocalCoordinateSystem(FloatMatrix &answer);

    /**
     * Updates the temporary triad at the centre to the state identified by given solution step.
     * The attribute tempQ is changed to reflect new state and tempQCounter is set to
     * solution step counter to avoid multiple updates.
     * @param tStep Solution step identifying reached state.
     */
    void updateTempQuaternion(TimeStep *tStep);
    /**
     * Compute the temporary curvature at the centre to the state identified by given solution step.
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
