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

#ifndef libeam3dnl_h
#define libeam3dnl_h

#include "nlstructuralelement.h"
#include "gausspoint.h"

///@name Input fields for LIBeam3dNL
//@{
#define _IFT_LIBeam3dNL_Name "libeam3dNL"
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
    LIBeam3dNL(int n, Domain *d);
    virtual ~LIBeam3dNL() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    //int computeGtoLRotationMatrix(FloatMatrix &answer);
    //void computeInitialStressMatrix(FloatMatrix& answer, TimeStep* tStep);

    virtual int computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    // definition & identification
    virtual const char *giveClassName() const { return "LIBeam3dNL"; }
    virtual classType giveClassID() const { return LIBeam3dNLClass; }
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
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)
    { _error("computeBmatrixAt: not implemented"); }
    //int computeGtoLRotationMatrix(FloatMatrix& answer);

    // nonlinear part of geometrical eqs. for i-th component of strain vector.
    //void computeNLBMatrixAt (FloatMatrix& answer, GaussPoint*, int ) ;
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();
    double giveLength();
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
