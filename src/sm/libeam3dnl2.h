/* $Header: /home/cvs/bp/oofem/sm/src/libeam3dnl2.h,v 1.6 2003/04/06 14:08:30 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ************************
//   *** CLASS LIBeam3dNL ***
//   ************************

#ifndef libeam3dnl2_h
#define libeam3dnl2_h

#include "nlstructuralelement.h"
#include "gausspnt.h"

class LIBeam3dNL2 : public NLStructuralElement
{
    /*
     * This class implements a 3-dimensional Linear Isoparametric
     * Mindlin theory beam element, with reduced integration.
     * Geometric nonlinearities are taken into account.
     * Based on Element due to Simo and Vu-Quoc, description taken from
     * Crisfield monograph.
     * Similar to Libeam3dNL, but rotational update is done using quaternions.
     */

private:

    /// initial length
    double l0;
    /// quaternion at the center (last equilibrated)
    FloatArray q;
    /// temporary quaternion at the center
    FloatArray tempQ;
    // curvature at the centre
    // FloatArray  kappa;
    /// time stamp of temporary centre quaternion
    StateCounterType tempQCounter;
    /// reference node
    int referenceNode;
public:

    LIBeam3dNL2(int, Domain *);                       // constructor
    ~LIBeam3dNL2()  { }                               // destructor

    // FloatMatrix*  ComputeConstitutiveMatrixAt (GaussPoint*) ;
    // FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;
    void          computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void          computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
    { computeLumpedMassMatrix(answer, tStep); }
    void          computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    //  int           computeGtoLRotationMatrix (FloatMatrix&);  // giveRotationMatrix () ;
    //  void          computeInitialStressMatrix (FloatMatrix& answer, TimeStep* tStep) ;

    virtual int            computeNumberOfDofs(EquationID ut) { return 12; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    double        computeVolumeAround(GaussPoint *);
    /**
     * Computes the global coordinates from given element's local coordinates.
     * @returns nonzero if successful
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    //
    // definition & identification
    //
    const char *giveClassName() const { return "LIBeam3dNL"; }
    classType             giveClassID()          const { return LIBeam3dNLClass; }
    IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
#endif


    virtual void          computeStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *, int useUpdatedGpRecord = 0);

protected:
    // edge load support
    void  computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *);
    void  giveEdgeDofMapping(IntArray &answer, int) const;
    double        computeEdgeVolumeAround(GaussPoint *, int);
    void          computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
    { computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) ); }
    int   computeLoadLEToLRotationMatrix(FloatMatrix &, int, GaussPoint *);
    int  computeLoadGToLRotationMtrx(FloatMatrix &answer);
    void  updateYourself(TimeStep *tStep);
    void  initForNewStep();
    //void          computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int, int)
    { _error("computeBmatrixAt: not implemented"); }
    //int           computeGtoLRotationMatrix (FloatMatrix& answer);
    contextIOResultType           saveContext(DataStream *stream, ContextMode mode, void *obj);
    contextIOResultType           restoreContext(DataStream *stream, ContextMode mode, void *obj);
    // nonlinear part of geometrical eqs. for i-th component of strain vector.
    // void          computeNLBMatrixAt (FloatMatrix& answer, GaussPoint*, int ) ;
    void          computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void          computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Line; }
    double        giveLength();
    //  double        givePitch () ;
    int           giveLocalCoordinateSystem(FloatMatrix &answer);
    /**
     * Updates the temporary triad at the centre to the state identified by given solution step.
     * The attribute tempQ is changed to reflect new state and tempQCounter is set to
     * solution step conter to avoid multiple updates.
     * @param answer returned triad
     * @param tStep solution step identifying reached state
     */
    void          updateTempQuaternion(TimeStep *tStep);
    /**
     * compute the temporary curvature at the centre to the state identified by given solution step.
     * @param tStep solution step identifying reached state
     */
    void          computeTempCurv(FloatArray &answer, TimeStep *tStep);

    /**
     * Evaluates the S matrix from given vector vec.
     * @param answer assembled result
     * @param vec source vector
     */
    void          computeSMtrx(FloatMatrix &answer, FloatArray &vec);
    /**
     * Evaluates the rotation matrix for large rotations according to Rodrigues formula for given
     * pseudovector psi.
     * @param answer result
     * @param psi pseudovector
     */
    void          computeRotMtrx(FloatMatrix &answer, FloatArray &psi);
    /**
     * Computes X mtrx at given solution state.
     * @param answer returned x matrix.
     * @param tStep determines solution state.
     */
    void          computeXMtrx(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes rotation matrix from given quaternion.
     * @param answer returned rotation matrix.
     * @param q input quaternion
     */

    void          computeRotMtrxFromQuaternion(FloatMatrix &answer, FloatArray &q);
    /**
     * Computes the normalized quatrnion from the given rotation matrix
     * @param answer computed quaternion
     * @param R source rotation matrix
     */
    void          computeQuaternionFromRotMtrx(FloatArray &answer, FloatMatrix &R);
    /**
     * Computes x_21' vector for given solution state.
     * @param answer returned x_21'
     * @param tStep determines solution state.
     */
    void          computeXdVector(FloatArray &answer, TimeStep *tStep);
};

#endif // libeam3dnl2_h
