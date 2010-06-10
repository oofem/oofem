/* $Header: /home/cvs/bp/oofem/sm/src/rershell.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   **************************
//   *** CLASS   RerShell   ***
//   **************************
#ifndef rershell_h
#define rershell_h

#include "cct.h"
#include "layeredcrosssection.h"

namespace oofem {

#ifndef __CHARTENSOR // termitovo
#define __CHARTENSOR
enum CharTensor {
    LocalStrainTensor,
    GlobalStrainTensor,
    LocalCurvatureTensor,
    GlobalCurvatureTensor,

    LocalForceTensor,
    GlobalForceTensor,
    LocalMomentumTensor,
    GlobalMomentumTensor
};
#endif

class RerShell : public CCTPlate
{
    /*
     * This class implements an triangular three-node  shell (CCT+linear plan stress)
     * curved finite element. Each node has 5 degrees of freedom.
     * DESCRIPTION :
     *
     * TASKS :
     *
     * - calculating its B,D,N matrices and dV.
     */

protected:

    double Rx, Ry, Rxy;
    FloatMatrix *GtoLRotationMatrix;

    // Transformation Matrix form GtoL(3,3) is stored
    // at the element level for computation
    // efficiency

public:
    RerShell(int, Domain *);                            // constructor
    ~RerShell()  { delete GtoLRotationMatrix; }         // destructor

    // FloatMatrix*       ComputeConstitutiveMatrixAt (GaussPoint *) ;
    void               computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void               computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    //    void               printOutputAt (TimeStep*) ;
    FloatMatrix *computeGtoLRotationMatrix();
    int                giveLocalCoordinateSystem(FloatMatrix &answer);
    void               giveLocalCoordinates(FloatArray &answer, FloatArray &);
    /**
     * Computes the element local (iso) coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    //
    void giveCharacteristicTensor(FloatMatrix & answer, CharTensor, GaussPoint *, TimeStep *);
    void               printOutputAt(FILE *, TimeStep *);
    //
    // layered cross section support functions
    //
    void               computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                                  GaussPoint *slaveGp, TimeStep *tStep);

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    virtual int            computeNumberOfDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);
    /**
     * @name The element interface required by ZZNodalRecoveryModel
     */
    //@{
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    /**
     * Returns the corresponding element to interface
     */
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    /**
     * Evaluates N matrix (interpolation estimated stress matrix).
     */
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);
    //@}
    /**
     * @name The element interface required by NodalAveragingRecoveryModel
     */
    //@{
    /**
     * Computes the element value in given node.
     * @param answer contains the result
     * @param node element node number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    /**
     * Computes the element value in given side.
     * @param answer contains the result
     * @param node element side number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    //@}

    /**
     * The element interface required by SPRNodalRecoveryModelInterface not implemented
     * because the leas square fit must be made in local shell coordinate system-not implemented
     */

    //
    // io routines
    //
#ifdef __OOFEG
    //   void          drawRawGeometry (oofegGraphicContext&);
    //   void          drawDeformedGeometry(oofegGraphicContext&);
    //     virtual void  drawScalar   (oofegGraphicContext& context);
    // void          drawInternalState (oofegGraphicContext&);
#endif
    //
    //
    // definition & identification
    //
    const char *giveClassName() const { return "RerShell"; }
    classType             giveClassID()          const { return RerShellClass; }
    IRResultType initializeFrom(InputRecord *ir);

protected:
    void               computeBodyLoadVectorAt(FloatArray &answer, Load *, TimeStep *, ValueModeType mode);
    //void               computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
    void               computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void               computeNmatrixAt(GaussPoint *, FloatMatrix &);
    int                computeGtoLRotationMatrix(FloatMatrix &); // giveRotationMatrix () ;
    //  int                computeGtoNRotationMatrix (FloatMatrix&);
    void               computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Triangle; }

    double             giveArea();
};

} // end namespace oofem
#endif // rershell_h
