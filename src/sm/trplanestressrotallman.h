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

#ifndef trplanestressrotallman_h
#define trplanestressrotallman_h

#include "trplanstrss.h"
#include "fei2dtrquad.h"

///@name Input fields for TrPlaneStrRot
//@{
#define _IFT_TrPlanestressRotAllman_Name "trplanestressrotallman"
//@}

namespace oofem {
/**
 * Class implements an triangular three-node  plane-
 * stress elasticity finite element with independentvertex rotations.
 * Each node has 3 degrees of freedom.
 * For reference, see: Allman, D.J. A compatible triangular element including vertex rotations for plane elas-
ticity analysis. Computers & Structures Vol. 19, No. 1-2, pp. 1-8, 1984.
 */
class TrPlanestressRotAllman : public TrPlaneStress2d
{
protected:
  static FEI2dTrQuad qinterpolation; // quadratic interpolation for constructing shape functons

public:
    TrPlanestressRotAllman(int, Domain *);
    virtual ~TrPlanestressRotAllman() { }

protected:
    virtual void computeGaussPoints();
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual double giveArea();
    virtual int giveApproxOrder() { return 2; }
    void computeLocalCoordinates(FloatArray lxy[6]);
    /** 
     * Computes the stiffness matrix stabilization of zero energy mode (equal rotations)
     * 
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    void computeStiffnessMatrixZeroEnergyStabilization(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
public:
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TrPlanestressRotAllman_Name; }
    virtual const char *giveClassName() const { return "TrPlanestressRotAllman"; }
    virtual classType giveClassID() const { return TrPlanestressRotAllmanClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    virtual integrationDomain giveIntegrationDomain() { return _Triangle; }
    /** Computes the stiffness matrix of receiver. Overloaded to add stabilization of zero-energy mode (equal rotations) */
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual int computeNumberOfDofs(EquationID ut) { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    Interface* giveInterface(InterfaceType interface);
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 2; };
    
    void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    int SPRNodalRecoveryMI_giveNumberOfIP();

};
} // end namespace oofem
#endif //  trplanestressrotallman_h
