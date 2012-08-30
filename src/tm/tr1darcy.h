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

#ifndef Q1TRBASE_H_
#define Q1TRBASE_H_
#include "transportelement.h"
#include "fei2dtrlin.h"
#include "domain.h"
#include "nodalaveragingrecoverymodel.h"
#include "transportmaterial.h"
   
namespace oofem {
/** Element class for the DarcyFlow engineering model. Linear description of seepage
 * @author Carl Sandström
 */
class Tr1Darcy : public TransportElement, public NodalAveragingRecoveryModelInterface
{
protected:
    int numberOfGaussPoints;
    static FEI2dTrLin interpolation_lin;

public:

    Tr1Darcy(int, Domain *);
    virtual ~Tr1Darcy();
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    virtual FEInterpolation *giveInterpolation() { return &interpolation_lin; }

    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLoadVector(FloatArray &answer, TimeStep *atTime);
    virtual void computeGaussPoints();
    virtual int computeNumberOfDofs(EquationID ut);

    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode) { }; // Empty function. Temporary workaround.

    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep);
    void computeInternalForcesVector(FloatArray &answer, TimeStep *atTime);

    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);

    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    // From NodalAveragingRecoveryModelInterface
    virtual const char *giveClassName() const { return "Tr1Darcy"; };
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type) { return 2; }
    virtual Interface *giveInterface(InterfaceType interface);
};
}

#endif /* Q1TRBASE_H_ */
