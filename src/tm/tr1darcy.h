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

#ifndef tr1darcy_h_
#define tr1darcy_h_

#include "transportelement.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Tr1Darcy_Name "tr1darcy"

namespace oofem {
class FEI2dTrLin;

/**
 * Element class for the DarcyFlow engineering model. Linear description of seepage
 * @author Carl Sandström
 * @author Mikael Öhman
 */
class Tr1Darcy : public TransportElement, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dTrLin interpolation_lin;

public:
    Tr1Darcy(int, Domain *);
    virtual ~Tr1Darcy();
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual FEInterpolation *giveInterpolation() const;

    virtual MaterialMode giveMaterialMode() { return _2dHeat; } ///@todo This isn't actually correct.

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep);

    virtual void computeGaussPoints();
    virtual int computeNumberOfDofs();

    void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode, int indx);
    void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);

    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual double giveThicknessAt(const FloatArray &gcoords);

    // From NodalAveragingRecoveryModelInterface
    virtual const char *giveInputRecordName() const { return _IFT_Tr1Darcy_Name; }
    virtual const char *giveClassName() const { return "Tr1Darcy"; }
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual Interface *giveInterface(InterfaceType interface);
};
}

#endif // tr1darcy_h_
