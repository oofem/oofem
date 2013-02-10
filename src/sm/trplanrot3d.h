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

#ifndef trplanrot3d_h
#define trplanrot3d_h

#include "trplanrot.h"

namespace oofem {

#ifndef __CHARTENSOR
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

/**
 * This class represent triangular plane stress element with rotational degree of freedom around normal
 * that can be arbitrary oriented in space, in contrast to base TrPlaneStrRot element that is
 * defined in xy-plane.
 *
 * Each node has 3 degrees of freedom.
 *
 * @author Ladislav Svoboda
 * @date 2010-5-25
 */
class TrPlaneStrRot3d : public TrPlaneStrRot
{
protected:
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;

public:
    TrPlaneStrRot3d(int n, Domain *d);
    virtual ~TrPlaneStrRot3d() { delete GtoLRotationMatrix; }

protected:
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    virtual void giveNodeCoordinates(FloatArray &x, FloatArray &y);

    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode);

    friend class TR_SHELL01;

public:
    // definition & identification
    virtual const char *giveClassName() const { return "TrPlaneStrRot3d"; }
    virtual classType giveClassID() const { return TrPlaneStrRot3dClass; }

    virtual int computeNumberOfDofs(EquationID ut) { return 9; }
    virtual int computeNumberOfGlobalDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    const FloatMatrix *computeGtoLRotationMatrix();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
};
} // end namespace oofem
#endif //  trplanrot3d_h
