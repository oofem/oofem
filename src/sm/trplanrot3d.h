/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

/* Author: L. Svoboda */

//   ****************************************************************
//   *** CLASS PLANE STRAIN WITH INDEPENDENT ROTATION FIELD in 3d ***
//   ****************************************************************
//   25.5.2010

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
 * This class represent triangular plane stress element with rotational degree of freedon around normal
 * that can be arbitrary orinted in space, in contract to base TrPlaneStrRot element that is
 * defined in xy plane
 */
class TrPlaneStrRot3d : public TrPlaneStrRot
{
    /*
     * This class implements an triangular three-node plane-stress elasticity finite element
     * with independent rotation field.
     * Each node has 3 degrees of freedom.
     * DESCRIPTION :
     *
     * TASKS :
     *
     * - calculating its B,D,N matrices and dV.
     */

protected:
    // Transformation Matrix form GtoL(3,3) is stored
    // at the element level for computation efficiency
    FloatMatrix *GtoLRotationMatrix;

public:
    TrPlaneStrRot3d(int, Domain *);                     // constructor
    ~TrPlaneStrRot3d() { delete GtoLRotationMatrix; }   // destructor

protected:
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    virtual void giveNodeCoordinates(FloatArray &x, FloatArray &y);

    GaussPoint *giveMiddleGaussPoint();

    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);
    int  giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    int  giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode);

    friend  class TR_SHELL01;

public:
    //
    // definition & identification
    //
    const char *giveClassName() const { return "TrPlaneStrRot3d"; }
    classType    giveClassID()   const { return TrPlaneStrRot3dClass; }

    virtual int  computeNumberOfDofs(EquationID ut) { return 9; }
    virtual int computeNumberOfL2GDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    const FloatMatrix *computeGtoLRotationMatrix();
    int computeGtoLRotationMatrix(FloatMatrix &answer);

    void printOutputAt(FILE *file, TimeStep *tStep);
};
} // end namespace oofem
#endif //  trplanrot3d_h
