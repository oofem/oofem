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

//   ****************************
//   *** CLASS CCTPlate in 3d ***
//   ****************************
//   25.5.2010

#ifndef cct3d_h
#define cct3d_h

#include "cct.h"

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
 * This class represent CCT plate element that can be arbitrary
 * orinted in space, in contract to base CCT element (CCTPlate) that is
 * defined in xy plane
 */
class CCTPlate3d : public CCTPlate
{
    /*
     * This class implements an triangular three-node plate CCT finite element.
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
    CCTPlate3d(int, Domain *);                     // constructor
    ~CCTPlate3d() { delete GtoLRotationMatrix; }   // destructor

protected:
    void giveLocalCoordinates(FloatArray &answer, FloatArray &global);
    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3,
                                     double &y1, double &y2, double &y3,
                                     double *z = NULL);

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
    const char *giveClassName() const { return "CCTPlate3d"; }
    classType    giveClassID()   const { return CCTPlate3dClass; }

    virtual int  computeNumberOfDofs(EquationID ut) { return 9; }
    virtual int computeNumberOfL2GDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    const FloatMatrix *computeGtoLRotationMatrix();
    int computeGtoLRotationMatrix(FloatMatrix &answer);

    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    /**
     * Computes the element local (iso) coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    int giveLocalCoordinateSystem(FloatMatrix &answer)
    { _error("cct3d :: giveLocalCoordinateSystem: calling of this function id not allowed"); return 0;}

    void printOutputAt(FILE *file, TimeStep *tStep);
};
} // end namespace oofem
#endif // cct3d_h
