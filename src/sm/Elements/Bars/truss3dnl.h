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

#ifndef truss3dnl_h
#define truss3dnl_h

#include "../sm/Elements/Bars/truss3d.h"

#define _IFT_Truss3dnl_Name "truss3dnl"


#define _IFT_Truss3dnl_initialStretch "initstretch"

namespace oofem {
/**
 * This class implements a nonlinear two-node truss bar element for three-dimensional
 * analysis.
 */
class Truss3dnl : public Truss3d
{
protected:
  double initialStretch;
public:
    Truss3dnl(int n, Domain * d);
    virtual ~Truss3dnl() { }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Truss3dnl_Name; }
    virtual const char *giveClassName() const { return "Truss3dnl"; }

    IRResultType initializeFrom(InputRecord *ir);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

protected:

    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin = false);
    void computeBlMatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void computeBnlMatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin = false);
    void computeInitialStressStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);

};
} // end namespace oofem
#endif // truss3dnl_h
