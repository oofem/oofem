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

#ifndef truss3dnl2_h
#define truss3dnl2_h

#include "../sm/Elements/Bars/truss3d.h"
#include "fei3dlinelin.h"

#define _IFT_Truss3dnl2_Name "truss3dnl2"


namespace oofem {
/**
 * This class implements a nonlinear two-node truss bar element for three-dimensional
 * analysis.
 */
class Truss3dnl2 : public Truss3d
{
protected:
  FEICellGeometry*  cellGeometryWrapper;
  // vector of reference coordinates;
  FloatArray X;
  // undeformed length
  double L;
  // init stress factor
  double initStressFactor = 0;


public:
    Truss3dnl2(int n, Domain * d);
    virtual ~Truss3dnl2() { }

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Truss3dnl2_Name; }
    const char *giveClassName() const override { return "Truss3dnl2"; }

    void initializeFrom(InputRecord &ir) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    
protected:
    FloatArray computeStrainVector(GaussPoint *gp, const FloatArray &u);
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, const FloatArray &u);
    void computeInitialStressStiffness(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, const FloatMatrix &B, const FloatArray &d);

    std::pair<double,double> computeDeformedLengthAt(GaussPoint *gp, const FloatArray &d);
    virtual FEICellGeometry* giveCellGeometryWrapper();

    void computeGaussPoints() override;


    double computeLength() override;
    double computeDeformedLength(const FloatArray &d);
    FloatMatrixF<6,6> giveAmatrix();
    FloatMatrixF<3,6> givePmatrix();

};
} // end namespace oofem
#endif // truss3dnl_h
