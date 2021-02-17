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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#ifndef qtrplstrslip_h
#define qtrplstrslip_h

#include "sm/Elements/structural2delement.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "sm/Elements/PlaneStress/qtrplstr.h"

#define _IFT_QTrPlaneStress2dSlip_Name "qtrplstrslip"

namespace oofem {
class FEI2dTrQuad;

/**
 * This class implements a quadratic triangular 6-node plane-
 * stress for structural problems with independent macroscopic reinforcement slip field.
 * Each node has 4 degrees of freedom.
 * Currently works only in FE2 setting, i.e., with StructuralSlipFE2Material
 *
 * @author Adam Sciegaj
 */
class QTrPlaneStress2dSlip : public QTrPlaneStress2d
{
public:
    QTrPlaneStress2dSlip(int n, Domain * d);
    virtual ~QTrPlaneStress2dSlip() { }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QTrPlaneStress2dSlip_Name; }
    const char *giveClassName() const override { return "QTrPlaneStress2dSlip"; }
    void initializeFrom(InputRecord &ir) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

protected:
    IntArray aMask = {1,2,5,6,9,10,13,14,17,18,21,22}; //dof numbers of displacement field
    IntArray bMask = {3,4,7,8,11,12,15,16,19,20,23,24}; //dof numbers of slip field

    void giveHomogenizedFields(FloatArray &stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip, const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep);
    void giveSensitivities(FloatMatrix &dStressdEps, FloatMatrix &dStressdS, FloatMatrix &dStressdG, FloatMatrix &dBStressdEps, FloatMatrix &dBStressdS,
        FloatMatrix &dBStressdG, FloatMatrix &dRStressdEps, FloatMatrix &dRStressdS, FloatMatrix &dRStressdG, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

};
} // end namespace oofem
#endif // qtrplstrslip_h