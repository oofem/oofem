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

#ifndef qplanstrssslip_h
#define qplanstrssslip_h

#include "sm/Elements/structural2delement.h"
#include "sm/Elements/PlaneStress/qplanstrss.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_QPlaneStress2dSlip_Name "qplanestress2dslip"

namespace oofem {
class FEI2dQuadQuad;

/**
 * This class implements a quadratic isoparametric 8-node quadrilateral plane-
 * stress elasticity finite element with independent slip field. Each node has 4 degrees of freedom.
 * Currently works only in FE2 setting, i.e., with StructuralSlipFE2Material.
 *
 * @author Adam Sciegaj
 */
class QPlaneStress2dSlip : public QPlaneStress2d
{
protected:
    static FEI2dQuadQuad interpolation;

public:
    QPlaneStress2dSlip(int n, Domain * d);
    virtual ~QPlaneStress2dSlip() { }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_QPlaneStress2dSlip_Name; }
    const char *giveClassName() const override { return "QPlaneStress2dSlip"; }
    void initializeFrom(InputRecord &ir) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep) override;


protected:
    IntArray aMask = {1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30}; //dof numbers of displacement field
    IntArray bMask = {3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32}; //dof numbers of slip field

    void giveHomogenizedFields(FloatArray &stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip, const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep);
    void giveSensitivities(FloatMatrix &dStressdEps, FloatMatrix &dStressdS, FloatMatrix &dStressdG, FloatMatrix &dBStressdEps, FloatMatrix &dBStressdS,
                           FloatMatrix &dBStressdG, FloatMatrix &dRStressdEps, FloatMatrix &dRStressdS, FloatMatrix &dRStressdG, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

};
} // end namespace oofem
#endif // qplanstrssslip_h
