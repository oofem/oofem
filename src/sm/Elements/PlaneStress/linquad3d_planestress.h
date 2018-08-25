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

#ifndef linquad3d_planestress_h
#define linquad3d_planestress_h

#include "sm/Elements/PlaneStress/planstrss.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"
#include "sm/ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

#define _IFT_LinQuad3DPlaneStress_Name "linquad3dplanestress"

namespace oofem {
class FEI2dQuadLin;

/**
 * This class implements an isoparametric four-node quadrilateral plane-
 * stress elasticity finite element. Each node has 2 degrees of freedom.
 */
class LinQuad3DPlaneStress : public PlaneStress2d
{
protected:
    /// Local vertex coordinates
    std :: vector< FloatArray > lc; 
    /**
     * Transformation Matrix form GtoL(3,3) is stored
     * at the element level for computation efficiency
     */
    FloatMatrix *GtoLRotationMatrix;

    enum CharTensor {
        LocalStrainTensor,
        GlobalStrainTensor,
        LocalCurvatureTensor,
        GlobalCurvatureTensor,

        LocalForceTensor,
        GlobalForceTensor,
        LocalMomentTensor,
        GlobalMomentTensor
    };

public:
    LinQuad3DPlaneStress(int n, Domain * d);
    virtual ~LinQuad3DPlaneStress();

    Interface *giveInterface(InterfaceType it) override;
    FEICellGeometry* giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates(std::vector< FloatArray > &lxy);

    int computeNumberOfDofs() override { return 8; }
    int computeNumberOfGlobalDofs() override { return 12; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix(FloatMatrix &answer) override;
    int computeLoadGToLRotationMtrx(FloatMatrix &answer) override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LinQuad3DPlaneStress_Name; }
    const char *giveClassName() const override { return "LinQuad3DPlaneStress"; }
};
} // end namespace oofem
#endif // linquad3d_planestress_h
