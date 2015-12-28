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

#include "Elements/PlaneStress/planstrss.h"
#include "ErrorEstimators/directerrorindicatorrc.h"
#include "ErrorEstimators/huertaerrorestimator.h"
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
        LocalMomentumTensor,
        GlobalMomentumTensor
    };

public:
    LinQuad3DPlaneStress(int n, Domain * d);
    virtual ~LinQuad3DPlaneStress();

    virtual Interface *giveInterface(InterfaceType it);
    virtual FEICellGeometry* giveCellGeometryWrapper();
    void computeLocalNodalCoordinates(std::vector< FloatArray > &lxy);

    virtual int computeNumberOfDofs() { return 8; }
    virtual int computeNumberOfGlobalDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;


    const FloatMatrix *computeGtoLRotationMatrix();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    void giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LinQuad3DPlaneStress_Name; }
    virtual const char *giveClassName() const { return "LinQuad3DPlaneStress"; }
};
} // end namespace oofem
#endif // linquad3d_planestress_h
