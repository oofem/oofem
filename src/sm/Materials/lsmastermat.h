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

#ifndef lsmastermat_h
#define lsmastermat_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for LargeStrainMasterMaterial
//@{
#define _IFT_LargeStrainMasterMaterial_Name "lsmastermat"
#define _IFT_LargeStrainMasterMaterial_m "m"
#define _IFT_LargeStrainMasterMaterial_slaveMat "slavemat"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * Large strain master material.
 * Stress and stiffness are computed from small strain(slaveMat) material model
 * using a strain tensor from the Seth-Hill strain tensors family (depends on parameter m,
 * m = 0 logarithmic strain,m = 1 Green-Lagrange strain ...)
 * then stress and stiffness are transformed 2.PK stress and appropriate stiffness
 */
class LargeStrainMasterMaterial : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial = nullptr;

    /// 'slave' material model number.
    int slaveMat = 0;
    /// Specifies the strain tensor.
    double m = 0.;

public:
    LargeStrainMasterMaterial(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveInputRecordName() const override { return _IFT_LargeStrainMasterMaterial_Name; }
    const char *giveClassName() const override { return "LargeStrainMasterMaterial"; }

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    FloatMatrixF<9,9> give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                            GaussPoint * gp,
                                            TimeStep * tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &, GaussPoint *, TimeStep *) const override
    { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    FloatArrayF<9> giveFirstPKStressVector_3d(const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep) const override;

    /// transformation matrices
    FloatMatrixF<6,6> constructTransformationMatrix(const FloatMatrixF<3,3> &eigenVectors) const;
    std::pair<FloatMatrixF<6,6>, FloatMatrixF<6,6>> constructL1L2TransformationMatrices(const FloatArrayF<3> &eigenValues, const FloatArrayF<6> &stress, double E1, double E2, double E3) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};

//=============================================================================


class LargeStrainMasterMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatMatrixF<6,6> Pmatrix = eye<6>();
    FloatMatrixF<6,6> TLmatrix, transformationMatrix;
    Domain *domain = nullptr;
    int slaveMat = 0.;

public:
    LargeStrainMasterMaterialStatus(GaussPoint *g, Domain *d, int s);

    const FloatMatrixF<6,6> &givePmatrix() const { return Pmatrix; }
    const FloatMatrixF<6,6> &giveTLmatrix() const { return TLmatrix; }
    const FloatMatrixF<6,6> &giveTransformationMatrix() const { return transformationMatrix; }

    void setPmatrix(const FloatMatrixF<6,6> &values) { Pmatrix = values; }
    void setTLmatrix(const FloatMatrixF<6,6> &values) { TLmatrix = values; }
    void setTransformationMatrix(const FloatMatrixF<6,6> &values) { transformationMatrix = values; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void initTempStatus() override;
    void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "LargeStrainMasterMaterialStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
