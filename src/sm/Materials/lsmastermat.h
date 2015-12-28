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

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/linearelasticmaterial.h"
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
    LinearElasticMaterial *linearElasticMaterial;

    /// 'slave' material model number.
    int slaveMat;
    /// Specifies the strain tensor.
    double m;


public:
    LargeStrainMasterMaterial(int n, Domain * d);
    virtual ~LargeStrainMasterMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveInputRecordName() const { return _IFT_LargeStrainMasterMaterial_Name; }
    virtual const char *giveClassName() const { return "LargeStrainMasterMaterial"; }

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                               MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *, const FloatArray &, TimeStep *)
    { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

    /// transformation matrices
    void constructTransformationMatrix(FloatMatrix &answer, const FloatMatrix &eigenVectors);
    void constructL1L2TransformationMatrices(FloatMatrix &answer1, FloatMatrix &answer2, const FloatMatrix &eigenVectors, FloatArray &stress, double E1, double E2, double E3);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class LargeStrainMasterMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatMatrix Pmatrix, TLmatrix, transformationMatrix;
    int slaveMat;

public:
    LargeStrainMasterMaterialStatus(int n, Domain * d, GaussPoint * g, int s);
    virtual ~LargeStrainMasterMaterialStatus();


    const FloatMatrix &givePmatrix() { return Pmatrix; }
    const FloatMatrix &giveTLmatrix() { return TLmatrix; }
    const FloatMatrix &giveTransformationMatrix() { return transformationMatrix; }

    void setPmatrix(const FloatMatrix &values) { Pmatrix = values; }
    void setTLmatrix(const FloatMatrix &values) { TLmatrix = values; }
    void setTransformationMatrix(const FloatMatrix &values) { transformationMatrix = values; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "LargeStrainMasterMaterialStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
