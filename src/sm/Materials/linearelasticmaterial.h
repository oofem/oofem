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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef linearelasticmaterial_h
#define linearelasticmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"

///@name Input fields for LinearElasticMaterial
//@{
#define _IFT_LinearElasticMaterial_preCastStiffRed "precaststiffred"
//@}


namespace oofem {
/**
 * This class is a abstract base class for all linear elastic material models
 * in a finite element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 * Efficient implementation of services for obtaining characteristic matrices
 * for several material modes is provided depending on other abstract services.
 * Also general implementation of giveRealStressVector service is provided,
 * computing the stress increment vector from strain increment multiplied by
 * stiffness.
 *
 * Tasks:
 * - Returning standard material stiffness and flexibility matrices for 3d-case.
 *   according to current state determined by using data stored in Gausspoint.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class LinearElasticMaterial : public StructuralMaterial
{
protected:
    /// artificial isotropic damage to reflect reduction in stiffness for time < castingTime.
    double preCastStiffnessReduction = 0.;

    /// Preconstructed 3d tangent
    FloatMatrixF<6,6> tangent;
    FloatMatrixF<4,4> tangentPlaneStrain;
    FloatMatrixF<3,3> tangentPlaneStress;

    /// Thermal expansion
    FloatArrayF<6> alpha;

public:
    /// Constructor.
    LinearElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    const FloatMatrixF<6,6> &giveTangent() const { return tangent; }
    const FloatArrayF<6> &giveAlpha() const { return alpha; }
    
    /**
     * Computes the plane strain and plane stress tangents from the 3D tangent.
     */
    void computesSubTangents();

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    ///@todo Should this be virtual? It's never used. It's not part of the base class.
    virtual FloatArrayF<6> giveRealStressVector_3dDegeneratedShell(const FloatArrayF<6> &reducedF, GaussPoint *gp, TimeStep *tStep) const;
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const  override;
    FloatArrayF<2> giveRealStressVector_Warping(const FloatArrayF<2> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<2> giveRealStressVector_2dBeamLayer(const FloatArrayF<2> &reducedE, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<5> giveRealStressVector_PlateLayer(const FloatArrayF<5> &reducedE, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStressVector_Fiber(const FloatArrayF<3> &reducedE, GaussPoint *gp, TimeStep *tStep) const override;

    void giveEshelbyStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep) override;
    double giveEnergyDensity(GaussPoint *gp, TimeStep *tStep);
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    ///@todo This makes no sense in this  base class, it should belong to isotropiclinearelastic material.
    virtual double giveShearModulus() const { return 1.; }
    bool hasCastingTimeSupport() const override { return true; }
    const char *giveClassName() const override { return "LinearElasticMaterial"; }
};
} // end namespace oofem
#endif // linearelasticmaterial_h
