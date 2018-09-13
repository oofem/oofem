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
    double preCastStiffnessReduction;

    /// Preconstructed 3d tangent
    FloatMatrixF<6,6> tangent;
    FloatMatrixF<4,4> tangentPlaneStrain;
    FloatMatrixF<3,3> tangentPlaneStress;

    /// Thermal expansion
    FloatArrayF<6> alpha;

public:
    /// Constructor.
    LinearElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    /// Destructor.
    virtual ~LinearElasticMaterial() { }

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    /**
     * Computes the plane strain and plane stress tangents from the 3D tangent.
     */
    void computesSubTangents();

    void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseMode mode, GaussPoint *gp,
                                       TimeStep *tStep) override;
    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep) override;
    void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    ///@todo Should this be virtual? It's never used. It's not part of the base class.
    virtual void giveRealStressVector_3dDegeneratedShell(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep);
    void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStressVector_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveRealStressVector_2dBeamLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;
    void giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;
    void giveRealStressVector_Fiber(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;

    void giveEshelbyStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep) override;
    double giveEnergyDensity(GaussPoint *gp, TimeStep *tStep);
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    ///@todo This makes no sense in this  base class, it should belong to isotropiclinearelastic material.
    virtual double giveShearModulus() { return 1.; }
    int hasCastingTimeSupport() override { return 1.; }
    const char *giveClassName() const override { return "LinearElasticMaterial"; }
};
} // end namespace oofem
#endif // linearelasticmaterial_h
