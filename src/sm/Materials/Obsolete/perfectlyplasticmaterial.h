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

#ifndef perfectlyplasticmaterial_h
#define perfectlyplasticmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"

namespace oofem {
/**
 * This class implements associated Material Status to PerfectlyPlasticMaterial.
 *
 * This class contains only yield_flag to inform, whether
 * gp is yielding or is in elastic state.
 * temp_yield_flag is used to indicate yield_flag during a process when we try
 * to find a equilibrium state. After This State is reached we update the yield_Flag
 * accordingly.
 * When lastly reached state is equilibrium state we call updateYourself(), so
 * yield_flag containing last reached status of yielding is
 * copied into yield_flag indicating yield status at previously
 * reached equilibrium state.
 *
 * Tasks:
 * - Storing the yield_flag
 * - Storing the plastic strain + increment.
 */
class PerfectlyPlasticMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray plasticStrainVector; // reduced form to save memory
    FloatArray plasticStrainIncrementVector;
    int yield_flag;
    int temp_yield_flag;

public:
    PerfectlyPlasticMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    int setTempYieldFlag(int i) { return temp_yield_flag = i; }
    int giveTempYieldFlag() { return temp_yield_flag; }
    int giveYieldFlag() { return yield_flag; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    const FloatArray &givePlasticStrainVector() const { return plasticStrainVector; }
    const FloatArray &givePlasticStrainIncrementVector() const { return plasticStrainIncrementVector; }
    void letPlasticStrainVectorBe(FloatArray v) { plasticStrainVector = std :: move(v); }
    void letPlasticStrainIncrementVectorBe(FloatArray v) { plasticStrainIncrementVector = std :: move(v); }

    const char *giveClassName() const override { return "PerfectlyPlasticMaterialStatus"; }
};

/**
 * This class implements a perfectly plastic material in a finite element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 *
 * This class is mainly abstract, no new methods or data added. This class is
 * mainly intended to be a base class for perfectly plastic materials.
 * (for perfectly plastic material a linear elastic behaviour is assumed
 * before yielding surface is being reached, and also unloading is assumed to
 * be linear elastic).
 * This class takes care about possible failure and fracture and hardening,
 * even if it is not required
 * for perfectly-plastic material, but this class is assumed to be a parent
 * class for more sophisticated materials where these phenomena should appear.
 * Tasks:
 * - Returning standard material stiffness and flexibility matrices for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at Gauss point for 3d-case.
 */
class PerfectlyPlasticMaterial : public StructuralMaterial
{
protected:
    int yieldCriteria = 0;
    int loadingCriteria = 0;
    LinearElasticMaterial *linearElasticMaterial = nullptr;

public:

    PerfectlyPlasticMaterial(int n, Domain * d) : StructuralMaterial(n, d) {}

    virtual ~PerfectlyPlasticMaterial() {
        delete linearElasticMaterial;
    }

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PerfectlyPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PerfectlyPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PerfectlyPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PerfectlyPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<2> giveRealStressVector_2dBeamLayer(const FloatArrayF<2> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PerfectlyPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<5> giveRealStressVector_PlateLayer(const FloatArrayF<5> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PerfectlyPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    // can be moved to parent classes
    virtual void updateIfFailure(GaussPoint *gp,
                                 FloatArray *stressVector3d,
                                 FloatArray *PlasticStrainVector3d) { }

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "PerfectlyPlasticMaterial"; }
    void initializeFrom(InputRecord &ir) override;
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    double give(int aProperty, GaussPoint *gp) const override;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    // two functions used to initialize and updating temporary variables in
    // gps status. These variables are used to control process, when
    // we try to find equilibrium state

    virtual void giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode mode, GaussPoint *gp,
                                             TimeStep *tStep);
    virtual void giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseMode mode, GaussPoint *gp,
                                                      TimeStep *tStep);

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<2,2> give2dBeamLayerStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<5,5> givePlateLayerStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &strainIncrement, TimeStep *tStep);
    void computePlasticStiffnessAt(FloatMatrix &answer,
                                   GaussPoint *gp,
                                   FloatArray *currentStressVector,
                                   FloatArray *currentPlasticStrainVector,
                                   FloatArray *strainIncrement3d,
                                   TimeStep *tStep,
                                   double &lambda);

    FloatArray *GiveStressCorrectionBackToYieldSurface(GaussPoint *gp,
                                                       FloatArray *stressVector3d,
                                                       FloatArray *plasticVector3d);
    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual double computeYCValueAt(GaussPoint *gp, FloatArray *, FloatArray *) { return -1.0; }
    virtual FloatArray *GiveYCStressGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    virtual FloatArray *GiveLCStressGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    virtual FloatArray *GiveYCPlasticStrainGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    virtual FloatArray *GiveLCPlasticStrainGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    //virtual FloatArray* GiveHardeningGradient (GaussPoint *gp,FloatArray *,FloatArray *) {return NULL;}
    virtual void updateTempYC(GaussPoint *gp, FloatArray *, FloatArray *) { }
    virtual void updateTempLC(GaussPoint *gp, FloatArray *, FloatArray *) { }
    //virtual int updateYieldStatus(GaussPoint* gp, FloatArray* strainIncrementIn3d);
};
} // end namespace oofem
#endif // perfectlyplasticmaterial_h
