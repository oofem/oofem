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

#ifndef tutorialmaterial_h
#define tutorialmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"

///@name Input fields for TutorialMaterial
//@{
#define _IFT_TutorialMaterial_Name "tutorialmaterial"
#define _IFT_TutorialMaterial_yieldstress "sigy"
#define _IFT_TutorialMaterial_hardeningmoduli "h"
//@}

namespace oofem {
class Domain;

/**
 * This class implements a isotropic plastic linear material (J2 plasticity condition is used).
 * @author Jim Brozoulis
 */
class TutorialMaterial : public StructuralMaterial
{
protected:
    /// Hardening modulus.
    double H;

    /// Initial (uniaxial) yield stress.
    double sig0;

    IsotropicLinearElasticMaterial D;

public:
    TutorialMaterial(int n, Domain * d);
    virtual ~TutorialMaterial();

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_TutorialMaterial_Name; }
    const char *giveClassName() const override { return "TutorialMaterial"; }
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) override { return true; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;

    void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override;

protected:
    static void giveDeviatoricProjectionMatrix(FloatMatrix &answer);

    static void computeSphDevPartOf(const FloatArray &sigV, FloatArray &sigSph, FloatArray &sigDev); 
};


class TutorialMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Temporary plastic strain (the given iteration)
    FloatArray tempPlasticStrain;

    ///  Last equilibriated plastic strain (end of last time step)
    FloatArray plasticStrain;

    FloatArray tempDevTrialStress;

    double tempK;
    double k;

public:
    TutorialMaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~TutorialMaterialStatus() {}

    const FloatArray &givePlasticStrain() { return plasticStrain; }

    void letTempPlasticStrainBe(const FloatArray &values) { tempPlasticStrain = values; }

    double giveK() { return this->k; }

    void letTempKBe(double value) { tempK = value; }

    void letTempDevTrialStressBe(const FloatArray &values) { tempDevTrialStress = values; }
    const FloatArray &giveTempDevTrialStress() { return tempDevTrialStress; }

    const char *giveClassName() const override { return "TutorialMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    // semi optional methods
    //void printOutputAt(FILE *file, TimeStep *tStep) override;
    //void saveContext(DataStream &stream, ContextMode mode) override;
    //void restoreContext(DataStream &stream, ContextMode mode) override;
};

} // end namespace oofem
#endif // tutorialmaterial_h
