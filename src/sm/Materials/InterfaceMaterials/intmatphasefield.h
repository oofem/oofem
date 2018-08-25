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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
/*
 *
 *
 *  Author: Jim Brouzoulis
 */

#ifndef intmatphasefield
#define intmatphasefield

#include "structuralinterfacematerialphf.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatPhaseField
//@{
#define _IFT_IntMatPhaseField_Name "intmatphf"
#define _IFT_IntMatPhaseField_kn "k"
#define _IFT_IntMatPhaseField_gc "gc"
//@}


namespace oofem {

/**
 * Development cz-model using phase field
 */
class IntMatPhaseFieldStatus : public StructuralInterfaceMaterialStatus
{
public:
    IntMatPhaseFieldStatus(int n, Domain * d, GaussPoint * g);
    virtual ~IntMatPhaseFieldStatus() { }

    /// damage variable
    double tempDamage;
    double giveDamage() override { return tempDamage; }

    double tempDrivingEnergy;
    double drivingEnergy;
    double giveTempDrivingEnergy() { return tempDrivingEnergy; }
    double giveDrivingEnergy() { return drivingEnergy; }
    void letTempDrivingEnergyBe(double val) { this->tempDrivingEnergy = val; }

    const char *giveClassName() const override { return "IntMatPhaseFieldStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
};


class IntMatPhaseField : public StructuralInterfaceMaterialPhF
{
protected:
    double k;
    double Gc;

public:
    IntMatPhaseField(int n, Domain * d);
    virtual ~IntMatPhaseField();

    int hasMaterialModeCapability(MaterialMode mode) override;
    const char *giveClassName() const override { return "IntMatPhaseField"; }
    const char *giveInputRecordName() const override { return _IFT_IntMatPhaseField_Name; }

    void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep) override;
    void give3dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveTangents(FloatMatrix &jj, FloatMatrix &jd, FloatMatrix &dj, FloatMatrix &dd, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new IntMatPhaseFieldStatus(1, domain, gp); }; 

    bool hasAnalyticalTangentStiffness() const override { return true; };

    double giveDrivingForce(GaussPoint *gp) override;
    double giveDrivingForcePrime(GaussPoint *gp) override;
    //double compute_fPrime(const double d);
    double compute_g(const double d);
    double compute_gPrime(const double d);
    double compute_gBis(const double d);
};

} /* namespace oofem */
#endif 
