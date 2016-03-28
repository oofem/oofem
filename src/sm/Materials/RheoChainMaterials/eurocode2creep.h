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

#ifndef eurocode2creep_h
#define eurocode2creep_h

#include "kelvinChM.h"

///@name Input fields for B3Material
//@{
#define _IFT_Eurocode2CreepMaterial_Name "ec2creepmat"
#define _IFT_Eurocode2CreepMaterial_fcm28 "fcm28"
#define _IFT_Eurocode2CreepMaterial_stiffnessFactor "stiffnessfactor"
#define _IFT_Eurocode2CreepMaterial_t0 "t0"
#define _IFT_Eurocode2CreepMaterial_cemType "cemtype"
#define _IFT_Eurocode2CreepMaterial_henv "henv"
#define _IFT_Eurocode2CreepMaterial_h0 "h0"
#define _IFT_Eurocode2CreepMaterial_shType "shtype"
#define _IFT_Eurocode2CreepMaterial_spectrum "spectrum"
#define _IFT_Eurocode2CreepMaterial_temperatureDependent "temperaturedependent"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to Eurocode2CreepMaterial.
 */
class Eurocode2CreepMaterialStatus : public KelvinChainMaterialStatus
{
protected:
    double maturity;
    double tempMaturity;
    double temperature;
    double tempTemperature;

public:
    Eurocode2CreepMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~Eurocode2CreepMaterialStatus() { }

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    double giveConcreteMaturity() const { return maturity; }
    void setTempConcreteMaturity(double src) { tempMaturity = src; }

    double giveTemperature() const { return temperature; }
    void setTempTemperature(double src) { tempTemperature = src; }


    // definition
    virtual const char *giveClassName() const { return "Eurocode2CreepMaterialStatus"; }
};


/**
 * This class implements a model for concrete creep and shrinkage according to EuroCode 2
 * The creep is assumed to be linear (formula from section 3.7 is not considered here)
 */
class Eurocode2CreepMaterial : public KelvinChainMaterial
{
protected:

    /// fixed retardation time of the first unit
    // the reason it is introduced is the application of the correction factors
    // to achieve a better approximation of the compliance function by the retardation spectrum
    double tau1;

    /// stiffness of the zeroth Kelvin unit
    double EspringVal;

    // ELASTICITY + SHORT TERM + STRENGTH
    /// mean compressive strength at 28 days default - to be specified in units of the analysis (e.g. 30.e6 + stiffnessFacotr 1. or 30. + stiffnessFactor 1.e6)
    double fcm28;

    /// Young's modulus at 28 days default [MPa]
    double Ecm28;

    /// factor unifying stiffnesses (Ecm is predicted from fcm...)
    double stiffnessFactor;

    /// parameter determined by cement type
    double s;

    // CREEP
    /// onset of drying [day, sec, ...]
    double t0;

    /// drying creep coefficient
    double phi_RH;

    /// drying creep coefficient
    double beta_fcm;

    /// drying creep coefficient
    double beta_H;

    /// effective thickness [mm]
    double h0;

    /// influence of cement type on concrete equivalent age, B.9 in EC2
    double alpha_T_cement;


    // DRYING SHRINKAGE
    /// age (absolute) when concrete started drying, must be in days!
    double begOfDrying;

    /// drying shrinkage coefficient
    double kh;

    /// asymptotic value of drying shrinkage at zero relative humidity, B.11 in EC2
    double eps_cd_0;


    // AUTEGENOUS SHRINKAGE
    /// asymptotic value of autogenous shrinakge, 3.12 in EC2
    double eps_ca_infty;


    enum ec2ShrinkageType { EC2_NoShrinkage, EC2_TotalShrinkage, EC2_DryingShrinkage, EC2_AutogenousShrinkage } shType;

    /**
     * If true, analysis of retardation spectrum is used for evaluation of Kelvin units moduli
     * If false, least-squares method is used for evaluation of Kelvin units moduli (default)
     */
    bool retardationSpectrumApproximation;

    bool temperatureDependent;


public:
    Eurocode2CreepMaterial(int n, Domain *d) : KelvinChainMaterial(n, d) {
        shType = EC2_NoShrinkage;
    }
    virtual ~Eurocode2CreepMaterial() { }

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual const char *giveClassName() const { return "Eurocode2CreepMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_Eurocode2CreepMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual double computeConcreteStrengthAtAge(double age);

    virtual double computeMeanElasticModulusAtAge(double age);


    /// Evaluation of the compliance function
    virtual double computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep);

    /// Evaluation of the compliance function
    virtual double computeCreepCoefficient(double t, double t_prime, GaussPoint *gp, TimeStep *tStep);


protected:
    virtual int hasIncrementalShrinkageFormulation() { return 1; }

    //    virtual void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep);

    /// implements B.9
    virtual double computeEquivalentAge(GaussPoint *gp, TimeStep *tStep);

    /// implements B.10
    virtual double computeEquivalentMaturity(GaussPoint *gp, TimeStep *tStep);


    /// sets parameters for elasticity and strength according to formulas from EC
    void computeElasticityStrengthParams(int);

    /// sets parameters for shrinkage according to formulas from EC
    void computeShrinkageParams(int, double);

    /// sets parameters for creep according to formulas from EC
    void computeCreepParams(int, double);

    virtual void computeCharTimes();

    /// computes correction factor which multiplies the retardation times
    double computeRetardationTimeCorrection(int mu);

    /// evaluates retardation spectrum at given time (t-t')
    double evaluateSpectrumAt(double tau);

    /// Evaluation of characteristic moduli of the Kelvin chain.
    virtual void computeCharCoefficients(FloatArray &answer, double tPrime, GaussPoint *gp, TimeStep *tStep);

    void computeIncrementOfDryingShrinkageVector(FloatArray &answer, GaussPoint *gp, double tNow, double tThen);

    void computeIncrementOfAutogenousShrinkageVector(FloatArray &answer, GaussPoint *gp, double tNow, double tThen);
};
} // end namespace oofem
#endif // eurocode2creep_h
