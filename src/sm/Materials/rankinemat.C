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

#include "rankinemat.h"
#include "isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(RankineMat);

// constructor
RankineMat :: RankineMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    E = 0.;
    nu = 0.;
    H0 = 0.;
    sig0 = 0.;
    delSigY = 0.;
    ep = 0.;
    md = 0.;
    damlaw = 1;
    param1 = 0.;
    param2 = 0.;
    param3 = 0.;
    param4 = 0.;
    param5 = 0.;
}


// specifies whether a given material mode is supported by this model
int
RankineMat :: hasMaterialModeCapability(MaterialMode mode)
{
    return ( ( mode == _PlaneStress ) || ( mode == _1dMat ) );
}


// reads the model parameters from the input file
IRResultType
RankineMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    result = linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;

    E = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveYoungsModulus();
    nu = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->givePoissonsRatio();

    IR_GIVE_FIELD(ir, sig0, _IFT_RankineMat_sig0); // uniaxial yield stress

    H0 = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H0, _IFT_RankineMat_h); // hardening modulus

    plasthardtype = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, plasthardtype, _IFT_RankineMat_plasthardtype); // type of plastic hardening (0=linear, 1=exponential, 2= prepeak hardening + linear softening )

    delSigY = 0.;
    if ( plasthardtype == 0 ) {
        //no extra required variables
    } else if ( plasthardtype == 1 ) {
        IR_GIVE_FIELD(ir, delSigY, _IFT_RankineMat_delsigy); // final increment of yield stress (at infinite cumulative plastic strain)
    } else if ( plasthardtype == 2 ) {
        IR_GIVE_FIELD(ir, ep, _IFT_RankineMat_ep);
        ep = ep - sig0 / E; // user input is strain at peak stress sig0 and is converted to plastic strain at peak stress sig0
        md = 1. / log(50. * E * ep / sig0); // exponent used on the 1st plasticity branch
    } else {
        OOFEM_WARNING("Plasticity hardening type number  %d is unknown", plasthardtype);
        return IRRT_BAD_FORMAT;
    }

    yieldtol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldtol, _IFT_RankineMat_yieldtol); // relative tolerance in yield condition

    damlaw = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, damlaw, _IFT_RankineMat_damlaw); // type of damage law (0=exponential, 1=exponential and  damage starts after peak stress sig0)

    if  ( damlaw == 0 ) {
        a = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_RankineMat_a); // coefficient in damage law
    } else if ( damlaw == 1 ) {
        IR_GIVE_FIELD(ir, param1, _IFT_RankineMat_param1); // coefficient in damage law
        IR_GIVE_FIELD(ir, param2, _IFT_RankineMat_param2); // coefficient in damage law. If b<1 use only stiffmode=1
    } else if ( damlaw == 2 ) {
        IR_GIVE_FIELD(ir, param1, _IFT_RankineMat_param1); // coefficients in damage law
        IR_GIVE_FIELD(ir, param2, _IFT_RankineMat_param2);
        IR_GIVE_FIELD(ir, param3, _IFT_RankineMat_param3);
        IR_GIVE_FIELD(ir, param4, _IFT_RankineMat_param4);
        IR_GIVE_FIELD(ir, param5, _IFT_RankineMat_param5);
    } else {
        OOFEM_WARNING("Damage law number  %d is unknown", damlaw);
        return IRRT_BAD_FORMAT;
    }

    double gf = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, gf, _IFT_RankineMat_gf); // dissipated energy per unit VOLUME

    if ( ( a != 0. ) && ( gf != 0 ) ) {
        OOFEM_WARNING("parameters a and gf cannot be prescribed simultaneously");
        return IRRT_BAD_FORMAT;
    }

    if ( gf > 0. ) {
        // evaluate parameter "a" from given "gf"
        double A = H0 * ( 1. + H0 / E );
        double B = sig0 * ( 1. + H0 / E );
        double C = 0.5 * sig0 * sig0 / E - gf;
        if ( C >= 0. ) {
            OOFEM_ERROR("parameter gf is too low");
        }

        double kappaf = ( -B + sqrt(B * B - 4. * A * C) ) / ( 2. * A );
        a = 1. / kappaf;
    }

    return IRRT_OK;
}


// creates a new material status  corresponding to this class
MaterialStatus *
RankineMat :: CreateStatus(GaussPoint *gp) const
{
    RankineMatStatus *status;
    status = new RankineMatStatus(1, this->giveDomain(), gp);
    return status;
}


// computes the stress vector corresponding to given (final) strain
void
RankineMat :: giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );

    // initialization
    this->initTempStatus(gp);

    // elastoplasticity
    this->performPlasticityReturn(gp, totalStrain);

    // damage
    double omega = computeDamage(gp, tStep);
    answer.beScaled(1. - omega, status->giveTempEffectiveStress());

    // store variables in status
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    status->computeWork_1d(gp, gf);
#endif
}


void
RankineMat :: giveRealStressVector_PlaneStress(FloatArray &answer,
                                               GaussPoint *gp,
                                               const FloatArray &totalStrain,
                                               TimeStep *tStep)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );

    // initialization
    this->initTempStatus(gp);

    // elastoplasticity
    this->performPlasticityReturn(gp, totalStrain);

    // damage
    double omega = computeDamage(gp, tStep);
    answer.beScaled(1. - omega, status->giveTempEffectiveStress());

    // store variables in status
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    status->computeWork_PlaneStress(gp, gf);
#endif
}


double
RankineMat :: evalYieldFunction(const FloatArray &sigPrinc, const double kappa)
{
    return sigPrinc.at(1) - evalYieldStress(kappa);
}


double
RankineMat :: evalYieldStress(const double kappa)
{
    double yieldStress = 0.;
    if ( plasthardtype == 0 ) { // linear hardening
        yieldStress = sig0 + H0 * kappa;
    } else if ( plasthardtype == 1 ) { // exponential hardening
        if ( delSigY == 0. ) {
            yieldStress = sig0;
        } else {
            yieldStress = sig0 + delSigY * ( 1. - exp(-H0 * kappa / delSigY) );
        }
    } else if ( plasthardtype == 2 ) { // exponential hardening before the peak stress sig0 and linear after the peak stress sig0
        if ( kappa <= ep ) {
            //1st branch in the rankine variation 2 trying to match the 1st branch of the smooth extended damage law reported by Grassl and Jirasek (2010)
            yieldStress = 50. *E *kappa *exp(-pow ( kappa / ep, md ) / md);
        } else {  //linear hardening branch
            yieldStress = sig0 + H0 * kappa;
        }
    }
    return yieldStress;
}

double
RankineMat :: evalPlasticModulus(const double kappa)
{
    double plasticModulus = 0.;
    if ( plasthardtype == 0 ) { // linear hardening
        plasticModulus = H0;
    } else if ( plasthardtype == 1 ) { // exponential hardening
        plasticModulus = H0 * exp(-H0 * kappa / delSigY);
    } else if ( plasthardtype == 2 ) { // exponential hardening before the peak stress sig0 and linear after the peak stress sig0
        if ( kappa <= ep ) { //1st branch of yield stress
            double aux = pow(kappa / ep, md);
            plasticModulus = 50. *E *exp(-aux / md) * ( 1. - aux );
        } else {  //2nd branch of yield stress
            plasticModulus = H0;
        }
    }
    return plasticModulus;
}


// computes the stress according to elastoplasticity
// (return of trial stress to the yield surface)
void
RankineMat :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain)
{
    double kappa, tempKappa, H;
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    MaterialMode mode = gp->giveMaterialMode();
    // get the initial plastic strain and initial kappa from the status
    FloatArray tempPlasticStrain = status->givePlasticStrain();
    kappa = tempKappa = status->giveCumulativePlasticStrain();

    // elastic predictor
    StrainVector elStrain(totalStrain, mode);
    elStrain.subtract(tempPlasticStrain);
    StressVector finalStress(mode);
    elStrain.applyElasticStiffness(finalStress, E, nu);
    FloatArray sigPrinc;
    FloatMatrix nPrinc;
    // get principal trial stresses (ordered) and principal stress directions
    finalStress.computePrincipalValDir(sigPrinc, nPrinc);
    double ftrial = evalYieldFunction(sigPrinc, tempKappa);
    if ( mode == _1dMat ) { //1d case
        //// Plastic Corrector
        if ( ftrial > 0. ) {
            double f = ftrial;
            // calculate the increment of cumulative plastic strain
            int i = 1;
            do {
                i = i + 1;
                if ( i > 1000 ) {
                    printf("kappa, ftrial: %g %g\n", kappa, ftrial);
                    OOFEM_ERROR("no convergence of regular stress return algorithm");
                }

                double ddKappa = f / ( E + evalPlasticModulus(tempKappa) );
                finalStress.at(1) -= E * ddKappa;
                tempKappa += ddKappa;
                f = finalStress.at(1) - evalYieldStress(tempKappa);
            } while ( fabs(f) > yieldtol * sig0 );
        }

        tempPlasticStrain.at(1) = tempKappa;
        //End of Dimitris change
    } else {  //Plane stress case
        double difPrincTrialStresses = sigPrinc.at(1) - sigPrinc.at(2);
        double tanG = E / ( 2. * ( 1. + nu ) );

        // plastic corrector - regular case
        bool vertex_case = false;
        if ( ftrial > 0. ) {
            double f = ftrial;
            double Enu = E / ( 1. - nu * nu );
            // calculate the increment of cumulative plastic strain
            int i = 1;
            do {
                if ( i++ > 50 ) {
                    finalStress.computePrincipalValDir(sigPrinc, nPrinc);
                    sigPrinc.pY();
                    printf("kappa, ftrial: %g %g\n", kappa, ftrial);
                    OOFEM_ERROR("no convergence of regular stress return algorithm");
                }

                H = evalPlasticModulus(tempKappa);
                double ddKappa = f / ( Enu + H );
                sigPrinc.at(1) -= Enu * ddKappa;
                sigPrinc.at(2) -= nu * Enu * ddKappa;
                tempKappa += ddKappa;
                f = evalYieldFunction(sigPrinc, tempKappa);
            } while ( fabs(f) > yieldtol * sig0 );

            if ( sigPrinc.at(2) > sigPrinc.at(1) ) {
                // plastic corrector - vertex case
                // ---------------------------------
                vertex_case = true;
                // recompute trial principal stresses
                finalStress.computePrincipalValDir(sigPrinc, nPrinc);
                // calculate the increment of cumulative plastic strain
                double sigstar = ( sigPrinc.at(1) - nu * sigPrinc.at(2) ) / ( 1. - nu );
                double alpha = E / ( 1. - nu );
                double dkap0 = ( sigPrinc.at(1) - sigPrinc.at(2) ) * ( 1. + nu ) / E;
                tempKappa = kappa;
                f = sigstar - evalYieldStress(tempKappa);
                double dkap1 = 0.;
                H = evalPlasticModulus(tempKappa);
                double C = alpha +  H * ( 1. + sqrt(2.) ) / 2.;
                i = 1;
                do {
                    if ( i++ > 20 ) {
                        finalStress.computePrincipalValDir(sigPrinc, nPrinc);
                        sigPrinc.pY();
                        printf("kappa, ftrial: %g %g\n", kappa, ftrial);
                        OOFEM_ERROR("no convergence of vertex stress return algorithm");
                    }

                    dkap1 += f / C;
                    tempKappa = kappa + sqrt( dkap1 * dkap1 + ( dkap1 - dkap0 ) * ( dkap1 - dkap0 ) );
                    f = sigstar - evalYieldStress(tempKappa) - alpha * dkap1;
                    double aux = dkap1 * dkap1 + ( dkap1 - dkap0 ) * ( dkap1 - dkap0 );
                    H = evalPlasticModulus(tempKappa);
                    if ( aux > 0. ) {
                        C = alpha + H * ( 2. * dkap1 - dkap0 ) / sqrt(aux);
                    } else {
                        C = alpha + H *sqrt(2.);
                    }
                } while ( fabs(f) > yieldtol * sig0 );

                sigPrinc.at(1) = sigPrinc.at(2) = sigstar - alpha * dkap1;
                status->setDKappa(dkap1, dkap1 - dkap0);
            }

            // principal stresses
            double sig1 = sigPrinc.at(1);
            double sig2 = sigPrinc.at(2);
            // compose the stress in global coordinates
            //   the first subscript refers to coordinate
            //   the second subscript refers to eigenvalue
            double n11 = nPrinc.at(1, 1);
            double n12 = nPrinc.at(1, 2);
            double n21 = nPrinc.at(2, 1);
            double n22 = nPrinc.at(2, 2);
            finalStress.at(1) = sig1 * n11 * n11 + sig2 * n12 * n12;
            finalStress.at(2) = sig1 * n21 * n21 + sig2 * n22 * n22;
            finalStress.at(3) = sig1 * n11 * n21 + sig2 * n12 * n22;
            // add the increment of plastic strain
            if ( !vertex_case ) {
                tempPlasticStrain.at(1) += ( tempKappa - kappa ) * n11 * n11;
                tempPlasticStrain.at(2) += ( tempKappa - kappa ) * n21 * n21;
                tempPlasticStrain.at(3) += 2. * ( tempKappa - kappa ) * n11 * n21;
            } else {
                double dkap1 = status->giveDKappa(1);
                double dkap2 = status->giveDKappa(2);
                tempPlasticStrain.at(1) += dkap1 * n11 * n11 + dkap2 * n12 * n12;
                tempPlasticStrain.at(2) += dkap1 * n21 * n21 + dkap2 * n22 * n22;
                tempPlasticStrain.at(3) += 2. * ( dkap1 * n11 * n21 + dkap2 * n12 * n22 );
            }

            // evaluate the tangent shear stiffness
            if ( difPrincTrialStresses != 0. ) {
                double factor = ( sig1 - sig2 ) / difPrincTrialStresses;
                if ( factor > 0. && factor <= 1. ) {
                    tanG *= factor;
                }
            }
        }

        status->setTangentShearStiffness(tanG); // store shear stiffness.Used in 2d/3d cases
    }

    // store the effective stress, plastic strain and cumulative plastic strain
    status->letTempEffectiveStressBe(finalStress);
    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->setTempCumulativePlasticStrain(tempKappa);
}

double
RankineMat :: computeDamageParam(double tempKappa)
{
    double tempDam = 0.;
    if ( tempKappa > 0. ) {
        if ( damlaw == 0 ) {
            tempDam = 1.0 - exp(-a * tempKappa);
        } else if ( damlaw == 1 && tempKappa > ep ) {
            tempDam = 1.0 - exp( -param1 * pow( ( tempKappa - ep ) / ep, param2 ) );
        } else if ( damlaw == 2 && tempKappa > ep ) {
            tempDam = 1.0 - param5 *exp( -param1 *pow ( ( tempKappa - ep ) / ep, param2 ) ) - ( 1. - param5 ) * exp( -param3 * pow( ( tempKappa - ep ) / ep, param4 ) );
        }
    }

    return tempDam;
}

double
RankineMat :: computeDamageParamPrime(double tempKappa)
{
    double tempDam = 0.;
    if ( tempKappa >= 0. ) {
        if ( damlaw == 0 ) {
            tempDam = a * exp(-a * tempKappa);
        } else if ( damlaw == 1 && tempKappa >= ep ) {
            tempDam = param1 * param2 * pow( ( tempKappa - ep ) / ep, param2 - 1 ) / ep *exp( -param1 *pow ( ( tempKappa - ep ) / ep, param2 ) );
        } else if ( damlaw == 2 && tempKappa >= ep ) {
            tempDam = param5 * param1 * param2 * pow( ( tempKappa - ep ) / ep, param2 - 1 ) / ep *exp( -param1 *pow ( ( tempKappa - ep ) / ep, param2 ) ) + ( 1. - param5 ) * param3 * param4 * pow( ( tempKappa - ep ) / ep, param4 - 1 ) / ep *exp( -param3 *pow ( ( tempKappa - ep ) / ep, param4 ) );
        }
    }

    return tempDam;
}

double
RankineMat :: computeDamage(GaussPoint *gp,  TimeStep *tStep)
{
    double tempKappa, dam;
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    dam = status->giveDamage();
    computeCumPlastStrain(tempKappa, gp, tStep);
    double tempDam = computeDamageParam(tempKappa);
    if ( dam > tempDam ) {
        tempDam = dam;
    }

    return tempDam;
}

void RankineMat :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *tStep)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    tempKappa = status->giveTempCumulativePlasticStrain();
}

// returns the consistent (algorithmic) stiffness matrix
void
RankineMat :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double gprime = computeDamageParamPrime(tempKappa);
    evaluatePlaneStressStiffMtrx(answer, mode, gp, tStep, gprime);
}

void
RankineMat :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    answer.resize(1, 1);
    answer.at(1, 1) = this->E;
    if ( mode == ElasticStiffness ) {
    } else if ( mode == SecantStiffness ) {
        double om = status->giveTempDamage();
        answer.times(1.0 - om);
    } else {
        OOFEM_ERROR("unknown type of stiffness (secant stiffness not implemented for 1d)");
    }
}

// this method is also used by the gradient version,
// with gprime replaced by gprime*m and evaluated for kappa hat
void
RankineMat :: evaluatePlaneStressStiffMtrx(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep, double gprime)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    if ( mode == ElasticStiffness || mode == SecantStiffness ) {
        // start from the elastic stiffness
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
        if ( mode == SecantStiffness ) {
            // transform to secant stiffness
            double damage = status->giveTempDamage();
            answer.times(1. - damage);
        }

        return;
    }

    // check the unloading condition
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    if ( tempKappa <= kappa ) { // tangent matrix requested, but unloading takes place - use secant
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
        double damage = status->giveTempDamage();
        answer.times(1. - damage);
        return;
    }

    // tangent stiffness requested, loading

    // elastoplastic tangent matrix in principal stress coordinates
    answer.resize(3, 3);
    answer.zero(); // just to be sure (components 13,23,31,32 must be zero)

    double eta1, eta2, dkap2;
    double dkap1 = status->giveDKappa(1);
    double H = evalPlasticModulus(tempKappa);

    if ( dkap1 == 0. ) {
        // regular case

        double Enu = E / ( 1. - nu * nu );
        double aux = Enu / ( Enu + H );
        answer.at(1, 1) = aux * H;
        answer.at(1, 2) = answer.at(2, 1) = nu * aux * H;
        answer.at(2, 2) = aux * ( E + H );
        answer.at(3, 3) = status->giveTangentShearStiffness();
        eta1 = aux;
        eta2 = nu * aux;
    } else {
        // vertex case

        dkap2 = status->giveDKappa(2);
        double denom = E * dkap1 + H * ( 1. - nu ) * ( dkap1 + dkap2 );
        eta1 = E * dkap1 / denom;
        eta2 = E * dkap2 / denom;
        answer.at(1, 1) = answer.at(2, 1) = H * eta1;
        answer.at(1, 2) = answer.at(2, 2) = H * eta2;
        answer.at(3, 3) = 0.; //E*1.e-6; // small stiffness, to suppress singularity
    }

    // the elastoplastic stiffness has been computed
    // now add the effect of damage

    double damage = status->giveTempDamage();
    answer.times(1. - damage);

    FloatArray sigPrinc(2);
    FloatMatrix nPrinc(2, 2);
    StressVector effStress(status->giveTempEffectiveStress(), _PlaneStress);
    effStress.computePrincipalValDir(sigPrinc, nPrinc);
    // sometimes the method is called with gprime=0., then we can save some work
    if ( gprime != 0. ) {
        FloatMatrix correction(3, 3);
        correction.zero();
        correction.at(1, 1) = sigPrinc.at(1) * eta1;
        correction.at(1, 2) = sigPrinc.at(1) * eta2;
        correction.at(2, 1) = sigPrinc.at(2) * eta1;
        correction.at(2, 2) = sigPrinc.at(2) * eta2;
        correction.times(gprime); // input parameter gprime used here
        answer.subtract(correction);
    }

    // transform to global coordinates
    FloatMatrix T;
    givePlaneStressVectorTranformationMtrx(T, nPrinc, true);
    answer.rotatedWith(T, 't');
}

// derivatives of final kappa with respect to final strain
void
RankineMat :: computeEta(FloatArray &answer, RankineMatStatus *status)
{
    FloatArray eta(3);
    double dkap1 = status->giveDKappa(1);
    double kap = status->giveTempCumulativePlasticStrain();
    double H = evalPlasticModulus(kap);

    // evaluate in principal coordinates

    if ( dkap1 == 0. ) {
        // regular case
        double Estar = E / ( 1. - nu * nu );
        double aux = Estar / ( H + Estar );
        eta.at(1) = aux;
        eta.at(2) = nu * aux;
        eta.at(3) = 0.;
    } else {
        // vertex case
        double dkap2 = status->giveDKappa(2);
        double denom = E * dkap1 + H * ( 1. - nu ) * ( dkap1 + dkap2 );
        eta.at(1) = E * dkap1 / denom;
        eta.at(2) = E * dkap2 / denom;
        eta.at(3) = 0.;
    }

    // transform to global coordinates

    FloatArray sigPrinc(2);
    FloatMatrix nPrinc(2, 2);
    StressVector effStress(status->giveTempEffectiveStress(), _PlaneStress);
    effStress.computePrincipalValDir(sigPrinc, nPrinc);

    FloatMatrix T(3, 3);
    givePlaneStressVectorTranformationMtrx(T, nPrinc, true);
    answer.beProductOf(T, eta);
}

int
RankineMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    RankineMatStatus *status = static_cast< RankineMatStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        const FloatArray ep = status->givePlasDef();
        answer.resize(6);
        answer.at(1) = ep.at(1);
        answer.at(2) = ep.at(2);
        answer.at(3) = 0.; ///@todo Fix this value!
        answer.at(4) = 0.;
        answer.at(5) = 0.;
        answer.at(6) = ep.at(3);
        return 1;
    } else if ( type == IST_CumPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_StressWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
        return 1;
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


//=============================================================================
// RANKINE MATERIAL STATUS
//=============================================================================

RankineMatStatus :: RankineMatStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrain(), tempPlasticStrain()
{
    damage = tempDamage = 0.;
    kappa = tempKappa = 0.;
    dKappa1 = dKappa2 = 0.;
    tanG = 0.;
    //effStress.resize(3);
    //tempEffStress.resize(3);
#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}


RankineMatStatus :: ~RankineMatStatus()
{ }


void
RankineMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    //int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "status { ");
    /*
     * // this would not be correct, since the status is already updated and kappa is now the "final" value
     * if ( tempKappa > kappa ) {
     * fprintf(file, " Yielding, ");
     * } else {
     * fprintf(file, " Unloading, ");
     * }
     */
    fprintf(file, "kappa %g, damage %g ", kappa, tempDamage);
#ifdef keep_track_of_dissipated_energy
    fprintf(file, ", dissW %g, freeE %g, stressW %g ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif
#if 0
    // print the plastic strain
    fprintf(file, " plastic_strains ");
    for ( auto &val : plasticStrain ) {
        fprintf( file, " %.4e", val );
    }
#endif
    // print the cumulative plastic strain
    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void RankineMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrain.zero();
    }

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
#ifdef keep_track_of_dissipated_energy
    tempStressWork = stressWork;
    tempDissWork = dissWork;
#endif
}


// updates internal variables when equilibrium is reached
void
RankineMatStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    damage = tempDamage;
#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork;
    dissWork = tempDissWork;
#endif
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
RankineMatStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    // write plastic strain (vector)
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write cumulative plastic strain (scalar)
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write damage (scalar)
    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream.write(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK;
}


contextIOResultType
RankineMatStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read plastic strain (vector)
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read cumulative plastic strain (scalar)
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read damage (scalar)
    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream.read(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK; // return success
}


#ifdef keep_track_of_dissipated_energy
void
RankineMatStatus :: computeWork_PlaneStress(GaussPoint *gp, double gf)
{
    // int n = deps.giveSize(); // would not work for gradient version
    int n = 3;

    // strain increment
    FloatArray deps;
    deps.beDifferenceOf(tempStrainVector, strainVector, n);

    // increment of stress work density
    double dSW = ( tempStressVector.dotProduct(deps, n) + stressVector.dotProduct(deps, n) ) / 2.;
    tempStressWork = stressWork + dSW;

    // elastically stored energy density
    FloatArray tempElasticStrainVector;
    tempElasticStrainVector.beDifferenceOf(tempStrainVector, tempPlasticStrain, n);
    double We = tempStressVector.dotProduct(tempElasticStrainVector, n) / 2.;

    // dissipative work density
    tempDissWork = tempStressWork - We;
    // to avoid extremely small negative dissipation due to round-off error
    // (note: gf is the dissipation density at complete failure, per unit volume)
    if ( fabs(tempDissWork) < 1.e-12 * gf ) {
        tempDissWork = 0.;
    }
}

void
RankineMatStatus :: computeWork_1d(GaussPoint *gp, double gf)
{
    // int n = deps.giveSize(); // would not work for gradient version
    int n = 1;


    // strain increment
    FloatArray deps;
    deps.beDifferenceOf(tempStrainVector, strainVector, n);

    // increment of stress work density
    double dSW = ( tempStressVector.dotProduct(deps, n) + stressVector.dotProduct(deps, n) ) / 2.;
    tempStressWork = stressWork + dSW;

    // elastically stored energy density
    FloatArray tempElasticStrainVector;
    tempElasticStrainVector.beDifferenceOf(tempStrainVector, tempPlasticStrain, n);
    double We = tempStressVector.dotProduct(tempElasticStrainVector, n) / 2.;

    // dissipative work density
    tempDissWork = tempStressWork - We;
    // to avoid extremely small negative dissipation due to round-off error
    // (note: gf is the dissipation density at complete failure, per unit volume)
    if ( fabs(tempDissWork) < 1.e-12 * gf ) {
        tempDissWork = 0.;
    }
}

#endif
} // end namespace oofem
