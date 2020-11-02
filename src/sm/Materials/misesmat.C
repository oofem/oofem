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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "misesmat.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "function.h"
#include "intarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "fieldmanager.h"
#include "sm/Elements/structuralelement.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(MisesMat);


MisesMat::MisesMat(int n, Domain *d) : StructuralMaterial(n, d),
    linearElasticMaterial(n, d)
{}


void
MisesMat::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);
    linearElasticMaterial.initializeFrom(ir); // takes care of elastic constants

    G = linearElasticMaterial.giveShearModulus();
    K = linearElasticMaterial.giveBulkModulus();


    hType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hType, _IFT_MisesMat_htype); //hardening type

    if ( hType == 0 ) {
        IR_GIVE_FIELD(ir, sig0, _IFT_MisesMat_sig0); // uniaxial yield stress
        H = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_MisesMat_h); // hardening modulus
    } else if ( hType == 1 ) {     //user defined hardening function
        IR_GIVE_FIELD(ir, h_eps, _IFT_MisesMat_h_eps);
        IR_GIVE_FIELD(ir, h_function_eps, _IFT_MisesMat_h_function_eps);

        if ( h_eps.at(1) != 0. ) {
            throw ValueInputException(ir, _IFT_MisesMat_h_eps, "The first entry in h_eps must be 0.");
        }

        if ( h_eps.giveSize() != h_function_eps.giveSize() ) {
            throw ValueInputException(ir, _IFT_MisesMat_h_function_eps, "the size of 'h_eps' and 'h(eps)' must be the same");
        }
    } else {
        throw ValueInputException(ir, _IFT_MisesMat_htype, "Unknown htype. Should be either 0 or 1.\n");
    }

    omega_crit = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, omega_crit, _IFT_MisesMat_omega_crit); // critical damage

    a = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_MisesMat_a); // exponent in damage law

    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_MisesMat_yieldTol); // tolerance in the yield condition
}

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMat::CreateStatus(GaussPoint *gp) const
{
    return new MisesMatStatus(gp);
}

FloatArrayF< 1 >
MisesMat::giveRealStressVector_1d(const FloatArrayF< 1 > &totalStrain,
                                  GaussPoint *gp,
                                  TimeStep *tStep) const
{
    /// @note: One should obtain the same answer using the iterations in the default implementation (this is verified for this model).
#if 1
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );

    // subtract stress independent part
    FloatArray strainR;
    this->giveStressDependentPartOfStrainVector(strainR, gp, totalStrain, tStep, VM_Total);
    this->performPlasticityReturn(strainR, gp, tStep);
    double omega = computeDamage(gp, tStep);
    FloatArrayF< 6 >stress = status->giveTempEffectiveStress() * ( 1 - omega );

    // Compute the other components of the strain:
    double E = linearElasticMaterial.give('E', gp), nu = linearElasticMaterial.give('n', gp);

    auto strain = status->getTempPlasticStrain();
    strain [ 0 ] = totalStrain [ 0 ];
    strain [ 1 ] -= nu / E * status->giveTempEffectiveStress() [ 0 ];
    strain [ 2 ] -= nu / E * status->giveTempEffectiveStress() [ 0 ];

    status->letTempStrainVectorBe(strain);
    status->setTempDamage(omega);
    status->letTempStressVectorBe(stress);
    return stress [ { 0 } ];

#else
    return StructuralMaterial::giveRealStressVector_1d(totalStrain, gp, tStep);

#endif
}


FloatArrayF< 3 >
MisesMat::giveRealStressVector_PlaneStress(const FloatArrayF< 3 > &totalStrain,
                                           GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    // initialization
    const_cast< MisesMat * >( this )->initTempStatus(gp);
    this->performPlasticityReturn_PlaneStress(totalStrain, gp, tStep);
    double omega = computeDamage(gp, tStep);
    FloatArrayF< 3 >stress = status->giveTempEffectiveStress() * ( 1 - omega );
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    return stress;
}


FloatArrayF< 6 >
MisesMat::giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp,
                                  TimeStep *tStep) const
{
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    // subtract stress independent part
    FloatArray strainR(6);
    this->giveStressDependentPartOfStrainVector(strainR, gp, strain, tStep, VM_Total);

    this->performPlasticityReturn(strainR, gp, tStep);
    double omega = computeDamage(gp, tStep);
    auto stress = status->giveTempEffectiveStress() * ( 1 - omega );
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    return stress;
}


void
MisesMat::performPlasticityReturn(const FloatArray &totalStrain, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa;
    FloatArray plStrain;
    FloatArray fullStress;
    // get the initial plastic strain and initial kappa from the status
    plStrain = status->givePlasticStrain();
    kappa = status->giveCumulativePlasticStrain();

    double dKappa = 0.;
    // === radial return algorithm ===
    if ( totalStrain.giveSize() == 1 ) {
        double E = linearElasticMaterial.give('E', gp);
        /*trial stress*/
        fullStress.resize(6);
        fullStress.at(1) = E * ( totalStrain.at(1) - plStrain.at(1) );
        double trialS = fabs(fullStress.at(1) );
        /*yield function*/
        double yieldValue = trialS - computeYieldStress(kappa, gp, tStep);
        // === radial return algorithm ===
        if ( yieldValue > 0 ) {
            if ( hType == 0 ) {
                dKappa = yieldValue / ( computeYieldStressPrime(kappa) + E );
            }
            if ( hType == 1 ) {
                dKappa += yieldValue / ( computeYieldStressPrime(kappa) + E );
                yieldValue = trialS - checkYieldStress(dKappa, kappa, gp, tStep);
                yieldValue -= E * dKappa;
                if ( yieldValue > 1.e-10 ) {
                    dKappa += yieldValue / ( computeYieldStressPrime(kappa + dKappa) + E );
                }
            }

            kappa += dKappa;
            plStrain.at(1) += dKappa * signum(fullStress.at(1) );
            plStrain.at(2) -= 0.5 * dKappa * signum(fullStress.at(1) );
            plStrain.at(3) -= 0.5 * dKappa * signum(fullStress.at(1) );
            fullStress.at(1) -= dKappa * E * signum(fullStress.at(1) );
        }
    } else {
        // elastic predictor
        FloatArray elStrain = totalStrain;
        elStrain.subtract(plStrain);
        //auto [elStrainDev, elStrainVol] = computeDeviatoricVolumetricSplit(elStrain); // c++17
        auto tmp = computeDeviatoricVolumetricSplit(elStrain);
        auto elStrainDev = tmp.first;
        auto elStrainVol = tmp.second;
        FloatArray trialStressDev = applyDeviatoricElasticStiffness(elStrainDev, G);
        /**************************************************************/
        double trialStressVol = 3 * K * elStrainVol;
        /**************************************************************/

        // store the deviatoric and trial stress (reused by algorithmic stiffness)
        status->letTrialStressDevBe(trialStressDev);
        status->setTrialStressVol(trialStressVol);
        // check the yield condition at the trial state
        double trialS = computeStressNorm(trialStressDev);
        double yieldValue = sqrt(3. / 2.) * trialS - ( computeYieldStress(kappa, gp, tStep) );
        if ( yieldValue > 0. ) {
            // increment of cumulative plastic strain
            double dKappa = yieldValue / ( computeYieldStressPrime(kappa) + 3. * G );
            kappa += dKappa;
            // the following line is equivalent to multiplication by scaling matrix P
            FloatArray dPlStrain = applyDeviatoricElasticCompliance(trialStressDev, 0.5);
            // increment of plastic strain
            plStrain.add(sqrt(3. / 2.) * dKappa / trialS, dPlStrain);
            // scaling of deviatoric trial stress
            trialStressDev.times(1. - sqrt(6.) * G * dKappa / trialS);
        }

        // assemble the stress from the elastically computed volumetric part
        // and scaled deviatoric part

        fullStress = computeDeviatoricVolumetricSum(trialStressDev, trialStressVol);
    }

    // store the effective stress in status
    status->letTempEffectiveStressBe(fullStress);
    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
}

void
MisesMat::performPlasticityReturn_PlaneStress(const FloatArrayF< 3 > &totalStrain, GaussPoint *gp, TimeStep *tStep) const
{
    double E = linearElasticMaterial.give('E', gp);
    double nu = linearElasticMaterial.give('n', gp);
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa;
    FloatArray plStrain, redPlStrain;
    FloatArray fullStress;
    // get the initial plastic strain and initial kappa from the status
    plStrain = status->givePlasticStrain();
    StructuralMaterial::giveReducedSymVectorForm(redPlStrain, plStrain, _PlaneStress);
    kappa = status->giveCumulativePlasticStrain();
    FloatMatrix Ps, Pe;

    FloatArray elStrain = totalStrain;
    elStrain.subtract(redPlStrain);
    FloatArray elStrainDev;
    double elStrainVol;
    // Elastic predictor: Compute elastic trial state
    // Volumetric strain
    double factor = 2. * G / ( K + 4. / 3. * G );
    elStrainVol = ( elStrain.at(1) + elStrain.at(2) ) * factor;
    //Elastic trial deviatoric strain
    elStrainDev = elStrain;
    elStrainDev.at(1) -= elStrainVol / 3.;
    elStrainDev.at(2) -= elStrainVol / 3.;
    // Elastic trial stress components
    FloatArray trialStressDev;
    trialStressDev = elStrainDev;
    trialStressDev.times(2. * G);
    trialStressDev.at(3) /= 2.;
    //applyDeviatoricElasticStiffness(trialStressDev, elStrainDev, G);
    double trialStressVol = 3 * K * elStrainVol;
    FloatArray trialStress(trialStressDev);
    trialStress.at(1) += trialStressVol / 3.;
    trialStress.at(2) += trialStressVol / 3.;
    // Compute yield function value at trial state
    double a1 = ( trialStress.at(1) + trialStress.at(2) ) * ( trialStress.at(1) + trialStress.at(2) );
    double a2 = ( trialStress.at(2) - trialStress.at(1) ) * ( trialStress.at(2) - trialStress.at(1) );
    double a3 = trialStress.at(3) * trialStress.at(3);
    double sigmaY = this->computeYieldStress(kappa, gp, tStep);
    double xi =  1. / 6. * a1 + 0.5 * a2 + 2. * a3;
    //Yield function
    double yieldValue = 0.5 * xi - 1. / 3. * sigmaY * sigmaY;
    // double dGamma ; // not used
    // Check for plastic admissibility
    if ( yieldValue / sigmaY > yieldTol ) {
        // Plastic step: Apply return mapping - use Newton-Raphson algorithm
        //        to solve the plane stress-projected return mapping
        //               equation for the plastic multiplier (Box 9.5)
        double dKappa = 0;
        int iter = 0;
        double f = yieldValue;
        double denom1 = 1.;
        double denom2 = 1.;
        while ( true ) {
            double HiP = this->computeYieldStressPrime(kappa + dKappa * sqrt(2. * xi / 3.) );
            double dXi = -a1 * E / ( 1. - nu ) / 9. / denom1 / denom1 / denom1 - 2. * G * ( a2 + 4. * a3 ) / denom2 / denom2 / denom2;
            double Hbar = 2. * sigmaY * HiP * sqrt(2. / 3.) * ( sqrt(xi) + dKappa * dXi / ( 2. * sqrt(xi) ) );
            double df = 0.5 * dXi - 1. / 3. * Hbar;
            dKappa -= f / df;
            // Compute new residual (yield function value)
            denom1 =  ( 1. + E * dKappa / 3. / ( 1. - nu ) );
            denom2 = ( 1 + 2. * G * dKappa );
            xi = a1 / 6. / denom1 / denom1 +  ( 0.5 * a2 + 2. * a3 ) / denom2 / denom2;
            sigmaY = this->computeYieldStress(kappa + sqrt(2. * xi / 3.) * dKappa, gp, tStep);

            f = 1. / 2. * xi - 1. / 3. * sigmaY * sigmaY;

            if ( fabs(f / sigmaY) < yieldTol ) {
                break;
            }
            iter++;
            if ( iter > 400 ) {
                OOFEM_WARNING("No convergence of the stress return algorithm in MisesMat :: performPlasticityReturn_PlaneStress\n");
                break;
            }
        }
        // update accumulated plastic strain
        // dGamma = dKappa;
        kappa += sqrt(2. * xi / 3.) * dKappa;
        // update stress components:   sigma := A sigma^trial
        double As1 = 3. * ( 1. - nu ) / ( 3 * ( 1 - nu ) + E * dKappa );
        double As2 = 1. / ( 1. + 2. * G * dKappa );
        FloatMatrix A(3, 3);
        A.at(1, 1) = 0.5 * ( As1 + As2 );
        A.at(2, 2) = A.at(1, 1);
        A.at(1, 2) = 0.5 * ( As1 - As2 );
        A.at(2, 1) = A.at(1, 2);
        A.at(3, 3) = As2;
        fullStress.beProductOf(A, trialStress);
        elStrainVol = ( fullStress.at(1) + fullStress.at(2) ) / 3. / K;
        // compute corresponding elastic (engineering) strain components
        redPlStrain.at(1) = totalStrain.at(1) - ( 2. / 3. * fullStress.at(1) - 1. / 3. * fullStress.at(2) ) / 2. / G - elStrainVol / 3;
        redPlStrain.at(2) = totalStrain.at(2) - ( 2. / 3. * fullStress.at(2) - 1. / 3. * fullStress.at(1) ) / 2. / G  - elStrainVol / 3;
        redPlStrain.at(3) = totalStrain.at(3) - fullStress.at(3) / G;
        StructuralMaterial::giveFullSymVectorForm(plStrain, redPlStrain, _PlaneStress);
        // incompresibility condition
        plStrain.at(3) = -( plStrain.at(1) + plStrain.at(2) );
        // store the plastic strain and cumulative plastic strain
    } else {
        fullStress = trialStress;
    }


    // storing into the status
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
    status->letTempEffectiveStressBe(fullStress);
}


double
MisesMat::checkYieldStress(double &dKappa, double kappa, GaussPoint *gp, TimeStep *tStep) const
{
    double yieldStress = 0.;
    if ( hType == 1 ) {
        if ( kappa + dKappa > h_eps.at(h_eps.giveSize() ) ) {
            OOFEM_ERROR("kappa outside range of specified hardening law/n");
        }

        for ( int i = 1; i < h_eps.giveSize(); i++ ) {
            if ( kappa + dKappa >= h_eps.at(i) && kappa + dKappa < h_eps.at(i + 1) && kappa < h_eps.at(i) ) {
                yieldStress = h_function_eps.at(i);
                dKappa = h_eps.at(i) - kappa;
                return yieldStress;
            } else if ( kappa >= h_eps.at(i) && kappa < h_eps.at(i + 1) && kappa + dKappa >= h_eps.at(i) && kappa + dKappa < h_eps.at(i + 1) ) {
                yieldStress = h_function_eps.at(i) + ( kappa + dKappa - h_eps.at(i) ) / ( h_eps.at(i + 1) - h_eps.at(i) ) * ( h_function_eps.at(i + 1) - h_function_eps.at(i) );
                return yieldStress;
            }
        }
    } else {
        OOFEM_ERROR("MisesMat: Should not check yield stress for htype = 0\n");
    }
}

double
MisesMat::computeYieldStress(double kappa, GaussPoint *gp, TimeStep *tStep) const
{
    double yieldStress = 0.;
    if ( hType == 0 ) {
        return this->give('s', gp, tStep) + this->H * kappa; // + ( this->sigInf - this->sig0 ) * (1. - exp(-expD*kappa));
    } else {
        if ( kappa > h_eps.at(h_eps.giveSize() ) ) {
            OOFEM_ERROR("kappa outside range of specified hardening law/n");
        }

        for ( int i = 1; i < h_eps.giveSize(); i++ ) {
            if ( kappa >= h_eps.at(i) && kappa < h_eps.at(i + 1) ) {
                yieldStress = h_function_eps.at(i) + ( kappa - h_eps.at(i) ) / ( h_eps.at(i + 1) - h_eps.at(i) ) * ( h_function_eps.at(i + 1) - h_function_eps.at(i) );
                return yieldStress;
            }
        }
    }
}


double
MisesMat::computeYieldStressPrime(double kappa) const
{
    double yieldStressPrime;
    if ( hType == 0 ) {
        yieldStressPrime = this->H;
        return yieldStressPrime;
    } else {
        if ( kappa > h_eps.at(h_eps.giveSize() ) ) {
            OOFEM_ERROR("kappa outside range of specified hardening law/n");
        }


        for ( int i = 1; i < h_eps.giveSize(); i++ ) {
            if ( kappa >= h_eps.at(i) && kappa < h_eps.at(i + 1) ) {
                yieldStressPrime = ( h_function_eps.at(i + 1) - h_function_eps.at(i) ) / ( h_eps.at(i + 1) - h_eps.at(i) );
                return yieldStressPrime;
            }
        }
    }
}



double
MisesMat::computeDamageParam(double tempKappa) const
{
    if ( tempKappa > 0. ) {
        return omega_crit * ( 1.0 - exp(-a * tempKappa) );
    } else {
        return 0.;
    }
}

double
MisesMat::computeDamageParamPrime(double tempKappa) const
{
    if ( tempKappa >= 0. ) {
        return omega_crit * a * exp(-a * tempKappa);
    } else {
        return 0.;
    }
}


double
MisesMat::computeDamage(GaussPoint *gp,  TimeStep *tStep) const
{
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double dam = status->giveDamage();
    double tempKappa = computeCumPlastStrain(gp, tStep);
    double tempDam = computeDamageParam(tempKappa);
    if ( dam > tempDam ) {
        tempDam = dam;
    }

    return tempDam;
}


double MisesMat::computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    return status->giveTempCumulativePlasticStrain();
}



// returns the consistent (algorithmic) tangent stiffness matrix
FloatMatrixF< 6, 6 >
MisesMat::give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const
{
    // start from the elastic stiffness
    auto d = this->linearElasticMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);
    if ( mode != TangentStiffness ) {
        return d;
    }

    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = tempKappa - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return d;
    }

    // === plastic loading ===

    // yield stress at the beginning of the step
    double sigmaY = computeYieldStress(kappa, gp, tStep);

    // trial deviatoric stress and its norm
    const FloatArrayF< 6 >trialStressDev = status->giveTrialStressDev();
    //double trialStressVol = status->giveTrialStressVol();
    double trialS = computeStressNorm(trialStressDev);

    // one correction term
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( computeYieldStressPrime(kappa) + 3. * G ) * trialS * trialS );
    d += factor1 * dyad(trialStressDev, trialStressDev);

    // another correction term
    double factor2 = factor * dKappa;
    d += factor2 * I_dev6;

    //influence of damage
    //    double omega = computeDamageParam(tempKappa);
    double omega = status->giveTempDamage();
    d *= 1. - omega;
    const FloatArrayF< 6 >effStress = status->giveTempEffectiveStress();
    double omegaPrime = computeDamageParamPrime(tempKappa);
    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + computeYieldStressPrime(kappa) ) / trialS;
    d += scalar * dyad(effStress, trialStressDev);
    return d;
}


FloatMatrixF< 3, 3 >
MisesMat::givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    // start from the elastic stiffness
    auto d = linearElasticMaterial.givePlaneStressStiffMtrx(mmode, gp, tStep);
    if ( mmode != TangentStiffness ) {
        double omega = status->giveTempDamage();
        return d * ( 1. - omega );
    }

    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = tempKappa - kappa;
    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return d;
    }
    // Compute elastoplastic consistent tangent (Box 9.6)
    FloatArray stress, fullStress;
    fullStress = status->giveTempStressVector();
    StructuralMaterial::giveReducedSymVectorForm(stress, fullStress, _PlaneStress);
    // Compute xi
    double xi = 2. / 3. * ( stress.at(1) * stress.at(1) + stress.at(2) * stress.at(2) - stress.at(1) * stress.at(2) ) + 2. * stress.at(3) * stress.at(3);
    // compute dGamma
    double dGamma = dKappa * sqrt(3. / 2. / xi);
    // Hardening slope
    double HiP = this->computeYieldStressPrime(tempKappa);
    // Matrix E components

    double E = linearElasticMaterial.give('E', gp);
    double nu = linearElasticMaterial.give('n', gp);
    double Es1 = 3. * E / ( 3. * ( 1. - nu ) + E * dGamma );
    double Es2 = 2. * G / ( 1. + 2. * G * dGamma );
    double Es3 = Es2 / 2.;
    // Components of the matrix product EP
    double EPs1 = 1. / 3. * Es1;
    double EPs2 = Es2;
    double EPs3 = EPs2;
    FloatMatrix EP(3, 3);
    EP.at(1, 1) = 0.5 * ( EPs1 + EPs2 );
    EP.at(2, 2) = EP.at(1, 1);
    EP.at(1, 2) = 0.5 * ( EPs1 - EPs2 );
    EP.at(2, 1) = EP.at(1, 2);
    EP.at(3, 3) = EPs3;
    // Vector n
    FloatArray n(3);
    n.beProductOf(EP, stress);
    // Scalar alpha
    double denom1 = stress.at(1) * ( 2. / 3. * n.at(1) - 1. / 3. * n.at(2) ) + stress.at(2) * ( 2. / 3. * n.at(2) - 1. / 3. * n.at(1) ) + 2. * stress.at(3) * n.at(3);
    double denom2 = 2. * xi * HiP / ( 3. - 2. * HiP * dGamma );
    double alpha = 1. / ( denom1 + denom2 );
    FloatMatrix correction;
    correction.beDyadicProductOf(n, n);
    correction.times(alpha);

    FloatMatrixF< 3, 3 >answer;
    answer.at(1, 1) = 0.5 * ( Es1 + Es2 );
    answer.at(2, 2) = answer.at(1, 1);
    answer.at(1, 2) = 0.5 * ( Es1 - Es2 );
    answer.at(2, 1) = answer.at(1, 2);
    answer.at(3, 3) = Es3;
    answer -= FloatMatrixF< 3, 3 >(correction);
    //@todo: add damage part of the stiffness
    return answer;
}


FloatMatrixF< 1, 1 >
MisesMat::give1dStressStiffMtrx(MatResponseMode mode,
                                GaussPoint *gp,
                                TimeStep *tStep) const
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double omega = status->giveTempDamage();
    auto elastic = this->linearElasticMaterial.give1dStressStiffMtrx(mode, gp, tStep);
    double E = elastic.at(1, 1);
    if ( mode != TangentStiffness ) {
        return elastic;
    }

    if ( tempKappa <= kappa ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return elastic * ( 1 - omega );
    }

    // === plastic loading ===
    const FloatArray &stressVector = status->giveTempEffectiveStress();
    double stress = stressVector.at(1);
    return {
               ( 1 - omega ) * E * computeYieldStressPrime(kappa) / ( E + computeYieldStressPrime(kappa) ) - computeDamageParamPrime(tempKappa) * E / ( E + computeYieldStressPrime(kappa) ) * stress * signum(stress)
    };
}

#ifdef __OOFEG
#endif

int
MisesMat::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_YieldStrength ) {
        answer.resize(1);
        answer.at(1) = this->give('s', gp, tStep);
        return 1;
    } else {
        return StructuralMaterial::giveIPValue(answer, gp, type, tStep);
    }
}

//=============================================================================

MisesMatStatus::MisesMatStatus(GaussPoint *g) :
    StructuralMaterialStatus(g), plasticStrain(6), tempPlasticStrain(), trialStressD()
{
    stressVector.resize(6);
    strainVector.resize(6);
    effStress.resize(6);
    tempEffStress.resize(6);
}


void
MisesMatStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, "              plastic  ");
    for ( auto &val : this->plasticStrain ) {
        fprintf(file, "%.4e ", val);
    }
    fprintf(file, "\n");

    fprintf(file, "status { ");
    // print damage
    fprintf(file, "damage %.4e", tempDamage);
    // print the cumulative plastic strain
    fprintf(file, ", kappa ");
    fprintf(file, " %.4e", kappa);

    fprintf(file, "}\n");

    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatStatus::initTempStatus()
{
    StructuralMaterialStatus::initTempStatus();

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStressD.clear(); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
MisesMatStatus::updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus::updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    damage = tempDamage;
    trialStressD.clear(); // to indicate that it is not defined any more
}


double
MisesMat::give(int aProperty, GaussPoint *gp, TimeStep *tStep) const
{
    if ( aProperty == 's' ) {
        ///FIXME: const cast workaround, until all methods have been properly marked const properly:

        return sig0.eval({ { "te", giveTemperature(gp, tStep) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp, tStep) );
    }

    return Material::give(aProperty, gp);
}

double MisesMat::giveTemperature(GaussPoint *gp, TimeStep *tStep) const
{
    FieldManager *fm = this->domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr tf;
    int err;
    if ( ( tf = fm->giveField(FT_Temperature) ) ) {
        // temperature field registered
        FloatArray gcoords, answer;
        static_cast< StructuralElement * >( gp->giveElement() )->computeGlobalCoordinates(gcoords, gp->giveNaturalCoordinates() );
        if ( ( err = tf->evaluateAt(answer, gcoords, VM_Total, tStep) ) ) {
            OOFEM_ERROR("tf->evaluateAt failed, element %d, error code %d", gp->giveElement()->giveNumber(), err);
        }
        return answer.at(1);
    }
    return 0.;
}


void
MisesMatStatus::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus::saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
MisesMatStatus::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus::restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem
