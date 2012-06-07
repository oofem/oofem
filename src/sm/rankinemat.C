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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "rankinemat.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"

#include "structuralcrosssection.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
// constructor
RankineMat :: RankineMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    E = 0.;
    nu = 0.;
    H0 = 0.;
    sig0 = 0.;
    delSigY = 0.;
}


// specifies whether a given material mode is supported by this model
int
RankineMat :: hasMaterialModeCapability(MaterialMode mode)
{
    return ( mode == _PlaneStress );
}


// reads the model parameters from the input file
IRResultType
RankineMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    E = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveYoungsModulus();
    nu = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->givePoissonsRatio();

    IR_GIVE_FIELD(ir, sig0, IFT_RankineMat_sig0, "sig0"); // uniaxial yield stress

    H0 = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H0, IFT_RankineMat_h, "h"); // hardening modulus

    plasthardtype = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, plasthardtype, IFT_RankineMat_plasthardtype, "plasthardtype"); // type of plastic hardening (0=linear, 1=exponential)

    delSigY = 0.;
    if ( plasthardtype == 1 ){
      IR_GIVE_FIELD(ir, delSigY, IFT_RankineMat_delsigy, "delsigy"); // final increment of yield stress (at infinite cumulative plastic strain)
    }

    yieldtol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldtol, IFT_RankineMat_yieldtol, "yieldtol"); // relative tolerance in yield condition

    a = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, a, IFT_RankineMat_a, "a"); // coefficient in damage law

    double gf = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, gf, IFT_RankineMat_a, "gf"); // dissipated energy per unit VOLUME

    if ( ( a != 0. ) && ( gf != 0 ) ) {
        OOFEM_ERROR("RankineMat: parameters a and gf cannot be prescribed simultaneously");
    }

    if ( gf > 0. ) {
        // evaluate parameter "a" from given "gf"
        double A = H0 * ( 1. + H0 / E );
        double B = sig0 * ( 1. + H0 / E );
        double C = 0.5 * sig0 * sig0 / E - gf;
        if ( C >= 0. ) {
            OOFEM_ERROR("RankineMat: parameter gf is too low");
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
RankineMat :: giveRealStressVector(FloatArray &answer,
                                   MatResponseForm form,
                                   GaussPoint *gp,
                                   const FloatArray &totalStrain,
                                   TimeStep *atTime)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode != _PlaneStress ) {
        OOFEM_ERROR("RankineMat::giveRealStressVector : unknown material response mode");
    }

    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);

    // initialization
    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // elastoplasticity
    this->performPlasticityReturn(gp, totalStrain, mode);

    // damage
    double omega = computeDamage(gp, atTime);
    status->giveTempEffectiveStress(answer);
    answer.times(1. - omega);

    // store variables in status
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    double gf = sig0 * sig0 / E; // only estimated, but OK for this purpose
    status->computeWork(gp, mode, gf);
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
    if ( plasthardtype == 0 ) { // linear hardening
        return sig0 + H0 * kappa;
    } else { // exponential hardening
        if ( delSigY == 0. ) {
        return sig0;
        } else {
        return sig0 + delSigY * ( 1. - exp(-H0*kappa/delSigY) );
        }
    }
}

double
RankineMat :: evalPlasticModulus(const double kappa)
{
    if ( plasthardtype == 0 ) { // linear hardening
        return H0;
    } else { // exponential hardening
        return H0 * exp(-H0*kappa/delSigY);
    }
}


// computes the stress according to elastoplasticity
// (return of trial stress to the yield surface)
void
RankineMat :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, MaterialMode mode)
{
  double kappa, tempKappa, H;
    FloatArray reducedStress;
    FloatArray strain, tempPlasticStrain;
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);

    // get the initial plastic strain and initial kappa from the status
    status->givePlasticStrain(tempPlasticStrain);
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
                OOFEM_ERROR("RankineMat::giveRealStressVector : no convergence of regular stress return algorithm");
            }

            H = evalPlasticModulus(tempKappa);
            double ddKappa = f / ( Enu + H );
            sigPrinc.at(1) -= Enu * ddKappa;
            sigPrinc.at(2) -= nu * Enu * ddKappa;
            tempKappa += ddKappa;
            f = evalYieldFunction(sigPrinc, tempKappa);
        } while ( fabs(f) > yieldtol*sig0 );

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
                    OOFEM_ERROR("RankineMat::giveRealStressVector : no convergence of vertex stress return algorithm");
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
            } while ( fabs(f) > yieldtol*sig0 );

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

    // store the effective stress, plastic strain and cumulative plastic strain
    status->letTempEffectiveStressBe(finalStress);
    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->setTempCumulativePlasticStrain(tempKappa);
    status->setTangentShearStiffness(tanG);
}

double
RankineMat :: computeDamageParam(double tempKappa)
{
    double tempDam;
    if ( tempKappa > 0. ) {
        tempDam = 1.0 - exp(-a * tempKappa);
    } else {
        tempDam = 0.;
    }

    return tempDam;
}

double
RankineMat :: computeDamageParamPrime(double tempKappa)
{
    double tempDam;
    if ( tempKappa >= 0. ) {
        tempDam = a * exp(-a * tempKappa);
    } else {
        tempDam = 0.;
    }

    return tempDam;
}

double
RankineMat :: computeDamage(GaussPoint *gp,  TimeStep *atTime)
{
    double tempKappa, dam;
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);
    dam = status->giveDamage();
    computeCumPlastStrain(tempKappa, gp, atTime);
    double tempDam = computeDamageParam(tempKappa);
    if ( dam > tempDam ) {
        tempDam = dam;
    }

    return tempDam;
}

void RankineMat :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *atTime)
{
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);
    tempKappa = status->giveTempCumulativePlasticStrain();
}

// returns the consistent (algorithmic) stiffness matrix
void
RankineMat :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *atTime)
{
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double gprime = computeDamageParamPrime(tempKappa);
    evaluatePlaneStressStiffMtrx(answer, form, mode, gp, atTime, gprime);
}

// this method is also used by the gradient version,
// with gprime replaced by gprime*m and evaluated for kappa hat
void
RankineMat :: evaluatePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime, double gprime)
{
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(gp);
    if ( mode == ElasticStiffness || mode == SecantStiffness ) {
        // start from the elastic stiffness
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
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
        this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
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
    StressVector effStress(_PlaneStress);
    status->giveTempEffectiveStress(effStress);
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
    FloatMatrix T(3, 3), TT;
    givePlaneStressVectorTranformationMtrx(T, nPrinc, true);
    TT.beTranspositionOf(T);
    answer.rotatedWith(TT);
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
    StressVector effStress(_PlaneStress);
    status->giveTempEffectiveStress(effStress);
    effStress.computePrincipalValDir(sigPrinc, nPrinc);

    FloatMatrix T(3, 3);
    givePlaneStressVectorTranformationMtrx(T, nPrinc, true);
    answer.beProductOf(T, eta);
}

int
RankineMat :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    RankineMatStatus *status = ( RankineMatStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_PlasticStrainTensor ) {
        answer  = * status->givePlasDef();
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
        answer.at(1) = status->giveDamage();
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
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

InternalStateValueType
RankineMat :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_PlasticStrainTensor ) || ( type == IST_DamageTensor ) ) {
        return ISVT_TENSOR_S3;
    } else if ( ( type == IST_CumPlasticStrain ) || ( type == IST_DamageScalar ) ) {
        return ISVT_SCALAR;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_DissWorkDensity ) {
        return ISVT_SCALAR;
    } else if ( type == IST_StressWorkDensity ) {
        return ISVT_SCALAR;
    } else if ( type == IST_FreeEnergyDensity ) {
        return ISVT_SCALAR;

#endif
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
RankineMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_PlasticStrainTensor ) || ( type == IST_DamageTensor ) ) {
        if ( ( mmode == _3dMat ) || ( mmode == _3dMat_F ) ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(4) = 4;
            answer.at(5) = 5;
            answer.at(6) = 6;
            return 1;
        } else if ( mmode == _PlaneStress || mmode == _PlaneStressGrad ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
            return 1;
        } else if ( mmode == _1dMat ) {
            answer.resize(1);
            answer.at(1) = 1;
            return 1;
        }
    } else if ( type == IST_CumPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( ( type == IST_DissWorkDensity ) || ( type == IST_StressWorkDensity ) || ( type == IST_FreeEnergyDensity ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;

#endif
    }

    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}


int
RankineMat :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( ( type == IST_PlasticStrainTensor ) || ( type == IST_DamageTensor ) ) {
        MaterialMode mode = gp->giveMaterialMode();
        if ( mode == _3dMat || mode == _3dMat_F ) {
            return 6;
        } else if ( mode == _PlaneStrain ) {
            return 4;
        } else if ( mode == _PlaneStress || mode == _PlaneStressGrad ) {
            return 3;
        } else if ( mode == _1dMat ) {
            return 1;
        } else {
            return 0;
        }
    } else if ( type == IST_CumPlasticStrain ) {
        return 1;
    } else if ( type == IST_DamageScalar ) {
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( ( type == IST_StressWorkDensity ) ||
               ( type == IST_DissWorkDensity ) || ( type == IST_FreeEnergyDensity ) ) {
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValueSize(type, gp);
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
    /*
     * // print the plastic strain
     *  n = plasticStrain.giveSize();
     *  fprintf(file, " plastic_strains ");
     *  for ( i = 1; i <= n; i++ ) {
     *      fprintf( file, " % .4e", plasticStrain.at(i) );
     *  }
     */
    // print the cumulative plastic strain
    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void RankineMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
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
RankineMatStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

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
RankineMatStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    // write plastic strain (vector)
    if ( ( iores = plasticStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write cumulative plastic strain (scalar)
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write damage (scalar)
    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->write(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK;
}


contextIOResultType
RankineMatStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
    if ( ( iores = plasticStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read cumulative plastic strain (scalar)
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read damage (scalar)
    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream->read(& stressWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& dissWork, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif

    return CIO_OK; // return success
}


#ifdef keep_track_of_dissipated_energy
void
RankineMatStatus :: computeWork(GaussPoint *gp, MaterialMode mode, double gf)
{
    // int n = deps.giveSize(); // would not work for gradient version
    int n;
    if ( mode == _PlaneStress || mode == _PlaneStressGrad ) {
        n = 3;
    } else {
        _error("Inappropriate material mode in RankineMatStatus :: computeWork\n");
    }

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
