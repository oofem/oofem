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

// This code is based on the anisotropic damage model proposed by Desmorat, Gatuingt and Ragueneau in
// their paper "Nonlocal anisotropic damage model and related computational aspects for quasi-brittle material"
// published in Engineering Fracture Mechanics 74 (2007) 1539-1560.

#include "anisodamagemodel.h"
#include "../sm/Materials/structuralmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(AnisotropicDamageMaterial);

AnisotropicDamageMaterial :: AnisotropicDamageMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    E = 0.;
    nu = 0.;
    equivStrainType = EST_Unknown;
    damageLawType = DLT_Unknown;
    kappa0 = 0.;
    kappaf = 0.;
    aA = 0.;
}

AnisotropicDamageMaterial :: ~AnisotropicDamageMaterial()
//
// destructor
//
{
    delete linearElasticMaterial;
}

int
AnisotropicDamageMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports the given mode
//
{
    return mode == _3dMat || mode == _PlaneStress;
    //return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}


//********************************************************
// plane stress implementation by Milan Jirasek
//********************************************************

void
AnisotropicDamageMaterial :: giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp,
                                                              const FloatArray &totalStrain, TimeStep *atTime)
//
// special version of the stress-evaluation algorithm applicable under plane stress
// based on the report by Jirasek & Suarez, 25 April 2014
//
// uses the following special methods:
//   computePrincValDir2D  ...  evaluation of eigenvalues and eigenvector of a symmetric 2x2 matrix
//   computeDamage  ...   application of the damage law
//     computeTraceD(double)  ...  function "g" evaluating trace(D) from kappa
//     checkPrincVal2D  ...  checking whether principal values are <= 1
//   computeOutOfPlaneStrain  ...  evaluation of epsZ for given in-plane strain and damage
//   computeDimensionlessOutOfPlaneStress  ...  evaluation of scaled sigZ for given strain and damage
//   computeInplaneStress  ...  evaluation of in-plane stress for given strain and damage
//
// Note: Only Mazars' equivalent strain is implemented for this case,
//       but the damage law can be varied and there are 4 cases implemented in computeTraceD(double).
//       There is another method computeTraceD(...), which is used by the 3D algorithm.
{
#define AD_TOLERANCE 1.e-10 // convergence tolerance for the internal iteration used under plane stress
#define AD_MAXITER 20       // maximum number of internal iterations used under plane stress
    //this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    // subtract the stress-independent part of strains (e.g. due to temperature)
    FloatArray eps;
    this->giveStressDependentPartOfStrainVector(eps, gp, totalStrain, atTime, VM_Total);

    // compute in-plane principal strains: epsilon1 and epsilon2
    // and the components of the first principal strain direction: ceps and seps
    double epsilon1, epsilon2, ceps, seps;
    this->computePrincValDir2D(epsilon1, epsilon2, ceps, seps, eps.at(1), eps.at(2), eps.at(3) / 2.);

    // compute the equivalent strain resulting from the in-plane strains only: equivStrainInPlane
    // only Mazars' strain implemented so far
    double equivStrainInPlane = 0.;
    if ( epsilon1 > 0. ) {
        equivStrainInPlane += epsilon1 * epsilon1;
        if ( epsilon2 > 0. ) { // note that we always have epsilon2<=epsilon1
            equivStrainInPlane += epsilon2 * epsilon2;
        }
    }
    equivStrainInPlane = sqrt(equivStrainInPlane);

    // compute ez0 = maximum value of the out-of-plane strain that still has no influence on damage
    // formula (85)
    double ez0 = 0.;
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double kappa = status->giveKappa();
    if ( equivStrainInPlane < kappa ) {
        ez0 = sqrt(kappa * kappa - equivStrainInPlane * equivStrainInPlane);
    }

    // compute the trial out-of-plane strain, assuming that the equivalent strain is equal to equivStrainInPlane
    FloatMatrix Dn = status->giveDamage();
    FloatMatrix tempDamage = Dn;
    // update damage, if needed
    if ( equivStrainInPlane > kappa ) {
        // formula (76)
        this->computeDamage(tempDamage, Dn, kappa, epsilon1, epsilon2, ceps, seps, 0.);
    }
    // compute the first trial value (with volumetric compression assumed, i.e., bulk damage deactivated)
    // formula (79)
    double ez1 = this->computeOutOfPlaneStrain(eps, tempDamage, false);

    double strainTraceInPlane = eps.at(1) + eps.at(2);
    if ( ( ez1 + strainTraceInPlane > 0. ) || ( ez1 > ez0 ) ) {     // first trial value is not admissible
        bool iteration_needed = true;
        if ( ez0 + strainTraceInPlane > 0. ) {
            // compute the second trial value (with volumetric tension assumed, i.e., bulk damage activated)
            ez1 = this->computeOutOfPlaneStrain(eps, tempDamage, true);
            iteration_needed = ( ( ez1 > ez0 ) || ( ez1 + strainTraceInPlane < 0. ) );
        }
        if ( iteration_needed ) {     // the second trial value is not admissible
            // iterative procedure - out-of-plane strain does have an effect on damage
            // scaled formula (78)
            double s0 = this->computeDimensionlessOutOfPlaneStress(eps, ez0, tempDamage);
            this->computeDamage(tempDamage, Dn, kappa, epsilon1, epsilon2, ceps, seps, ez1);
            double s1 = this->computeDimensionlessOutOfPlaneStress(eps, ez1, tempDamage);
            int count = 0;
            while ( fabs(s1) > AD_TOLERANCE * this->kappa0 ) {
                if ( ++count > AD_MAXITER ) {
                    OOFEM_ERROR("No convergence in AnisotropicDamageMaterial :: giveRealStressVector for the plane stress case\n");
                }
                double ez2 = ( ez1 * s0 - ez0 * s1 ) / ( s0 - s1 );
                this->computeDamage(tempDamage, Dn, kappa, epsilon1, epsilon2, ceps, seps, ez2);
                double s2 = this->computeDimensionlessOutOfPlaneStress(eps, ez2, tempDamage);
                ez0 = ez1;
                s0 = s1;
                ez1 = ez2;
                s1 = s2;
            }
        }
    }
    status->setTempStrainZ(ez1);
    status->setTempDamage(tempDamage);
    double equivStrain = sqrt( equivStrainInPlane * equivStrainInPlane + macbra(ez1) * macbra(ez1) );
    if ( equivStrain > kappa ) {
        status->setTempKappa(equivStrain);
    } else {
        status->setTempKappa(kappa);
    }
    // formulae (93)-(100)
    computeInplaneStress(answer, eps, ez1, tempDamage);
    //this->correctBigValues(stressTensor); // ???
    status->setTempDamage(tempDamage);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
}

void
AnisotropicDamageMaterial :: computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy)
//
// computes the principal values and directions of a symmetric second-order tensor in 2D
// input: Dx, Dy, Dxy ... components of the tensor wrt global coordinates
// output: D1, D2 ... ordered principal values, D1>=D2
// output: c, s ... components of the unit principal vector associated with D1
//                  (cosine and sine of the angle between the major principal direction and the global x-axis)
{
    // based on formulae (88)-(89) from the report by Jirasek & Suarez, 25 April 2014
    double aux1 = ( Dx + Dy ) / 2.;
    double aux2 = ( Dx - Dy ) / 2.;
    double aux3 = sqrt(aux2 * aux2 + Dxy * Dxy);
    D1 = aux1 + aux3;
    D2 = aux1 - aux3;
    // formulae (90)-(92) and the two cases preceding them
    c = 1.;
    s = 0.;         // cases 1 and 2a
    if ( Dxy != 0. ) {  // case 3
        double t = ( D1 - Dx ) / Dxy;
        c = 1. / sqrt(1. + t * t);
        s = c * t;
    } else if ( Dx < Dy ) { // case 2b
        c = 0.;
        s = 1.;
    }
    return;
}

bool
AnisotropicDamageMaterial :: checkPrincVal2D(double Dx, double Dy, double Dxy)
//
// checks whether both eigenvalues of a symmetric second-order tensor in 2D are <= 1
//
{
    if ( Dx + Dy > 2. ) {
        return false;
    }
    if ( Dx + Dy > 1. + Dx * Dy - Dxy * Dxy ) {
        return false;
    }
    return true;
}

void
AnisotropicDamageMaterial :: computeDamage(FloatMatrix &tempDamage, const FloatMatrix &damage, double kappa, double eps1, double eps2, double ceps, double seps, double epsZ)
//
// evaluates the final damage "tempDamage" from the initial damage "damage", initial history variable "kappa"
//    and the final in-plane strain
//    (given by principal values eps1>=eps2, components of first principal direction ceps and seps)
//    and out-of-plane strain
//
{
    // set final damage to initial damage
    tempDamage = damage;

    // evaluate final equivalent strain using Mazars' definition
    double tempEpsEq = 0.;
    if ( eps1 > 0. ) {
        tempEpsEq += eps1 * eps1;
        if ( eps2 > 0. ) {
            tempEpsEq += eps2 * eps2;
        }
    }
    if ( epsZ > 0. ) {
        tempEpsEq += epsZ * epsZ;
    }
    tempEpsEq = sqrt(tempEpsEq);

    // check for damage growth
    if ( tempEpsEq <= kappa ) { // no damage growth
        return;
    }

    // apply incremental damage evolution law
    // formula (75)
    double deltaLambda = ( computeTraceD(tempEpsEq) - computeTraceD(kappa) ) / tempEpsEq / tempEpsEq;
    double eps1p = macbra(eps1);
    double eps2p = macbra(eps2);
    double epsZp = macbra(epsZ);
    double aux1 = deltaLambda * eps1p * eps1p;
    double aux2 = deltaLambda * eps2p * eps2p;
    // formula (76), with the square of positive strain expressed using the spectral decomposition
    tempDamage.at(1, 1) += aux1 * ceps * ceps + aux2 * seps * seps;
    tempDamage.at(1, 2) += ( aux1 - aux2 ) * ceps * seps;
    tempDamage.at(2, 1) = tempDamage.at(1, 2);
    tempDamage.at(2, 2) += aux1 * seps * seps + aux2 * ceps * ceps;
    tempDamage.at(3, 3) += deltaLambda * epsZp * epsZp;

    // treat the case when the out-of-plane damage exceeds 1
    if ( tempDamage.at(3, 3) > 1. ) {
        tempDamage.at(3, 3) = 1.;
    }
    // check whether in-plane principal damages do not exceed 1
    if ( this->checkPrincVal2D( tempDamage.at(1, 1), tempDamage.at(2, 2), tempDamage.at(1, 2) ) ) {
        return;
    }
    // treat the case when the in-plane damage exceeds 1
    double D1, D2, cdam, sdam;
    // compute principal values and direction of the initial damage
    this->computePrincValDir2D( D1, D2, cdam, sdam, damage.at(1, 1), damage.at(2, 2), damage.at(1, 2) );
    // compute damage increment components 11 and 22 in the principal damage coordinates
    double dDx = tempDamage.at(1, 1) - damage.at(1, 1);
    double dDy = tempDamage.at(2, 2) - damage.at(2, 2);
    double dDxy = tempDamage.at(1, 2) - damage.at(1, 2);
    double dD11 = cdam * cdam * dDx + sdam * sdam * dDy + 2. * cdam * sdam * dDxy;
    double dD22 = sdam * sdam * dDx + cdam * cdam * dDy - 2. * cdam * sdam * dDxy;
    // compute new principal damages, truncated to 1
    D1 += dD11;
    D2 += dD22;
    if ( D1 > 1. ) {
        D1 = 1.;
    }
    if ( D2 > 1. ) {
        D2 = 1.;
    }
    // evaluate global components of the new damage tensor
    tempDamage.at(1, 1) = cdam * cdam * D1 + sdam * sdam * D2;
    tempDamage.at(2, 2) = sdam * sdam * D1 + cdam * cdam * D2;
    tempDamage.at(1, 2) = cdam * sdam * ( D1 - D2 );
    tempDamage.at(2, 1) = tempDamage.at(1, 2);
}

double
AnisotropicDamageMaterial :: computeTraceD(double equivStrain)
{
    double knu, aux;
    double answer = 0.;
    if ( equivStrain > this->kappa0 ) {
        switch ( this->damageLawType ) {
        case DLT_Desmorat1:
            answer = ( equivStrain - this->kappa0 ) / ( this->kappaf - this->kappa0 );
            break;
        case DLT_Desmorat2:
            answer = this->aA * ( atan(equivStrain / this->kappaf) - atan(this->kappa0 / this->kappaf) );
            break;
        case DLT_Linear:
            knu = 2. * ( 1. + this->nu ) / 9.;
            answer =  this->kappaf * ( equivStrain - this->kappa0 ) / ( equivStrain * ( this->kappaf - this->kappa0 ) - knu * this->kappa0 * ( this->kappaf - equivStrain ) );
            break;
        case DLT_Exponential:
            knu = 2. * ( 1. + this->nu ) / 9.;
            aux = exp( -( equivStrain - this->kappa0 ) / ( this->kappaf - this->kappa0 ) );
            answer = ( equivStrain - aux * this->kappa0 ) / ( equivStrain - aux * knu * this->kappa0 );
            break;
        default:
            OOFEM_ERROR("Unknown type of damage law.\n");
        }
    }
    return answer;
}

double
AnisotropicDamageMaterial :: computeOutOfPlaneStrain(const FloatArray &inplaneStrain, const FloatMatrix &dam, bool tens_flag)
//
// evaluate the out-of-plane strain from the condition of zero out-of-plane stress for the given damage
// based on formula (79) from the report by Jirasek & Suarez, 25 April 2014
//
{
    // evaluate principal damage values and directions (in-plane)
    double D1, D2, cdam, sdam;
    this->computePrincValDir2D( D1, D2, cdam, sdam, dam.at(1, 1), dam.at(2, 2), dam.at(1, 2) );

    // evaluate normal strain components in directions of principal damage (in-plane)
    // formulae (93)-(94)
    double eps11 = cdam * cdam * inplaneStrain.at(1) + cdam *sdam *inplaneStrain.at(3) + sdam *sdam *inplaneStrain.at(2);
    double eps22 = sdam * sdam * inplaneStrain.at(1) - cdam *sdam *inplaneStrain.at(3) + cdam *cdam *inplaneStrain.at(2);

    // out-of-plane damage
    double Dz = dam.at(3, 3);
    // formula (4)
    double Bv = 1.;
    if ( tens_flag ) {
        Bv = macbra(1. - D1 - D2 - Dz);
    }
    // evaluate auxiliary constants
    // Q corresponds to K*B*Bv, Q1 corresponds to 2*G*Bz*B1 and Q2 corresponds to 2*G*Bz*B2
    // but all of them divided by E and multiplied by 3*(1+nu)*(1-2*nu)
    double Q = ( 3. - D1 - D2 - Dz ) * Bv * ( 1. + nu );
    double Q1 = 3. * ( 1. - D1 ) * ( 1. - Dz ) * ( 1. - 2. * nu );
    double Q2 = 3. * ( 1. - D2 ) * ( 1. - Dz ) * ( 1. - 2. * nu );

    // evaluate out-of-plane strain which would give zero out-of-plane stress (at the given damage)
    // formula (79)
    double answer = ( ( Q1 - Q ) * eps11 + ( Q2 - Q ) * eps22 ) / ( Q + Q1 + Q2 );
    return answer;
}

double
AnisotropicDamageMaterial :: computeDimensionlessOutOfPlaneStress(const FloatArray &inplaneStrain, double epsZ, const FloatMatrix &dam)
//
// evaluate the dimensionless out-of-plane stress (i.e., stress divided by E and multiplied by 3*B*(1+nu)*(1-2*nu))
// for the given final in-plane strain, out-of-plane strain and damage (which is already updated for this strain)
// based on formula (78) from the report by Jirasek & Suarez, 25 April 2014
//
{
    // evaluate principal damage values and directions (in-plane)
    double D1, D2, cdam, sdam;
    this->computePrincValDir2D( D1, D2, cdam, sdam, dam.at(1, 1), dam.at(2, 2), dam.at(1, 2) );

    // evaluate normal strain components in directions of principal damage (in-plane)
    // formulae (93)-(94)
    double eps11 = cdam * cdam * inplaneStrain.at(1) + cdam *sdam *inplaneStrain.at(3) + sdam *sdam *inplaneStrain.at(2);
    double eps22 = sdam * sdam * inplaneStrain.at(1) - cdam *sdam *inplaneStrain.at(3) + cdam *cdam *inplaneStrain.at(2);

    // out-of-plane damage
    double Dz = dam.at(3, 3);
    // formula (4)
    double Bv = 1.;
    double strainTrace = inplaneStrain.at(1) + inplaneStrain.at(2) + epsZ;
    if ( strainTrace > 0. ) {
        Bv = macbra(1. - D1 - D2 - Dz);
    }
    // evaluate auxiliary constants
    // Q corresponds to K*B*Bv, Q1 corresponds to 2*G*Bz*B1 and Q2 corresponds to 2*G*Bz*B2
    // but all of them divided by E and multiplied by 3*(1+nu)*(1-2*nu)
    double Q = ( 3. - D1 - D2 - Dz ) * Bv * ( 1. + nu );
    double Q1 = 3. * ( 1. - D1 ) * ( 1. - Dz ) * ( 1. - 2. * nu );
    double Q2 = 3. * ( 1. - D2 ) * ( 1. - Dz ) * ( 1. - 2. * nu );
    // formula (78) divided by E and multiplied by 3*B*(1+nu)*(1-2*nu)
    double answer = ( Q - Q1 ) * eps11 + ( Q - Q2 ) * eps22 + ( Q + Q1 + Q2 ) * epsZ;
    //printf("%g %g %g %g %g %g %g\n",epsZ,eps11,eps22,Q,Q1,Q2,answer);
    return answer;
}

void
AnisotropicDamageMaterial :: computeInplaneStress(FloatArray &inplaneStress, const FloatArray &inplaneStrain, double epsZ, const FloatMatrix &dam)
//
// evaluate the in-plane stress components
// for the given final in-plane strain, out-of-plane strain and damage (which is already updated for this strain)
//
{
    // evaluate principal damage values and directions (in-plane)
    double D1, D2, cdam, sdam;
    this->computePrincValDir2D( D1, D2, cdam, sdam, dam.at(1, 1), dam.at(2, 2), dam.at(1, 2) );

    // evaluate strain components with respect to the principal damage coordinates (in-plane)
    // formulae (93)-(95), with both sides of (95) multiplied by factor 2
    double eps11 = cdam * cdam * inplaneStrain.at(1) + cdam *sdam *inplaneStrain.at(3) + sdam *sdam *inplaneStrain.at(2);
    double eps22 = sdam * sdam * inplaneStrain.at(1) - cdam *sdam *inplaneStrain.at(3) + cdam *cdam *inplaneStrain.at(2);
    double gam12 = 2. * cdam * sdam * ( inplaneStrain.at(2) - inplaneStrain.at(1) ) + ( cdam * cdam - sdam * sdam ) * inplaneStrain.at(3);

    /* OLD STYLE
     * // evaluate in-plane effective stress
     * double Eaux = E / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
     * double sig11eff = Eaux * ( ( 1. - nu ) * eps11 + nu * ( eps22 + epsZ ) );
     * double sig22eff = Eaux * ( ( 1. - nu ) * eps22 + nu * ( eps11 + epsZ ) );
     * double sigZeff  = Eaux * ( ( 1. - nu ) * epsZ  + nu * ( eps11 + eps22 ) );
     * double G = E / ( 2. * ( 1 + nu ) );
     * double sig12eff =  G * gam12;
     *
     * // evaluate auxiliary constants
     * double Dz = dam.at(3, 3);
     * double Bv = 1.;
     * double strainTrace = eps11 + eps22 + epsZ;
     * double damTrace = D1 + D2 + Dz;
     * if ( strainTrace > 0. ) {
     *  Bv = macbra(1. - damTrace);
     * }
     * double K = E / ( 3. * ( 1 - 2 * nu ) );
     * double sigm = Bv * K * strainTrace;
     * double Bsig = ( 1. - D1 ) * sig11eff + ( 1. - D2 ) * sig22eff + ( 1. - Dz ) * sigZeff;
     * Bsig /= ( 3. - damTrace );
     *
     * // evaluate nominal stress from effective stress (in principal damage coordinate system)
     * double sig11 = ( 1. - D1 ) * ( sig11eff - Bsig ) + sigm;
     * double sig22 = ( 1. - D2 ) * ( sig22eff - Bsig ) + sigm;
     * double sig12 = sqrt( ( 1. - D1 ) * ( 1. - D2 ) ) * sig12eff;
     */

    // evaluate auxiliary constants
    double K = E / ( 3. * ( 1 - 2 * nu ) );
    double G = E / ( 2. * ( 1 + nu ) );
    double Bv = 1.;
    double strainTrace = eps11 + eps22 + epsZ;
    double damTrace = D1 + D2 + dam.at(3, 3);
    if ( strainTrace > 0. ) {
        Bv = macbra(1. - damTrace);
    }
    double B = 3. - damTrace;
    double B1 = 1. - D1;
    double B2 = 1. - D2;
    double Bz = 1. - dam.at(3, 3);

    // evaluate components of nominal in-plane stress wrt principal damage coordinates
    // formulae (96)-(97)
    double sigm = Bv * K * strainTrace;
    double sig11 = 2. * G * B1 * ( B2 * ( eps11 - eps22 ) + Bz * ( eps11 - epsZ ) ) / B + sigm;
    double sig22 = 2. * G * B2 * ( B1 * ( eps22 - eps11 ) + Bz * ( eps22 - epsZ ) ) / B + sigm;
    double sig12 = G * sqrt(B1 * B2) * gam12;

    // evaluate global components of nominal in-plane stress
    // formulae (98)-(100)
    inplaneStress.resize(3);
    inplaneStress.at(1) = cdam * cdam * sig11 - 2. * cdam * sdam * sig12 + sdam * sdam * sig22;
    inplaneStress.at(2) = sdam * sdam * sig11 + 2. * cdam * sdam * sig12 + cdam * cdam * sig22;
    inplaneStress.at(3) = cdam * sdam * ( sig11 - sig22 ) + ( cdam * cdam - sdam * sdam ) * sig12;
}

//********************************************************
// end of the plane stress implementation by Milan Jirasek
//********************************************************

void
AnisotropicDamageMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver
// computed from the state at the beginning of the step and strain at the end of the step
//
{
    //this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    MaterialMode mode = gp->giveMaterialMode();
    // subtract the stress-independent part of strains (e.g. due to temperature)
    FloatArray reducedTotalStrainVector;
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);

    // evaluate stress under general triaxial stress conditions
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    IsotropicLinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatMatrix de, tempDamage;
    double equivStrain, kappa = 0.0, tempKappa = 0.0, traceTempD;
    FloatArray eVals, effectiveStressVector, fullEffectiveStressVector, stressVector;
    FloatMatrix eVecs, effectiveStressTensor, stressTensor, strainTensor;
    FloatMatrix Dn = status->giveDamage();

    this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);
    kappa = this->computeKappa(Dn);
    if ( equivStrain <= kappa ) { // damage does not grow
        tempDamage = Dn;
        tempKappa = kappa;
    } else {
        this->computeDamageTensor(tempDamage, gp, reducedTotalStrainVector, equivStrain, atTime);
        tempKappa = computeKappa(tempDamage);
    }

    if ( ( equivStrain <= kappa ) && ( kappa <= this->kappa0 ) ) {      // elastic behavior
        lmat->giveStiffnessMatrix(de, ElasticStiffness, gp, atTime);
        effectiveStressVector.beProductOf(de, reducedTotalStrainVector);
        answer = effectiveStressVector;
    } else {
        //      StructuralMaterial :: giveFullSymVectorForm(fullEffectiveStressVector, effectiveStressVector, mode);
        //      effectiveStressTensor.beMatrixForm(fullEffectiveStressVector);
        strainTensor.resize(3, 3);
        strainTensor.zero();
        if ( mode == _PlaneStress ) {
            /*
             * this->computePlaneStressStrain(strainTensor, tempDamage, reducedTotalStrainVector, gp, atTime);
             * effectiveStressTensor.resize(3, 3);
             * effectiveStressTensor.zero();
             * double aux;
             * aux = E / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
             * effectiveStressTensor.at(1, 1) = aux * ( ( 1 - nu ) * strainTensor.at(1, 1)     +       nu * ( strainTensor.at(2, 2) + strainTensor.at(3, 3) ) );
             * effectiveStressTensor.at(2, 2) = aux * ( ( 1 - nu ) * strainTensor.at(2, 2)     +       nu * ( strainTensor.at(1, 1) + strainTensor.at(3, 3) ) );
             * effectiveStressTensor.at(3, 3) = aux * ( ( 1 - nu ) * strainTensor.at(3, 3)     +       nu * ( strainTensor.at(1, 1) + strainTensor.at(2, 2) ) );
             * effectiveStressTensor.at(1, 2) = ( E / ( 1 + nu ) ) * strainTensor.at(1, 2);
             * effectiveStressTensor.at(2, 1) = effectiveStressTensor.at(1, 2);
             */
        } else {
            strainTensor.at(1, 1) = reducedTotalStrainVector.at(1);
            strainTensor.at(2, 2) = reducedTotalStrainVector.at(2);
            strainTensor.at(3, 3) = reducedTotalStrainVector.at(3);
            strainTensor.at(2, 3) =  strainTensor.at(3, 2) = reducedTotalStrainVector.at(4) / 2.0;
            strainTensor.at(1, 3) = strainTensor.at(3, 1) = reducedTotalStrainVector.at(5) / 2.0;
            strainTensor.at(1, 2) = strainTensor.at(2, 1) = reducedTotalStrainVector.at(6) / 2.0;
            lmat->giveStiffnessMatrix(de, ElasticStiffness, gp, atTime);
            effectiveStressVector.beProductOf(de, reducedTotalStrainVector);
            StructuralMaterial :: giveFullSymVectorForm(fullEffectiveStressVector, effectiveStressVector, mode);
            effectiveStressTensor.beMatrixForm(fullEffectiveStressVector);
        }
        /*        lmat->giveStiffnessMatrix(de, ElasticStiffness, gp, atTime);
         *      effectiveStressVector.beProductOf(de, reducedTotalStrainVector);
         *      StructuralMaterial :: giveFullSymVectorForm(fullEffectiveStressVector, effectiveStressVector, mode);
         *      effectiveStressTensor.beMatrixForm(fullEffectiveStressVector);*/
        //traceTempD=tempDamage.at(1,1)+tempDamage.at(2,2)+tempDamage.at(3,3);
        double effectiveStressTrace = effectiveStressTensor.at(1, 1) + effectiveStressTensor.at(2, 2) + effectiveStressTensor.at(3, 3);
        FloatMatrix Part1, Part2, Part3;
        // First term of the equation 53 of the reference paper****************************************************************************************
        FloatMatrix AuxMatrix;
        // Compute (1-D) (called ImD here)
        FloatMatrix ImD, sqrtImD;
        ImD.resize(3, 3);
        ImD.zero();
        ImD.at(1, 1) = ImD.at(2, 2) = ImD.at(3, 3) = 1;
        ImD.subtract(tempDamage);
        // Compute the square root of (1-D), needed in the equation 53 of the reference paper
        //int checker1 = this->checkSymmetry(ImD);
        ImD.jaco_(eVals, eVecs, 40);
        sqrtImD.resize(3, 3);
        sqrtImD.zero();
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                for ( int k = 1; k <= 3; k++ ) {
                    if ( eVals.at(k) < 0.0 ) {
                        eVals.at(k) = 0.0;
                    }
                }
                sqrtImD.at(i, j) = sqrt( eVals.at(1) ) * eVecs.at(i, 1) * eVecs.at(j, 1) + sqrt( eVals.at(2) ) * eVecs.at(i, 2) * eVecs.at(j, 2) + sqrt( eVals.at(3) ) * eVecs.at(i, 3) * eVecs.at(j, 3);
            }
        }

        AuxMatrix.beProductOf(effectiveStressTensor, sqrtImD);
        //stressTensor.beProductOf(sqrtImD,AuxMatrix);
        Part1.beProductOf(sqrtImD, AuxMatrix);

        // Second term of the equation 53 of the reference paper*****************************************************************************************
        /// @todo check if computeTraceD is necessary and check its implementation
        // Correct the trace of the damage tensor if necessary (see section 8.1 of the reference paper)
        traceTempD = computeTraceD(tempDamage, strainTensor, gp);
        double scalar = 0;
        AuxMatrix.zero();
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                scalar += ImD.at(i, j) * effectiveStressTensor.at(i, j);
            }
        }
        if ( ( 3. - traceTempD ) < 0.0001 ) {
            scalar = 0.0;
        } else {
            scalar = scalar / ( 3. - traceTempD );
        }
        //scalar /= 3.-traceTempD;
        AuxMatrix = ImD;
        AuxMatrix.times(scalar);
        Part2 = ImD;
        Part2.times(scalar);


        //stressTensor.subtract(AuxMatrix);
        // Third term of the equation 53 of the reference paper********************************************************************************************
        AuxMatrix.zero();
        AuxMatrix.at(1, 1) = AuxMatrix.at(2, 2) = AuxMatrix.at(3, 3) = 1. / 3.;
        if ( effectiveStressTrace > 0 ) {
            AuxMatrix.times( ( 1 - traceTempD ) * effectiveStressTrace );
        } else {
            AuxMatrix.times(effectiveStressTrace);
        }
        Part3 = AuxMatrix;
        //stressTensor.add(AuxMatrix);
        stressTensor = Part1;
        stressTensor.subtract(Part2);
        stressTensor.add(Part3);
        double factor;
        factor = computeCorrectionFactor(tempDamage, strainTensor, gp);
        stressTensor.times(factor);
        stressVector.beSymVectorForm(stressTensor);
        StructuralMaterial :: giveReducedSymVectorForm(answer, stressVector, mode);
    }

    // update gp
    this->correctBigValues(stressTensor);
    //    int checker20 = this->checkSymmetry(tempDamage);
    //    int checker21 = this->checkSymmetry(stressTensor);
    //    int checker22 = this->checkSymmetry(strainTensor);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempDamage(tempDamage);
    status->setTempKappa(tempKappa);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
}

void
AnisotropicDamageMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }


    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains, fullstrain;

        StructuralMaterial :: giveFullSymVectorForm( fullstrain, strain, gp->giveMaterialMode() );

        // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            //                fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
            fullstrain.at(3) = status->giveTempStrainZ();
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            fullstrain.at(2) = -nu *fullstrain.at(1);
            fullstrain.at(3) = -nu *fullstrain.at(1);
        }

        this->computePrincipalValues(principalStrains, fullstrain, principal_strain);

        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        kappa = sqrt(posNorm);
    } else {
        OOFEM_ERROR("computeEquivalentStrain: unknown EquivStrainType");
    }
    /*
     *  if ( this->equivStrainType == EST_Mazars ) {
     *      double posNorm = 0.0;
     *      FloatArray principalStrains, fullstrain;
     *
     *      StructuralMaterial :: giveFullSymVectorForm( fullstrain, strain, gp->giveMaterialMode() );
     *
     *      // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
     *      if ( gp->giveMaterialMode() == _PlaneStress ) {
     *          fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
     *      } else if ( gp->giveMaterialMode() == _1dMat ) {
     *          fullstrain.at(2) = -nu *fullstrain.at(1);
     *          fullstrain.at(3) = -nu *fullstrain.at(1);
     *      }
     *
     *      this->computePrincipalValues(principalStrains, fullstrain, principal_strain);
     *
     *      for ( int i = 1; i <= 3; i++ ) {
     *          if ( principalStrains.at(i) > 0.0 ) {
     *              posNorm += principalStrains.at(i) * principalStrains.at(i);
     *          }
     *      }
     *
     *      kappa = sqrt(posNorm);
     *  } else if ( ( this->equivStrainType == EST_Rankine_Smooth ) || ( this->equivStrainType == EST_Rankine_Standard ) ) {
     *      // EST_Rankine equiv strain measure
     *      double sum = 0.;
     *      FloatArray stress, fullStress, principalStress;
     *      FloatMatrix de;
     *
     *      lmat->giveStiffnessMatrix(de, SecantStiffness, gp, atTime);
     *      stress.beProductOf(de, strain);
     *      StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
     *      this->computePrincipalValues(principalStress, fullStress, principal_stress);
     *      for ( int i = 1; i <= 3; i++ ) {
     *          if ( principalStress.at(i) > 0.0 ) {
     *              if ( this->equivStrainType == EST_Rankine_Smooth ) {
     *                  sum += principalStress.at(i) * principalStress.at(i);
     *              } else if ( sum < principalStress.at(i) ) {
     *                  sum = principalStress.at(i);
     *              }
     *          } else if ( sum < principalStress.at(i) ) {
     *              sum = principalStress.at(i);
     *          }
     *      }
     *
     *      if ( this->equivStrainType == EST_Rankine_Smooth ) {
     *          sum = sqrt(sum);
     *      }
     *
     *      kappa = sum / lmat->give('E', gp);
     *  } else if ( ( this->equivStrainType == EST_ElasticEnergy ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) || ( this->equivStrainType == EST_ElasticEnergyPositiveStrain ) ) {
     *      // equivalent strain expressions based on elastic energy
     *      FloatMatrix de;
     *      FloatArray stress;
     *      double sum;
     *
     *      lmat->giveStiffnessMatrix(de, SecantStiffness, gp, atTime);
     *      if ( this->equivStrainType == EST_ElasticEnergy ) {
     *          // standard elastic energy
     *          stress.beProductOf(de, strain);
     *          sum = strain.dotProduct(stress);
     *      } else if ( this->equivStrainType == EST_ElasticEnergyPositiveStress ) {
     *          // elastic energy corresponding to positive part of stress
     *          FloatArray fullStress, principalStress;
     *          StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
     *          this->computePrincipalValues(principalStress, fullStress, principal_stress);
     *          // TO BE FINISHED
     *          sum = 0.;
     *          OOFEM_ERROR("Elastic energy corresponding to positive part of stress not finished");
     *      } else {
     *          // elastic energy corresponding to positive part of strain
     *          // TO BE DONE
     *          sum = 0.;
     *          OOFEM_ERROR("Elastic energy corresponding to positive part of strain not finished");
     *      }
     *
     *      kappa = sqrt( sum / lmat->give('E', gp) );
     *  } else if ( this->equivStrainType == EST_Griffith ) {
     *      double sum = 0.;
     *      FloatArray stress, fullStress, principalStress;
     *      FloatMatrix de;
     *
     *      lmat->giveStiffnessMatrix(de, SecantStiffness, gp, atTime);
     *      stress.beProductOf(de, strain);
     *      StructuralMaterial :: giveFullSymVectorForm( fullStress, stress, gp->giveMaterialMode() );
     *      this->computePrincipalValues(principalStress, fullStress, principal_stress);
     *      for ( int i = 1; i <= 3; i++ ) {
     *          if ( principalStress.at(i) > 0.0 && sum < principalStress.at(i) ) {
     *              sum = principalStress.at(i);
     *          }
     *      }
     *
     *      //Use Griffith criterion if Rankine not applied
     *      if (sum == 0.){
     *          sum = -pow(principalStress.at(1)-principalStress.at(3),2.)/8./(principalStress.at(1)+principalStress.at(3));
     *      }
     *      sum = max(sum,0.);
     *      kappa = sum / lmat->give('E', gp);
     *  } else {
     *      OOFEM_ERROR("computeEquivalentStrain: unknown EquivStrainType");
     *  }
     */
}

// Computes Kappa according to the first damage law proposed in reference paper.
double
AnisotropicDamageMaterial :: computeKappa(FloatMatrix damageTensor)
{
    double trace = damageTensor.giveTrace();
    double answer = ( this->kappaf - this->kappa0 ) * trace + this->kappa0;
    return answer;
}


// Computes the percentage of deltaD that needs to be added so the first eigenvalue of (tempDamageTensor + alpha*deltaD) reaches Dc
// To do this, the middle point algorithm is used.
// @TODO: this algorithm is not particularly efficient and another algorithm could be implemented.
double
AnisotropicDamageMaterial :: obtainAlpha1(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, double damageThreshold)
{
    double alpha_a, alpha_b, newAlpha, eps, maxDamage, size;
    FloatMatrix deltaD, positiveStrainTensorSquared, resultingDamageTensor, eVecs;
    FloatArray eVals;
    int cont;
    cont = 1;
    alpha_a = 0;
    alpha_b = 1;
    newAlpha = ( alpha_a + alpha_b ) / 2;
    positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
    deltaD = positiveStrainTensorSquared;
    deltaD.times(newAlpha * deltaLambda);
    this->correctBigValues(deltaD);
    resultingDamageTensor = tempDamageTensor;
    resultingDamageTensor.add(deltaD);
    //int checker2 = this->checkSymmetry(resultingDamageTensor);
    resultingDamageTensor.jaco_(eVals, eVecs, 20);
    size = eVals.giveSize();
    maxDamage = eVals.at(1);
    for ( int i = 2; i <= size; i++ ) {
        if ( eVals.at(i) > maxDamage ) {
            maxDamage = eVals.at(i);
        }
    }
    eps = maxDamage - damageThreshold;
    do {
        if ( eps > 0.0 ) {
            alpha_b = newAlpha;
            newAlpha = ( alpha_a + alpha_b ) / 2;
            deltaD = positiveStrainTensorSquared;
            deltaD.times(newAlpha * deltaLambda);
            this->correctBigValues(deltaD);
            resultingDamageTensor = tempDamageTensor;
            resultingDamageTensor.add(deltaD);
            //int checker3 = this->checkSymmetry(resultingDamageTensor);
            resultingDamageTensor.jaco_(eVals, eVecs, 20);
            size = eVals.giveSize();
            maxDamage = eVals.at(1);
            for ( int i = 2; i <= size; i++ ) {
                if ( eVals.at(i) > maxDamage ) {
                    maxDamage = eVals.at(i);
                }
            }
            eps = maxDamage - damageThreshold;
            cont = cont + 1;
            if ( cont == 100 ) {
                return newAlpha;
            }
        } else {
            alpha_a = newAlpha;
            newAlpha = ( alpha_a + alpha_b ) / 2;
            deltaD = positiveStrainTensorSquared;
            deltaD.times(newAlpha * deltaLambda);
            this->correctBigValues(deltaD);
            resultingDamageTensor = tempDamageTensor;
            resultingDamageTensor.add(deltaD);
            //int checker4 = this->checkSymmetry(resultingDamageTensor);
            resultingDamageTensor.jaco_(eVals, eVecs, 20);
            size = eVals.giveSize();
            maxDamage = eVals.at(1);
            for ( int i = 2; i <= size; i++ ) {
                if ( eVals.at(i) > maxDamage ) {
                    maxDamage = eVals.at(i);
                }
            }
            eps = maxDamage - damageThreshold;
            cont = cont + 1;
            if ( cont == 100 ) {
                return newAlpha;
            }
        }
    } while ( fabs(eps) > 1.0e-15 );
    return newAlpha;
}

// Computes the percentage of deltaD that needs to be added so the second eigenvalue of (tempDamageTensor + alpha*deltaD) reaches Dc
// To do this, the middle point algorithm is used.
// @TODO: this algorithm is not particularly efficient and another algorithm could be implemented.
double
AnisotropicDamageMaterial :: obtainAlpha2(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, FloatMatrix projPosStrainTensor, double damageThreshold)
{
    double alpha_a, alpha_b, newAlpha, eps, maxDamage, size, minVal;
    FloatMatrix deltaD, resultingDamageTensor, eVecs;
    FloatArray eVals;
    int cont;
	maxDamage = 0.0;
    cont = 1;
    alpha_a = 0;
    alpha_b = 1;
    newAlpha = ( alpha_a + alpha_b ) / 2;
    deltaD = projPosStrainTensor;
    deltaD.times(newAlpha * deltaLambda);
    this->correctBigValues(deltaD);
    resultingDamageTensor = tempDamageTensor;
    resultingDamageTensor.add(deltaD);
    //int checker5 = this->checkSymmetry(resultingDamageTensor);
    resultingDamageTensor.jaco_(eVals, eVecs, 40);
    size = eVals.giveSize();
    minVal = eVals.at(1);
    for ( int i = 2; i <= size; i++ ) {
        if ( eVals.at(i) < minVal ) {
            minVal = eVals.at(i);
            maxDamage = ( ( eVals.at(1) + eVals.at(2) + eVals.at(3) ) - eVals.at(i) ) / 2;
        }
    }
    eps = maxDamage - damageThreshold;
    do {
        if ( eps > 0.0 ) {
            alpha_b = newAlpha;
            newAlpha = ( alpha_a + alpha_b ) / 2;
            deltaD = projPosStrainTensor;
            deltaD.times(newAlpha * deltaLambda);
            this->correctBigValues(deltaD);
            resultingDamageTensor = tempDamageTensor;
            resultingDamageTensor.add(deltaD);
            //int checker6 = this->checkSymmetry(resultingDamageTensor);
            resultingDamageTensor.jaco_(eVals, eVecs, 40);
            size = eVals.giveSize();
            minVal = eVals.at(1);
            for ( int i = 2; i <= size; i++ ) {
                if ( eVals.at(i) < minVal ) {
                    minVal = eVals.at(i);
                    maxDamage = ( ( eVals.at(1) + eVals.at(2) + eVals.at(3) ) - eVals.at(i) ) / 2;
                }
            }
            eps = maxDamage - damageThreshold;
            cont = cont + 1;
            if ( cont == 100 ) {
                return newAlpha;
            }
        } else {
            alpha_a = newAlpha;
            newAlpha = ( alpha_a + alpha_b ) / 2;
            deltaD = projPosStrainTensor;
            deltaD.times(newAlpha * deltaLambda);
            this->correctBigValues(deltaD);
            resultingDamageTensor = tempDamageTensor;
            resultingDamageTensor.add(deltaD);
            //int checker7 = this->checkSymmetry(resultingDamageTensor);
            resultingDamageTensor.jaco_(eVals, eVecs, 40);
            size = eVals.giveSize();
            minVal = eVals.at(1);
            for ( int i = 2; i <= size; i++ ) {
                if ( eVals.at(i) < minVal ) {
                    minVal = eVals.at(i);
                    maxDamage = ( ( eVals.at(1) + eVals.at(2) + eVals.at(3) ) - eVals.at(i) ) / 2;
                }
            }
            eps = maxDamage - damageThreshold;
            cont = cont + 1;
            if ( cont == 100 ) {
                return newAlpha;
            }
        }
    } while ( fabs(eps) > 1.0e-15 );
    return newAlpha;
}

// Computes the percentage of deltaD that needs to be added so the third eigenvalue of (tempDamageTensor + alpha*deltaD) reaches Dc
// To do this, the middle point algorithm is used.
// @TODO: this algorithm is not particularly efficient and another algorithm could be implemented.
double
AnisotropicDamageMaterial :: obtainAlpha3(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, FloatArray vec3, double damageThreshold)
{
    double alpha_a, alpha_b, newAlpha, eps, aux = 0;
    FloatMatrix deltaD, positiveStrainTensorSquared, resultingDamageTensor, eVecs;
    FloatArray eVals;
    int cont;
    cont = 1;
    alpha_a = 0;
    alpha_b = 1;
    newAlpha = ( alpha_a + alpha_b ) / 2;
    positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
    for ( int i = 1; i <= 3; i++ ) {
        aux = aux + vec3.at(i) * ( positiveStrainTensorSquared.at(i, 1) * vec3.at(1) + positiveStrainTensorSquared.at(i, 2) * vec3.at(2) + positiveStrainTensorSquared.at(i, 3) * vec3.at(3) );
    }
    deltaD.beDyadicProductOf(vec3, vec3);
    deltaD.times(newAlpha * deltaLambda * aux);
    this->correctBigValues(deltaD);
    resultingDamageTensor = tempDamageTensor;
    resultingDamageTensor.add(deltaD);
    //int checker8 = this->checkSymmetry(resultingDamageTensor);
    resultingDamageTensor.jaco_(eVals, eVecs, 20);
    eps = eVals.at(1) + eVals.at(2) + eVals.at(3) - 3 * damageThreshold;
    do {
        if ( eps > 0.0 ) {
            alpha_b = newAlpha;
            newAlpha = ( alpha_a + alpha_b ) / 2;
            deltaD.beDyadicProductOf(vec3, vec3);
            deltaD.times(newAlpha * deltaLambda * aux);
            this->correctBigValues(deltaD);
            resultingDamageTensor = tempDamageTensor;
            resultingDamageTensor.add(deltaD);
            //int checker9 = this->checkSymmetry(resultingDamageTensor);
            resultingDamageTensor.jaco_(eVals, eVecs, 20);
            eps = eVals.at(1) + eVals.at(2) + eVals.at(3) - 3 * damageThreshold;
            cont = cont + 1;
            if ( cont == 100 ) {
                return newAlpha;
            }
        } else {
            alpha_a = newAlpha;
            newAlpha = ( alpha_a + alpha_b ) / 2;
            deltaD.beDyadicProductOf(vec3, vec3);
            deltaD.times(newAlpha * deltaLambda * aux);
            this->correctBigValues(deltaD);
            resultingDamageTensor = tempDamageTensor;
            resultingDamageTensor.add(deltaD);
            //int checker10 = this->checkSymmetry(resultingDamageTensor);
            resultingDamageTensor.jaco_(eVals, eVecs, 20);
            eps = eVals.at(1) + eVals.at(2) + eVals.at(3) - 3 * damageThreshold;
            cont = cont + 1;
            if ( cont == 100 ) {
                return newAlpha;
            }
        }
    } while ( fabs(eps) > 1.0e-15 );
    return newAlpha;
}
//To check symmetry: delete this function when everything works fine
double
AnisotropicDamageMaterial :: checkSymmetry(FloatMatrix matrix)
{
    int a = 0;
    int nRows = matrix.giveNumberOfRows();
    for ( int i = 1; i <= nRows; i++ ) {
        for ( int j = 1; j <= nRows; j++ ) {
            if ( fabs( matrix.at(i, j) - matrix.at(j, i) ) < 1.e-6 ) {
                ;
            } else {
                a = 1;
            }
        }
    }
    if ( a == 1 ) {
        a = 1;
    }
    return a;
}

void
AnisotropicDamageMaterial :: correctBigValues(FloatMatrix &matrix)
{
    int nRows = matrix.giveNumberOfRows();
    for ( int i = 1; i <= nRows; i++ ) {
        for ( int j = 1; j <= nRows; j++ ) {
            if ( matrix.at(i, j) != matrix.at(j, i) ) {
                double Aux = ( matrix.at(i, j) + matrix.at(j, i) ) / 2.0;
                matrix.at(i, j) = Aux;
                matrix.at(j, i) = Aux;
            }
        }
    }
}

double
AnisotropicDamageMaterial :: computeTraceD(FloatMatrix tempDamageTensor, FloatMatrix strainTensor, GaussPoint *gp)
{
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    int flag = status->giveFlag();
    //int tempFlag=status->giveTempFlag();
    double Dc = 1.00, trD = 0;
    // If flag = 0, the trace of the damage tensor has never been greater than 1 before
    if ( flag == 0 ) {
        if ( ( strainTensor.at(1, 1) + strainTensor.at(2, 2) + strainTensor.at(3, 3) ) < 0 ) { // Compression
            trD = tempDamageTensor.at(1, 1) + tempDamageTensor.at(2, 2) + tempDamageTensor.at(3, 3);
            if ( trD >= 1 ) {
                status->setTempFlag(1);
            }                                                   // The trace of the damage tensor is greater than 1 for the first time, then, flag turns into 1
        } else {                                                                                                                                                // Tension
            if ( ( tempDamageTensor.at(1, 1) + tempDamageTensor.at(2, 2) + tempDamageTensor.at(3, 3) ) >= 1 ) {
                trD = Dc;
            } else {
                trD = tempDamageTensor.at(1, 1) + tempDamageTensor.at(2, 2) + tempDamageTensor.at(3, 3);
            }
        }
    }
    // If flag = 1, the trace of the damage tensor has become greater than 1 before
    if ( flag == 1 ) {
        if ( ( strainTensor.at(1, 1) + strainTensor.at(2, 2) + strainTensor.at(3, 3) ) < 0 ) { // Compression
            trD = tempDamageTensor.at(1, 1) + tempDamageTensor.at(2, 2) + tempDamageTensor.at(3, 3);
        } else {
            trD = Dc;
        }                                                                                                                                                   // Tension
    }
    return trD;
}

double
AnisotropicDamageMaterial :: computeCorrectionFactor(FloatMatrix tempDamageTensor, FloatMatrix strainTensor, GaussPoint *gp)
{
    // In the case that the material has experimented some damaged under compression, this must affect the material behaviour when it is
    // under tension in the future
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    int tempFlag = status->giveFlag();
    double tempStoredFactor = status->giveStoredFactor();
    double Dc = 1.00, trD = 0;
    double factor = 0.0;
    trD = tempDamageTensor.at(1, 1) + tempDamageTensor.at(2, 2) + tempDamageTensor.at(3, 3);
    // If flag = 0, the trace of the damage tensor has never been greater than Dc under compression before
    if ( tempFlag == 0 ) {
        if ( ( strainTensor.at(1, 1) + strainTensor.at(2, 2) + strainTensor.at(3, 3) ) < 0 ) { // Compression
            factor = 1.0;
            if ( ( 1.0 - trD ) < tempStoredFactor ) {
                status->setTempStoredFactor(1.0 - trD);
            }
            if ( ( 1.0 - trD ) <= ( 1.0 - Dc ) ) {
                tempFlag = 1;
                status->setTempFlag(1);                 // The trace of the damage tensor is greater than Dc for the first time, then, flag turns into 1
                status->setTempStoredFactor(1. - Dc);
            }
        } else {                                                                                                                                                // Tension
            factor = tempStoredFactor;
        }
    }

    // If flag = 1, the trace of the damage tensor has become greater than 1 before
    if ( tempFlag == 1 ) {
        if ( ( strainTensor.at(1, 1) + strainTensor.at(2, 2) + strainTensor.at(3, 3) ) < 0 ) { // Compression
            factor = 1.0;
        } else {
            //{factor=status->giveStoredFactor();}    // Tension
            factor = 1. - Dc;
        }
    }
    return factor;
}

void
AnisotropicDamageMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
//
// Implementation of the 3D stiffness matrix, according to the equations 56 and 57 of the reference paper.
{
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, atTime);
    } else {
        const FloatArray &totalStrain = status->giveTempStrainVector();
        FloatArray reducedTotalStrainVector;
        this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);
        //FloatArray totalStrain;
        FloatMatrix damageTensor, strainTensor;
        // The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
        // must be divided by 2
        strainTensor.resize(3, 3);
        strainTensor.zero();
        strainTensor.at(1, 1) = reducedTotalStrainVector.at(1);
        strainTensor.at(2, 2) = reducedTotalStrainVector.at(2);
        strainTensor.at(3, 3) = reducedTotalStrainVector.at(3);
        strainTensor.at(2, 3) = reducedTotalStrainVector.at(4) / 2.0;
        strainTensor.at(3, 2) = reducedTotalStrainVector.at(4) / 2.0;
        strainTensor.at(1, 3) = reducedTotalStrainVector.at(5) / 2.0;
        strainTensor.at(3, 1) = reducedTotalStrainVector.at(5) / 2.0;
        strainTensor.at(1, 2) = reducedTotalStrainVector.at(6) / 2.0;
        strainTensor.at(2, 1) = reducedTotalStrainVector.at(6) / 2.0;
        // The damage tensor is read
        damageTensor = status->giveTempDamage();
        AnisotropicDamageMaterial :: computeSecantOperator(answer, strainTensor, damageTensor, gp);
        for ( int j = 4; j <= 6; j++ ) {
            for ( int i = 1; i <= 6; i++ ) {
                answer.at(i, j) = answer.at(i, j) / 2.0;
            }
        }
    }
}


void AnisotropicDamageMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *atTime)
{
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, atTime);
    } else {
        FloatArray totalStrain = status->giveTempStrainVector();
        FloatArray reducedTotalStrainVector;
        this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);
        FloatMatrix damageTensor, strainTensor;
        // The damage tensor is read
        damageTensor.resize(3, 3);
        damageTensor.zero();
        damageTensor = status->giveTempDamage();
        this->computePlaneStressStrain(strainTensor, damageTensor, reducedTotalStrainVector, gp, atTime);
        // The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
        // must be divided by 2
        FloatMatrix secantOperator;
        secantOperator.resize(6, 6);
        AnisotropicDamageMaterial :: computeSecantOperator(secantOperator, strainTensor, damageTensor, gp);
        double C11, C12, C13, C16, C21, C22, C23, C26, C61, C62, C63, C66, q, r, s;
        C11 = secantOperator.at(1, 1);
        C12 = secantOperator.at(1, 2);
        C13 = secantOperator.at(1, 3);
        C16 = secantOperator.at(1, 6);
        C21 = secantOperator.at(2, 1);
        C22 = secantOperator.at(2, 2);
        C23 = secantOperator.at(2, 3);
        C26 = secantOperator.at(2, 6);
        C61 = secantOperator.at(6, 1);
        C62 = secantOperator.at(6, 2);
        C63 = secantOperator.at(6, 3);
        C66 = secantOperator.at(6, 6);
        //Change 14-march-2014
        FloatArray stressVector;
        FloatMatrix Dn;
        stressVector = status->giveTempStressVector();
        Dn = status->giveTempDamage();
        this->computePlaneStressStrain(strainTensor, Dn, reducedTotalStrainVector, gp, atTime);
        if ( ( stressVector.at(1) + stressVector.at(2) ) < 1.0e-5 ) {
            q = -nu / E;
        } else {
            q = strainTensor.at(3, 3) / ( stressVector.at(1) + stressVector.at(2) );
        }
        //q = strainTensor.at(3,3)/(stressVector.at(1)+stressVector.at(2));
        //q = -nu/E;
        r = 1. / ( 1. - C13 * q );
        s = 1. / ( 1. - q * C23 - C23 * q * q * r * C13 );
        answer.resize(3, 3);
        answer.at(2, 1) = s * ( C21 + C11 * C23 * q * r );
        answer.at(2, 2) = s * ( C22 + C12 * C23 * q * r );
        answer.at(2, 3) = s * ( C26 + C16 * C23 * q * r ) * 1. / 2.;
        answer.at(1, 1) = r * ( C11 + C13 * q * answer.at(2, 1) );
        answer.at(1, 2) = r * ( C12 + C13 * q * answer.at(2, 2) );
        answer.at(1, 3) = r * ( C16 + C13 * q * answer.at(2, 3) ) * 1 / 2;
        answer.at(3, 1) = C61 + C63 * q * ( answer.at(1, 1) + answer.at(2, 1) );
        answer.at(3, 2) = C62 + C63 * q * ( answer.at(1, 2) + answer.at(2, 2) );
        answer.at(3, 3) = ( C66 + C63 * q * ( answer.at(1, 3) + answer.at(2, 3) ) ) * 1. / 2.;
    }
}

void AnisotropicDamageMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *atTime)
{}

void AnisotropicDamageMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                        GaussPoint *gp, TimeStep *atTime)
{
    // Implementation of the 3D stiffness matrix
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( mode == ElasticStiffness ) {
        this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, atTime);
    } else {
        FloatArray strain = status->giveTempStrainVector();
        if ( ( strain.at(1) + strain.at(2) + strain.at(3) )  > 0 ) {
            //@todo eq 56
        } else {
            //@todo eq 57
        }
    }
}

void
AnisotropicDamageMaterial :: computePlaneStressStrain(FloatMatrix &answer, FloatMatrix damageTensor, FloatArray reducedTotalStrainVector, GaussPoint *gp,
                                                      TimeStep *atTime)
//
{
    FloatMatrix inPlaneStrain, B, eVecs;
    FloatArray eVals;
    double outOfPlaneStrain;
    B.resize(2, 2);
    B.at(1, 1) = 1. - damageTensor.at(1, 1);
    B.at(1, 2) = 0. - damageTensor.at(1, 2);
    B.at(2, 1) = 0. - damageTensor.at(2, 1);
    B.at(2, 2) = 1. - damageTensor.at(2, 2);
    //The eigenVectors for the change of base must be computed using the damageTensor, NOT the B matrix!!!
    //B.jaco_(eVals, eVecs, 40);
    FloatMatrix Auxiliar;
    Auxiliar.resize(2, 2);
    Auxiliar.at(1, 1) = damageTensor.at(1, 1);
    Auxiliar.at(2, 2) = damageTensor.at(2, 2);
    Auxiliar.at(1, 2) = damageTensor.at(1, 2);
    Auxiliar.at(2, 1) = damageTensor.at(2, 1);

    Auxiliar.jaco_(eVals, eVecs, 40);

    // Change of base of reducedTotalStrainVector from the canonical to the base formed by the eigenvectors of the damageTensor
    inPlaneStrain.resize(2, 2);
    inPlaneStrain.at(1, 1) = reducedTotalStrainVector.at(1);
    inPlaneStrain.at(2, 2) = reducedTotalStrainVector.at(2);
    inPlaneStrain.at(1, 2) = reducedTotalStrainVector.at(3) / 2.0;
    inPlaneStrain.at(2, 1) = reducedTotalStrainVector.at(3) / 2.0;

    double term1, term2, term3, B1, B2, Bz, trD, h, epsilon11, epsilon22;
    B1 = 1.0 - eVals.at(1);
    B2 = 1.0 - eVals.at(2);
    Bz = 1. - damageTensor.at(3, 3);
    FloatArray vector1, vector2, auxVector;
    vector1.resize(2);
    vector2.resize(2);
    vector1.at(1) = eVecs.at(1, 1);
    vector1.at(2) = eVecs.at(2, 1);
    vector2.at(1) = eVecs.at(1, 2);
    vector2.at(2) = eVecs.at(2, 2);
    auxVector.beProductOf(inPlaneStrain, vector1);
    epsilon11 = vector1.at(1) * auxVector.at(1) + vector1.at(2) * auxVector.at(2);
    auxVector.beProductOf(inPlaneStrain, vector2);
    epsilon22 = vector2.at(1) * auxVector.at(1) + vector2.at(2) * auxVector.at(2);
    trD = damageTensor.at(1, 1) + damageTensor.at(2, 2) + damageTensor.at(3, 3);
    // Assuming a Tension state --> h = 1.0
    h = 1.0;
    term1 = 3. * Bz * B1 * ( 1. - 2. * nu ) - ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu );
    term2 = 3. * Bz * B2 * ( 1. - 2. * nu ) - ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu );
    term3 = 3. * Bz * ( 1. - 2. * nu ) * ( B1 + B2 ) + ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu );
    if ( Bz < 0.00001 ) {
        outOfPlaneStrain = -( epsilon11 + epsilon22 );
    } else {
        outOfPlaneStrain = ( term1 * epsilon11 + term2 * epsilon22 ) / term3;
    }
    /*    if (outOfPlaneStrain != outOfPlaneStrain){
     *      outOfPlaneStrain=-(epsilon11+epsilon22);
     *  }*/
    double trStrain;
    trStrain = inPlaneStrain.at(1, 1) + inPlaneStrain.at(2, 2) + outOfPlaneStrain;
    // Check if actually under Tension, if not, recalculate term1, term2 and term3 with h = 0.0
    if ( trStrain < 0.0 ) {
        h = 0.0;
        term1 = 3. * Bz * B1 * ( 1. - 2. * nu ) - ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu );
        term2 = 3. * Bz * B2 * ( 1. - 2. * nu ) - ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu );
        term3 = 3. * Bz * ( 1. - 2. * nu ) * ( B1 + B2 ) + ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu );
        if ( Bz < 0.00001 ) {
            outOfPlaneStrain = -( epsilon11 + epsilon22 );
        } else {
            outOfPlaneStrain = ( term1 * epsilon11 + term2 * epsilon22 ) / term3;
        }
    }
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = inPlaneStrain.at(1, 1);
    answer.at(1, 2) = inPlaneStrain.at(1, 2);
    answer.at(2, 1) = inPlaneStrain.at(2, 1);
    answer.at(2, 2) = inPlaneStrain.at(2, 2);
    answer.at(3, 3) = outOfPlaneStrain;
}

void
AnisotropicDamageMaterial :: computePlaneStressSigmaZ(double &answer, FloatMatrix damageTensor, FloatArray reducedTotalStrainVector,
                                                      double epsilonZ, GaussPoint *gp, TimeStep *atTime)
//
{
    FloatMatrix Auxiliar, inPlaneStrain;
    FloatArray eVals;
    FloatMatrix eVecs;
    Auxiliar.resize(2, 2);
    Auxiliar.at(1, 1) = damageTensor.at(1, 1);
    Auxiliar.at(2, 2) = damageTensor.at(2, 2);
    Auxiliar.at(1, 2) = damageTensor.at(1, 2);
    Auxiliar.at(2, 1) = damageTensor.at(2, 1);
    Auxiliar.jaco_(eVals, eVecs, 40);
    inPlaneStrain.resize(2, 2);
    inPlaneStrain.at(1, 1) = reducedTotalStrainVector.at(1);
    inPlaneStrain.at(2, 2) = reducedTotalStrainVector.at(2);
    inPlaneStrain.at(1, 2) = reducedTotalStrainVector.at(3) / 2.0;
    inPlaneStrain.at(2, 1) = reducedTotalStrainVector.at(3) / 2.0;
    double term1, term2, termZ, B1, B2, Bz, trD, h, epsilon11, epsilon22;
    B1 = 1.0 - eVals.at(1);
    B2 = 1.0 - eVals.at(2);
    Bz = 1. - damageTensor.at(3, 3);
    FloatArray vector1, vector2, auxVector;
    vector1.resize(2);
    vector2.resize(2);
    vector1.at(1) = eVecs.at(1, 1);
    vector1.at(2) = eVecs.at(2, 1);
    vector2.at(1) = eVecs.at(1, 2);
    vector2.at(2) = eVecs.at(2, 2);
    auxVector.beProductOf(inPlaneStrain, vector1);
    epsilon11 = vector1.at(1) * auxVector.at(1) + vector1.at(2) * auxVector.at(2);
    auxVector.beProductOf(inPlaneStrain, vector2);
    epsilon22 = vector2.at(1) * auxVector.at(1) + vector2.at(2) * auxVector.at(2);
    trD = damageTensor.giveTrace();
    double Estar;
    Estar = E / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
    if ( ( epsilon11 + epsilon22 + epsilonZ ) >= 0. ) {
        h = 1.;
    } else {
        h = 0.;
    }
    term1 = ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu ) - 3. * Bz * B1 * ( 1. - 2. * nu );
    term2 = ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu ) - 3. * Bz * B2 * ( 1. - 2. * nu );
    termZ = ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu ) + 3. * Bz * ( 1. - 2. * nu ) * ( B1 + B2 );
    // Finally the expression of sigmaZ is composed
    answer = ( Estar / ( 3. * ( B1 + B2 + Bz ) ) ) * ( epsilon11 * term1 + epsilon22 * term2 + epsilonZ * termZ );

#if 0
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double B1, B2, Bz, eps11, eps22, Estar, term1, term2, termZ, trD, h;
    Estar=E / ((1. + nu)*(1. - 2. * nu));
    // Compute the eigenvalues of the in-plane damage tensor
    double eVal1, eVal2, aux1, aux2;
    aux1 = (damageTensor.at(1,1) + damageTensor.at(2,2))/2.0;
    aux2 = sqrt(pow((damageTensor.at(1,1) - damageTensor.at(2,2)) / 2. , 2.) + damageTensor.at(1,2) * damageTensor.at(2,1));
    eVal1 = aux1 + aux2 ;
    eVal2 = aux1 - aux2 ;
    B1 = 1. - eVal1;
    B2 = 1. - eVal2;
    Bz = 1. - damageTensor.at(3, 3);
    eps11 = reducedTotalStrainVector.at(1);
    eps22 = reducedTotalStrainVector.at(2);
    if ((eps11 + eps22 + epsZ)>=0.) {
        h = 1.;
    } else {
        h = 0.;
    }
    trD = damageTensor.giveTrace();
    term1 = ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu ) - 3. * Bz * B1 * ( 1. - 2. * nu ) ;
    term2 = ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu ) - 3. * Bz * B2 * ( 1. - 2. * nu ) ;
    termZ = ( 3. - trD ) * ( 1. - h * trD ) * ( 1. + nu ) + 3. * Bz * ( 1. - 2. * nu ) * ( B1 + B2 ) ;
    // Finally the expression of sigmaZ is composed
    answer = (Estar / (3. * (B1 + B2 + Bz))) * (eps11 * term1 + eps22 * term2 + epsZ * termZ);
    Estar = Estar + 0.0;
#endif
}

void
AnisotropicDamageMaterial :: computeDamageTensor(FloatMatrix &answer, GaussPoint *gp,
                                                 const FloatArray &reducedTotalStrainVector, double equivStrain,
                                                 TimeStep *atTime)
//
{
    //
    // returns real stress vector in 3d stress space of receiver according to
    // previous level of stress and current strain increment, the only way,
    // how to correctly update gp records
    //
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double Dc = 1.00;
    double Kappa;
    FloatMatrix tempDamageTensor;
    //this->computeEquivalentStrain(equivStrain, reducedTotalStrainVector, gp, atTime);
    FloatMatrix Dn = status->giveDamage();
    Kappa = this->computeKappa(Dn);
    //damageTensor = status->giveDamage();

    if ( equivStrain <= Kappa ) {                       // damage does not grow. Elastic behaviour
        answer.resize(3, 3);
        answer.zero();
        answer = Dn;
    } else {                                                            // damage grows
        Kappa = equivStrain;
        double deltaLambda;
        FloatArray eVals, fullStrainVector;
        FloatMatrix eVecs, strainTensor, positiveStrainTensor, positiveStrainTensorSquared, tempDamageTensor0;
        MaterialMode mode = gp->giveMaterialMode();
        // Compute square of positive part of strain tensor
        //1.- converts strain vector to full form;
        StructuralMaterial :: giveFullSymVectorForm(fullStrainVector, reducedTotalStrainVector, mode);
        // The strain vector is turned into a tensor; for that, the elements that are out of the diagonal
        // must be divided by 2
        strainTensor.resize(3, 3);
        strainTensor.zero();
        if ( mode == _PlaneStress ) {
            this->computePlaneStressStrain(strainTensor, Dn, reducedTotalStrainVector, gp, atTime);
            strainTensor.at(3, 3) = status->giveTempStrainZ();
        } else {
            strainTensor.at(1, 1) = reducedTotalStrainVector.at(1);
            strainTensor.at(2, 2) = reducedTotalStrainVector.at(2);
            strainTensor.at(3, 3) = reducedTotalStrainVector.at(3);
            strainTensor.at(2, 3) = reducedTotalStrainVector.at(4) / 2.0;
            strainTensor.at(3, 2) = reducedTotalStrainVector.at(4) / 2.0;
            strainTensor.at(1, 3) = reducedTotalStrainVector.at(5) / 2.0;
            strainTensor.at(3, 1) = reducedTotalStrainVector.at(5) / 2.0;
            strainTensor.at(1, 2) = reducedTotalStrainVector.at(6) / 2.0;
            strainTensor.at(2, 1) = reducedTotalStrainVector.at(6) / 2.0;
        }
        // computes polar decomposition and negative eigenvalues set to zero
        //int checker14 = this->checkSymmetry(strainTensor);
        strainTensor.jaco_(eVals, eVecs, 40);
        for ( int i = 1; i <= 3; i++ ) {
            if ( eVals.at(i) < 0 ) {
                eVals.at(i) = 0;
            }
        }
        // computes the positive part of the strain tensor
        positiveStrainTensor.resize(3, 3);
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                positiveStrainTensor.at(i, j) = eVals.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eVals.at(2) * eVecs.at(i, 2) * eVecs.at(j, 2) +  eVals.at(3) * eVecs.at(i, 3) * eVecs.at(j, 3);
            }
        }
        // computes the square of positiveStrainTensor
        positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
        //compute delta Lambda
        //double traceD = damageTensor.at(1,1) + damageTensor.at(2,2) + damageTensor.at(3,3);
        double traceD = Dn.at(1, 1) + Dn.at(2, 2) + Dn.at(3, 3);

        double traceTempD = this->computeTraceD(equivStrain);
        // equation 50 of the reference paper
        deltaLambda =  ( traceTempD - traceD ) / equivStrain / equivStrain;
        //compute delta D: equation 48 of the reference paper
        FloatMatrix deltaD;
        deltaD = positiveStrainTensorSquared;
        deltaD.times(deltaLambda);
        this->correctBigValues(deltaD);
        // compute new damage tensor
        tempDamageTensor0 = Dn;
        tempDamageTensor0.add(deltaD);                          //damage tensor is equal to tempDamageTensor+deltaD
        // The following loop implements the algorithm for a case in which the maximum damage threshold is reached,
        // in such a case, the remaining damage is projected on the other two directions available and, finally,
        // if the damage threshold is reached in two of the three possible directions, the remaining damage is projected
        // on the remaining third direction. If the threshold is reached in the three possible directions, the damage
        // tensor remains unchanged in the future and with all their eigenvalues equal to the damage threshold Dc.
        // This part of the code is based on the section 8.2 of the reference paper.
        //int checker15 = this->checkSymmetry(tempDamageTensor0);
        tempDamageTensor0.jaco_(eVals, eVecs, 20);
        if ( ( eVals.at(1) > ( Dc * 1.001 ) ) || ( eVals.at(2) > ( Dc * 1.001 ) ) || ( eVals.at(3) > ( Dc * 1.001 ) ) ) {
            double alpha = 0, deltaLambda1 = 0, Aux1 = 0, Aux2 = 0, Aux3 = 0;
            FloatMatrix deltaD1(3, 3), positiveStrainTensorSquared(3, 3), tempDamageTensor1(3, 3), deltaD2(3, 3), N11(3, 3), N12(3, 3), N13(3, 3), N12sym(3, 3), N13sym(3, 3), projPosStrainTensor(3, 3);
            FloatArray auxVals(3), auxVec1(3), auxVec2(3), auxVec3(3);

            // the percentage alpha of deltaD that needs to be added so the first eigenvalue reaches Dc is obtained
            alpha = obtainAlpha1(Dn, deltaLambda, positiveStrainTensor, Dc);

            // deltaD1 is obtained --> deltaD1=alpha*deltaD
            positiveStrainTensorSquared.beProductOf(positiveStrainTensor, positiveStrainTensor);
            deltaD1 = positiveStrainTensorSquared;
            deltaD1.times(alpha * deltaLambda);
            this->correctBigValues(deltaD1);
            tempDamageTensor1 = Dn;
            tempDamageTensor1.add(deltaD1);
            // The following lines describe the process to apply the equation 64 of the reference paper. First, the
            // eigenvalues and eigenvectors of the damage tensor resulting at the moment when the threshold is reached
            // are obtained
            // (note: the equation 64 is not correctly written in the paper, it should be as implemented here:
            // D_dot = lambda_dot * [ <e>^2-(nI·<e>^2 nI)(nI x nI) - 2(nII·<e>^2 nI)(nI x nII)_sym - 2(nIII·<e>^2 nI)(nI x nIII)_sym ]
            //int checker16 = this->checkSymmetry(tempDamageTensor1);

            tempDamageTensor1.jaco_(eVals, eVecs, 40);
            // The eigenvalues and eigenvectors are ordered, with the maximum eigenvalue being I, as its corresponding
            // eigenvector, and the other two being II and III. This is necessary so the equation 64 of the reference
            // paper can be applied
            if ( eVals.at(1) >= eVals.at(2) &&  eVals.at(1) >= eVals.at(3) ) {
                auxVals.at(1) = eVals.at(1);
                auxVec1.at(1) = eVecs.at(1, 1);
                auxVec1.at(2) = eVecs.at(2, 1);
                auxVec1.at(3) = eVecs.at(3, 1);
                auxVals.at(2) = eVals.at(2);
                auxVec2.at(1) = eVecs.at(1, 2);
                auxVec2.at(2) = eVecs.at(2, 2);
                auxVec2.at(3) = eVecs.at(3, 2);
                auxVals.at(3) = eVals.at(3);
                auxVec3.at(1) = eVecs.at(1, 3);
                auxVec3.at(2) = eVecs.at(2, 3);
                auxVec3.at(3) = eVecs.at(3, 3);
            } else if ( eVals.at(2) >= eVals.at(1) &&  eVals.at(2) >= eVals.at(3) ) {
                auxVals.at(1) = eVals.at(2);
                auxVec1.at(1) = eVecs.at(1, 2);
                auxVec1.at(2) = eVecs.at(2, 2);
                auxVec1.at(3) = eVecs.at(3, 2);
                auxVals.at(2) = eVals.at(1);
                auxVec2.at(1) = eVecs.at(1, 1);
                auxVec2.at(2) = eVecs.at(2, 1);
                auxVec2.at(3) = eVecs.at(3, 1);
                auxVals.at(3) = eVals.at(3);
                auxVec3.at(1) = eVecs.at(1, 3);
                auxVec3.at(2) = eVecs.at(2, 3);
                auxVec3.at(3) = eVecs.at(3, 3);
            } else {
                auxVals.at(1) = eVals.at(3);
                auxVec1.at(1) = eVecs.at(1, 3);
                auxVec1.at(2) = eVecs.at(2, 3);
                auxVec1.at(3) = eVecs.at(3, 3);
                auxVals.at(2) = eVals.at(2);
                auxVec2.at(1) = eVecs.at(1, 2);
                auxVec2.at(2) = eVecs.at(2, 2);
                auxVec2.at(3) = eVecs.at(3, 2);
                auxVals.at(3) = eVals.at(1);
                auxVec3.at(1) = eVecs.at(1, 1);
                auxVec3.at(2) = eVecs.at(2, 1);
                auxVec3.at(3) = eVecs.at(3, 1);
            }

            // The symmetric part of the dyadic product of eigenvectors n1 and n2 is obtained
            N11.beDyadicProductOf(auxVec1, auxVec1);
            N12.beDyadicProductOf(auxVec1, auxVec2);
            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    N12sym.at(i, j) = 0.5 * ( N12.at(i, j) + N12.at(j, i) );
                }
            }

            // The symmetric part of the dyadic product of eigenvectors n1 and n3 is obtained
            N13.beDyadicProductOf(auxVec1, auxVec3);
            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    N13sym.at(i, j) = 0.5 * ( N13.at(i, j) + N13.at(j, i) );
                }
            }

            //The projected positive strain tensor is obtained
            for ( int i = 1; i <= 3; i++ ) {
                Aux1 = Aux1 + auxVec1.at(i) * ( positiveStrainTensorSquared.at(i, 1) * auxVec1.at(1) + positiveStrainTensorSquared.at(i, 2) * auxVec1.at(2) + positiveStrainTensorSquared.at(i, 3) * auxVec1.at(3) );                    //==(n1*(<e>^2*n1))     (eq. 64 of the reference paper)
                Aux2 = Aux2 + auxVec2.at(i) * ( positiveStrainTensorSquared.at(i, 1) * auxVec1.at(1) + positiveStrainTensorSquared.at(i, 2) * auxVec1.at(2) + positiveStrainTensorSquared.at(i, 3) * auxVec1.at(3) );                    //==(n2*(<e>^2*n1))_sym (eq. 64 of the reference paper)
                Aux3 = Aux3 + auxVec3.at(i) * ( positiveStrainTensorSquared.at(i, 1) * auxVec1.at(1) + positiveStrainTensorSquared.at(i, 2) * auxVec1.at(2) + positiveStrainTensorSquared.at(i, 3) * auxVec1.at(3) );                    //==(n3*(<e>^2*n1))_sym (eq. 64 of the reference paper)
            }
            N11.times(Aux1);
            N12sym.times(2 * Aux2);
            N13sym.times(2 * Aux3);

            // Finally, the expression between brackets in equation 64 is built up and called projPosStrainTensor
            projPosStrainTensor = positiveStrainTensorSquared;
            projPosStrainTensor.subtract(N11);
            projPosStrainTensor.subtract(N12sym);
            projPosStrainTensor.subtract(N13sym);

            // The following loop avoids numerical problems in the case that the trace of projPosStrainTensor is very small
            if ( ( projPosStrainTensor.at(1, 1) + projPosStrainTensor.at(2, 2) + projPosStrainTensor.at(3, 3) ) < traceTempD * 1e-10 ) {
                deltaLambda1 = 0.;
            } else {
                deltaLambda1 = ( traceTempD - ( tempDamageTensor1.at(1, 1) + tempDamageTensor1.at(2, 2) + tempDamageTensor1.at(3, 3) ) ) / ( projPosStrainTensor.at(1, 1) + projPosStrainTensor.at(2, 2) + projPosStrainTensor.at(3, 3) );
            }
            //projPosStrainTensor.symmetrized();
            deltaD2 = projPosStrainTensor;
            deltaD2.times(deltaLambda1);
            this->correctBigValues(deltaD2);
            tempDamageTensor1.add(deltaD2);                             //damage tensor is equal to tempDamageTensor+deltaD1+deltaD2
            // The following loop checks if after the addition of D2, any other eigenvalue of the damage tensor
            // has reached the threshold. If it has, it repeats the process, but this time projecting the
            // remaining damage on the direction of the remaining eigenvector
            //int checker17 = this->checkSymmetry(tempDamageTensor1);
            tempDamageTensor1.jaco_(eVals, eVecs, 40);
            if ( ( eVals.at(1) > ( Dc * 1.001 ) ) || ( eVals.at(2) > ( Dc * 1.001 ) ) || ( eVals.at(3) > ( Dc * 1.001 ) ) ) {
                FloatMatrix deltaD3(3, 3), projPosStrainTensor_new(3, 3), tempDamageTensor2(3, 3), deltaD4(3, 3);
                FloatArray vec3(3);
                double alpha2 = 0, deltaLambda2 = 0, Aux4 = 0;
                // double val3=0;
                // Restoring the value of tempDamageTensor1 = tempDamageTensor + deltaD1
                tempDamageTensor1 = Dn;
                tempDamageTensor1.add(deltaD1);
                // the percentage alpha2 of deltaD2 that needs to be added to (tempDamageTensor+deltaD1) so the second eigenvalue
                // reaches Dc is obtained
                alpha2 = obtainAlpha2(tempDamageTensor1, deltaLambda1, positiveStrainTensor, projPosStrainTensor, Dc);
                deltaD3 = deltaD2;
                deltaD3.times(alpha2);
                tempDamageTensor2 = tempDamageTensor1;
                tempDamageTensor2.add(deltaD3);
                // The smallest eigenvalue is detected and its eigenvector is used to build the new projPosStrainTensor
                //int checker18 = this->checkSymmetry(tempDamageTensor2);
                tempDamageTensor2.jaco_(eVals, eVecs, 40);
                if ( eVals.at(1) <= eVals.at(2) && eVals.at(1) <= eVals.at(3) ) {
                    //val3=eVals.at(1);
                    vec3.at(1) = eVecs.at(1, 1);
                    vec3.at(2) = eVecs.at(2, 1);
                    vec3.at(3) = eVecs.at(3, 1);
                } else if ( eVals.at(2) <= eVals.at(1) && eVals.at(2) <= eVals.at(3) ) {
                    //val3=eVals.at(2);
                    vec3.at(1) = eVecs.at(1, 2);
                    vec3.at(2) = eVecs.at(2, 2);
                    vec3.at(3) = eVecs.at(3, 2);
                } else {
                    //val3=eVals.at(3);
                    vec3.at(1) = eVecs.at(1, 3);
                    vec3.at(2) = eVecs.at(2, 3);
                    vec3.at(3) = eVecs.at(3, 3);
                }

                // The following loop computes nIII·<e>^2 nIII
                for ( int i = 1; i <= 3; i++ ) {
                    Aux4 = Aux4 + vec3.at(i) * ( positiveStrainTensorSquared.at(i, 1) * vec3.at(1) + positiveStrainTensorSquared.at(i, 2) * vec3.at(2) + positiveStrainTensorSquared.at(i, 3) * vec3.at(3) );
                }
                projPosStrainTensor_new.beDyadicProductOf(vec3, vec3);
                //
                projPosStrainTensor_new.times(Aux4);

                // The following loop avoids numerical problems in the case that the trace of projPosStrainTensor is very small
                if ( ( projPosStrainTensor_new.at(1, 1) + projPosStrainTensor_new.at(2, 2) + projPosStrainTensor_new.at(3, 3) ) < traceTempD * 1e-10 ) {
                    deltaLambda2 = 0;
                } else {
                    deltaLambda2 = ( traceTempD - ( tempDamageTensor2.at(1, 1) + tempDamageTensor2.at(2, 2) + tempDamageTensor2.at(3, 3) ) ) / ( projPosStrainTensor_new.at(1, 1) + projPosStrainTensor_new.at(2, 2) + projPosStrainTensor_new.at(3, 3) );
                }
                deltaD4 = projPosStrainTensor_new;
                deltaD4.times(deltaLambda2);
                tempDamageTensor2.add(deltaD4);                                 //damage tensor is equal to tempDamageTensor+deltaD1+deltaD3+deltaD4

                // The following loop checks if after the addition of D4, the remaining eigenvalue of the damage tensor
                // has reached the threshold. If it has, it computes a damage tensor with all its eigenvalues equal
                // to the damage threshold Dc
                //int checker19 = this->checkSymmetry(tempDamageTensor2);
                tempDamageTensor2.jaco_(eVals, eVecs, 40);
                if ( ( eVals.at(1) > ( Dc * 1.001 ) ) || ( eVals.at(2) > ( Dc * 1.001 ) ) || ( eVals.at(3) > ( Dc * 1.001 ) ) ) {
                    double alpha3 = 0;
                    FloatMatrix deltaD5, tempDamageTensor3;
                    tempDamageTensor2 = Dn;
                    tempDamageTensor2.add(deltaD1);
                    tempDamageTensor2.add(deltaD3);
                    alpha3 = obtainAlpha3(tempDamageTensor2, deltaLambda2, positiveStrainTensor, vec3, Dc);
                    deltaD5 = deltaD4;
                    deltaD5.times(alpha3);
                    tempDamageTensor3 = tempDamageTensor2;
                    tempDamageTensor3.add(deltaD5);
                    tempDamageTensor = tempDamageTensor3;
                    tempDamageTensor3.jaco_(eVals, eVecs, 40);
                    if ( ( eVals.at(1) > ( Dc * 1.001 ) ) || ( eVals.at(2) > ( Dc * 1.001 ) ) || ( eVals.at(3) > ( Dc * 1.001 ) ) ) {
                        tempDamageTensor3.zero();
                        for ( int i = 1; i <= 3; i++ ) {
                            if ( eVals.at(i) > Dc * 1.001 ) {
                                eVals.at(i) = Dc;
                            }
                        }
                        for ( int i = 1; i <= 3; i++ ) {
                            for ( int j = 1; j <= 3; j++ ) {
                                tempDamageTensor3.at(i, j) = eVals.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1)       +       eVals.at(2) * eVecs.at(i, 2) * eVecs.at(j, 2) +       eVals.at(3) * eVecs.at(i, 3) * eVecs.at(j, 3);
                            }
                        }

                        /*
                         *                                        tempDamageTensor3.zero();
                         *                                        tempDamageTensor3.at(1,1)=Dc;
                         *                                        tempDamageTensor3.at(2,2)=Dc;
                         *                                        tempDamageTensor3.at(3,3)=Dc;*/
                    }
                    tempDamageTensor = tempDamageTensor3;
                } else {
                    tempDamageTensor = tempDamageTensor2;
                }
            } else {
                tempDamageTensor = tempDamageTensor1;
            }
        } else {
            tempDamageTensor = tempDamageTensor0;
        }
        answer = tempDamageTensor;
    }
}

void
AnisotropicDamageMaterial :: computeSecantOperator(FloatMatrix &answer, FloatMatrix strainTensor, FloatMatrix damageTensor, GaussPoint *gp)
//
// Implementation of the 3D stiffness matrix, according to the equations 56 and 57 of the reference paper.
{
    //    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    double G, K;
    double traceD, Aux;
    FloatMatrix ImD, sqrtImD, eVecs, Imatrix;
    FloatArray eVals;
    //MaterialMode mode = gp->giveMaterialMode();
    G = E / ( 2.0 * ( 1.0 + nu ) );
    K = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
    //Compute the trace of the damage tensor, correcting it if necessary (see section 8.1 of the reference paper)
    traceD = computeTraceD(damageTensor, strainTensor, gp);

    if ( fabs(3. - traceD) < 0.001 ) {
        Aux = 0.0;
    } else {
        Aux = ( 1. / ( 3. - traceD ) );
    }

    // compute square root of (1-D)
    ImD.resize(3, 3);
    ImD.zero();
    ImD.at(1, 1) = ImD.at(2, 2) = ImD.at(3, 3) = 1.0;
    ImD.subtract(damageTensor);

    // computes square of positive part of strain tensor
    //int checker11=this->checkSymmetry(ImD);
    ImD.jaco_(eVals, eVecs, 40);
    sqrtImD.resize(3, 3);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            for ( int k = 1; k <= 3; k++ ) {
                if ( eVals.at(k) < 0.0 ) {
                    eVals.at(k) = 0.0;
                }
            }
            sqrtImD.at(i, j) = sqrt( eVals.at(1) ) * eVecs.at(i, 1) * eVecs.at(j, 1) + sqrt( eVals.at(2) ) * eVecs.at(i, 2) * eVecs.at(j, 2) + sqrt( eVals.at(3) ) * eVecs.at(i, 3) * eVecs.at(j, 3);
        }
    }
    // To compute the expresions 56 and 57 of the reference paper, we need to work with fourth order tensors. To do this,
    // a structured called fourthOrderTensor is defined. This structure is composed by nine 3x3 FloatMatrix objects
    struct fourthOrderTensor
    {
        FloatMatrix Matrix_11kl;
        FloatMatrix Matrix_12kl;
        FloatMatrix Matrix_13kl;
        FloatMatrix Matrix_21kl;
        FloatMatrix Matrix_22kl;
        FloatMatrix Matrix_23kl;
        FloatMatrix Matrix_31kl;
        FloatMatrix Matrix_32kl;
        FloatMatrix Matrix_33kl;
    };
    // Four fourthOrderTensor structures are defined
    fourthOrderTensor Block1, Block2, Block3, secantOperator;
    Imatrix.resize(3, 3);
    Imatrix.zero();
    Imatrix.at(1, 1) = Imatrix.at(2, 2) = Imatrix.at(3, 3) = 1.0;

    // The fourthOrderTensor structures are initialised
    Block1.Matrix_11kl.resize(3, 3);
    Block1.Matrix_12kl.resize(3, 3);
    Block1.Matrix_13kl.resize(3, 3);
    Block1.Matrix_21kl.resize(3, 3);
    Block1.Matrix_22kl.resize(3, 3);
    Block1.Matrix_23kl.resize(3, 3);
    Block1.Matrix_31kl.resize(3, 3);
    Block1.Matrix_32kl.resize(3, 3);
    Block1.Matrix_33kl.resize(3, 3);
    Block2.Matrix_11kl.resize(3, 3);
    Block2.Matrix_12kl.resize(3, 3);
    Block2.Matrix_13kl.resize(3, 3);
    Block2.Matrix_21kl.resize(3, 3);
    Block2.Matrix_22kl.resize(3, 3);
    Block2.Matrix_23kl.resize(3, 3);
    Block2.Matrix_31kl.resize(3, 3);
    Block2.Matrix_32kl.resize(3, 3);
    Block2.Matrix_33kl.resize(3, 3);
    Block3.Matrix_11kl.resize(3, 3);
    Block3.Matrix_12kl.resize(3, 3);
    Block3.Matrix_13kl.resize(3, 3);
    Block3.Matrix_21kl.resize(3, 3);
    Block3.Matrix_22kl.resize(3, 3);
    Block3.Matrix_23kl.resize(3, 3);
    Block3.Matrix_31kl.resize(3, 3);
    Block3.Matrix_32kl.resize(3, 3);
    Block3.Matrix_33kl.resize(3, 3);
    secantOperator.Matrix_11kl.resize(3, 3);
    secantOperator.Matrix_12kl.resize(3, 3);
    secantOperator.Matrix_13kl.resize(3, 3);
    secantOperator.Matrix_21kl.resize(3, 3);
    secantOperator.Matrix_22kl.resize(3, 3);
    secantOperator.Matrix_23kl.resize(3, 3);
    secantOperator.Matrix_31kl.resize(3, 3);
    secantOperator.Matrix_32kl.resize(3, 3);
    secantOperator.Matrix_33kl.resize(3, 3);

    for ( int k = 1; k <= 3; k++ ) {
        for ( int l = 1; l <= 3; l++ ) {
            //The first block inside the brackets is obtained --> (1-D)^(1/2) _x_ (1-D)^(1/2)
            Block1.Matrix_11kl.at(k, l) = sqrtImD.at(1, k) * sqrtImD.at(1, l);
            Block1.Matrix_12kl.at(k, l) = sqrtImD.at(1, k) * sqrtImD.at(2, l);
            Block1.Matrix_13kl.at(k, l) = sqrtImD.at(1, k) * sqrtImD.at(3, l);
            Block1.Matrix_21kl.at(k, l) = sqrtImD.at(2, k) * sqrtImD.at(1, l);
            Block1.Matrix_22kl.at(k, l) = sqrtImD.at(2, k) * sqrtImD.at(2, l);
            Block1.Matrix_23kl.at(k, l) = sqrtImD.at(2, k) * sqrtImD.at(3, l);
            Block1.Matrix_31kl.at(k, l) = sqrtImD.at(3, k) * sqrtImD.at(1, l);
            Block1.Matrix_32kl.at(k, l) = sqrtImD.at(3, k) * sqrtImD.at(2, l);
            Block1.Matrix_33kl.at(k, l) = sqrtImD.at(3, k) * sqrtImD.at(3, l);
            //The second block inside the brackets is obtained --> ((1-D) x (1-D))/(3-trD)
            Block2.Matrix_11kl.at(k, l) = ImD.at(1, 1) * ImD.at(k, l) * Aux;
            Block2.Matrix_12kl.at(k, l) = ImD.at(1, 2) * ImD.at(k, l) * Aux;
            Block2.Matrix_13kl.at(k, l) = ImD.at(1, 3) * ImD.at(k, l) * Aux;
            Block2.Matrix_21kl.at(k, l) = ImD.at(2, 1) * ImD.at(k, l) * Aux;
            Block2.Matrix_22kl.at(k, l) = ImD.at(2, 2) * ImD.at(k, l) * Aux;
            Block2.Matrix_23kl.at(k, l) = ImD.at(2, 3) * ImD.at(k, l) * Aux;
            Block2.Matrix_31kl.at(k, l) = ImD.at(3, 1) * ImD.at(k, l) * Aux;
            Block2.Matrix_32kl.at(k, l) = ImD.at(3, 2) * ImD.at(k, l) * Aux;
            Block2.Matrix_33kl.at(k, l) = ImD.at(3, 3) * ImD.at(k, l) * Aux;
            //The crossed-product of two identity tensors is obtained --> 1 x 1
            Block3.Matrix_11kl.at(k, l) = Imatrix.at(1, 1) * Imatrix.at(k, l);
            Block3.Matrix_12kl.at(k, l) = Imatrix.at(1, 2) * Imatrix.at(k, l);
            Block3.Matrix_13kl.at(k, l) = Imatrix.at(1, 3) * Imatrix.at(k, l);
            Block3.Matrix_21kl.at(k, l) = Imatrix.at(2, 1) * Imatrix.at(k, l);
            Block3.Matrix_22kl.at(k, l) = Imatrix.at(2, 2) * Imatrix.at(k, l);
            Block3.Matrix_23kl.at(k, l) = Imatrix.at(2, 3) * Imatrix.at(k, l);
            Block3.Matrix_31kl.at(k, l) = Imatrix.at(3, 1) * Imatrix.at(k, l);
            Block3.Matrix_32kl.at(k, l) = Imatrix.at(3, 2) * Imatrix.at(k, l);
            Block3.Matrix_33kl.at(k, l) = Imatrix.at(3, 3) * Imatrix.at(k, l);
        }
    }
    // equation 56 of the reference paper
    if ( ( strainTensor.at(1, 1) + strainTensor.at(2, 2) + strainTensor.at(3, 3) )  > 0.0 ) {
        for ( int k = 1; k <= 3; k++ ) {
            for ( int l = 1; l <= 3; l++ ) {
                secantOperator.Matrix_11kl.at(k, l) = 2 * G * ( Block1.Matrix_11kl.at(k, l) - Block2.Matrix_11kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_11kl.at(k, l);
                secantOperator.Matrix_12kl.at(k, l) = 2 * G * ( Block1.Matrix_12kl.at(k, l) - Block2.Matrix_12kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_12kl.at(k, l);
                secantOperator.Matrix_13kl.at(k, l) = 2 * G * ( Block1.Matrix_13kl.at(k, l) - Block2.Matrix_13kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_13kl.at(k, l);
                secantOperator.Matrix_21kl.at(k, l) = 2 * G * ( Block1.Matrix_21kl.at(k, l) - Block2.Matrix_21kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_21kl.at(k, l);
                secantOperator.Matrix_22kl.at(k, l) = 2 * G * ( Block1.Matrix_22kl.at(k, l) - Block2.Matrix_22kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_22kl.at(k, l);
                secantOperator.Matrix_23kl.at(k, l) = 2 * G * ( Block1.Matrix_23kl.at(k, l) - Block2.Matrix_23kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_23kl.at(k, l);
                secantOperator.Matrix_31kl.at(k, l) = 2 * G * ( Block1.Matrix_31kl.at(k, l) - Block2.Matrix_31kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_31kl.at(k, l);
                secantOperator.Matrix_32kl.at(k, l) = 2 * G * ( Block1.Matrix_32kl.at(k, l) - Block2.Matrix_32kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_32kl.at(k, l);
                secantOperator.Matrix_33kl.at(k, l) = 2 * G * ( Block1.Matrix_33kl.at(k, l) - Block2.Matrix_33kl.at(k, l) ) + K * ( 1 - traceD ) * Block3.Matrix_33kl.at(k, l);
            }
        }
        // equation 57 of the reference paper
    } else {
        for ( int k = 1; k <= 3; k++ ) {
            for ( int l = 1; l <= 3; l++ ) {
                secantOperator.Matrix_11kl.at(k, l) = 2 * G * ( Block1.Matrix_11kl.at(k, l) - Block2.Matrix_11kl.at(k, l) ) + K *Block3.Matrix_11kl.at(k, l);
                secantOperator.Matrix_12kl.at(k, l) = 2 * G * ( Block1.Matrix_12kl.at(k, l) - Block2.Matrix_12kl.at(k, l) ) + K *Block3.Matrix_12kl.at(k, l);
                secantOperator.Matrix_13kl.at(k, l) = 2 * G * ( Block1.Matrix_13kl.at(k, l) - Block2.Matrix_13kl.at(k, l) ) + K *Block3.Matrix_13kl.at(k, l);
                secantOperator.Matrix_21kl.at(k, l) = 2 * G * ( Block1.Matrix_21kl.at(k, l) - Block2.Matrix_21kl.at(k, l) ) + K *Block3.Matrix_21kl.at(k, l);
                secantOperator.Matrix_22kl.at(k, l) = 2 * G * ( Block1.Matrix_22kl.at(k, l) - Block2.Matrix_22kl.at(k, l) ) + K *Block3.Matrix_22kl.at(k, l);
                secantOperator.Matrix_23kl.at(k, l) = 2 * G * ( Block1.Matrix_23kl.at(k, l) - Block2.Matrix_23kl.at(k, l) ) + K *Block3.Matrix_23kl.at(k, l);
                secantOperator.Matrix_31kl.at(k, l) = 2 * G * ( Block1.Matrix_31kl.at(k, l) - Block2.Matrix_31kl.at(k, l) ) + K *Block3.Matrix_31kl.at(k, l);
                secantOperator.Matrix_32kl.at(k, l) = 2 * G * ( Block1.Matrix_32kl.at(k, l) - Block2.Matrix_32kl.at(k, l) ) + K *Block3.Matrix_32kl.at(k, l);
                secantOperator.Matrix_33kl.at(k, l) = 2 * G * ( Block1.Matrix_33kl.at(k, l) - Block2.Matrix_33kl.at(k, l) ) + K *Block3.Matrix_33kl.at(k, l);
            }
        }
    }
    // The resulting material stiffness matrix is built
    answer.resize(6, 6);
    answer.at(1, 1) = secantOperator.Matrix_11kl.at(1, 1);
    answer.at(1, 2) = secantOperator.Matrix_11kl.at(2, 2);
    answer.at(1, 3) = secantOperator.Matrix_11kl.at(3, 3);
    answer.at(1, 4) = secantOperator.Matrix_11kl.at(2, 3);
    answer.at(1, 5) = secantOperator.Matrix_11kl.at(3, 1);
    answer.at(1, 6) = secantOperator.Matrix_11kl.at(1, 2);
    answer.at(2, 1) = secantOperator.Matrix_22kl.at(1, 1);
    answer.at(2, 2) = secantOperator.Matrix_22kl.at(2, 2);
    answer.at(2, 3) = secantOperator.Matrix_22kl.at(3, 3);
    answer.at(2, 4) = secantOperator.Matrix_22kl.at(2, 3);
    answer.at(2, 5) = secantOperator.Matrix_22kl.at(3, 1);
    answer.at(2, 6) = secantOperator.Matrix_22kl.at(1, 2);
    answer.at(3, 1) = secantOperator.Matrix_33kl.at(1, 1);
    answer.at(3, 2) = secantOperator.Matrix_33kl.at(2, 2);
    answer.at(3, 3) = secantOperator.Matrix_33kl.at(3, 3);
    answer.at(3, 4) = secantOperator.Matrix_33kl.at(2, 3);
    answer.at(3, 5) = secantOperator.Matrix_33kl.at(3, 1);
    answer.at(3, 6) = secantOperator.Matrix_33kl.at(1, 2);
    answer.at(4, 1) = secantOperator.Matrix_23kl.at(1, 1);
    answer.at(4, 2) = secantOperator.Matrix_23kl.at(2, 2);
    answer.at(4, 3) = secantOperator.Matrix_23kl.at(3, 3);
    answer.at(4, 4) = secantOperator.Matrix_23kl.at(2, 3);
    answer.at(4, 5) = secantOperator.Matrix_23kl.at(3, 1);
    answer.at(4, 6) = secantOperator.Matrix_23kl.at(1, 2);
    answer.at(5, 1) = secantOperator.Matrix_31kl.at(1, 1);
    answer.at(5, 2) = secantOperator.Matrix_31kl.at(2, 2);
    answer.at(5, 3) = secantOperator.Matrix_31kl.at(3, 3);
    answer.at(5, 4) = secantOperator.Matrix_31kl.at(2, 3);
    answer.at(5, 5) = secantOperator.Matrix_31kl.at(3, 1);
    answer.at(5, 6) = secantOperator.Matrix_31kl.at(1, 2);
    answer.at(6, 1) = secantOperator.Matrix_12kl.at(1, 1);
    answer.at(6, 2) = secantOperator.Matrix_12kl.at(2, 2);
    answer.at(6, 3) = secantOperator.Matrix_12kl.at(3, 3);
    answer.at(6, 4) = secantOperator.Matrix_12kl.at(2, 3);
    answer.at(6, 5) = secantOperator.Matrix_12kl.at(3, 1);
    answer.at(6, 6) = secantOperator.Matrix_12kl.at(1, 2);
}

int
AnisotropicDamageMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    AnisotropicDamageMaterialStatus *status = static_cast< AnisotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) { // returning the trace of the damage tensor
        answer.resize(1);
        answer.at(1) = status->giveDamage().at(1, 1) + status->giveDamage().at(2, 2) + status->giveDamage().at(3, 3);
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveDamage().at(1, 1);
        answer.at(2) = status->giveDamage().at(2, 2);
        answer.at(3) = status->giveDamage().at(3, 3);
        answer.at(4) = status->giveDamage().at(2, 3);
        answer.at(5) = status->giveDamage().at(1, 3);
        answer.at(6) = status->giveDamage().at(1, 2);
        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        //int checker12=this->checkSymmetry(status->giveDamage());
        FloatMatrix dam = status->giveDamage();
        FloatMatrix eVecs;
        dam.jaco_(answer, eVecs, 20);
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveTempDamage().at(1, 1);
        answer.at(2) = status->giveTempDamage().at(2, 2);
        answer.at(3) = status->giveTempDamage().at(3, 3);
        answer.at(4) = status->giveTempDamage().at(2, 3);
        answer.at(5) = status->giveTempDamage().at(1, 3);
        answer.at(6) = status->giveTempDamage().at(1, 2);
        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        //int checker13=this->checkSymmetry(status->giveTempDamage());
        FloatMatrix dam = status->giveTempDamage();
        FloatMatrix eVecs;
        dam.jaco_(answer, eVecs, 20);
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_StressWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }

    return 1; // to make the compiler happy
}


IRResultType
AnisotropicDamageMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    linearElasticMaterial->initializeFrom(ir);
    E = linearElasticMaterial->giveYoungsModulus();
    nu = linearElasticMaterial->givePoissonsRatio();
    int eqStrain = 0;
    // specify the type of formula for equivalent strain
    // currently only the Mazars formula is allowed !!! (this should be generalized later)
    //
    // IR_GIVE_OPTIONAL_FIELD(ir, eqStrain, _IFT_AnisotropicDamageMaterial_equivStrainType);
    //
    switch ( eqStrain ) {
    case 1: this->equivStrainType = EST_Rankine_Smooth;
        break;
    case 2: this->equivStrainType = EST_ElasticEnergy;
        break;
    case 3: this->equivStrainType = EST_Mises;
        // IR_GIVE_FIELD(ir, k, _IFT_IsotropicDamageMaterial1_k);
        break;
    case 4: this->equivStrainType = EST_Rankine_Standard;
        break;
    case 5: this->equivStrainType = EST_ElasticEnergyPositiveStress;
        break;
    case 6: this->equivStrainType = EST_ElasticEnergyPositiveStrain;
        break;
    case 7: this->equivStrainType = EST_Griffith;
        break;
    default: this->equivStrainType = EST_Mazars;
    }

    int damlaw = 0;
    // specify the type of damage law (which affects the shape of the stress-strain curve)
    IR_GIVE_OPTIONAL_FIELD(ir, damlaw, _IFT_AnisotropicDamageMaterial_damageLawType);
    switch ( damlaw ) {
    case 1: this->damageLawType = DLT_Desmorat1;
        break;
    case 2: this->damageLawType = DLT_Desmorat2;
        break;
    case 3: this->damageLawType = DLT_Linear;
        break;
    case 4: this->damageLawType = DLT_Exponential;
        break;
    default: this->damageLawType = DLT_Desmorat2;
    }

    IR_GIVE_FIELD(ir, kappa0, _IFT_AnisotropicDamageMaterial_kappa0);
    IR_GIVE_FIELD(ir, kappaf, _IFT_AnisotropicDamageMaterial_kappaf);
    if ( damageLawType == DLT_Desmorat2 ) {
        IR_GIVE_FIELD(ir, aA, _IFT_AnisotropicDamageMaterial_aA);
    }

    return StructuralMaterial :: initializeFrom(ir);
}

void
AnisotropicDamageMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->kappa0, _IFT_AnisotropicDamageMaterial_kappa0);
}


AnisotropicDamageMaterialStatus :: AnisotropicDamageMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    damage.resize(3, 3);
    damage.zero();
    tempDamage.resize(3, 3);
    tempDamage.zero();
    strainZ = tempStrainZ = 0.0;
    flag = tempFlag = 0;
    storedFactor = 1.0;
    tempStoredFactor = 1.0;

#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}

AnisotropicDamageMaterialStatus :: ~AnisotropicDamageMaterialStatus()
{ }

void
AnisotropicDamageMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _PlaneStress ) { // special treatment of the out-of-plane strain
        FloatArray helpVec;
        MaterialStatus :: printOutputAt(file, tStep);
        fprintf(file, "  strains ");
        StructuralMaterial :: giveFullSymVectorForm(helpVec, strainVector, mode);
        helpVec.at(3) = this->strainZ;
        for ( auto &v : helpVec ) {
            fprintf( file, " %.4e", v );
        }
        fprintf(file, "\n              stresses");
        StructuralMaterial :: giveFullSymVectorForm(helpVec, stressVector, mode);
        for ( auto &v : helpVec ) {
            fprintf( file, " %.4e", v );
        }
        fprintf(file, "\n");
    } else {
        StructuralMaterialStatus :: printOutputAt(file, tStep); // standard treatment of strains and stresses
    }

    fprintf(file, "status { ");
    fprintf(file, "kappa %g", this->kappa);
    double damtrace = tempDamage.giveTrace();
    if ( damtrace > 0.0 ) {
        fprintf(file, ", damage");
        int n = tempDamage.giveNumberOfRows();
        for ( int i = 1; i <= n; i++ ) {
            for ( int j = 1; j <= n; j++ ) {
                fprintf( file, " %g", tempDamage.at(i, j) );
            }
        }
    }

#ifdef keep_track_of_dissipated_energy
    fprintf(file, ", dissW %f, freeE %f, stressW %f ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif

    fprintf(file, "}\n");
}



void
AnisotropicDamageMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempDamage = this->damage;
    this->tempStrainZ = this->strainZ;
    this->tempStoredFactor = this->storedFactor;
#ifdef keep_track_of_dissipated_energy
    this->tempStressWork = this->stressWork;
    this->tempDissWork = this->dissWork;
#endif
}


void
AnisotropicDamageMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
    this->damage = this->tempDamage;
    this->strainZ = this->tempStrainZ;
    this->flag = this->tempFlag;
    this->storedFactor = this->tempStoredFactor;
#ifdef keep_track_of_dissipated_energy
    this->stressWork = this->tempStressWork;
    this->dissWork = this->tempDissWork;
#endif
}


contextIOResultType
AnisotropicDamageMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    /*    contextIOResultType iores;
     *
     *  // save parent class status
     *  if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
     *      THROW_CIOERR(iores);
     *  }
     *
     *  // write raw data
     *  if ( !stream.write(kappa) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     *  if ( !stream.write(damage) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     * #ifdef keep_track_of_dissipated_energy
     *  if ( !stream.write(stressWork) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     *  if ( !stream.write(dissWork) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     * #endif
     */
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // write damage (vector)
    if ( ( iores = damage.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
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
AnisotropicDamageMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    /*    contextIOResultType iores;
     *
     *  // read parent class status
     *  if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
     *      THROW_CIOERR(iores);
     *  }
     *
     *  // read raw data
     *  if ( !stream.read(kappa) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     *  if ( !stream.read(damage) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     * #ifdef keep_track_of_dissipated_energy
     *  if ( !stream.read(stressWork) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     *  if ( !stream.read(dissWork) ) {
     *      THROW_CIOERR(CIO_IOERR);
     *  }
     *
     * #endif
     */

    contextIOResultType iores;


    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // read damage (vector)
    if ( ( iores = damage.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream.read(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif


    return CIO_OK;
}

#ifdef keep_track_of_dissipated_energy
void
AnisotropicDamageMaterialStatus :: computeWork(GaussPoint *gp)
{
    /*
     * // strain increment
     * FloatArray deps;
     * deps.beDifferenceOf(tempStrainVector, strainVector);
     *
     * // increment of stress work density
     * double dSW = ( tempStressVector.dotProduct(deps) + stressVector.dotProduct(deps) ) / 2.;
     * tempStressWork = stressWork + dSW;
     *
     * // elastically stored energy density
     * double We = tempStressVector.dotProduct(tempStrainVector) / 2.;
     *
     * // dissipative work density
     * tempDissWork = tempStressWork - We;
     */
}
#endif
} // end namespace oofem
