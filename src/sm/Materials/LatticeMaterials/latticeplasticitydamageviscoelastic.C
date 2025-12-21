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

#include "latticeplasticitydamageviscoelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(LatticePlasticityDamageViscoelastic);

LatticePlasticityDamageViscoelastic::LatticePlasticityDamageViscoelastic(int n, Domain *d) : LatticePlasticityDamage(n, d)
{}

void
LatticePlasticityDamageViscoelastic::initializeFrom(InputRecord &ir)
{
    LatticePlasticityDamage::initializeFrom(ir);

    IR_GIVE_FIELD(ir, viscoMat, _IFT_LatticePlasticityDamageViscoelastic_viscoMat); // number of slave material

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    // override elastic modulus by the one given by the compliance function
    timeFactor = 0.;
    IR_GIVE_FIELD(ir, timeFactor, _IFT_LatticePlasticityDamageViscoelastic_timeFactor); // timeConversion equal to 1 day in current time units of the analysis

    double E28 = 1. / rheoMat->computeCreepFunction(timeFactor * 28.01, timeFactor * 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min

    this->eNormalMean = E28; // swap elastic modulus/stiffness

    if ( ir.hasField(_IFT_LatticePlasticityDamageViscoelastic_timedepfracturing) ) {
        this->fib = true;
        IR_GIVE_FIELD(ir, fib_fcm28, _IFT_LatticePlasticityDamageViscoelastic_fcm28);
        IR_GIVE_FIELD(ir, fib_s, _IFT_LatticePlasticityDamageViscoelastic_fib_s);
        IR_GIVE_FIELD(ir, stiffnessFactor, _IFT_LatticePlasticityDamageViscoelastic_stiffnessFactor);
    }
}


std::unique_ptr<MaterialStatus> 
LatticePlasticityDamageViscoelastic::CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<LatticePlasticityDamageViscoelasticStatus>(1, LatticePlasticityDamageViscoelastic::domain, gp);
}


FloatArrayF< 6 >
LatticePlasticityDamageViscoelastic::giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep)
{
    double tol = 1.e-12; // error in order of approx. Pascals

    auto status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();
    // just a dummy yet essential call to create a status of viscomaterial. Otherwise initTempStatus() would fail.
    rheoMat->giveStatus(rChGP);

    status->initTempStatus();

    FloatArrayF< 6 >reducedStrainForViscoMat;
    FloatArrayF< 6 >reducedStrain;
    FloatArrayF< 6 >quasiReducedStrain;

    FloatArray tempStressVE;
    FloatArrayF< 6 >viscoStress;
    FloatArrayF< 6 >plastDamStress;

    int itercount = 1;

    FloatArray indepStrain;
    indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);

    auto elasticStiffnessMatrix = LatticeLinearElastic::give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    double tolerance = 1.;

    FloatArrayF< 6 >inelasticTrialStrain;

    double sigmaResid;
    double sigmaResidPlus = this->eNormalMean;
    double sigmaResidMinus = -sigmaResidPlus;


    FloatArrayF< 6 >inelastPlus;
    FloatArrayF< 6 >inelastMinus;

    bool plusFlag = false;
    bool minusFlag = false;

    int iterBisection = 10;

    do {
        if ( itercount > 100 ) {
            OOFEM_WARNING("Algorithm not converging");

            printf("Algorithm not converging, giving up\n");
            printf("Unable to reach equilibrium between viscoelastic and CDPM materials, Element %d\n", gp->giveElement()->giveNumber() );
            printf("tolerance = %e > %e!, stress error = %e, itercount = %d\n", tolerance, tol, sigmaResid, itercount);

            if ( plusFlag ) {
                printf("inelastPlus  = ");
                for ( auto &val : inelastPlus ) {
                    printf(" %.10e", val);
                }
                printf(" sigmaResidPlus  = %e \n", sigmaResidPlus);
            }

            if ( minusFlag ) {
                printf("inelastMinus = ");
                for ( auto &val : inelastMinus ) {
                    printf(" %.10e", val);
                }
                printf(" sigmaResidMinus = %e \n", sigmaResidMinus);
            }

            break;
        }

        reducedStrainForViscoMat = totalStrain;

        if ( indepStrain.giveSize() > 0 ) {
            reducedStrainForViscoMat -= FloatArrayF< 6 >(indepStrain);
        }

        if ( itercount < iterBisection || ( !plusFlag || !minusFlag ) ) {
            inelasticTrialStrain = status->giveTempPlasticLatticeStrain() + FloatArrayF< 6 >( status->giveTempDamageLatticeStrain() );
        }

        reducedStrainForViscoMat -= inelasticTrialStrain;

        rheoMat->giveRealStressVector(tempStressVE, rChGP, reducedStrainForViscoMat, tStep);
        viscoStress = FloatArrayF< 6 >(tempStressVE);

        for ( int i = 1; i <= 6; i++ ) {
            quasiReducedStrain.at(i) =  viscoStress.at(i) / elasticStiffnessMatrix.at(i, i) + inelasticTrialStrain.at(i);
        }

        plastDamStress = this->performPlasticityReturn(gp, quasiReducedStrain, tStep);
        this->performDamageEvaluation(gp, quasiReducedStrain, tStep);
        double tempDamage = status->giveTempDamage();
        plastDamStress *= ( 1. - tempDamage );

        inelasticTrialStrain = status->giveTempPlasticLatticeStrain();
        inelasticTrialStrain += FloatArrayF< 6 >(status->giveTempDamageLatticeStrain() );

        tolerance = norm(plastDamStress - viscoStress) / this->eNormalMean;

        // try bisection method - prepare bounds, might be needed
        sigmaResid = viscoStress [ 0 ] - plastDamStress [ 0 ];

        if ( sigmaResid > 0. ) {
            if ( sigmaResid < sigmaResidPlus ) {
                plusFlag = true;
                sigmaResidPlus = sigmaResid;
                inelastPlus = inelasticTrialStrain;
            }
        } else {
            if ( sigmaResid > sigmaResidMinus ) {
                minusFlag = true;
                sigmaResidMinus = sigmaResid;
                inelastMinus = inelasticTrialStrain;
            }
        }

        if ( itercount >= iterBisection && plusFlag && minusFlag ) {
            inelasticTrialStrain = ( sigmaResidPlus * inelastMinus - sigmaResidMinus * inelastPlus ) / ( sigmaResidPlus - sigmaResidMinus );
        }

        itercount++;
    } while ( tolerance >= tol );

    status->letTempLatticeStrainBe(totalStrain);
    status->letTempLatticeStressBe(viscoStress);
    status->letTempReducedLatticeStrainBe(quasiReducedStrain);

    return viscoStress;
}

double
LatticePlasticityDamageViscoelastic::giveCompressiveStrength(GaussPoint *gp, TimeStep *tStep) const
{
    if ( fib ) {
        double equivalentTime = this->giveEquivalentTime(gp, tStep);
        double fcm = exp(this->fib_s * ( 1. - sqrt(28. * this->timeFactor / equivalentTime) ) ) * this->fib_fcm28;
        return fcm / this->fib_fcm28 * LatticePlasticityDamage::giveCompressiveStrength(gp, tStep);
    } else {
        return LatticePlasticityDamage::giveCompressiveStrength(gp, tStep);
    }
}


double
LatticePlasticityDamageViscoelastic::giveTensileStrength(GaussPoint *gp, TimeStep *tStep) const
{
    // check if the fracture properties are constant or time-dependent
    double ftm = 1.;
    double ftm28 = 1.;

    if ( fib ) {
        double equivalentTime = this->giveEquivalentTime(gp, tStep);
        double fcm = exp(this->fib_s * ( 1. - sqrt(28. * this->timeFactor / equivalentTime) ) ) * this->fib_fcm28;

        //Calculate the aged tensile strength
        if ( fcm >= 58. ) {
            ftm = 2.12 * log(1. + 0.1 * fcm) * 1.e6 / this->stiffnessFactor;
        } else if ( fcm <= 20. ) {
            ftm = 0.07862 * fcm * 1.e6 / this->stiffnessFactor; // 12^(2/3) * 0.3 / 20 = 0.07862
        } else {
            ftm = 0.3 * pow(fcm - 8., 2. / 3.) * 1.e6 / this->stiffnessFactor; //5.1-3a
        }

        //Do the same for the 28 day compressive strength
        if ( this->fib_fcm28 >= 58. ) {
            ftm28 = 2.12 * log(1. + 0.1 * this->fib_fcm28) * 1.e6 / this->stiffnessFactor;
        } else if ( this->fib_fcm28 <= 20. ) {
            ftm28 = 0.07862 * this->fib_fcm28 * 1.e6 / this->stiffnessFactor; // 12^(2/3) * 0.3 / 20 = 0.07862
        } else {
            ftm28 = 0.3 * pow(this->fib_fcm28 - 8., 2. / 3.) * 1.e6 / this->stiffnessFactor; //5.1-3a
        }
        return ftm / ftm28 * LatticePlasticityDamage::giveTensileStrength(gp, tStep);
    } else {
        return LatticePlasticityDamage::giveTensileStrength(gp, tStep);
    }
}


FloatMatrixF< 6, 6 >
LatticePlasticityDamageViscoelastic::give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();

    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    auto answer = LatticeLinearElastic::give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    answer *= ( Eincr / this->eNormalMean );

    if ( rmode == ElasticStiffness ) {
        return answer;
    } else if ( ( rmode == SecantStiffness ) || ( rmode == TangentStiffness ) ) {
        double omega = min(status->giveTempDamage(), 0.99999);
        return answer * ( 1. - omega );
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
    }
}

int
LatticePlasticityDamageViscoelastic::giveIPValue(FloatArray &answer,
                                                 GaussPoint *gp,
                                                 InternalStateType type,
                                                 TimeStep *tStep)
{
    if ( ( type == IST_DryingShrinkageTensor ) ||
         ( type == IST_AutogenousShrinkageTensor ) ||
         ( type == IST_TotalShrinkageTensor ) ||
         ( type == IST_CreepStrainTensor ) ||
         ( type == IST_DryingShrinkageTensor ) ||
         ( type == IST_ThermalStrainTensor ) ) {
        LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

        RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
        return rheoMat->giveIPValue(answer, status->giveSlaveGaussPointVisco(), type, tStep);
    }

    return LatticePlasticityDamage::giveIPValue(answer, gp, type, tStep);
}


int LatticePlasticityDamageViscoelastic::checkConsistency()
{
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    if ( rheoMat->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }

    if ( rheoMat->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }

    GaussPoint *noGP = NULL;
    if ( rheoMat->give(tAlpha, noGP) != 0. ) {
        OOFEM_ERROR("tAlpha must be set to 0. in slave viscoelastic material");
    }


    return FEMComponent::checkConsistency();
}

double
LatticePlasticityDamageViscoelastic::giveEquivalentTime(GaussPoint *gp, TimeStep *tStep) const
{
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

    return rheoMat->giveEquivalentTime(status->giveSlaveGaussPointVisco(), tStep);
}


LatticePlasticityDamageViscoelasticStatus::LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *gp) :
    LatticePlasticityDamageStatus(n, d, gp), slaveGpVisco(std::make_unique< GaussPoint >(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
{}

void
LatticePlasticityDamageViscoelasticStatus::initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
{
    LatticePlasticityDamageStatus::initTempStatus();

    // at this point rheomat :: giveStatus (viscoGP) has to be called first
    RheoChainMaterialStatus *rheoStatus = static_cast< RheoChainMaterialStatus * >( this->giveSlaveGaussPointVisco()->giveMaterialStatus() );

    rheoStatus->initTempStatus();
}

void
LatticePlasticityDamageViscoelasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticePlasticityDamageStatus::printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->printOutputAt(file, tStep);

    fprintf(file, "\n");
}


void
LatticePlasticityDamageViscoelasticStatus::updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->updateYourself(tStep);

    LatticePlasticityDamageStatus::updateYourself(tStep);
}

void
LatticePlasticityDamageViscoelasticStatus::saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    LatticePlasticityDamageStatus::saveContext(stream, mode);

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->saveContext(stream, mode);
}

void
LatticePlasticityDamageViscoelasticStatus::restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticePlasticityDamageStatus::restoreContext(stream, mode);

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
