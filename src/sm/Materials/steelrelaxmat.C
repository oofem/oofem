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

#include "steelrelaxmat.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(SteelRelaxMat);

SteelRelaxMat::SteelRelaxMat(int n, Domain *d) : StructuralMaterial(n, d)
{}


bool
SteelRelaxMat::hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _1dMat;
}

// reads the model parameters from the input file
void
SteelRelaxMat::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->E, _IFT_SteelRelaxMat_E); // Young's modulus

    int reinfClass;
    IR_GIVE_FIELD(ir, reinfClass, _IFT_SteelRelaxMat_reinfClass); // class of reinforcement according to MC2010 or EC2

    // Attention!!! rho1000 must be provided in % and not dimensionless
    if ( reinfClass == 1 ) { // class 1 (wires or cables with normal relaxation)
        this->k1 = 5.39;
        this->k2 = 6.7;
        this->rho1000 = 8.;
    } else if ( reinfClass == 2 ) { // class 2 (wires or cables with lowered relaxation)
        this->k1 = 0.66;
        this->k2 = 9.1;
        this->rho1000 = 2.5;
    } else if ( reinfClass == 3 ) { // class 3 (hot-rolled or modified rods)
        this->k1 = 1.98;
        this->k2 = 8.0;
        this->rho1000 = 4.;
    } else {
        OOFEM_ERROR("unsupported value of reinfClass");
    }

    // possibility to overwrite
    IR_GIVE_OPTIONAL_FIELD(ir, this->k1, _IFT_SteelRelaxMat_k1);
    IR_GIVE_OPTIONAL_FIELD(ir, this->k2, _IFT_SteelRelaxMat_k2);
    IR_GIVE_OPTIONAL_FIELD(ir, this->rho1000, _IFT_SteelRelaxMat_rho1000);

    IR_GIVE_FIELD(ir, this->timeFactor, _IFT_SteelRelaxMat_timeFactor);

    IR_GIVE_FIELD(ir, this->charStrength, _IFT_SteelRelaxMat_charStrength);

    this->relRelaxBound = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->relRelaxBound, _IFT_SteelRelaxMat_relRelaxBound);

    int approach;
    IR_GIVE_FIELD(ir, approach, _IFT_SteelRelaxMat_approach);

    this->Approach = ( approachType ) approach;

    // in incremental linear statics it seems that it is necessary to iterate equilibrium at the material point level (for the Bazant approach only) because otherwise the strain relaxation would be calculated from too high stress
    // tolerance = 1 Pa
    // MPa -> SF = 1e6 -> tol = 1
    // Pa -> SF = 1 -> tol = 1
    this->tolerance =  1. / 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->tolerance, _IFT_SteelRelaxMat_tolerance);
}

// creates a new material status  corresponding to this class
std::unique_ptr<MaterialStatus> 
SteelRelaxMat::CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<SteelRelaxMatStatus>(gp);
}

void
SteelRelaxMat::giveRealStressVector(FloatArray &answer,
                                    GaussPoint *gp,
                                    const FloatArray &totalStrain,
                                    TimeStep *tStep) const
{
    FloatArray reducedStrain, strainIncrement, stressVector;
    double stressIncrement;

    SteelRelaxMatStatus *status = static_cast< SteelRelaxMatStatus * >( this->giveStatus(gp) );

    if ( !this->isActivated(tStep) ) {
        stressVector.resize(1);
        stressVector.zero();

        status->letTempStrainVectorBe(totalStrain);
        status->letTempStressVectorBe(stressVector);

        answer.resize(1);
        answer.zero();
        return;
    }

    if ( this->Approach == EquivTime_EC2 ) {
        double lossIncrement;

        if ( status->giveStressVector().giveSize() ) {
            stressVector = status->giveStressVector();
        } else {
            stressVector.resize(1);
            stressVector.zero();
        }

        StructuralMaterial::giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Incremental);

        strainIncrement.beDifferenceOf(reducedStrain, status->giveStrainVector() );
        stressIncrement = strainIncrement.at(1) * this->E;

        stressVector.at(1) += stressIncrement;

        if ( stressVector.at(1) > 0. ) {
            this->computeIncrOfPrestressLossAtVarStrain(lossIncrement, gp, tStep, stressVector.at(1) );
            stressVector.at(1) -= lossIncrement;
        }

        status->letTempStrainVectorBe(totalStrain);
        status->letTempStressVectorBe(stressVector);
    } else if ( this->Approach == Bazant_EC2 ) {
        double prevIterTempStress;
        int i = 0;

        do {
            prevIterTempStress = status->giveTempStressVector().at(1);

            if ( status->giveStressVector().giveSize() ) {
                stressVector = status->giveStressVector();
            } else {
                stressVector.resize(1);
                stressVector.zero();
            }

            // get strain increment without strain due to cable relaxation and thermal strain
            this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Incremental);

            strainIncrement.beDifferenceOf(reducedStrain, status->giveStrainVector() );
            stressIncrement = strainIncrement.at(1) * this->E;

            stressVector.at(1) += stressIncrement;


            // store the final strain and stress in the status
            status->letTempStrainVectorBe(totalStrain);
            status->letTempStressVectorBe(stressVector);

            i++;

            if ( i > 1000 ) {
                OOFEM_ERROR("Algorithm not converging");
            }
        } while ( fabs(prevIterTempStress - status->giveTempStressVector().at(1) ) >= this->tolerance );


        if ( i > 50 ) {
            OOFEM_WARNING("Criterion of the algorithm reached in %d iterations, consider increasing tolerance", i);
        }
    }

    if ( stressVector.at(1) > this->charStrength ) {
        OOFEM_ERROR("Stress %f exeeds the characteristic strength of the material!", stressVector.at(1) );
    }

    answer.resize(1);
    answer.at(1) = stressVector.at(1);
}


FloatMatrixF< 1, 1 >
SteelRelaxMat::give1dStressStiffMtrx(MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *tStep) const
{
    if ( this->isActivated(tStep) ) {
        return { this->E };
    } else {
        return { 0. };
    }
}


void
SteelRelaxMat::giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep, ValueModeType mode) const
{
    FloatArray temperatureFreeStrain;
    FloatArray relaxationStrain;

    // subtracts temperature strain
    StructuralMaterial::giveStressDependentPartOfStrainVector(temperatureFreeStrain, gp, totalStrain, tStep, mode);

    // strains due to relaxation
    this->computeStressRelaxationStrainVector(relaxationStrain, gp, totalStrain, tStep, mode);

    answer = temperatureFreeStrain;
    answer.subtract(relaxationStrain);
}


void
SteelRelaxMat::evalStressRelaxationAtConstStrain(double &answer, GaussPoint *gp, double dt) const
{
    double rho;
    double k;
    double lambda;
    double mu;
    double prestress;

    SteelRelaxMatStatus *status = static_cast< SteelRelaxMatStatus * >( this->giveStatus(gp) );
    // prestress loss is not calculated in the time step when prestressing is applied
    prestress = status->givePrestress();

    answer = 0.;

    if ( prestress >  this->relRelaxBound * this->charStrength ) {
        mu = prestress / this->charStrength;

        rho = this->k1 * this->rho1000 * exp(this->k2 * mu) * 1.e-5;
        k = 0.75 * ( 1. - mu );
        lambda = 1000. * this->timeFactor / 24.;

        answer = prestress * rho * pow( ( dt / lambda ), k);
    }
}

// this method is active only in the case of EC2 approach
void
SteelRelaxMat::computeIncrOfPrestressLossAtVarStrain(double &answer, GaussPoint *gp, TimeStep *tStep, double stress) const
{
    double rho;
    double k;
    double lambda;
    double mu;
    double prestress;
    double lossesUpTillNow;
    double loss;

    SteelRelaxMatStatus *status = static_cast< SteelRelaxMatStatus * >( this->giveStatus(gp) );
    lossesUpTillNow =  status->giveRelaxIntVariable();

    // "initial" value of prestress is the sum of the current stress plus the cumulative subsequent relaxation
    prestress = stress + lossesUpTillNow;
    status->setTempPrestress(prestress);

    answer = 0.;

    if ( prestress > this->relRelaxBound * this->charStrength ) {
        mu = prestress / this->charStrength;
        rho = this->k1 * this->rho1000 * exp(this->k2 * mu) * 1.e-5;
        k = 0.75 * ( 1. - mu );
        lambda = 1000. * this->timeFactor / 24.;

        // compute total loss for updated prestress and time equiv
        double t_equiv;
        t_equiv = pow( ( lossesUpTillNow / ( prestress * rho ) ), ( 1. / k ) ) * lambda;
        this->evalStressRelaxationAtConstStrain(loss, gp, t_equiv + tStep->giveTimeIncrement() );

        // set temporary sum of losses
        status->setTempRelaxIntVariable(loss);

        // subtract the preceding sum of losses to get the increment
        answer = loss - lossesUpTillNow;
    }
}



void
SteelRelaxMat::computeStressRelaxationStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep, ValueModeType mode) const
{
    SteelRelaxMatStatus *status = static_cast< SteelRelaxMatStatus * >( this->giveStatus(gp) );
    double averageStress = 0.;
    double mu = 0.;
    double k, rho, lambda;
    double dt = tStep->giveTimeIncrement();

    // get average strain due to stress relaxation
    double averageRelaxationStrain = 0.5 * ( status->giveTempRelaxIntVariable() + status->giveRelaxIntVariable() );

    // get average stress
    int n = 0;
    if ( status->giveStressVector().at(1) > 0. ) {
        averageStress += status->giveStressVector().at(1);
        n++;
    }

    if ( status->giveTempStressVector().at(1) > 0 ) {
        averageStress += status->giveTempStressVector().at(1);
        n++;
    }

    if ( n > 0 ) {
        averageStress /= n;
    }

    if ( !status->givePrestress() ) {
        status->setTempPrestress(averageStress);
    }


    double averageMechStrain = averageRelaxationStrain + averageStress / this->E;

    double F = averageMechStrain * this->E;

    double relaxStrainIncrement;
    // different approach for the first step after prestressing and for the following steps
    // when the prestressing is applied, the losses are assumed as zero
    if (  status->giveRelaxIntVariable() == 0. ) {
        // use Eurocode function for the first prestress loss
        this->evalStressRelaxationAtConstStrain(relaxStrainIncrement, gp, dt);
        // convert into strain
        relaxStrainIncrement /= this->E;
    } else {
        // governing equations 3.28-3.30 rewritten to the following form
        //
        // original:
        // mu = sigma_0 / fk
        // t in hours
        // delta sigma_p / sigma_0 = k1 * rho1000 * exp( k2 * mu) * (t/1000)^(0.75*(1-mu) ) * 1e-5;

        // new: (according to Bazant & Yu 2013)
        // delta sigma_p / F(epsilon) = rho * (t/lambda)^k;
        // where ... rho = k1 * rho1000 * exp( k2 * mu) * 1e-5
        //       ... k = 0.75*(1-mu)
        //       ... lambda = 1000 * timeFactor / 24
        //       ... timeFactor = 1 for day / 24 for hour / 86400 for sec, etc.

        relaxStrainIncrement = 0.;
        double prestress = status->giveTempPrestress();


        if ( averageStress > this->relRelaxBound * this->charStrength ) {
            mu = prestress / this->charStrength;
            rho = this->k1 * this->rho1000 * exp(this->k2 * mu) * 1.e-5;
            k = 0.75 * ( 1. - mu );
            lambda = 1000. * this->timeFactor / 24.;
            relaxStrainIncrement = k * pow(rho, 1. / k) * F * dt;
            relaxStrainIncrement /= this->E * lambda * pow( ( 1. - averageStress / F ), 1. / k - 1.);
        }
    }

    status->setTempRelaxIntVariable(relaxStrainIncrement + status->giveRelaxIntVariable() );

    answer.resize(1);

    if ( mode == VM_Incremental ) {
        answer.at(1) = relaxStrainIncrement;
    } else {
        answer.at(1) = relaxStrainIncrement + status->giveRelaxIntVariable();
    }

    return;
}


int
SteelRelaxMat::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return StructuralMaterial::giveIPValue(answer, gp, type, tStep);
}

//=============================================================================

SteelRelaxMatStatus::SteelRelaxMatStatus(GaussPoint *g) : StructuralMaterialStatus(g)
{}


void
SteelRelaxMatStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, " relaxationInternalVariable  ");
    fprintf(file, "%.4e ", relaxIntVariable);
    fprintf(file, "\n");

    fprintf(file, " prestress  ");
    fprintf(file, "%.4e ", prestress);
    fprintf(file, "\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void SteelRelaxMatStatus::initTempStatus()
{
    StructuralMaterialStatus::initTempStatus();

    tempRelaxIntVariable = relaxIntVariable;
    tempPrestress = prestress;
}


// updates internal variables when equilibrium is reached
void
SteelRelaxMatStatus::updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus::updateYourself(tStep);

    relaxIntVariable = tempRelaxIntVariable;

    prestress = tempPrestress;

    //    if ( this->prestress > 0.) {
    // prestressFlag = true
}


void
SteelRelaxMatStatus::saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus::saveContext(stream, mode);

    if ( !stream.write(relaxIntVariable) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(prestress) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
SteelRelaxMatStatus::restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus::restoreContext(stream, mode);

    if ( !stream.read(relaxIntVariable) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(prestress) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem
