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

// constructor
SteelRelaxMat :: SteelRelaxMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    E = 0.;
    k1 = 0.;
    k2 = 0.;
    rho1000 = 0.;
    mu = 0.;
    timeFactor = 0.;
    //stiffnessFactor = 1.e6;
    //    prestress = 0.;
}

// destructor
SteelRelaxMat :: ~SteelRelaxMat()
{}


int
SteelRelaxMat :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether the receiver supports the given mode
//
{
    return mode == _1dMat;
}



// reads the model parameters from the input file
IRResultType
SteelRelaxMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }


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

    //IR_GIVE_OPTIONAL_FIELD(ir, this->stiffnessFactor, _IFT_SteelRelaxMat_stiffnessFactor);

    IR_GIVE_FIELD(ir, this->timeFactor, _IFT_SteelRelaxMat_timeFactor);

    //    IR_GIVE_OPTIONAL_FIELD(ir, this->prestress, _IFT_SteelRelaxMat_prestress);

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
    //this->tolerance =  1. / this->stiffnessFactor;
    this->tolerance =  1. / 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->tolerance, _IFT_SteelRelaxMat_tolerance);


    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
SteelRelaxMat :: CreateStatus(GaussPoint *gp) const
{
    return new SteelRelaxMatStatus(1, this->giveDomain(), gp);
}

void
SteelRelaxMat :: giveRealStressVector(FloatArray &answer,
                                      GaussPoint *gp,
                                      const FloatArray &totalStrain,
                                      TimeStep *tStep)
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

        StructuralMaterial :: giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Incremental);

        strainIncrement.beDifferenceOf( reducedStrain, status->giveStrainVector() );
        stressIncrement = strainIncrement.at(1) * this->E;

        //   subtract stress increment due to prestress losses in the current time step


        if ( status->giveStressVector().giveSize() ) {
            stressVector = status->giveStressVector();
        } else {
            stressVector.resize(1);
            stressVector.zero();
        }

        stressVector.at(1) += stressIncrement;

        if ( stressVector.at(1) > 0. ) {
            this->computeIncrOfPrestressLossAtVarStrain( lossIncrement, gp, tStep, stressVector.at(1) );
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

            // subtracts both thermal strain increment and strain due to cable relaxation
            this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Incremental);

            strainIncrement.beDifferenceOf( reducedStrain, status->giveStrainVector() );
            stressIncrement = strainIncrement.at(1) * this->E;

            stressVector.at(1) += stressIncrement;


            // store the final strain and stress in the status
            status->letTempStrainVectorBe(totalStrain);
            status->letTempStressVectorBe(stressVector);

            i++;

            if ( i > 1000 ) {
                OOFEM_ERROR("Algorithm not converging");
            }
        } while ( fabs( prevIterTempStress - status->giveTempStressVector().at(1) ) >= this->tolerance );


        if ( i > 30 ) {
            OOFEM_WARNING("Criterion of the algorithm reached in %d iterations, consider increasing tolerance", i);
        }
    }

    if ( stressVector.at(1) > this->charStrength ) {
        OOFEM_ERROR( "Stress %f exeeds the characteristic strength of the material!", stressVector.at(1) );
    }

    answer.resize(1);
    answer.at(1) = stressVector.at(1);
}


void
SteelRelaxMat :: give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.zero();
    if ( this->isActivated(tStep) ) {
        answer.at(1, 1) = this->E;
    }
}




void
SteelRelaxMat :: giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep, ValueModeType mode)
{
    FloatArray temperatureFreeStrain;
    FloatArray relaxationStrain;

    // subtracts temperature strain
    StructuralMaterial :: giveStressDependentPartOfStrainVector(temperatureFreeStrain, gp, totalStrain, tStep, mode);

    // strains due to relaxation
    this->computeStressRelaxationStrainVector(relaxationStrain, gp, totalStrain, tStep, mode);

    answer = temperatureFreeStrain;
    answer.subtract(relaxationStrain);

    return;
}


void
SteelRelaxMat :: evalStressRelaxationAtConstStrain(double &answer, GaussPoint *gp, double dt)
{
    double rho;
    double k;
    double lambda;
    double mu;
    double prestress;

    SteelRelaxMatStatus *status = static_cast< SteelRelaxMatStatus * >( this->giveStatus(gp) );
    prestress = status->givePrestress();

    answer = 0.;

    if ( prestress >  this->relRelaxBound * this->charStrength ) {
        mu = prestress / this->charStrength;

        rho = this->k1 * this->rho1000 * exp(this->k2 * mu) * 1.e-5;
        k = 0.75 * ( 1. - mu );
        lambda = 1000. * this->timeFactor / 24.;

        answer = prestress * rho * pow( ( dt / lambda ), k );
    }
}


void
SteelRelaxMat :: computeIncrOfPrestressLossAtVarStrain(double &answer, GaussPoint *gp, TimeStep *tStep, double stress)
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
    status->setPrestress(prestress);

    answer = 0.;

    if ( prestress > this->relRelaxBound * this->charStrength ) {
        mu = prestress / this->charStrength;
        rho = this->k1 * this->rho1000 * exp(this->k2 * mu) * 1.e-5;
        k = 0.75 * ( 1. - mu );
        lambda = 1000. * this->timeFactor / 24.;


        // compute total loss for updated prestress and time equiv
        double t_equiv;
        t_equiv = pow( ( lossesUpTillNow / ( prestress * rho ) ), ( 1. / k ) ) * lambda;
        this->evalStressRelaxationAtConstStrain( loss, gp, t_equiv + tStep->giveTimeIncrement() );


        // set temporary sum of losses
        status->setTempRelaxIntVariable(loss);

        // subtract the preceding sum of losses to get the increment
        loss -= lossesUpTillNow;

        answer = loss;
    }
}







void
SteelRelaxMat :: computeStressRelaxationStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep, ValueModeType mode)
{
    // here we deal with total strain vector with subtracted tempreature

    SteelRelaxMatStatus *status = static_cast< SteelRelaxMatStatus * >( this->giveStatus(gp) );

    double averageStress = 0.;
    int n = 0;
    double mu = 0.;

    double k, rho, lambda;

    double dt = tStep->giveTimeIncrement();

    FloatArray temperFreeStrain, temperFreeStrainIncrement;
    StructuralMaterial :: giveStressDependentPartOfStrainVector(temperFreeStrain, gp, totalStrain, tStep, VM_Total);

    double averageMechStrain = temperFreeStrain.at(1);

    double prestress;
    prestress = status->givePrestress();

    if ( prestress > 0. ) {
        temperFreeStrainIncrement = totalStrain; // result = eps i+1 tot
        temperFreeStrainIncrement.subtract( status->giveStrainVector() ); // result = delta eps tot

        // epsilon temperature increment
        FloatArray deltaEpsTemperature;
        this->computeStressIndependentStrainVector(deltaEpsTemperature, gp, tStep, VM_Incremental);

        temperFreeStrainIncrement.subtract(deltaEpsTemperature); // results = delta epsilon stress

        averageMechStrain -= 0.5 * temperFreeStrainIncrement.at(1);
    }

    double F = averageMechStrain * this->E;


    if ( prestress == 0. ) {
        prestress = F;
        status->setPrestress(prestress);
    }

    double relaxStrainIncrement;

    // different approach for the first step and for the following steps
    if (  status->giveRelaxIntVariable() == 0. ) {
        // assume that the strain is constant and equal to zero - eval from Eurocode function
        this->evalStressRelaxationAtConstStrain(relaxStrainIncrement, gp, dt);
        // convert into strain
        relaxStrainIncrement /= this->E;
    } else {
        if ( status->giveStressVector().giveSize() ) {
            averageStress += status->giveStressVector().at(1);
            n++;
        }

        if ( status->giveTempStressVector().giveSize() ) {
            averageStress += status->giveTempStressVector().at(1);
            n++;
        }

        if ( n > 0 ) {
            averageStress /= n;
        }

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

        if ( averageStress > this->relRelaxBound * this->charStrength ) {
            mu = prestress / this->charStrength;
            rho = this->k1 * this->rho1000 * exp(this->k2 * mu) * 1.e-5;
            k = 0.75 * ( 1. - mu );
            lambda = 1000. * this->timeFactor / 24.;
            relaxStrainIncrement = k * pow(rho, 1. / k) * F * dt;
            relaxStrainIncrement /= this->E * lambda * pow( ( 1. - averageStress / F ), 1. / k - 1. );
        }
    }

    status->setTempRelaxIntVariable( relaxStrainIncrement + status->giveRelaxIntVariable() );

    answer.resize(1);

    if ( mode == VM_Incremental ) {
        answer.at(1) = relaxStrainIncrement;
    } else {
        answer.at(1) = relaxStrainIncrement + status->giveRelaxIntVariable();
    }

    return;
}


int
SteelRelaxMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
}

//=============================================================================

SteelRelaxMatStatus :: SteelRelaxMatStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    relaxIntVariable = tempRelaxIntVariable = 0.;
    prestress = 0.;
}

SteelRelaxMatStatus :: ~SteelRelaxMatStatus()
{ }

void
SteelRelaxMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, " relaxationInternalVariable  ");
    fprintf(file, "%.4e ", relaxIntVariable);
    fprintf(file, "\n");

    fprintf(file, " initialPrestress  ");
    fprintf(file, "%.4e ", prestress);
    fprintf(file, "\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void SteelRelaxMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    tempRelaxIntVariable = relaxIntVariable;
}


// updates internal variables when equilibrium is reached
void
SteelRelaxMatStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    relaxIntVariable = tempRelaxIntVariable;
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
SteelRelaxMatStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( !stream.write(relaxIntVariable) ) {
        return CIO_IOERR;
    }

    if ( !stream.write(prestress) ) {
        return CIO_IOERR;
    }

    return CIO_OK;
}



contextIOResultType
SteelRelaxMatStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(relaxIntVariable) ) {
        return CIO_IOERR;
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
