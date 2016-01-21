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

#include "hydratingconcretemat.h"
#include "gausspoint.h"
#include "timestep.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HydratingConcreteMat);

HydratingConcreteMat :: HydratingConcreteMat(int n, Domain *d) : IsotropicHeatTransferMaterial(n, d)
{
    // constructor
    maxModelIntegrationTime = 0.;
    minModelTimeStepIntegrations = 0;
    P1 = 0.;
}


HydratingConcreteMat :: ~HydratingConcreteMat()
{
    // destructor
}


IRResultType
HydratingConcreteMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    // set conductivity k and capacity c
    result = IsotropicHeatTransferMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    referenceTemperature = 25.;
    IR_GIVE_OPTIONAL_FIELD(ir, referenceTemperature, _IFT_HydratingConcreteMat_referenceTemperature);

    IR_GIVE_FIELD(ir, hydrationModelType, _IFT_HydratingConcreteMat_hydrationModelType);

    /*hydrationModelType==1:exponential hydration model, summarized in A.K. Schindler and K.J. Folliard: Heat of Hydration Models for Cementitious
     *  Materials, ACI Materials Journal, 2005
     * hydrationModelType==2: affinity hydration model inspired by Miguel Cervera and Javier Oliver and Tomas Prato:Thermo-chemo-mechanical model
     *  for concrete. I: Hydration and aging, Journal of Engineering Mechanics ASCE, 125(9), 1018-1027, 1999
     */
    if ( hydrationModelType == 1 ) {
        IR_GIVE_FIELD(ir, tau, _IFT_HydratingConcreteMat_tau); // [s]
        IR_GIVE_FIELD(ir, beta, _IFT_HydratingConcreteMat_beta); // [-]
        IR_GIVE_FIELD(ir, DoHInf, _IFT_HydratingConcreteMat_DoHInf);
    } else if ( hydrationModelType == 2 ) {
        IR_GIVE_FIELD(ir, B1, _IFT_HydratingConcreteMat_B1); // [1/s]
        IR_GIVE_FIELD(ir, B2, _IFT_HydratingConcreteMat_B2);
        IR_GIVE_FIELD(ir, eta, _IFT_HydratingConcreteMat_eta);
        IR_GIVE_FIELD(ir, DoHInf, _IFT_HydratingConcreteMat_DoHInf);
        IR_GIVE_OPTIONAL_FIELD(ir, DoH1, _IFT_HydratingConcreteMat_DoH1);
        IR_GIVE_OPTIONAL_FIELD(ir, P1, _IFT_HydratingConcreteMat_P1);
    } else {
        OOFEM_WARNING("Unknown hdyration model type %d", hydrationModelType);
        return IRRT_BAD_FORMAT;
    }

    IR_GIVE_FIELD(ir, Qpot, _IFT_HydratingConcreteMat_qpot); // [1/s]

    IR_GIVE_FIELD(ir, massCement, _IFT_HydratingConcreteMat_massCement);

    maxModelIntegrationTime = 36000;
    IR_GIVE_OPTIONAL_FIELD(ir, maxModelIntegrationTime, _IFT_HydratingConcreteMat_maxModelIntegrationTime);

    minModelTimeStepIntegrations = 30.;
    IR_GIVE_OPTIONAL_FIELD(ir, minModelTimeStepIntegrations, _IFT_HydratingConcreteMat_minModelTimeStepIntegrations);


    conductivityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, conductivityType, _IFT_HydratingConcreteMat_conductivitytype);
    capacityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, capacityType, _IFT_HydratingConcreteMat_capacitytype);
    densityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, densityType, _IFT_HydratingConcreteMat_densitytype);

    activationEnergy = 38400; //J/mol/K
    IR_GIVE_OPTIONAL_FIELD(ir, activationEnergy, _IFT_HydratingConcreteMat_activationEnergy);

    reinforcementDegree = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, reinforcementDegree, _IFT_HydratingConcreteMat_reinforcementDegree);

    return IRRT_OK;
}

// returns hydration power [W/m3 of concrete]
void
HydratingConcreteMat :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
{
    val.resize(1);
    if ( mode == VM_Total ) {
        val.at(1) = this->GivePower(tStep, gp, mode);
    } else {
        OOFEM_ERROR("Undefined mode %s\n", __ValueModeTypeToString(mode) );
    }
}


double
HydratingConcreteMat :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( mode == Capacity ) {
        return ( giveConcreteCapacity(gp) * giveConcreteDensity(gp) );
    } else if ( mode == IntSource ) { //for nonlinear solver, return dHeat/dTemperature
        HydratingConcreteMatStatus *ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
        //it suffices to compute derivative of scaling Arrhenius equation with respect to temporary temperature
        double stateVec = ms->giveField().at(1) + 273.15;
        double tempStateVec = ms->giveTempField().at(1) + 273.15;
        return this->activationEnergy / ( 8.314 * tempStateVec * tempStateVec ) * exp(1. / stateVec -  1. / tempStateVec);
    } else {
        OOFEM_ERROR("unknown mode (%s)\n", __MatResponseModeToString(mode) );
    }

    return 0.;
}


double HydratingConcreteMat :: giveIsotropicConductivity(GaussPoint *gp)
{
    HydratingConcreteMatStatus *ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
    double conduct;

    if ( conductivityType == 0 ) { //given from input file
        conduct = IsotropicHeatTransferMaterial :: give('k', gp);
    } else if ( conductivityType == 1 ) { //compute according to Ruiz, Schindler, Rasmussen. Kim, Chang: Concrete temperature modeling and strength prediction using maturity concepts in the FHWA HIPERPAV software, 7th international conference on concrete pavements, Orlando (FL), USA, 2001
        conduct = IsotropicHeatTransferMaterial :: give('k', gp) * ( 1.0 - 0.33 / 1.33 * ms->giveDoHActual() );
    } else {
        OOFEM_ERROR("Unknown conductivityType %d\n", conductivityType);
        conduct = 0.;
    }

    //Parallel Voigt model, 20 W/m/K for steel
    conduct = conduct * ( 1. - this->reinforcementDegree ) + 20. * this->reinforcementDegree;

    if ( conduct < 0.3 || conduct > 5 ) {
        OOFEM_WARNING("Weird concrete thermal conductivity %f W/m/K\n", conduct);
    }

    return conduct;
}

//normally it returns J/kg/K of concrete
double HydratingConcreteMat :: giveConcreteCapacity(GaussPoint *gp)
{
    double capacityConcrete;

    if ( capacityType == 0 ) { //given from OOFEM input file
        capacityConcrete = IsotropicHeatTransferMaterial :: give('c', gp);
    } else if ( capacityType == 1 ) { ////calculate from 5-component model
        OOFEM_ERROR("Calculate from 5-component model, not implemented in capacityType %d\n", capacityType);
        capacityConcrete = 0.;
    } else {
        OOFEM_ERROR("Unknown capacityType %d\n", capacityType);
        capacityConcrete = 0.;
    }

    //Parallel Voigt model, 500 J/kg/K for steel
    capacityConcrete = capacityConcrete * ( 1. - this->reinforcementDegree ) + 500. * this->reinforcementDegree;

    if ( capacityConcrete < 500 || capacityConcrete > 2000 ) {
        OOFEM_WARNING("Weird concrete heat capacity %f J/kg/K\n", capacityConcrete);
    }

    return capacityConcrete;
}


double HydratingConcreteMat :: giveConcreteDensity(GaussPoint *gp)
{
    double concreteBulkDensity;

    if ( densityType == 0 ) { //get from input file
        concreteBulkDensity = IsotropicHeatTransferMaterial :: give('d', gp);
    } else if ( densityType == 1 ) { //calculate from 5-component model - not implemented
        OOFEM_ERROR("Calculate from 5-component model, not implemented in densityType %d\n", densityType);
        concreteBulkDensity = 0.;
    } else {
        OOFEM_ERROR("Unknown densityType %d\n", densityType);
        concreteBulkDensity = 0.;
    }

    //Parallel Voigt model, 7850 kg/m3 for steel
    concreteBulkDensity = concreteBulkDensity * ( 1. - this->reinforcementDegree ) + 7850. * this->reinforcementDegree;

    if ( concreteBulkDensity < 1000 || concreteBulkDensity > 4000 ) {
        OOFEM_WARNING("Weird concrete density %f kg/m3\n", concreteBulkDensity);
    }

    return concreteBulkDensity;
}


int
HydratingConcreteMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    // printf ("IP %d::giveIPValue, IST %d", giveNumber(), type);
    if ( type == IST_HydrationDegree ) {
        HydratingConcreteMatStatus *status = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
        answer.resize(1);
        answer.at(1) = status->giveDoHActual();
        //else answer.at(1) = 0;
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


MaterialStatus *
HydratingConcreteMat :: CreateStatus(GaussPoint *gp) const
{
    return new HydratingConcreteMatStatus(1, domain, gp);
}


HydratingConcreteMatStatus :: HydratingConcreteMatStatus(int n, Domain *d, GaussPoint *g) : TransportMaterialStatus(n, d, g)
{
    //constructor
    power = 0.;
    lastEvalTime = -1.e20; //start from begining (set to -1.e6 s)
    degreeOfHydration = 0.;
    lastDegreeOfHydration = 0.;
    lastEquivalentTime = 0.;
    equivalentTime = 0.;
}


HydratingConcreteMatStatus :: ~HydratingConcreteMatStatus()
{
    //destructor
}

//linear solver (NonStationaryTransportProblem) IntrinsicTime = TargetTime
//nonlinear solver (NLTransientTransportProblem) IntrinsicTime depends on alpha
double HydratingConcreteMat :: GivePower(TimeStep *tStep, GaussPoint *gp, ValueModeType mode)
{
    HydratingConcreteMatStatus *ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
    double castingTime = this->giveCastingTime();
    double evalTime = tStep->giveIntrinsicTime();
    double targTime = tStep->giveTargetTime();

    if ( tStep->giveNumber() == 0 ) {
        return 0;
    }

    //do not calculate anything before casting time
    if ( targTime - castingTime <= 0 ) {
        ms->power = 0.;
        return 0.;
    }
    
    //return computed power (homexportmodule)
    if ( evalTime == ms->lastEvalTime ) {
        return ms->power;
    }

    if ( this->hydrationModelType == 1 ) { //exponential affinity hydration model, need to keep equivalent time
        ms->equivalentTime = ms->lastEquivalentTime + ( evalTime - ms->lastEvalTime ) * scaleTemperature(gp);
        if ( ms->equivalentTime != 0. ) {
            ms->degreeOfHydration = this->DoHInf * exp( -pow(this->tau / ms->equivalentTime, this->beta) );
            //printf("%f %f %f %f\n", equivalentTime, this->lastEquivalentTime, evalTime, lastEvalTime);
        }
    } else if ( this->hydrationModelType == 2 ) { //affinity hydration model inspired by Miguel Cervera et al.
        //determine timeStep for integration
        double alphaTrialOld, alphaTrialNew = 0.0;
        double time = ms->lastEvalTime;
        double timeStep = ( evalTime - time ) / this->minModelTimeStepIntegrations;
        if ( timeStep > this->maxModelIntegrationTime ) {
            timeStep = this->maxModelIntegrationTime;
        }
        ms->degreeOfHydration = ms->lastDegreeOfHydration;
        //integration loop through hydration model at a given TimeStep
        while ( time < evalTime ) {
            if ( time + timeStep > evalTime ) {
                timeStep = evalTime - time;
            } else {
                time += timeStep;
            }
            //printf("%f %f %f %f\n", time, affinity, scaleTemperature(), degreeOfHydration);
            alphaTrialOld = ms->degreeOfHydration + scaleTemperature(gp) * affinity25(ms->degreeOfHydration) * timeStep; //predictor
            //http://en.wikipedia.org/wiki/Predictor%E2%80%93corrector_method
            //corrector - integration through trapezoidal rule
            //3 loops normally suffices
            for ( int i = 0; i < 4; i++ ) {
                alphaTrialNew = ms->degreeOfHydration + scaleTemperature(gp) * timeStep / 2. * ( affinity25(ms->degreeOfHydration) + affinity25(alphaTrialOld) );
                alphaTrialOld = alphaTrialNew;
            }
            ms->degreeOfHydration = alphaTrialNew;
        }
    } else {
        OOFEM_ERROR("Unknown hydration model type %d", this->hydrationModelType);
    }

    ms->power = this->Qpot * ( ms->degreeOfHydration - ms->lastDegreeOfHydration ) / ( evalTime - ms->lastEvalTime );
    ms->power *= 1000 * this->massCement; // W/m3 of concrete

    //internal variables are updated in HydratingConcreteMatStatus :: updateYourself()
    return ms->power;
}


double HydratingConcreteMat :: scaleTemperature(GaussPoint *gp)
{
    HydratingConcreteMatStatus *ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
    return exp( this->activationEnergy / 8.314 * ( 1. / ( 273.15 + this->referenceTemperature ) - 1. / ( 273.15 + ms->giveTempField().at(1) ) ) );
}

double HydratingConcreteMat :: affinity25(double DoH)
{
    double result =  this->B1 * ( this->B2 / this->DoHInf + DoH ) * ( this->DoHInf - DoH ) * exp(-this->eta * DoH / this->DoHInf);
    if ( result < 0. ) { //numerical instabilities
        return 0.;
    }

    //add slag reaction
    if ( this->P1 != 0. && DoH >= this->DoH1 ) {
        result *= 1. + this->P1 * ( DoH - this->DoH1 );
    }
    return result;
}

double HydratingConcreteMatStatus :: giveDoHActual()
{
    return degreeOfHydration;
}


void
HydratingConcreteMatStatus :: updateYourself(TimeStep *tStep)
{
    HydratingConcreteMat *mat = static_cast< HydratingConcreteMat * >( this->gp->giveMaterial() );
    this->lastEvalTime = tStep->giveIntrinsicTime(); //where heat power was evaluated in the last equilibrium
    this->lastEquivalentTime = this->equivalentTime;
    this->lastDegreeOfHydration = this->degreeOfHydration;
    //average from last and current temperatures, in C*hour
    if ( !tStep->isIcApply() && mat->giveCastingTime() < tStep->giveIntrinsicTime() ) {
        this->maturity += ( ( this->giveField().at(1) + this->giveTempField().at(1) ) / 2. - mat->giveMaturityT0() ) * tStep->giveTimeIncrement() / 3600.;
    }
    TransportMaterialStatus :: updateYourself(tStep);
}


void
HydratingConcreteMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    HydratingConcreteMat *mat = static_cast< HydratingConcreteMat * >( this->gp->giveMaterial() );
    TransportMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "   status {");
    fprintf( file, "EvaluatingTime %e  DoH %f HeatPower %f [W/m3 of concrete] Temperature %f conductivity %f  capacity %f  density %f", tStep->giveIntrinsicTime(), this->giveDoHActual(), this->power, this->giveTempField().at(1), mat->giveIsotropicConductivity(this->gp), mat->giveConcreteCapacity(this->gp), mat->giveConcreteDensity(this->gp) );
    fprintf(file, "}\n");
}
} // end namespace oofem
