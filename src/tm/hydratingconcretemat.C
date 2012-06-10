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

#include "hydratingconcretemat.h"
#include "gausspnt.h"
#include "timestep.h"
#include "mathfem.h"

namespace oofem {
HydratingConcreteMat :: HydratingConcreteMat(int n, Domain *d) : IsotropicHeatTransferMaterial(n, d)
{
    // constructor
    maxModelIntegrationTime = 0.;
    minModelTimeStepIntegrations = 0;
}


HydratingConcreteMat :: ~HydratingConcreteMat()
{
    // destructor
}


IRResultType
HydratingConcreteMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    // set conductivity k and capacity c
    IsotropicHeatTransferMaterial :: initializeFrom(ir);

    referenceTemperature = 25.;
    IR_GIVE_OPTIONAL_FIELD(ir, referenceTemperature, IFT_HydratingConcreteMat_referenceTemperature, "referencetemperature"); // Macro

    IR_GIVE_FIELD(ir, hydrationModelType, IFT_HydratingConcreteMat_hydrationModelType, "hydrationmodeltype");

    /*hydrationModelType==1:exponential hydration model, summarized in A.K. Schindler and K.J. Folliard: Heat of Hydration Models for Cementitious
     *  Materials, ACI Materials Journal, 2005
     * hydrationModelType==2: affinity hydration model inspired by Miguel Cervera and Javier Oliver and Tomas Prato:Thermo-chemo-mechanical model
     *  for concrete. I: Hydration and aging, Journal of Engineering Mechanics ASCE, 125(9), 1018-1027, 1999
     */
    if ( hydrationModelType == 1 ) {
        IR_GIVE_FIELD(ir, tau, IFT_HydratingConcreteMat_tau, "tau"); // [s]
        IR_GIVE_FIELD(ir, beta, IFT_HydratingConcreteMat_beta, "beta"); // [-]
        IR_GIVE_FIELD(ir, DoHInf, IFT_HydratingConcreteMat_DoHInf, "dohinf");
    } else if ( hydrationModelType == 2 ) {
        IR_GIVE_FIELD(ir, B1, IFT_HydratingConcreteMat_B1, "b1"); // [1/s]
        IR_GIVE_FIELD(ir, B2, IFT_HydratingConcreteMat_B2, "b2");
        IR_GIVE_FIELD(ir, eta, IFT_HydratingConcreteMat_eta, "eta");
        IR_GIVE_FIELD(ir, DoHInf, IFT_HydratingConcreteMat_DoHInf, "dohinf");
    } else {
        OOFEM_ERROR2("Unknown hdyration model type %d", hydrationModelType);
    }

    IR_GIVE_FIELD(ir, Qpot, IFT_HydratingConcreteMat_qpot, "qpot"); // [1/s]

    IR_GIVE_FIELD(ir, massCement, IFT_HydratingConcreteMat_massCement, "masscement");

    maxModelIntegrationTime = 36000;
    IR_GIVE_OPTIONAL_FIELD(ir, maxModelIntegrationTime, IFT_maxModelIntegrationTime, "maxmodelintegrationtime"); // Macro

    minModelTimeStepIntegrations = 30.;
    IR_GIVE_OPTIONAL_FIELD(ir, minModelTimeStepIntegrations, IFT_minModelTimeStepIntegrations, "minmodeltimestepintegrations"); // Macro


    conductivityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, conductivityType, IFT_HydratingConcreteMat_conductivitytype, "conductivitytype"); // Macro
    capacityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, capacityType, IFT_HydratingConcreteMat_capacitytype, "capacitytype"); // Macro
    densityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, densityType, IFT_HydratingConcreteMat_densitytype, "densitytype"); // Macro

    activationEnergy = 38400; //J/mol/K
    IR_GIVE_OPTIONAL_FIELD(ir, activationEnergy, IFT_HydratingConcreteMat_activationEnergy, "activationenergy"); // Macro

    reinforcementDegree = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, reinforcementDegree, IFT_HydratingConcreteMat_reinforcementDegree, "reinforcementdegree"); // Macro

    return IRRT_OK;
}

// returns hydration power [W/m3 of concrete]
void
HydratingConcreteMat :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    HydratingConcreteMatStatus *ms = ( HydratingConcreteMatStatus * ) this->giveStatus(gp);
    val.resize(1);
    if ( mode == VM_Total ) {
        val.at(1) = ms->GivePower(atTime);
    } else {
        OOFEM_ERROR2( "Undefined mode %s\n", __ValueModeTypeToString(mode) );
    }
}


void
HydratingConcreteMat :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *atTime)
{
    HydratingConcreteMatStatus *ms = ( HydratingConcreteMatStatus * ) this->giveStatus(gp);
    ms->letTempStateVectorBe(vec);
}


double
HydratingConcreteMat :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    if ( mode == Capacity ) {
        return ( giveConcreteCapacity(gp) * giveConcreteDensity(gp) );
    } else if ( mode == IntSource ) { //for nonlinear solver, return dHeat/dTemperature
        HydratingConcreteMatStatus *ms = ( HydratingConcreteMatStatus * ) this->giveStatus(gp);
        //it suffices to compute derivative of scaling Arrhenius equation with respect to temporary temperature
        double stateVec = ms->giveStateVector().at(1) + 273.15;
        double tempStateVec = ms->giveTempStateVector().at(1) + 273.15;
        return this->activationEnergy / ( 8.314 * tempStateVec * tempStateVec) * exp(1. / stateVec -  1. / tempStateVec ) ;
    } else {
        OOFEM_ERROR2( "giveCharacteristicValue : unknown mode (%s)\n", __MatResponseModeToString(mode) );
    }

    return 0.;
}


double HydratingConcreteMat :: giveConcreteConductivity(GaussPoint *gp)
{
    HydratingConcreteMatStatus *ms = ( HydratingConcreteMatStatus * ) this->giveStatus(gp);
    double conduct;

    if ( conductivityType == 0 ) { //given from input file
        conduct = IsotropicHeatTransferMaterial :: give('k', gp);
    } else if ( conductivityType == 1 ) { //compute according to Ruiz, Schindler, Rasmussen. Kim, Chang: Concrete temperature modeling and strength prediction using maturity concepts in the FHWA HIPERPAV software, 7th international conference on concrete pavements, Orlando (FL), USA, 2001
        conduct = IsotropicHeatTransferMaterial :: give('k', gp) * ( 1.0 - 0.33 / 1.33 * ms->giveDoHActual() );
    } else {
        OOFEM_ERROR2("Unknown conductivityType %d\n", conductivityType);
    }

    //Parallel Voigt model, 20 W/m/K for steel
    conduct = conduct * ( 1. - this->reinforcementDegree ) + 20. * this->reinforcementDegree;

    if ( conduct < 0.3 || conduct > 5 ) {
        OOFEM_WARNING2("Weird concrete thermal conductivity %f W/m/K\n", conduct);
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
        //not yet implemented
    } else {
        OOFEM_ERROR2("Unknown capacityType %d\n", capacityType);
    }

    //Parallel Voigt model, 500 J/kg/K for steel
    capacityConcrete = capacityConcrete * ( 1. - this->reinforcementDegree ) + 500. * this->reinforcementDegree;

    if ( capacityConcrete < 500 || capacityConcrete > 2000 ) {
        OOFEM_WARNING2("Weird concrete heat capacity %f J/kg/K\n", capacityConcrete);
    }

    return capacityConcrete;
}


double HydratingConcreteMat :: giveConcreteDensity(GaussPoint *gp)
{
    double concreteBulkDensity;

    if ( densityType == 0 ) { //get from input file
        concreteBulkDensity = IsotropicHeatTransferMaterial :: give('d', gp);
    } else if ( densityType == 1 ) { //calculate from 5-component model - not implemented
        //not yet implemented
    } else {
        OOFEM_ERROR2("Unknown densityType %d\n", densityType);
    }

    //Parallel Voigt model, 7850 kg/m3 for steel
    concreteBulkDensity = concreteBulkDensity * ( 1. - this->reinforcementDegree ) + 7850. * this->reinforcementDegree;

    if ( concreteBulkDensity < 1000 || concreteBulkDensity > 4000 ) {
        OOFEM_WARNING2("Weird concrete density %f kg/m3\n", concreteBulkDensity);
    }

    return concreteBulkDensity;
}


int
HydratingConcreteMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    // printf ("IP %d::giveIPValue, IST %d", giveNumber(), type);
    if ( type == IST_HydrationDegree ) {
        HydratingConcreteMatStatus *status = ( HydratingConcreteMatStatus * ) this->giveStatus(gp);
        answer.resize(1);
        answer.at(1) = status->giveDoHActual();
        //else answer.at(1) = 0;
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, atTime);
    }
}


InternalStateValueType
HydratingConcreteMat :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_HydrationDegree ) {
        return ISVT_SCALAR;
    } else {
        return TransportMaterial :: giveIPValueType(type);
    }
}


int
HydratingConcreteMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return TransportMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
HydratingConcreteMat :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_HydrationDegree ) {
        return 1;
    } else {
        return TransportMaterial :: giveIPValueSize(type, aGaussPoint);
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
    lastIntrinsicTime = -1.e6; //start from begining (set to -1.e6 s)
    degreeOfHydration = 0.;
    lastDegreeOfHydration = 0.;
    lastEquivalentTime = 0.;
    equivalentTime = 0.;
}


HydratingConcreteMatStatus :: ~HydratingConcreteMatStatus()
{
    //destructor
}


double HydratingConcreteMatStatus :: GivePower(TimeStep *atTime)
{
    double castingTime = this->gp->giveMaterial()->giveCastingTime();
    HydratingConcreteMat *mat = ( HydratingConcreteMat * ) this->gp->giveMaterial();
    double intrinsicTime = atTime->giveIntrinsicTime();
    double targTime = atTime->giveTargetTime();

    if ( atTime->giveNumber() == 0 ) {
        return 0;
    }

    //do not calculate anything before casting time
    if ( targTime - castingTime <= 0 ) {
        power = 0.;
        return 0.;
    }

    if ( mat->hydrationModelType == 1 ) { //exponential affinity hydration model, need to keep equivalent time
        equivalentTime = this->lastEquivalentTime + ( intrinsicTime - lastIntrinsicTime ) * scaleTemperature();
        if ( equivalentTime != 0. ) {
            degreeOfHydration = mat->DoHInf * exp( -pow(mat->tau / equivalentTime, mat->beta) );
            //printf("%f %f %f %f\n", equivalentTime, this->lastEquivalentTime, intrinsicTime, lastIntrinsicTime);
        }
    } else if ( mat->hydrationModelType == 2 ) { //affinity hydration model inspired by Miguel Cervera et al.
        //determine timeStep for integration
        double alphaTrialOld, alphaTrialNew;
        double time = lastIntrinsicTime;
        double timeStep = ( intrinsicTime - time ) / mat->minModelTimeStepIntegrations;
        if ( timeStep > mat->maxModelIntegrationTime ) {
            timeStep = mat->maxModelIntegrationTime;
        }
        degreeOfHydration = lastDegreeOfHydration;
        //integration loop through hydration model at a given TimeStep
        while ( time < intrinsicTime ) {
            if ( time + timeStep > intrinsicTime ) {
                timeStep = intrinsicTime - time;
            } else {
                time += timeStep;
            }
            //printf("%f %f %f %f\n", time, affinity, scaleTemperature(), degreeOfHydration);
            alphaTrialOld = degreeOfHydration + scaleTemperature() * affinity25(degreeOfHydration) * timeStep;//predictor
            //http://en.wikipedia.org/wiki/Predictor%E2%80%93corrector_method
            //corrector - integration through trapezoidal rule
            //3 loops normally suffices
            for(int i=0; i<4; i++){
                alphaTrialNew = degreeOfHydration + scaleTemperature()*timeStep/2.*(affinity25(degreeOfHydration)+affinity25(alphaTrialOld) );
                alphaTrialOld = alphaTrialNew;
            }
            degreeOfHydration= alphaTrialNew;
        }
    } else {
        OOFEM_ERROR2("Unknown hydration model type %d", mat->hydrationModelType);
    }

    power = mat->Qpot * ( degreeOfHydration - lastDegreeOfHydration ) / ( intrinsicTime - lastIntrinsicTime );
    power *= 1000 * mat->massCement; // W/m3 of concrete

    //internal variables are updated in HydratingConcreteMatStatus :: updateYourself()
    return power;
}


double HydratingConcreteMatStatus :: scaleTemperature(void)
{
    HydratingConcreteMat *mat = ( HydratingConcreteMat * ) this->gp->giveMaterial();
    return exp( mat->activationEnergy / 8.314 * ( 1. / ( 273.15 + mat->referenceTemperature ) - 1. / ( 273.15 + this->giveTempStateVector().at(1) ) ) );
}

double HydratingConcreteMatStatus :: affinity25(double DoH)
{
    HydratingConcreteMat *mat = ( HydratingConcreteMat * ) this->gp->giveMaterial();
    double result =  mat->B1 * ( mat->B2 / mat->DoHInf + DoH ) * ( mat->DoHInf - DoH ) * exp(-mat->eta * DoH / mat->DoHInf);
    if (result<0.){//numerical instabilities
        return 0.;
    }
    return result;
}

double HydratingConcreteMatStatus :: giveDoHActual(void)
{
    return degreeOfHydration;
}


void
HydratingConcreteMatStatus :: updateYourself(TimeStep *atTime)
{
    this->lastIntrinsicTime = atTime->giveIntrinsicTime(); //where heat power was evaluated in last equilibrium
    this->lastEquivalentTime = this->equivalentTime;
    this->lastDegreeOfHydration = this->degreeOfHydration;
    TransportMaterialStatus :: updateYourself(atTime);
}


void
HydratingConcreteMatStatus :: printOutputAt(FILE *file, TimeStep *atTime)
{
    HydratingConcreteMat *mat = ( HydratingConcreteMat * ) this->gp->giveMaterial();
    TransportMaterialStatus :: printOutputAt(file, atTime);
    fprintf(file, "   status {");
    fprintf( file, "IntrinsicTime %e  DoH %f HeatPower %f [W/m3 of concrete] conductivity %f  capacity %f  density %f", atTime->giveIntrinsicTime(), this->giveDoHActual(), this->power, mat->giveConcreteConductivity(this->gp), mat->giveConcreteCapacity(this->gp), mat->giveConcreteDensity(this->gp) );
    fprintf(file, "}\n");
}

} // end namespace oofem
