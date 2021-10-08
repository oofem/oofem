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

#include "tm/Materials/hydratingconcretemat.h"
#include "gausspoint.h"
#include "timestep.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HydratingConcreteMat);

HydratingConcreteMat :: HydratingConcreteMat(int n, Domain *d) : IsotropicHeatTransferMaterial(n, d){ }
        
void
HydratingConcreteMat :: initializeFrom(InputRecord &ir)
{
    // set conductivity k and capacity c
    IsotropicHeatTransferMaterial :: initializeFrom(ir);

    activationEnergy = 38400; //J/mol/K
    referenceTemperature = 25.;//C
    IR_GIVE_OPTIONAL_FIELD(ir, activationEnergy, _IFT_HydratingConcreteMat_activationEnergy);

    IR_GIVE_FIELD(ir, hydrationModelType, _IFT_HydratingConcreteMat_hydrationModelType);

    /*hydrationModelType==1: exponential hydration model, summarized in A.K. Schindler and K.J. Folliard: Heat of Hydration Models for Cementitious
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
    } else if ( hydrationModelType == 3 ) { //Saeed Rahimi-Aghdam, Zdeněk P. Bažant, Gianluca Cusatis: Extended Microprestress-Solidification Theory (XMPS) for Long-Term Creep and Diffusion Size Effect in Concrete at Variable Environment, JEM-ASCE, 2019. Appendix A.
        referenceTemperature = 20.;//according to the authors
        IR_GIVE_FIELD(ir, wc, _IFT_HydratingConcreteMat_wc);
        IR_GIVE_FIELD(ir, ac, _IFT_HydratingConcreteMat_ac);
        rhoCem = 3150.; //kg/m3
        IR_GIVE_OPTIONAL_FIELD(ir, rhoCem, _IFT_HydratingConcreteMat_rhoCem);
        rhoAgg = 1600.; //kg/m3
        IR_GIVE_OPTIONAL_FIELD(ir, rhoAgg, _IFT_HydratingConcreteMat_rhoAgg);
        Vc0 = rhoAgg*1000./(rhoAgg*1000.+rhoCem*1000.*ac+rhoCem*rhoAgg*wc);
        Vw0 = rhoAgg*rhoCem*wc/(rhoAgg*1000.+rhoCem*1000.*ac+rhoCem*rhoAgg*wc);
        
        double Blaine=350.;//m2/kg
        IR_GIVE_OPTIONAL_FIELD(ir, Blaine, _IFT_HydratingConcreteMat_Blaine);
        this->a0 = 6.5e-6*350/Blaine;
        this->ng = Vc0/(4./3*M_PI*a0*a0*a0);
        
        double alphaSet0 = 0.05;
        IR_GIVE_OPTIONAL_FIELD(ir, alphaSet0, _IFT_HydratingConcreteMat_alphaSet0);
        
        timeSet = 3*3600; //s
        IR_GIVE_OPTIONAL_FIELD(ir, timeSet, _IFT_HydratingConcreteMat_timeSet);
        
        
        double alphaCrit0 = 0.36;
        IR_GIVE_OPTIONAL_FIELD(ir, alphaCrit0, _IFT_HydratingConcreteMat_alphaCrit0);
        
        this->B0 = 1.1e-11/3600./24.; //m2/s
        IR_GIVE_OPTIONAL_FIELD(ir, this->B0, _IFT_HydratingConcreteMat_B0);
        
        double f1 = 6.5e-6/a0;
        double f2 = 1+2.5*(wc-0.3);
        double f3 = exp( this->activationEnergy / 8.314 * ( 1. / ( 273.15 + 20. ) - 1. / ( 273.15 + this->referenceTemperature ) ) );
        alphaCrit = alphaCrit0 * f1 * f2 * f3;
        if (alphaCrit>0.65) { alphaCrit = 0.65; }
        alphaSet = alphaSet0 * alphaCrit / alphaCrit0;
        
        VCemSet = (1-alphaSet)*Vc0;
        VCHSet = 0.59*alphaSet*Vc0; //CH vol. per unit vol. of cement 
        VGelSet = 1.52*alphaSet*Vc0; //Gel vol. per unit vol. of cement 
        this->aSet = pow( VCemSet/(4./3.*M_PI*ng), 1./3.);
        this->zSet = pow( (VCemSet+VGelSet)/(4./3.*M_PI*ng), 1./3.);
    } else {
        throw ValueInputException(ir, _IFT_HydratingConcreteMat_hydrationModelType, "Unknown hdyration model");
    }

    IR_GIVE_FIELD(ir, Qpot, _IFT_HydratingConcreteMat_qpot); // [1/s]

    IR_GIVE_FIELD(ir, massCement, _IFT_HydratingConcreteMat_massCement);

    maxModelIntegrationTime = 36000;
    IR_GIVE_OPTIONAL_FIELD(ir, maxModelIntegrationTime, _IFT_HydratingConcreteMat_maxModelIntegrationTime);

    minModelTimeStepIntegrations = 30.;
    IR_GIVE_OPTIONAL_FIELD(ir, minModelTimeStepIntegrations, _IFT_HydratingConcreteMat_minModelTimeStepIntegrations);
    
    IR_GIVE_OPTIONAL_FIELD(ir, referenceTemperature, _IFT_HydratingConcreteMat_referenceTemperature);
    
    conductivityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, conductivityType, _IFT_HydratingConcreteMat_conductivitytype);
    capacityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, capacityType, _IFT_HydratingConcreteMat_capacitytype);
    densityType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, densityType, _IFT_HydratingConcreteMat_densitytype);

    reinforcementDegree = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, reinforcementDegree, _IFT_HydratingConcreteMat_reinforcementDegree);
    
    timeToSeconds = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, timeToSeconds, _IFT_HydratingConcreteMat_timeToSeconds);
    
}

// returns hydration power [W/m3 of concrete]
void
HydratingConcreteMat :: computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
{
    val.resize(1);
    if ( mode == VM_Total || mode == VM_TotalIntrinsic ) {
        val.at(1) = this->GivePower(tStep, gp, mode);
    } else {
        OOFEM_ERROR("Undefined mode %s\n", __ValueModeTypeToString(mode) );
    }
}


double
HydratingConcreteMat :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    if ( mode == Capacity ) {
        return ( giveConcreteCapacity(gp, tStep) * giveConcreteDensity(gp, tStep) );
    } else if ( mode == IntSource ) { //for nonlinear solver, return dHeat/dTemperature
        HydratingConcreteMatStatus *ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
        //it suffices to compute derivative of scaling Arrhenius equation with respect to temporary temperature
        double stateVec = ms->giveField() + 273.15;
        double tempStateVec = ms->giveTempField() + 273.15;
        return this->activationEnergy / ( 8.314 * tempStateVec * tempStateVec ) * exp(1. / stateVec -  1. / tempStateVec);
    } else {
        OOFEM_ERROR("unknown mode (%s)\n", __MatResponseModeToString(mode) );
    }

    return 0.;
}


double HydratingConcreteMat :: giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
    double conduct;

    if ( conductivityType == 0 ) { //given from input file
        conduct = IsotropicHeatTransferMaterial :: give('k', gp, tStep);
    } else if ( conductivityType == 1 ) { //compute according to Ruiz, Schindler, Rasmussen. Kim, Chang: Concrete temperature modeling and strength prediction using maturity concepts in the FHWA HIPERPAV software, 7th international conference on concrete pavements, Orlando (FL), USA, 2001
        conduct = IsotropicHeatTransferMaterial :: give('k', gp, tStep) * ( 1.0 - 0.33 / 1.33 * ms->giveDoHActual() );
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
double HydratingConcreteMat :: giveConcreteCapacity(GaussPoint *gp, TimeStep *tStep) const
{
    double capacityConcrete;

    if ( capacityType == 0 ) { //given from OOFEM input file
        capacityConcrete = IsotropicHeatTransferMaterial :: give('c', gp, tStep);
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


double HydratingConcreteMat :: giveConcreteDensity(GaussPoint *gp, TimeStep *tStep) const
{
    double concreteBulkDensity;

    if ( densityType == 0 ) { //get from input file
        concreteBulkDensity = IsotropicHeatTransferMaterial :: give('d', gp, tStep);
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
        auto status = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
        answer.resize(1);
        answer.at(1) = status->giveDoHActual();
        //else answer.at(1) = 0;
        return 1;
    } else if (type == IST_EquivalentTime) {
        auto status = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
        answer.resize(1);
        answer.at(1) = status->equivalentTime;
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


MaterialStatus *
HydratingConcreteMat :: CreateStatus(GaussPoint *gp) const
{
    return new HydratingConcreteMatStatus(gp);
}


HydratingConcreteMatStatus :: HydratingConcreteMatStatus(GaussPoint *g) : TransportMaterialStatus(g) { }


//linear solver (NonStationaryTransportProblem) IntrinsicTime = TargetTime
//nonlinear solver (NLTransientTransportProblem) IntrinsicTime depends on alpha
double HydratingConcreteMat :: GivePower(TimeStep *tStep, GaussPoint *gp, ValueModeType mode) const
{
    auto ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
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
        //determine dTime for integration
        double alphaTrialOld, alphaTrialNew = 0.0;
        double time = ms->lastEvalTime;
        double dTime = ( evalTime - time ) / this->minModelTimeStepIntegrations;
        if ( dTime > this->maxModelIntegrationTime ) {
            dTime = this->maxModelIntegrationTime;
        }
        ms->degreeOfHydration = ms->lastDegreeOfHydration;
        //integration loop through hydration model at a given TimeStep
        while ( time < evalTime ) {
            if ( time + dTime > evalTime ) {
                dTime = evalTime - time;
            } else {
                time += dTime;
            }
            //printf("%f %f %f %f\n", time, affinity, scaleTemperature(), degreeOfHydration);
            alphaTrialOld = ms->degreeOfHydration + scaleTemperature(gp) * affinity25(ms->degreeOfHydration) * dTime; //predictor
            //http://en.wikipedia.org/wiki/Predictor%E2%80%93corrector_method
            //corrector - integration through trapezoidal rule
            //3 loops normally suffices
            for ( int i = 0; i < 4; i++ ) {
                alphaTrialNew = ms->degreeOfHydration + scaleTemperature(gp) * dTime / 2. * ( affinity25(ms->degreeOfHydration) + affinity25(alphaTrialOld) );
                alphaTrialOld = alphaTrialNew;
            }
            ms->degreeOfHydration = alphaTrialNew;
        }
    
    } else if ( this->hydrationModelType == 3 ) { //Rahimi-Aghdam's model, still unfinished  
        double RH=0.99;//ToDo - relative humidity in pore checking from registered fields
        double hStar = 0.88;
        double cf = 0.;
        double nh = 8.;
        
        double time = ms->lastEvalTime;
        double dTime = ( evalTime - time ) / this->minModelTimeStepIntegrations;
        
        if ( dTime > this->maxModelIntegrationTime ) {
            dTime = this->maxModelIntegrationTime;
        }
        
        ms->degreeOfHydration = ms->lastDegreeOfHydration;
        ms->zShell = ms->lastZShell;
        ms->aCement = ms->lastACement;
        ms->VCem = ms->lastVCem;
        ms->VGel = ms->lastVGel;
        ms->VCH = ms->lastVCH;
    
        while ( time < evalTime ) {
            if ( time + dTime > evalTime ) {
                dTime = evalTime - time;
            } else {
                time += dTime;
            }
        
            if ( ms->degreeOfHydration < alphaSet ){//treat first dormant period
                //assume constant rate    
                ms->degreeOfHydration += scaleTemperature(gp) * (alphaSet/timeSet ) * dTime;
                ms->VCem = (1.-ms->degreeOfHydration)*Vc0;
                ms->VCH = VCHSet * ms->degreeOfHydration/alphaSet;
                ms->VGel = VGelSet * ms->degreeOfHydration/alphaSet;
            } else {
                double f0, f4;
                
                if (ms->zShell==0.){//setting time
                    ms->zShell = this->zSet;
                    ms->aCement = this->aSet;
                    ms->VCem=VCemSet;
                    ms->VGel=VGelSet;
                    ms->VCH=VCHSet;
                }
                
                f0 = cf + (1.-cf)/(1. + pow((1.-RH)/(1.-hStar) , nh) );
                if ( ms->degreeOfHydration > 0.75*alphaCrit ) {
                    double beta =  ms->degreeOfHydration - 0.75*alphaCrit + 0.75*alphaCrit*0.30/(alphaCrit/2.); 
                    double p = pow(beta/0.30, 1.8);
                    f4 = p * exp(-p);//error in the original article?
                } else {
                    double gamma = pow(ms->degreeOfHydration/(alphaCrit/2.), 1.8);
                    f4 = gamma * exp (-gamma);
                }
                
                double Beff = B0 * f0 * f4;
                // radius of the equivalent contact-free C-S-H shells
                double zShellWedge = ms->zShell / (1.+pow((ms->zShell-this->a0)/(this->a0/6.4),5.));
                double alphaU = 0.46 + 0.95*pow(wc-0.17,0.6);
                alphaU = min(alphaU,1.0);
//                 double hc= 0.77 + 0.22*pow(wc-0.17,0.5) + 0.15*(alphaU/ms->degreeOfHydration-1.);
//                 hc = min(hc,0.99); //wrong numbers
                double hc=0.78;
                
                double Qt1 = 4*M_PI*ms->aCement*ms->zShell*Beff*(RH-hc)/(ms->zShell - ms->aCement) * zShellWedge*zShellWedge / (ms->zShell*ms->zShell);
                
                double xsi_gc=1.21, xsi_CHc=0.59, xsi_wc=(1.21+1.13)/2.;//assuming half
                double xsi_cw=1/xsi_wc;
                double xsi_gw=xsi_gc*xsi_cw;
                double xsi_CHw = xsi_CHc*xsi_cw;
                
                double dVCem = -ng*Qt1*xsi_cw * dTime;
                ms->VCem += dVCem;
                double dVGel = ng*Qt1*xsi_gw * dTime;
                ms->VGel += dVGel;
                ms->VCH += ng*Qt1*xsi_CHw * dTime;
                
                ms->aCement += 1./(4.*M_PI*a0*a0*ng)*dVCem;
                ms->degreeOfHydration -= 3./(4.*M_PI*a0*a0*a0*ng)*dVCem;
                
//                 if (ms->degreeOfHydration > alphaCrit) {
                ms->zShell += (dVGel+dVCem)/(4.*M_PI*ng*zShellWedge*zShellWedge);//error in article, missing ng?
//                 }
                //ToDo relative humidity decrement, saturation degree
            }
        }
        
        
        
    } else {
        OOFEM_ERROR("Unknown hydration model type %d", this->hydrationModelType);
    }

    ms->power = this->Qpot * ( ms->degreeOfHydration - ms->lastDegreeOfHydration ) / ( evalTime - ms->lastEvalTime );
    ms->power *= 1000 * this->massCement; // W/m3 of concrete

    //internal variables are updated in HydratingConcreteMatStatus :: updateYourself()
    return ms->power;
}


double HydratingConcreteMat :: scaleTemperature(GaussPoint *gp) const
{
    auto ms = static_cast< HydratingConcreteMatStatus * >( this->giveStatus(gp) );
    return exp( this->activationEnergy / 8.314 * ( 1. / ( 273.15 + this->referenceTemperature ) - 1. / ( 273.15 + ms->giveTempField() ) ) );
}

double HydratingConcreteMat :: affinity25(double DoH) const
{
    double result =  this->B1 *timeToSeconds* ( this->B2 / this->DoHInf + DoH ) * ( this->DoHInf - DoH ) * exp(-this->eta * DoH / this->DoHInf);
    if ( result < 0. ) { //numerical instabilities
        return 0.;
    }

    //add slag reaction
    if ( this->P1 != 0. && DoH >= this->DoH1 ) {
        result *= 1. + this->P1 * ( DoH - this->DoH1 );
    }
    return result;
}

double HydratingConcreteMatStatus :: giveDoHActual() const
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
    this->lastZShell = this->zShell;
    this->lastACement = this->aCement;
    this->lastVCem = this->VCem;
    this->lastVGel = this->VGel;
    this->lastVCH = this->VCH; 
    
    //average from last and current temperatures, in C*hour
    if ( !tStep->isIcApply() && mat->giveCastingTime() < tStep->giveIntrinsicTime() ) {
        this->maturity += ( ( this->giveField() + this->giveTempField() ) / 2. - mat->giveMaturityT0() ) * tStep->giveTimeIncrement() / 3600.;
    }
    TransportMaterialStatus :: updateYourself(tStep);
}


void
HydratingConcreteMatStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    HydratingConcreteMat *mat = static_cast< HydratingConcreteMat * >( this->gp->giveMaterial() );
    TransportMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "   status {");
    fprintf( file, "EvaluatingTime %e  DoH %f HeatPower %f [W/m3 of concrete] Temperature %f conductivity %f  capacity %f  density %f", tStep->giveIntrinsicTime(), this->giveDoHActual(), this->power, this->giveTempField(), mat->giveIsotropicConductivity(this->gp, tStep), mat->giveConcreteCapacity(this->gp, tStep), mat->giveConcreteDensity(this->gp, tStep) );
    fprintf(file, "}\n");
}
} // end namespace oofem
