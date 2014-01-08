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

#include "latticetransmat.h"
#include "domain.h"
#include "gausspoint.h"
#include "latticetransportelement.h"
#include "math.h"
#include "staggeredproblem.h"
#include "classfactory.h"
#ifdef __SM_MODULE
 #include "latticestructuralelement.h"
#endif

namespace oofem {
REGISTER_Material(LatticeTransportMaterial);

IRResultType
LatticeTransportMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->viscosity, _IFT_LatticeTransportMaterial_vis);

    IR_GIVE_FIELD(ir, this->permeability, _IFT_LatticeTransportMaterial_k);

    IR_GIVE_FIELD(ir, this->thetaS, _IFT_LatticeTransportMaterial_thetas);

    IR_GIVE_FIELD(ir, this->thetaR, _IFT_LatticeTransportMaterial_thetar);

    conType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->conType, _IFT_LatticeTransportMaterial_contype);
    //Options are 0 = constant conductivity and capacity and 1 = van Genuchten conductivity and capactity

    if ( conType == 1 ) {
        IR_GIVE_FIELD(ir, this->paramM, _IFT_LatticeTransportMaterial_m);

        IR_GIVE_FIELD(ir, this->paramA, _IFT_LatticeTransportMaterial_a);

        //Assume that original van Genuchten is used.
        this->thetaM = this->thetaS;
        IR_GIVE_OPTIONAL_FIELD(ir, this->thetaM, _IFT_LatticeTransportMaterial_thetam);

        this->suctionAirEntry = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, this->suctionAirEntry, _IFT_LatticeTransportMaterial_paev);

        if ( thetaM < thetaS ) {
            _warning("initializeFrom : thetaM cannot be smaller than thetaS. Choose thetaM=thetaS.");
            thetaM = thetaS;
        }

        //Define value of either air entry pressure or thetaM for later use. Only based on input parameters.
        if ( suctionAirEntry == 0. ) {
            suctionAirEntry = paramA * pow( ( pow( ( thetaS - thetaR ) / ( thetaM - thetaR ), -1. / paramM ) - 1. ), ( 1 - paramM ) );
        } else {
            thetaM = ( thetaS - thetaR ) * pow(1. + pow( suctionAirEntry / paramA, 1. / ( 1. - paramM ) ), paramM) + thetaR;
        }
    } //end of contype condition
    else if ( conType != 0 && conType != 1 ) {
        _error("initializeFrom : unknown conType mode");
    }

    crackTortuosity = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->crackTortuosity, _IFT_LatticeTransportMaterial_ctor);



    return IRRT_OK;
}


void
LatticeTransportMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    ///@todo Is there anything meaningful to do here for lattice models?
}


double
LatticeTransportMaterial :: give(int aProperty, GaussPoint *gp)
{
    if ( aProperty == 'k' ) {
        return permeability;
    } else if ( ( aProperty == HeatCapaCoeff ) || ( aProperty == 'c' ) ) {
        return ( this->give('d', gp) );
    }

    return this->Material :: give(aProperty, gp);
}


double
LatticeTransportMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                    GaussPoint *gp,
                                                    TimeStep *tStep)
{
    LatticeTransportMaterialStatus *status = static_cast< LatticeTransportMaterialStatus * >( this->giveStatus(gp) );
    double suction = status->giveTempField().at(1);

    if ( mode == Capacity ) {
        return computeCapacity(suction, gp);
    } else if ( mode == Conductivity ) {
        return computeConductivity(suction, gp, tStep);
    } else {
        _error("giveCharacteristicValue : unknown mode");
    }

    return 0; // to make compiler happy
}


double
LatticeTransportMaterial :: computeConductivity(double suction,
                                                GaussPoint *gp,
                                                TimeStep *tStep)
{
    LatticeTransportMaterialStatus *status = static_cast< LatticeTransportMaterialStatus * >( this->giveStatus(gp) );

    this->density = this->give('d', gp);

    double relativePermeability;
    double conductivity;

    relativePermeability = 0.;
    conductivity = 0.;

    double crackWidth = 0.;

#ifdef __SM_MODULE
    IntArray coupledModels;
    if ( domain->giveEngngModel()->giveMasterEngngModel() ) {
        ( static_cast< StaggeredProblem * >( domain->giveEngngModel()->giveMasterEngngModel() ) )->giveCoupledModels(coupledModels);
        int couplingFlag = ( static_cast< LatticeTransportElement * >( gp->giveElement() ) )->giveCouplingFlag();

        if ( couplingFlag == 1 && coupledModels.at(1) != 0 && !tStep->isTheFirstStep() ) {
            int couplingNumber;
            couplingNumber = ( static_cast< LatticeTransportElement * >( gp->giveElement() ) )->giveCouplingNumber();
            LatticeStructuralElement *coupledElement;
            coupledElement  = static_cast< LatticeStructuralElement * >( domain->giveEngngModel()->giveMasterEngngModel()->giveSlaveProblem( coupledModels.at(1) )->giveDomain(1)->giveElement(couplingNumber) );
            crackWidth = coupledElement->giveCrackWidth();
        }
    }
#endif

    if ( !domain->giveEngngModel()->giveMasterEngngModel() ) {
        crackWidth = ( static_cast< LatticeTransportElement * >( gp->giveElement() ) )->giveCrackWidth();
    }


    double saturation, partOne, partTwo, numerator, denominator;
    if ( suction < this->suctionAirEntry || conType == 0 ) {
        saturation = 1.;
        relativePermeability = 1.;
    } else {
        partOne = pow( suction / this->paramA, 1. / ( 1. - this->paramM ) );

        saturation = ( this->thetaM - this->thetaR ) / ( this->thetaS - this->thetaR ) * pow(1. + partOne, -this->paramM);

        partTwo = ( this->thetaS - this->thetaR ) / ( this->thetaM - this->thetaR );

        numerator = ( 1. - pow(1. - pow(partTwo * saturation, 1 / this->paramM), this->paramM) );
        denominator = ( 1. - pow(1. - pow(partTwo, 1 / this->paramM), this->paramM) );

        relativePermeability = sqrt(saturation) * pow( ( numerator ) / ( denominator ), 2.0 );
    }

    //Calculate mass for postprocessing
    double mass = ( saturation * ( this->thetaS - this->thetaR ) + this->thetaR ) * this->density;

    status->setMass(mass);

    conductivity = this->permeability / this->viscosity * relativePermeability;

    //add crack contribution;
    double crackContribution = 0.;
    double width = ( static_cast< LatticeTransportElement * >( gp->giveElement() ) )->giveWidth();

    crackContribution = pow(crackWidth, 3.) / ( width * 12. * this->viscosity ) * this->crackTortuosity * relativePermeability;

    conductivity += crackContribution;

    return this->density * conductivity;
}



double
LatticeTransportMaterial :: computeCapacity(double suction, GaussPoint *gp)
{
    double capacity = 0.;

    this->density = this->give('d', gp);

    if ( suction < this->suctionAirEntry || conType == 0 ) {
        capacity = 0.;
    } else {
        double partOne = this->paramM / ( this->paramA * ( 1. - this->paramM ) );
        double partTwo = pow( suction / this->paramA, this->paramM / ( 1. - this->paramM ) );
        double partThree = pow(1. + pow( suction / this->paramA, 1. / ( 1. - this->paramM ) ), -this->paramM - 1.);
        capacity = ( this->thetaM - this->thetaR ) * partOne * partTwo * partThree;
    }

    return this->density * capacity;
}


MaterialStatus *
LatticeTransportMaterial :: CreateStatus(GaussPoint *gp) const
{
    LatticeTransportMaterialStatus *answer = new LatticeTransportMaterialStatus(1, LatticeTransportMaterial :: domain, gp);
    return answer;
}


void
LatticeTransportMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
{
    MaterialStatus :: printOutputAt(File, tNow);

    fprintf(File, "  state");

    for ( int i = 1; i <= field.giveSize(); i++ ) {
        fprintf( File, " % .4e", field.at(i) );
    }

    fprintf(File, "  mass %.8e", mass);
    fprintf(File, "\n");
}

void
LatticeTransportMaterialStatus :: updateYourself(TimeStep *tStep)
{
    TransportMaterialStatus :: updateYourself(tStep);
}


void
LatticeTransportMaterialStatus :: initTempStatus()
{
    TransportMaterialStatus :: initTempStatus();
}


LatticeTransportMaterialStatus :: LatticeTransportMaterialStatus(int n, Domain *d, GaussPoint *g) :
    TransportMaterialStatus(n, d, g)
{
    // stateVector.resize(1);
    // stateVector.zero();
    // tempStateVector = stateVector;
    mass = 0.;
}
} // end namespace oofem
