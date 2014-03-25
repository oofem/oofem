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

#include "linearfracturemechanics.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"
#include "errorestimator.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "element.h"
#include "gausspoint.h"


namespace oofem {
REGISTER_EngngModel(LinearFractureMechanics);

void
LinearFractureMechanics :: updateYourself(TimeStep *tStep)
{
    LinearStatic :: updateYourself(tStep);
    // evaluate fracture criterion
	fprintf(outputStream,"\n ::::kontrola 00 zacatek");
	bool remeshing = evaluateFractureCriterion(fractureCriterionType, tStep);
	fprintf(outputStream,"\n ::::kontrola 01 konec");

	/*
	//something like
	if( remeshing) {
		MesherInterface *mesher = classFactory.createMesherInterface( meshPackage, this->giveDomain(1) );
        Domain *newDomain;
        MesherInterface :: returnCode result = mesher->createMesh(tStep, 1, this->giveDomain(1)->giveSerialNumber() + 1, & newDomain);
	}
	*/
	

    
}

void
LinearFractureMechanics :: terminate(TimeStep *tStep)
{
    LinearStatic :: terminate(tStep);
    //
    // print estimated error
    //
    fprintf(outputStream, "\nRelative error estimate: %5.2f%%\n", this->defaultErrEstimator->giveValue(relativeErrorEstimateEEV, tStep) * 100.0);
}



contextIOResultType LinearFractureMechanics :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = LinearStatic :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

IRResultType
LinearFractureMechanics :: initializeFrom(InputRecord *ir)
// input from inputString
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LinearStatic :: initializeFrom(ir);

	int meshPackageId = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, meshPackageId, _IFT_LinearFractureMechanics_Meshpackage);
	// This should be changed, probably only one mesh package could be used within this class
	
    /*

    if ( meshPackageId == 1 ) {
        meshPackage = MPT_TARGE2;
    } else if ( meshPackageId == 2 ) {
        meshPackage = MPT_FREEM;
    } else if ( meshPackageId == 4 ) {
		meshPackage = MPT_VIENNAGRID;
	} else {
        meshPackage = MPT_T3D;
    }
	*/

	int val = 0;
    IR_GIVE_FIELD(ir, val, _IFT_LinearFractureMechanics_FractureCriterion);
	fractureCriterionType = ( FractureCriterion ) val;


	IR_GIVE_FIELD(ir, crackTipCoordinates, _IFT_LinearFractureMechanics_CrackTipCoordinates);



	if(fractureCriterionType == FC_SIF) {
		IR_GIVE_FIELD(ir, fractureTougness, _IFT_LinearFractureMechanics_FractureThougness);
		IR_GIVE_FIELD(ir, rmax, _IFT_LinearFractureMechanics_rmax);

	} else if(fractureCriterionType = FC_JIntegral) {
		IR_GIVE_FIELD(ir, fractureEnergy, _IFT_LinearFractureMechanics_FractureEnergy);	
	}

    return IRRT_OK;
}


bool LinearFractureMechanics :: evaluateFractureCriterion(FractureCriterion fractureCriterionType, TimeStep *tStep)
{
	if ( fractureCriterionType == FC_SIF ) {
		double SIF;
		this->giveStressIntensityFactor(SIF, tStep);
		if( SIF > fractureTougness )
			return 1;
		else 
			return 0;
	} else if( fractureCriterionType == FC_JIntegral) {
		return 1;
	} else {
		return 0;
	}
}


void
LinearFractureMechanics :: giveStressIntensityFactor(double &answer, TimeStep *tStep)
{
	    int nelems = (this->giveDomain(1))->giveNumberOfElements();
		fprintf(outputStream,"\n ::::kontrola - nelems=%d",nelems);
		for (int i = 1; i <= nelems; i++)
		{
			fprintf(outputStream,"\n ::::kontrola - prvek c.: %d\n",i);
			Element *elem = this->giveDomain(1)->giveElement(i);
			FloatArray elementSIF;
			this->giveElementStressIntensityFactor(elementSIF,elem, tStep);
		}
		answer = 0;
}


void
LinearFractureMechanics :: giveElementStressIntensityFactor(FloatArray &answer, Element *elem, TimeStep *tStep)
{
	fprintf(outputStream,"\n   ::::vlastnosti prvku");
    FloatArray val;
	answer.resize(3);
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();
    int result = 1, nip = iRule->giveNumberOfIntegrationPoints();	   
	FloatArray *gaussPointCoordinates;
    for ( int i = 0; i < nip; i++ ) {
		gaussPointCoordinates = iRule -> getIntegrationPoint(i)->giveCoordinates();
		double x = gaussPointCoordinates->at(1);
		double y = gaussPointCoordinates->at(2);
		double crackTipDistance, crackTipAngle;
		fprintf(outputStream,"\n   :::: souradnice x: %g  y: %g",x,y);
		crackTipDistance = sqrt((x-crackTipCoordinates.at(1))*(x-crackTipCoordinates.at(1)) + (y-crackTipCoordinates.at(2))*(y-crackTipCoordinates.at(2)));
		double x0 = crackTipCoordinates.at(1);
		double y0 = crackTipCoordinates.at(2);
		fprintf(outputStream,"\n   :::: prvek x0: %g  y0: %g",x0,y0);
		if ( (x<x0) & (y>y0) ) { crackTipAngle= asin(abs(x-x0)/crackTipDistance);}
		if ( (x<x0) & (y<y0) ) { crackTipAngle=PI - asin(abs(x-x0)/crackTipDistance);}
		if ( (x>x0) & (y<y0) ) { crackTipAngle= -(PI - asin(abs(x-x0)/crackTipDistance));}
		if ( (x>x0) & (y>y0) ) { crackTipAngle=-(asin(abs(x-x0)/crackTipDistance));}
		if (x=x0) { crackTipAngle= 0;}
		if (y=y0) { if (x<x0) {crackTipAngle=PI/2;} else crackTipAngle=3*PI/2;}
		if( crackTipDistance*crackTipDistance < rmax*rmax) {
			result = elem->giveIPValue(val, iRule->getIntegrationPoint(i), IST_StressTensor, tStep);   
			double sigmaX = val.at(1);
			double sigmaY = val.at(2); // ?
			double tauXY = val.at(6); // ?
			fprintf(outputStream, "\n   :::: napeti \n 1: %f\n 2: %f\n 3: %f\n", val.at(1),val.at(2),val.at(6));
			answer.at(1) = sigmaX*sqrt(2*PI*crackTipDistance)/(1-sin(x/2)*sin(3*x/2));
			answer.at(2) = sigmaY*sqrt(2*PI*crackTipDistance)/(1+sin(x/2)*sin(3*x/2));
			if ( (sin(x/2)*cos(x/2)*cos(3*x/2)) != 0) {
				answer.at(3) = tauXY*sqrt(2*PI*crackTipDistance)/(sin(x/2)*cos(x/2)*cos(3*x/2));
			}
			//fprintf(outputStream, "\n:: VYSTUP: %5.2f\n", answer.at(1));
		}
	}

}




void
LinearFractureMechanics :: updateDomainLinks()
{
    LinearStatic :: updateDomainLinks();
}
} // end namespace oofem
