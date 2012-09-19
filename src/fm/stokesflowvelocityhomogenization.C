/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
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

#include "stokesflow.h"
#include "stokesflowvelocityhomogenization.h"
#include "primaryfield.h"

#include "deadwght.h"
#include "tr21stokes.h"

namespace oofem
{
StokesFlowVelocityHomogenization :: StokesFlowVelocityHomogenization(int i, EngngModel *_master) : StokesFlow(i, _master)
{
	areaOfDomain = -1;
	areaOfRVE = -1;
}

StokesFlowVelocityHomogenization ::  ~StokesFlowVelocityHomogenization()
{}

IRResultType StokesFlowVelocityHomogenization ::  initializeFrom(InputRecord *ir)
{
	StokesFlow :: initializeFrom(ir);

	return IRRT_OK;
}

double
StokesFlowVelocityHomogenization :: giveAreaOfDomain()
{
	int i;

	if ( areaOfDomain >= 0 ) {
		return areaOfDomain;
	}

	areaOfDomain = 0;

	for ( i = 1; i <= this->giveDomain(1)->giveNumberOfElements(); i++ ) {
		areaOfDomain = areaOfDomain + giveDomain(1)->giveElement(i)->computeArea();
	}

	return areaOfDomain;
}

double
StokesFlowVelocityHomogenization :: giveAreaOfRVE()
{
	double xmin = 0, xmax = 0, ymin = 0, ymax = 0, x, y;

	if ( areaOfRVE > 0 ) {
		return areaOfRVE;
	}

	for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
		x = this->giveDomain(1)->giveDofManager(i)->giveCoordinate(1);
		y = this->giveDomain(1)->giveDofManager(i)->giveCoordinate(1);
		if ( x > xmax ) {
			xmax = x;
		}

		if ( x < xmin ) {
			xmin = x;
		}

		if ( y > ymax ) {
			ymax = y;
		}

		if ( y < ymin ) {
			ymin = y;
		}
	}

	areaOfRVE = ( xmax - xmin ) * ( ymax - ymin );
	return areaOfRVE;
}

void
StokesFlowVelocityHomogenization :: handlePrescribedValues()
{
	this->giveAreaOfDomain();

	for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfElements(); i++ ) {
		//		if (CarlTr *ThisElement = dynamic_cast<CarlTr *> ( this->giveDomain(1)->giveElement(i) ) ) {
		//
		//			ThisElement = (CarlTr *)this->giveDomain(1)->giveElement(i);
		//			ThisElement->numberOfElementsOnDomain = this->giveDomain(1)->giveNumberOfElements();
		//			ThisElement->totalAreaOfDomain = this->giveAreaOfDomain();
		//			ThisElement->specialUnknowns = & (this->SpecialUnknowns);
		//			DofMan = (Node *) ThisElement->giveDofManager(7);
		//			this->prescribedType = ThisElement->prescribedType;
		//
		//		}
	}
}

void
StokesFlowVelocityHomogenization :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
	DofIDItem type = iDof->giveDofID();
	if ( ( type == V_u ) ) {
		iDof->printSingleOutputAt(stream, atTime, 'u', EID_MomentumBalance, VM_Total, 1);
	} else if ( ( type == V_v ) ) {
		iDof->printSingleOutputAt(stream, atTime, 'v', EID_MomentumBalance, VM_Total, 1);
	} else if ( ( type == P_f ) ) {
		iDof->printSingleOutputAt(stream, atTime, 'p', EID_ConservationEquation, VM_Total, 1);
	} else {
		_error("printDofOutputAt: unsupported dof type");
	}
}


void
StokesFlowVelocityHomogenization :: solveYourselfAt(TimeStep *tStep)
{
	handlePrescribedValues();
	currentStep = tStep;
	StokesFlow :: solveYourselfAt(tStep);
}


void
StokesFlowVelocityHomogenization :: rveSetBoundaryConditions(int BCType, FloatArray eps)
{}

void
StokesFlowVelocityHomogenization :: getMeans(FloatArray &gradP, FloatArray &v, TimeStep *atTime)
{
	FloatMatrix gradPTemp, v_hatTemp;
	double Area = 0, AreaFull = 0; //(xmax-xmin)*(ymax-ymin);
	double xmax = 0, xmin = 0, ymax = 0, ymin = 0;

	gradP.resize(2);
	gradP.zero();

	v.resize(2);
	v.zero();

	for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfElements(); i++ ) {
		if ( Tr21Stokes * T = dynamic_cast< Tr21Stokes * >( this->giveDomain(1)->giveElement(i) ) ) {
			// The following only works for a rectangle

			for ( int j = 1; j <= T->giveNumberOfDofManagers(); j++ ) {
				if ( T->giveDofManager(j)->giveCoordinate(1) > xmax ) {
					xmax = T->giveDofManager(j)->giveCoordinate(1);
				}

				if ( T->giveDofManager(j)->giveCoordinate(1) < xmin ) {
					xmin = T->giveDofManager(j)->giveCoordinate(1);
				}

				if ( T->giveDofManager(j)->giveCoordinate(2) > ymax ) {
					ymax = T->giveDofManager(j)->giveCoordinate(2);
				}

				if ( T->giveDofManager(j)->giveCoordinate(2) < ymin ) {
					ymin = T->giveDofManager(j)->giveCoordinate(2);
				}
			}

			T->giveGradP(gradPTemp, atTime);
			T->giveIntegratedVelocity(v_hatTemp, atTime);

			gradP.at(1) = gradP.at(1) + gradPTemp.at(1, 1);
			gradP.at(2) = gradP.at(2) + gradPTemp.at(2, 1);

			v.at(1) = v.at(1) + v_hatTemp.at(1, 1);
			v.at(2) = v.at(2) + v_hatTemp.at(2, 1);

			Area = Area + T->computeArea();
		}
	}

	AreaFull = ( xmax - xmin ) * ( ymax - ymin );
	gradP.times(1 / Area);
	v.times(1 / AreaFull);
}

void
StokesFlowVelocityHomogenization :: updateC()
{
	OOFEM_LOG_ERROR("Uses StokesFlowVelocityHomogenization :: updateC()");
}

void
StokesFlowVelocityHomogenization :: rveGiveCharacteristicData(int DataType, void *input, void *answer, TimeStep *atTime)
{
	/*
	 * Datatype:
	 *  1 : Get velocity. *answer is a pointer to a FloatArray.
	 *  2 : Get tangent matrix. Linearization using the last tangent matrix using 1) but with other gradP
	 */

	switch ( DataType ) {
	case 1: {
		FloatArray *gradP;
		FloatArray thisGradP, v_hat;
		gradP = ( FloatArray * ) input;

		for ( int i = 1; i < this->giveDomain(1)->giveNumberOfElements(); i++ ) {
			if ( Tr21Stokes * T = dynamic_cast< Tr21Stokes * >( this->giveDomain(1)->giveElement(i) ) ) {

				IntArray *bodyLoad = T->giveBodyLoadArray();
				DeadWeight *load;
				load = dynamic_cast< DeadWeight * >( this->giveDomain(1)->giveLoad( bodyLoad->at(1) ) );

				FloatArray Components;
				Components.resize(2);
				Components.at(1) = gradP->at(1) * -1;
				Components.at(2) = gradP->at(2) * -1;
				load->setDeadWeighComponents(Components);

				break;
			}
		}

		solveYourselfAt(atTime);

		getMeans(thisGradP, v_hat, atTime);

		( ( FloatArray * ) answer )->resize(2);
		( ( FloatArray * ) answer )->at(1) = v_hat.at(1);
		( ( FloatArray * ) answer )->at(2) = v_hat.at(2);

		break;
	}
	case 2:	{
		FloatMatrix K;

		this->computeTangent(K, atTime);

		( ( FloatMatrix * ) answer )->resize(2, 2);

		( ( FloatMatrix * ) answer )->at(1, 1) = K.at(1, 1);
		( ( FloatMatrix * ) answer )->at(1, 2) = K.at(1, 2);
		( ( FloatMatrix * ) answer )->at(2, 1) = K.at(2, 1);
		( ( FloatMatrix * ) answer )->at(2, 2) = K.at(2, 2);

		return;
	}
	}
}
void
StokesFlowVelocityHomogenization :: computeTangent(FloatMatrix &answer, TimeStep *atTime)
{
	int ndof = 0;

	IntArray loc, col;
	FloatArray averagev;

	Domain *domain = this->giveDomain(1);
	ndof = this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);

	// Build F matrix
	FloatMatrix F, Fe;

	F.resize(ndof, 2);
	F.zero();
	col.resize(2);
	col.at(1) = 1;
	col.at(2) = 2;

	for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
		if ( Tr21Stokes * T = dynamic_cast< Tr21Stokes * >( this->giveDomain(1)->giveElement(i) ) ) {
			T->giveElementFMatrix(Fe);
			T->giveLocationArray( loc, EID_MomentumBalance_ConservationEquation, EModelDefaultEquationNumbering() );
			Fe.resizeWithData(15, 2);

			F.assemble(Fe, loc, col);
		}
	}

	FloatMatrix H, L;

	SparseLinearSystemNM *linMethod = CreateUsrDefSparseLinSolver(ST_Petsc, 1, this->giveDomain(1), this);

	FloatArray Hcol, Fcol;

	H.resize( F.giveNumberOfRows(), F.giveNumberOfColumns() );
	H.zero();
	Hcol.resize( F.giveNumberOfRows() );
	Hcol.zero();
	Fcol.resize( F.giveNumberOfRows() );
	Fcol.zero();

	for ( int i = 1; i <= F.giveNumberOfColumns(); i++ ) {
		for ( int j = 1; j <= F.giveNumberOfRows(); j++ ) {
			Fcol.at(j) = F.at(j, i);
		}

		linMethod->solve(stiffnessMatrix, & Fcol, & Hcol);

		for ( int j = 1; j <= H.giveNumberOfRows(); j++ ) {
			H.at(j, i) = Hcol.at(j);
		}
	}

	answer.beTProductOf(H, F);
	answer.times( 1 / this->giveAreaOfRVE() );

	delete(linMethod);
}

#ifdef __PETSC_MODULE
void
StokesFlowVelocityHomogenization :: initPetscContexts()
{
	PetscContext *petscContext;

	int i;
	petscContextList->growTo(ndomains);
	for ( i = 0; i < this->ndomains; i++ ) {
		petscContext =  new PetscContext(this, EID_MomentumBalance_ConservationEquation);
		petscContextList->put(i + 1, petscContext);
	}
}
#endif

double
StokesFlowVelocityHomogenization :: giveUnknownComponent(EquationID chc, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
	if ( ( chc == EID_ConservationEquation ) || ( chc == EID_MomentumBalance ) || ( chc == EID_MomentumBalance_ConservationEquation ) ) {
		return velocityPressureField->giveUnknownValue(dof, mode, tStep);
	} else {
		OOFEM_ERROR("giveUnknownComponent: Unknown is of undefined equation id for this problem");
	}

	return 0;
}
}
