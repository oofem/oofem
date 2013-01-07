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

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "activebc.h"
#include "weakperiodicbc.h"
#include "inputrecord.h"
#include "element.h"
#include "elementside.h"
#include "node.h"
#include "masterdof.h"
#include "sparsemtrx.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"

namespace oofem {

WeakPeriodicbc :: WeakPeriodicbc (int n, Domain *d) : ActiveBoundaryCondition (n, d)
{
	useBasisType = trigonometric;
	doUpdateSminmax = true;
}

IRResultType
WeakPeriodicbc :: initializeFrom(InputRecord *ir)
{

	const char *__proc = "initializeFrom";
	IRResultType result;

	ActiveBoundaryCondition :: initializeFrom(ir);

	orderOfPolygon=2;
	IR_GIVE_OPTIONAL_FIELD(ir, orderOfPolygon, IFT_WeakPeriodicBoundaryCondition_order, "order");

	int t= (int) trigonometric;
	IR_GIVE_OPTIONAL_FIELD(ir, t, IFT_WeakPeriodicBoundaryCondition_descritization, "descritizationtype");
	useBasisType= (basisType) t; // Fourierseries by default

	dofid=11;  // Pressure as default
	IR_GIVE_OPTIONAL_FIELD(ir, dofid, IFT_WeakPeriodicBoundaryCondition_order, "dofid");

	IntArray temp;
	IR_GIVE_OPTIONAL_FIELD(ir, temp, IFT_ActiveBoundaryCondition_elementSides, "elementsidespositive");
	for (int i=0; i<temp.giveSize()/2; i++) {
		//		printf("Add positive edge, element %u, side %u\n", temp.at(2*i+2), temp.at(2*i+1));
		side[0].push_back(temp.at(2*i+1));
		element[0].push_back(temp.at(2*i+2));
	}

	IR_GIVE_OPTIONAL_FIELD(ir, temp, IFT_ActiveBoundaryCondition_elementSides, "elementsidesnegative");
	for (int i=0; i<temp.giveSize()/2; i++) {
		//		printf("Add negative edge, element %u, side %u\n", temp.at(2*i+2), temp.at(2*i+1));
		side[1].push_back(temp.at(2*i+1));
		element[1].push_back(temp.at(2*i+2));
	}

	// Create dofs for coefficients
	bcID=this->giveNumber();
	gammaDman = new Node(100, this->domain);

	DofIDList.resize(orderOfPolygon+1);
	for (int i=0; i<=orderOfPolygon; i++) {
		gammaDman->appendDof(new MasterDof(i, gammaDman, (DofIDItem) (bcID*100+i)));
		DofIDList.at(i+1)=(bcID*100+i);
	}

	return IRRT_OK;
}

void WeakPeriodicbc :: giveEdgeNormal(FloatArray &answer, int element, int side)
{
	FloatArray Tangent;

	answer.resize(2);
	Tangent.resize(2);

	Element *thisElement;
	thisElement=this->domain->giveElement(element);

	if (side==1) {
		Tangent.at(1)=thisElement->giveNode(2)->giveCoordinate(1) - thisElement->giveNode(1)->giveCoordinate(1);	// x
		Tangent.at(2)=thisElement->giveNode(2)->giveCoordinate(2) - thisElement->giveNode(1)->giveCoordinate(2);	// y
	} else if (side==2) {
		Tangent.at(1)=thisElement->giveNode(3)->giveCoordinate(1) - thisElement->giveNode(2)->giveCoordinate(1);	// x
		Tangent.at(2)=thisElement->giveNode(3)->giveCoordinate(2) - thisElement->giveNode(2)->giveCoordinate(2);	// y
	} else if (side==3) {
		Tangent.at(1)=thisElement->giveNode(1)->giveCoordinate(1) - thisElement->giveNode(3)->giveCoordinate(1);	// x
		Tangent.at(2)=thisElement->giveNode(1)->giveCoordinate(2) - thisElement->giveNode(3)->giveCoordinate(2);	// y
	}

	// Normalize
	double l = sqrt(Tangent.at(1)*Tangent.at(1)+Tangent.at(2)*Tangent.at(2));
	Tangent.at(1)=Tangent.at(1)/( l );
	Tangent.at(2)=Tangent.at(2)/( l );

	answer.at(1)=Tangent.at(2);
	answer.at(2)=-Tangent.at(1);
}

void WeakPeriodicbc :: updateDirection()
{
	// Check orientation for s
	FloatArray normal;

	giveEdgeNormal(normal, element[0].at(0), side[0].at(0));
	if ( abs ( abs ( normal.at(1) ) -1) < 0.0001 ) {	// Normal point in x direction, thus set direction to y
		direction=2;
	} else {
		direction=1;
	}
}

void WeakPeriodicbc :: updateSminmax()
{

	if (doUpdateSminmax) {
		updateDirection();

		smin=this->domain->giveDofManager(1)->giveCoordinate(direction);
		smax=this->domain->giveDofManager(1)->giveCoordinate(direction);

		for (int i=1; i<=this->domain->giveNumberOfDofManagers(); i++) {
			double sValue = this->domain->giveDofManager(i)->giveCoordinate(direction);
			smin=std::min(smin, sValue);
			smax=std::max(smax,sValue);
		}
		//printf("smin=%f\tsmax=%f\n", smin, smax);
		doUpdateSminmax=false;
	}

}

void WeakPeriodicbc :: addElementSide(int newElement, int newSide)
{
	//	printf ("Add element %u, side %u\n", newElement, newSide);

	FloatArray normalNew, normal0;
	int addToList=0;

	if (element[0].size()>0) {
		// If there are elements in the list, compare normals in order to determine which list to store them in
		giveEdgeNormal(normalNew, newElement, newSide);
		//		normalNew.printYourself();
		giveEdgeNormal(normal0, element[0].at(0), side[0].at(0));
		double d=sqrt(pow(normalNew.at(1)-normal0.at(1), 2) + pow(normalNew.at(2)-normal0.at(2), 2));
		if (abs(d)<0.001) {
			addToList=0;
		} else {
			addToList=1;
		}

	} else {
		// Otherwise, check the normal in order to decide upon which direction the parameter runs (x or y)
		giveEdgeNormal(normalNew, newElement, newSide);
		//		normalNew.printYourself();
		if ( abs ( abs ( normalNew.at(1) ) -1) < 0.0001 ) {	// Normal point in x direction, thus set direction to y
			direction=2;
		} else {
			direction=1;
		}
	}
//	printf(" to list %u\n", addToList);
	element[addToList].push_back(newElement);
	side[addToList].push_back(newSide);
}

void WeakPeriodicbc :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain)
{

	GaussIntegrationRule *iRule;

	GaussPoint *gp;
	FloatArray *lcoords, gcoords, normal;
	int normalSign, dofCountOnBoundary;

	Element *thisElement;

	FEI2dTrLin interpolation_lin(1, 2);
	FEI2dTrQuad interpolation_quad(1, 2);
	FEInterpolation2d *interpolation;

	updateSminmax();

	// Assemble each side
	for (int thisSide=0; thisSide<=1; thisSide++) {
		giveEdgeNormal(normal, element[thisSide].at(0), side[thisSide].at(0));
		if ( (normal.at(1)+normal.at(2)) <= 0.0001 ) { // This is a south or west edge
			normalSign = -1;
		} else {
			normalSign = 1;
		}

		for (size_t ielement=0; ielement < element[thisSide].size(); ielement++) {	// Loop over each element on this edge

			FloatMatrix B, BT;

			thisElement = this->domain->giveElement(element[thisSide].at(ielement));

			iRule = new GaussIntegrationRule(1, thisElement, 1, 1);
			iRule->setUpIntegrationPoints(_Line, 3, _Unknown);

			// Find dofs for this element side
			IntArray tempSideLocation, sideLocation;

			thisElement->giveBoundaryLocationArray(tempSideLocation, side[thisSide].at(ielement), EID_MomentumBalance, r_s);

			// Find dofs for this element which should be periodic
			IntArray bNodes, nodeDofIDMask, periodicDofIDMask, nodalArray;
			periodicDofIDMask.resize(1);
			periodicDofIDMask.at(1)=dofid;

			thisElement->giveInterpolation()->boundaryGiveNodes(bNodes, side[thisSide].at(ielement));

			sideLocation.resize(0);
			dofCountOnBoundary=0;
			for (int i=1; i<=bNodes.giveSize(); i++) {
				thisElement->giveDofManDofIDMask(bNodes.at(i), EID_MomentumBalance_ConservationEquation, nodeDofIDMask);

				for (int j=1; j<=nodeDofIDMask.giveSize(); j++) {
					if (nodeDofIDMask.at(j)==dofid) {
						thisElement->giveDofManager(bNodes.at(i))->giveLocationArray(periodicDofIDMask, nodalArray, r_s);
						sideLocation.followedBy(nodalArray);
						dofCountOnBoundary++;
					}
				}
			}

			B.resize(dofCountOnBoundary,orderOfPolygon+1);

			// Use linear or quadratic interpolation?
			if (dofCountOnBoundary==2) {
				interpolation = &interpolation_lin;
			} else {
				interpolation = &interpolation_quad;
			}

			IntArray cloc;
			gammaDman->giveLocationArray(DofIDList, cloc, EModelDefaultEquationNumbering());

			B.zero();
			for (int i=0; i < iRule->getNumberOfIntegrationPoints(); i++) {

				gp = iRule->getIntegrationPoint(i);
				lcoords = gp->giveCoordinates();

				FloatArray N;

				// Find the value of parameter s which is the vert/horiz distance to 0
				interpolation->edgeLocal2global(gcoords, side[thisSide].at(ielement), *lcoords, FEIElementGeometryWrapper(thisElement) );
				// Compute base function values
				interpolation->boundaryEvalN(N, *lcoords, FEIElementGeometryWrapper(thisElement) );
				// Compute Jacobian
				double detJ = fabs(interpolation->edgeGiveTransformationJacobian(side[thisSide].at(ielement), * lcoords, FEIElementGeometryWrapper(thisElement)));
				double s=gcoords.at(direction);

				for (int k=0; k<dofCountOnBoundary; k++) {
					for (int j=0; j<=orderOfPolygon; j++) {
						double fVal=computeBaseFunctionValue(j, s);
						B.at(k+1, j+1)=B.at(k+1, j+1)+N.at(k+1)*fVal*detJ*normalSign*gp->giveWeight();
					}
				}
			}

			answer->assemble(sideLocation, cloc, B);
			BT.beTranspositionOf(B);
			answer->assemble(cloc, sideLocation, BT);

			delete iRule;

		}
	}
}

double WeakPeriodicbc :: computeBaseFunctionValue(int baseID, double coordinate)
{
	double fVal;
	double sideLength=smax-smin;

	if (useBasisType==monomial) {
		fVal=pow(coordinate,baseID);
	} else if(useBasisType==trigonometric) {
		if (baseID%2==0) { // Even
			fVal=cos(((double) baseID)/2*(coordinate*2*3.141593/sideLength));
		} else {
			fVal=sin(((double) baseID+1)/2*(coordinate*2*3.141593/sideLength));
		}
	}

	return fVal;
}

double WeakPeriodicbc :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
		CharType type, ValueModeType mode,
		const UnknownNumberingScheme &s, Domain *domain)
{
	double norm = 0.0;

	// Fetch unknowns of this boundary condition
	FloatArray gamma;
	gammaDman->giveUnknownVector(gamma, DofIDList, eid, mode, tStep);

	if (type==InternalForcesVector) {

		GaussIntegrationRule *iRule;

		GaussPoint *gp;
		FloatArray *lcoords, gcoords, normal;
		int normalSign, dofCountOnBoundary;

		Element *thisElement;

		FEI2dTrLin interpolation_lin(1, 2);
		FEI2dTrQuad interpolation_quad(1, 2);
		FEInterpolation2d *interpolation;

		updateSminmax();

		// Assemble each side
		for (int thisSide=0; thisSide<=1; thisSide++) {
			giveEdgeNormal(normal, element[thisSide].at(0), side[thisSide].at(0));
			if ( (normal.at(1)+normal.at(2)) <= 0.0001 ) { // This is a south or west edge
				normalSign = -1;
			} else {
				normalSign = 1;
			}

			for (size_t ielement=0; ielement < element[thisSide].size(); ielement++) {	// Loop over each element on this edge
				FloatMatrix B, BT;

				thisElement = this->domain->giveElement(element[thisSide].at(ielement));

				// Values from solution
				FloatArray a;

				iRule = new GaussIntegrationRule(1, thisElement, 1, 1);
				iRule->setUpIntegrationPoints(_Line, 3, _Unknown);

				// Find dofs for this element side
				IntArray tempSideLocation, sideLocation;

				thisElement->giveBoundaryLocationArray(tempSideLocation, side[thisSide].at(ielement), EID_MomentumBalance, s);

				// Find dofs for this element which should be periodic
				IntArray bNodes, nodeDofIDMask, periodicDofIDMask, nodalArray;
				periodicDofIDMask.resize(1);
				periodicDofIDMask.at(1)=dofid;

				thisElement->giveInterpolation()->boundaryGiveNodes(bNodes, side[thisSide].at(ielement));

				sideLocation.resize(0);
				a.resize(0);
				dofCountOnBoundary=0;
				for (int i=1; i<=bNodes.giveSize(); i++) {
					thisElement->giveDofManDofIDMask(bNodes.at(i), EID_MomentumBalance_ConservationEquation, nodeDofIDMask);

					for (int j=1; j<=nodeDofIDMask.giveSize(); j++) {
						if (nodeDofIDMask.at(j)==dofid) {
							thisElement->giveDofManager(bNodes.at(i))->giveLocationArray(periodicDofIDMask, nodalArray, s);
							double value=thisElement->giveDofManager(bNodes.at(i))->giveDof(j)->giveUnknown(eid, mode, tStep);
							sideLocation.followedBy(nodalArray);
							a.resize(sideLocation.giveSize());
							a.at(sideLocation.giveSize())=value;
							dofCountOnBoundary++;
						}
					}
				}

				B.resize(dofCountOnBoundary,orderOfPolygon+1);

				// Use linear or quadratic interpolation?
				if (dofCountOnBoundary==2) {
					interpolation = &interpolation_lin;
				} else {
					interpolation = &interpolation_quad;
				}

				IntArray cloc;
				gammaDman->giveLocationArray(DofIDList, cloc, EModelDefaultEquationNumbering());

				B.zero();
				for (int i=0; i < iRule->getNumberOfIntegrationPoints(); i++) {

					gp = iRule->getIntegrationPoint(i);
					lcoords = gp->giveCoordinates();

					FloatArray N;

					// Find the value of parameter s which is the vert/horiz distance to 0
					interpolation->edgeLocal2global(gcoords, side[thisSide].at(ielement), *lcoords, FEIElementGeometryWrapper(thisElement) );
					// Compute base function values
					interpolation->boundaryEvalN(N, *lcoords, FEIElementGeometryWrapper(thisElement) );
					// Compute Jacobian
					double detJ = fabs(interpolation->edgeGiveTransformationJacobian(side[thisSide].at(ielement), * lcoords, FEIElementGeometryWrapper(thisElement)));
					double s=gcoords.at(direction);

					for (int k=0; k<dofCountOnBoundary; k++) {
						for (int j=0; j<=orderOfPolygon; j++) {
							double fVal=computeBaseFunctionValue(j, s);
							B.at(k+1, j+1)=B.at(k+1, j+1)+N.at(k+1)*fVal*detJ*normalSign*gp->giveWeight();
						}
					}
				}

				FloatArray myProd, myProdGamma;
				myProd.beTProductOf(B, a);
				myProdGamma.beProductOf(B, gamma);

				norm += myProd.computeSquaredNorm();
				norm += myProdGamma.computeSquaredNorm();

				answer.assemble(myProd, cloc);
				answer.assemble(myProdGamma, sideLocation);

				delete iRule;
			}

		}
	}

	return norm;
}

int WeakPeriodicbc :: giveNumberOfInternalDofManagers()
{
	return 1;
}

DofManager * WeakPeriodicbc :: giveInternalDofManager(int i)
{
	if (i==1) {
		return gammaDman;
	} else {
		return NULL;
	}

}

}

