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

#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "activebc.h"
#include "weakperiodicbc.h"
#include "inputrecord.h"
#include "element.h"
#include "elementside.h"
#include "node.h"
#include "masterdof.h"
#include "sparsemtrx.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"
#include "classfactory.h"
#include "set.h"

namespace oofem {

REGISTER_BoundaryCondition( WeakPeriodicBoundaryCondition );

WeakPeriodicBoundaryCondition :: WeakPeriodicBoundaryCondition(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
	useBasisType = trigonometric;
	doUpdateSminmax = true;
}

IRResultType
WeakPeriodicBoundaryCondition :: initializeFrom(InputRecord *ir)
{
	const char *__proc = "initializeFrom";
	IRResultType result;

	ActiveBoundaryCondition :: initializeFrom(ir); ///@todo Carl, remove this line and use elementsidespositive/negative instead.

	orderOfPolygon = 2;
	IR_GIVE_OPTIONAL_FIELD(ir, orderOfPolygon, _IFT_WeakPeriodicBoundaryCondition_order);

	int t = ( int ) trigonometric;
	IR_GIVE_OPTIONAL_FIELD(ir, t, _IFT_WeakPeriodicBoundaryCondition_descritizationType);
	useBasisType = ( basisType ) t;  // Fourierseries by default

	dofid = P_f;    // Pressure as default
	IR_GIVE_OPTIONAL_FIELD(ir, dofid, _IFT_WeakPeriodicBoundaryCondition_dofid);

	ngp = -1;    // Pressure as default
	IR_GIVE_OPTIONAL_FIELD(ir, ngp, _IFT_WeakPeriodicBoundaryCondition_ngp);

	IntArray temp;

	posSet=-1;
	negSet=-1;

	IR_GIVE_OPTIONAL_FIELD(ir, posSet, _IFT_WeakPeriodicBoundaryCondition_elementSidesPositiveSet);
	IR_GIVE_OPTIONAL_FIELD(ir, negSet, _IFT_WeakPeriodicBoundaryCondition_elementSidesNegativeSet);

	if (posSet==-1) {
		IR_GIVE_OPTIONAL_FIELD(ir, temp, _IFT_WeakPeriodicBoundaryCondition_elementSidesPositive);
		for ( int i = 0; i < temp.giveSize() / 2; i++ ) {
			side [ 0 ].push_back( temp.at(2 * i + 1) );
			element [ 0 ].push_back( temp.at(2 * i + 2) );
		}
	}

	if (negSet==-1) {
		IR_GIVE_OPTIONAL_FIELD(ir, temp, _IFT_WeakPeriodicBoundaryCondition_elementSidesNegative);
		for ( int i = 0; i < temp.giveSize() / 2; i++ ) {
			side [ 1 ].push_back( temp.at(2 * i + 1) );
			element [ 1 ].push_back( temp.at(2 * i + 2) );
		}
	}

	if (this->domain->giveNumberOfSpatialDimensions()==2)
		ndof=orderOfPolygon+1;
	else if (this->domain->giveNumberOfSpatialDimensions()==3) {
		ndof=1;
		for (int i=1; i<=orderOfPolygon; i++) ndof=ndof + (i+1);
	}

	// Create dofs for coefficients
	bcID = this->giveNumber();
	gammaDman = new Node(0, this->domain);

	for ( int i = 0; i < ndof; i++ ) {
		gammaDman->appendDof( new MasterDof( i, gammaDman, (DofIDItem)this->domain->giveNextFreeDofID() ) );
	}

	return IRRT_OK;
}

void WeakPeriodicBoundaryCondition :: giveEdgeNormal(FloatArray &answer, int element, int side)
{

	FloatArray xi;

	if (this->domain->giveNumberOfSpatialDimensions()==3) {
		xi.resize(2);
		xi(0)=0.25;
		xi(1)=0.25;
	} else {
		xi.resize(1);
		xi(0)=0.5;
	}

	Element *thisElement = this->domain->giveElement( element );
	FEInterpolation *interpolation = thisElement->giveInterpolation( (DofIDItem)dofid );

	interpolation->boundaryEvalNormal(answer, side, xi, FEIElementGeometryWrapper(thisElement));

}

void WeakPeriodicBoundaryCondition :: updateDirection()
{
	// Check orientation for s
	FloatArray normal;

	if (this->domain->giveNumberOfSpatialDimensions()==2) {
		surfaceIndexes.resize(1);
		smin.resize(1);
		smax.resize(1);
		sideGeom = _Line;
	} else {
		surfaceIndexes.resize(2);
		smin.resize(2);
		smax.resize(2);
		sideGeom = _Triangle;
	}

	giveEdgeNormal( normal, element [ 0 ].at(0), side [ 0 ].at(0) );

	if ( fabs(normal.at(1))>0.99999 ) {              // Normal points in X direction
		direction = 1;
		if (this->domain->giveNumberOfSpatialDimensions()==2) {
			surfaceIndexes.at(1)=2;
		} else {
			surfaceIndexes.at(1)=2;
			surfaceIndexes.at(2)=3;
		}

	} else if ( fabs(normal.at(2))>0.99999) {         // Normal points in Y direction
		direction = 2;
		if (this->domain->giveNumberOfSpatialDimensions()==2) {
			surfaceIndexes.at(1)=1;
		} else {
			surfaceIndexes.at(1)=1;
			surfaceIndexes.at(2)=3;
		}
	} else if ( fabs(normal.at(3))>0.99 ) {         // Normal points in Z direction
		direction = 3;
		if (this->domain->giveNumberOfSpatialDimensions()==2) {
			_error1("3 dimensioal normal in a 2 dimensional problem.\n");
		} else {
			surfaceIndexes.at(1)=1;
			surfaceIndexes.at(2)=2;
		}
	} else {
		normal.printYourself();
		Element *thisElement=this->giveDomain()->giveElement(element[0].at(0));
		_error3("Only surfaces with normal in x, y or z direction supported. (el=%d, side=%d) \n", thisElement->giveLabel(), side [ 0 ].at(0));
	}
}

void WeakPeriodicBoundaryCondition :: updateSminmax()
{
	if ( doUpdateSminmax ) {

		// If sets are used, now is the time to update lists of elements and sides since the sets are unknown in initializeFrom
		if (posSet!=-1) {
			IntArray posBoundary, negBoundary;

			posBoundary=this->giveDomain()->giveSet(posSet)->giveBoundaryList();
			for (int i=0; i<posBoundary.giveSize() / 2; i++) {
				side [ 0 ].push_back( posBoundary.at(2 * i + 2) );
				element [ 0 ].push_back( posBoundary.at(2 * i + 1) );
			}

			negBoundary=this->giveDomain()->giveSet(negSet)->giveBoundaryList();
			for (int i=0; i<negBoundary.giveSize() / 2; i++) {
				side [ 1 ].push_back( negBoundary.at(2 * i + 2) );
				element [ 1 ].push_back( negBoundary.at(2 * i + 1) );
			}

		}

		// Determine which is the positive and which is the negative side
		FloatArray normal;
		giveEdgeNormal( normal, element [ 0 ].at(0), side [ 0 ].at(0) );

		double normalSum=-1;
		(this->giveDomain()->giveNumberOfSpatialDimensions()<=2) ? normalSum=normal.at(1)+normal.at(2) : normalSum=normal.at(1)+normal.at(2)+normal.at(3);

		if (normalSum > - 0.000001) { // No support for 3D yet
			sideSign [ 0 ] = 1;
			sideSign [ 1 ] = -1;
		} else {
			sideSign [ 0 ] = -1;
			sideSign [ 1 ] = 1;
		}

		updateDirection();

		for (int i=1; i<=surfaceIndexes.giveSize(); i++) {
			smin.at(i) = this->domain->giveDofManager(1)->giveCoordinate(surfaceIndexes.at(i));
			smax.at(i) = this->domain->giveDofManager(1)->giveCoordinate(surfaceIndexes.at(i));

			for ( int j = 1; j <= this->domain->giveNumberOfDofManagers(); j++ ) {
				double sValue = this->domain->giveDofManager(j)->giveCoordinate(surfaceIndexes.at(i));
				smin.at(i) = std :: min(smin.at(i), sValue);
				smax.at(i) = std :: max(smax.at(i), sValue);
			}

		}
		doUpdateSminmax = false;
	}
}

void WeakPeriodicBoundaryCondition :: addElementSide(int newElement, int newSide)
{
	//printf ("Add element %u, side %u\n", newElement, newSide);
	//	_error1("Not supported");

	FloatArray normalNew, normal0;
	int addToList = 0;

	if ( element [ 0 ].size() > 0 ) {
		// If there are elements in the list, compare normals in order to determine which list to store them in
		giveEdgeNormal(normalNew, newElement, newSide);
		//normalNew.printYourself();
		giveEdgeNormal( normal0, element [ 0 ].at(0), side [ 0 ].at(0) );
		double d = sqrt( pow(normalNew.at(1) - normal0.at(1), 2) + pow(normalNew.at(2) - normal0.at(2), 2) );
		if ( abs(d) < 0.001 ) {
			addToList = 0;
		} else {
			addToList = 1;
		}
	} else {
		// Otherwise, check the normal in order to decide upon which direction the parameter runs (x or y)
		giveEdgeNormal(normalNew, newElement, newSide);
		if ( abs(abs( normalNew.at(1) ) - 1) < 0.0001 ) {               // Normal point in x direction, thus set direction to y
			direction = 2;
		} else {
			direction = 1;
		}
	}

	element [ addToList ].push_back(newElement);
	side [ addToList ].push_back(newSide);
}

void WeakPeriodicBoundaryCondition :: computeElementTangent(FloatMatrix &B, Element *e, int boundary)
{
	FloatArray gcoords;
	IntArray bnodes;
	FEInterpolation *geoInterpolation = e->giveInterpolation();
	FEInterpolation *interpolation = e->giveInterpolation( (DofIDItem)dofid );

	interpolation->boundaryGiveNodes(bnodes, boundary);

	B.resize(bnodes.giveSize(), ndof);
	B.zero();
	///@todo Add this to all interpolators, should return _Point, _Line, _Triangle, ... etc.
	//integrationDomain sideGeom = geoInterpolation->boundaryGiveGeometry(boundary);

	GaussIntegrationRule iRule(1, e);
	if ( ngp == -1 ) {
		ngp = iRule.getRequiredNumberOfIntegrationPoints(sideGeom, 2 + orderOfPolygon);
	}
	iRule.setUpIntegrationPoints(sideGeom, ngp, _Unknown);

	for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
		GaussPoint *gp = iRule.getIntegrationPoint(i);
		FloatArray *lcoords = gp->giveCoordinates();

		FloatArray N;

		// Find the value of parameter s which is the vert/horiz distance to 0
		geoInterpolation->boundaryLocal2Global( gcoords, boundary, * lcoords, FEIElementGeometryWrapper(e) );
		// Compute base function values
		interpolation->boundaryEvalN( N, boundary, * lcoords, FEIElementGeometryWrapper(e) );
		// Compute Jacobian
		double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( boundary, * lcoords, FEIElementGeometryWrapper(e) ) );
		double s = gcoords.at(surfaceIndexes.at(1));

		for ( int j = 0; j < B.giveNumberOfColumns(); j++ ) {

			double fVal;

			if (this->domain->giveNumberOfSpatialDimensions()==2 )
				fVal = computeBaseFunctionValue(j, s);
			else {
				int a, b;
				getExponents(j+1, a, b);
				fVal=pow(gcoords.at(surfaceIndexes.at(1)), a)*pow(gcoords.at(surfaceIndexes.at(2)), b);
			}

			for ( int k = 0; k < B.giveNumberOfRows(); k++ ) {
				B(k, j) += + N(k) * fVal * detJ * gp->giveWeight();
			}
		}
	}
}

void WeakPeriodicBoundaryCondition :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
	if ( type != StiffnessMatrix ) {
		return;
	}

	IntArray c_loc, r_loc;
	gammaDman->giveCompleteLocationArray( r_loc, r_s );
	gammaDman->giveCompleteLocationArray( c_loc, c_s );

	FloatMatrix B, BT;
	FloatArray gcoords, normal;
	int normalSign, dofCountOnBoundary;

	updateSminmax();

	// Assemble each side
	for ( int thisSide = 0; thisSide <= 1; thisSide++ ) {

		normalSign = sideSign[thisSide];

		for ( size_t ielement = 0; ielement < element [ thisSide ].size(); ielement++ ) {       // Loop over each element on this edge
			Element *thisElement = this->domain->giveElement( element [ thisSide ].at(ielement) );

			// Find dofs for this element side
			IntArray r_sideLoc, c_sideLoc;

			// Find dofs for this element which should be periodic
			IntArray bNodes, nodeDofIDMask, periodicDofIDMask, nodalArray;
			periodicDofIDMask.resize(1);
			periodicDofIDMask.at(1) = dofid;

			FEInterpolation *interpolation = thisElement->giveInterpolation( (DofIDItem)dofid );
			interpolation->boundaryGiveNodes( bNodes, side [ thisSide ].at(ielement) );

            ///@todo See todo on assembleVector
            //thisElement->giveBoundaryLocationArray(r_sideLoc, bNodes, dofids, eid, r_s);
            //thisElement->giveBoundaryLocationArray(c_sideLoc, bNodes, dofids, eid, c_s);

			r_sideLoc.resize(0);
			c_sideLoc.resize(0);
			dofCountOnBoundary = 0;
			for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
				thisElement->giveDofManDofIDMask(bNodes.at(i), EID_MomentumBalance_ConservationEquation, nodeDofIDMask);

				for ( int j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
					if ( nodeDofIDMask.at(j) == dofid ) {
						thisElement->giveDofManager( bNodes.at(i) )->giveLocationArray(periodicDofIDMask, nodalArray, r_s);
						r_sideLoc.followedBy(nodalArray);
						thisElement->giveDofManager( bNodes.at(i) )->giveLocationArray(periodicDofIDMask, nodalArray, c_s);
						c_sideLoc.followedBy(nodalArray);
						dofCountOnBoundary++;
						break;
					}
				}
			}

			this->computeElementTangent(B, thisElement, side [ thisSide ].at(ielement));
			B.times( normalSign );

			BT.beTranspositionOf(B);

			answer->assemble(r_sideLoc, c_loc, B);
			answer->assemble(r_loc, c_sideLoc, BT);
		}
	}

}

double WeakPeriodicBoundaryCondition :: computeBaseFunctionValue(int baseID, double coordinate)
{
	double fVal=0.0;
	FloatArray sideLength;

	// compute side lengths
	sideLength.resize(smax.giveSize());
	for (int i=1; i<=smax.giveSize(); i++) {
		sideLength.at(i) = smax.at(i) - smin.at(i);
	}

	if ( useBasisType == monomial ) {
		fVal = pow(coordinate, baseID);
	} else if ( useBasisType == trigonometric ) {
		if ( baseID % 2 == 0 ) {   // Even (does not yet work in 3D)
			fVal = cos( ( ( double ) baseID ) / 2. * ( coordinate * 2. * M_PI / sideLength.at(1) ) );
		} else {
			fVal = sin( ( ( double ) baseID + 1 ) / 2. * ( coordinate * 2. * M_PI / sideLength.at(1) ) );
		}
	} else if ( useBasisType == legendre ) {
		double n = (double) baseID;
		coordinate = 2.0*coordinate-1.0;
		for ( int k = 0.; k <= baseID; k++ ) {
			fVal = fVal + binomial(n,k) * binomial(-n-1.0, k) * pow((1.0-coordinate)/2.0, (double)k);
		}
	}

	return fVal;
}

void WeakPeriodicBoundaryCondition :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
		CharType type, ValueModeType mode,
		const UnknownNumberingScheme &s, FloatArray *eNorms)
{
	if ( type != InternalForcesVector ) {
		return;
	}

	// Fetch unknowns of this boundary condition
	IntArray gammaLoc, gammaDofIDs;
	FloatArray gamma;
	gammaDman->giveCompleteUnknownVector( gamma, mode, tStep);
	gammaDman->giveCompleteLocationArray( gammaLoc, s );
	gammaDman->giveCompleteMasterDofIDArray( gammaDofIDs );

	// Values from solution
	FloatArray a;
	// Find dofs for this element side
	IntArray sideLocation, masterDofIDs;

	FloatMatrix B;

	FloatArray gcoords, normal;
	int normalSign, dofCountOnBoundary;

	updateSminmax();

	// Assemble each side
	for ( int thisSide = 0; thisSide <= 1; thisSide++ ) {

		normalSign = sideSign[thisSide];

		for ( size_t ielement = 0; ielement < element [ thisSide ].size(); ielement++ ) {           // Loop over each element on this edge

			Element *thisElement = this->domain->giveElement( element [ thisSide ].at(ielement) );

			// Find dofs for this element which should be periodic
			IntArray bNodes, nodeDofIDMask, periodicDofIDMask, nodalArray;
			periodicDofIDMask.resize(1);
			periodicDofIDMask.at(1) = dofid;

			FEInterpolation *interpolation = thisElement->giveInterpolation( (DofIDItem)dofid );
			interpolation->boundaryGiveNodes( bNodes, side [ thisSide ].at(ielement) );

            ///@todo Carl, change to this:
            //thisElement->giveBoundaryLocationArray(sideLocation, bNodes, &dofids, eid, s, &masterDofIDs);
            //thisElement->computeBoundaryVectorOf(bNodes, eid, VM_Total, tStep, a);

			sideLocation.resize(0);
			masterDofIDs.resize(0);
			a.resize(0);
			dofCountOnBoundary = 0;
			for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
				thisElement->giveDofManDofIDMask(bNodes.at(i), EID_MomentumBalance_ConservationEquation, nodeDofIDMask);

				for ( int j = 1; j <= nodeDofIDMask.giveSize(); j++ ) {
					if ( nodeDofIDMask.at(j) == dofid ) {
						thisElement->giveDofManager( bNodes.at(i) )->giveLocationArray(periodicDofIDMask, nodalArray, s);
						sideLocation.followedBy(nodalArray);
						thisElement->giveDofManager( bNodes.at(i) )->giveMasterDofIDArray(periodicDofIDMask, nodalArray);
						masterDofIDs.followedBy(nodalArray);

						double value = thisElement->giveDofManager( bNodes.at(i) )->giveDof(j)->giveUnknown(mode, tStep);
						a.resizeWithValues( sideLocation.giveSize() );
						a.at( sideLocation.giveSize() ) = value;
						dofCountOnBoundary++;
						break;
					}
				}
			}


			this->computeElementTangent(B, thisElement, side [ thisSide ].at(ielement));
			B.times( normalSign );

			FloatArray myProd, myProdGamma;
			myProd.beTProductOf(B, a);
			myProdGamma.beProductOf(B, gamma);

			if ( eNorms ) {
				for ( int i = 1; i <= gammaLoc.giveSize(); ++i ) {
					if ( gammaLoc.at(i) )
						eNorms->at(gammaDofIDs.at(i)) += myProd.at(i) * myProd.at(i);
				}
				for ( int i = 1; i <= sideLocation.giveSize(); ++i ) {
					if ( sideLocation.at(i) )
						eNorms->at(masterDofIDs.at(i)) += myProdGamma.at(i) * myProdGamma.at(i);
				}
			}

			answer.assemble(myProd, gammaLoc);
			answer.assemble(myProdGamma, sideLocation);
		}
	}
}

int WeakPeriodicBoundaryCondition :: giveNumberOfInternalDofManagers()
{
	return 1;
}

DofManager *WeakPeriodicBoundaryCondition :: giveInternalDofManager(int i)
{
	if ( i == 1 ) {
		return gammaDman;
	} else {
		return NULL;
	}
}

double WeakPeriodicBoundaryCondition::factorial(int n)
{
	int x = 1.0;
	for ( int i = 1; i <= n; i++ ) {
		x = x * i;
	}
	return x;
}

double WeakPeriodicBoundaryCondition::binomial(double n, int k)
{
	double f = 1.0;
	for ( int i = 1; i <= k; i++ ) {
		f = f * ( n - ( k - i ) ) / i;
	}
	return f;
}

void  WeakPeriodicBoundaryCondition::getExponents(int term, int &i, int &j)
{

	bool doContinue=true;

	// c is the number of the current term
	int c=0;

	// n is the order of the current polynomial (row in Pascals triangle)
	int n=0;

	while (doContinue) {
		for (int t=0; t<=n; t++) {
			c++;
			if (c==term) {
				i=n-t;
				j=t;
				doContinue=false;
				break;
			}
		}
		n++;
	}


}

}

