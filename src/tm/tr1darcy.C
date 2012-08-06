/*
 * Q1TrBase.C
 *
 *  Created on: Feb 25, 2010
 *      Author: carl
 */

#include "tr1darcy.h"
// #include "./../fm/fmelement.h"	// Fixa
#include "node.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "bcgeomtype.h"
#include "generalbc.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "crosssection.h"
#include "matresponseform.h"
#include "matresponsemode.h"
#include "fei2dtrlin.h"
// #include "./../fm/fluiddynamicmaterial.h" // Fixa

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#endif


namespace oofem {

FEI2dTrLin Tr1Darcy :: interpolation_lin(1, 2);

Tr1Darcy :: Tr1Darcy(int n, Domain *aDomain) : TransportElement (n, aDomain)
{
	/*
	 * Contructor
	 */
	numberOfDofMans  = 3;
	numberOfGaussPoints = 1;
}

Tr1Darcy :: ~Tr1Darcy ()
{
	/*
	 * Destructor
	 */
}

IRResultType Tr1Darcy :: initializeFrom(InputRecord *ir)
{
	/*
	 * Initalize the element. I.e set up Gauss points, geometry etc. Additional information is stored in argument ir.
	 */
    const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    this->TransportElement :: initializeFrom(ir);

    this->computeGaussPoints();
    //this->initGeometry();

    //this->area=this->computeVolume();

    return IRRT_OK;
}

void Tr1Darcy :: initGeometry()
{
	/*
	 * Setup geometry. Even though there are quadratic shape functions on the element, since the pressure is
	 * approximated using linear shape functions, the area is calculated according to the linear ones.
	 */
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3;

    node1 = giveNode(1);
    node2 = giveNode(2);
    node3 = giveNode(3);

    // init geometry data
    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    this->area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

}

void Tr1Darcy :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{

	FloatArray *lcoords = gp->giveCoordinates();
	this->interpolation_lin.evaldNdx(answer, *lcoords, FEIElementGeometryWrapper(this));
	return;
};

void Tr1Darcy :: computeGaussPoints()
{
	/*
	 *	Set up gausspoints for element
	 */
	if (!integrationRulesArray) {
		numberOfIntegrationRules = 1;
		integrationRulesArray = new IntegrationRule * [ 1 ];
		integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
		integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _2dHeat);//_2dFlow);
	}
}

void Tr1Darcy :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *atTime)
{
	/*
	 * Return Ke = integrate(B^T K B)
	 */

	FloatMatrix Bxieta,B, BT, K, K2, Ke_Temp1, Ke_Temp2, J, JT, dNdxy;
	FloatArray devstress, eps, a;
	FloatArray *lcoords;
	GaussPoint *gp;

	double detJ;

	TransportMaterial *mat = ( TransportMaterial * ) this->domain->giveMaterial(this->material);

	IntegrationRule *iRule = integrationRulesArray [ 0 ];


	K2.resize(2,2);

	//K.printYourself();

	answer.resize(3,3);
	answer.zero();

	this->computeVectorOf(EID_ConservationEquation, VM_Total, atTime, a);
	//a.printYourself();

	for (int i=0; i<iRule->getNumberOfIntegrationPoints(); i++) {

		gp = iRule->getIntegrationPoint(i);
		lcoords = gp->giveCoordinates();

		double detJ = this->interpolation_lin.giveTransformationJacobian( *lcoords, FEIElementGeometryWrapper(this));
		this->interpolation_lin.evaldNdx(B, *lcoords, FEIElementGeometryWrapper(this));
		eps.beTProductOf(B,a);
		eps.resize(3,3);	// To avoid memory trouble...

		mat->giveFluxVector(devstress, gp, eps, atTime);
		mat->giveCharacteristicMatrix(K, FullForm, TangentStiffness, gp, atTime);

		//mat->computeDeviatoricStressVector(devstress, gp, eps, atTime);
		//mat->giveDeviatoricStiffnessMatrix(K, TangentStiffness, gp, atTime);

		K2.at(1,1)=K.at(1,1);
		K2.at(1,2)=K.at(1,2);
		K2.at(2,1)=K.at(2,1);
		K2.at(2,2)=K.at(2,2);

		Ke_Temp1.beProductOf(B, K2);
		Ke_Temp2.beProductTOf(Ke_Temp1, B);
		Ke_Temp2.times(detJ*gp->giveWeight());

		answer.add(Ke_Temp2);

	}

	return;

}

void Tr1Darcy ::  giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                         TimeStep *tStep)
{
	/*
	 * Compute characteristic vector for this element. I.e the load vector(s)
	 *
	 * 	TODO: Implement support for body forces
	 *
	 */

	if( mtrx == ExternalForcesVector )  {
		this->computeLoadVector(answer, tStep);
	} else if ( mtrx == InternalForcesVector )  {
		this->computeInternalForcesVector(answer, tStep);
	}
	else {
		_error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
	}

	return;
}

void Tr1Darcy :: computeInternalForcesVector (FloatArray &answer, TimeStep *atTime)
{

	double detJ;

	FloatArray *lcoords, w, a, gradP, I;
	FloatMatrix N, B, p;
	GaussPoint *gp;

	TransportMaterial *mat = ( TransportMaterial * ) this->domain->giveMaterial(this->material);
	IntegrationRule *iRule = integrationRulesArray [ 0 ];

	// a.resize(3);

	this->computeVectorOf(EID_ConservationEquation, VM_Total, atTime, a);

	answer.resize(3);
	answer.zero();

	for (int i=0; i<iRule->getNumberOfIntegrationPoints(); i++) {
		gp = iRule->getIntegrationPoint(i);
		lcoords = gp->giveCoordinates();

		double detJ = this->interpolation_lin.giveTransformationJacobian( *lcoords, FEIElementGeometryWrapper(this));
		this->interpolation_lin.evaldNdx(B, *lcoords, FEIElementGeometryWrapper(this));

		gradP.beTProductOf(B,a);
		gradP.resize(3,3);

		mat->giveFluxVector(w, gp, gradP, atTime);

		I.beProductOf(B, w);
		I.times(-1*gp->giveWeight()*detJ);

		answer.add(I);

	}


}

void Tr1Darcy :: computeLoadVector(FloatArray &answer, TimeStep *atTime)
{
	FloatArray p, vec;
	FloatMatrix Ke;

	answer.resize(3);
	answer.zero();

	// Compute characteristic vector for Neumann boundary conditions.
    int i, load_number, load_id;
    GeneralBoundaryCondition *load;
    bcGeomType ltype;

   	int nLoads = boundaryLoadArray.giveSize() / 2;

    for (i=1;i<=nLoads;i++ ) {	// For each Neumann boundary condition ....
    	load_number = boundaryLoadArray.at(2*i-1);
    	load_id = boundaryLoadArray.at(2*i);
    	load = ( GeneralBoundaryCondition * ) domain->giveLoad(load_number);
    	ltype = load->giveBCGeoType();

    	if (ltype==EdgeLoadBGT) {
			this->computeEdgeBCSubVectorAt(vec, (Load *) load, load_id, atTime);
    	}

    	answer.add(vec);
    }
    answer.times(-1.0);

}

void Tr1Darcy :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep)
{
	/*
	 * Given the load *load, return it's contribution.
	 *
	 */

    answer.resize(3);

	answer.zero();

	if (load->giveType()==TransmissionBC) {			// Neumann boundary conditions (traction)

		BoundaryLoad *boundaryLoad;
		boundaryLoad = (BoundaryLoad *) load;

		int numberOfEdgeIPs;
		numberOfEdgeIPs = (int) ceil ((boundaryLoad->giveApproxOrder()+1.)/2.)*2;

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatMatrix N;
        FloatArray loadValue, reducedAnswer;
		reducedAnswer.resize(3);
		reducedAnswer.zero();
        IntArray mask;

        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);

        for (int i=0;i<iRule.getNumberOfIntegrationPoints(); i++) {

        	gp = iRule.getIntegrationPoint(i);
    		this->giveNEdge_xieta(N, gp);
			double dV=this->computeEdgeVolumeAround(gp, iEdge);

			if (boundaryLoad->giveFormulationType()==BoundaryLoad :: BL_EntityFormulation) {	// Edge load in xi-eta system
				boundaryLoad->computeValueAt(loadValue, tStep, *(gp->giveCoordinates()), VM_Total);
        	}
        	else {	// Edge load in x-y system

        	}

			reducedAnswer.at(1)+=N.at(1,1)*loadValue.at(1)*dV;
			reducedAnswer.at(2)+=N.at(1,2)*loadValue.at(1)*dV;
			reducedAnswer.at(3)+=N.at(1,3)*loadValue.at(1)*dV;

        }

		this->giveEdgeDofMappingV(mask, iEdge);

		mask.resize(3);
		mask.at(3)=0;

        answer.assemble(reducedAnswer, mask);
	}

}

void Tr1Darcy :: giveNEdge_xieta(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
	/*
	 * Returns value of quadratic shape functions on interval [-1 1]. User when computing contributions
	 * to load vector from boundary terms thru Gauss integration.
	 */

	double xi=-1, xj=0, xk=1, x;

	answer.resize(1,3);
	answer.zero();

	x=aGaussPoint->giveCoordinate(1);

	answer.at(1,1) = ((x - xj)*(x - xk))/((xi - xj)*(xi - xk));
	answer.at(1,2) = -(((x - xi)*(x - xk))/((xi - xj)*(xj - xk)));
	answer.at(1,3) = ((x - xi)*(x - xj))/((xi - xk)*(xj - xk));

}

double Tr1Darcy :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
	/*
	 * Gives length of edge iEdge.
	 */
	Node *node1, *node2;

	if (iEdge==1) {
		node1 = this->giveNode(1);
		node2 = this->giveNode(2); // 4
	}
	else if (iEdge==2) {
		node1 = this->giveNode(2);
		node2 = this->giveNode(3); // 5
	}
	else if (iEdge==3) {
		node1 = this->giveNode(3);
		node2 = this->giveNode(1); // 6
	}

	double dx = node1->giveCoordinate(1) - node2->giveCoordinate(1);
	double dy = node1->giveCoordinate(2) - node2->giveCoordinate(2);
	double length = sqrt(dx*dx+dy*dy);
	double thickness;
	thickness = 1;//this->giveCrossSection()->give('t');
	return length*thickness*gp->giveWeight()/2;

}

void Tr1Darcy :: giveEdgeDofMappingV(IntArray &answer, int iEdge)
{
	/*
	 * Given an edge iEdge, return velocity dofs at that edge
	 */
	answer.resize(2);

	if (iEdge==1) {
		answer.at(1)=1;
		answer.at(2)=2;
	} else if (iEdge==2) {
		answer.at(1)=2;
		answer.at(2)=3;
	} else if (iEdge==3) {
		answer.at(1)=3;
		answer.at(2)=1;
	}
}


void Tr1Darcy ::  giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
{
	/*
	 * Compute characteristic matrix for this element. The only option is the stiffness matrix...
	 */

	if ( mtrx == StiffnessMatrix )  {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }

    return;
}

void Tr1Darcy :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
	/*
	 * Returns the mask for node number inode of this element. The mask tells what quantities
	 * are held by each node. Since this element holds velocities (both in x and y direction),
	 * in six nodes and pressure in three nodes the answer depends on which node is requested.
	 */

	if ((inode==1)||(inode==2)||(inode==3)) {
		if ( ut == EID_ConservationEquation ) {
			answer.resize(1);
			answer.at(1) = P_f;
		} else {
			_error("giveDofManDofIDMask: Unknown equation id encountered");
		}

	}

}

void
Tr1Darcy :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                          InternalStateType type, TimeStep *tStep)
{

	GaussPoint *gp;
	TransportMaterial *mat = ( TransportMaterial * ) this->domain->giveMaterial(this->material);
	// FluidDynamicMaterial *mat = ( FluidDynamicMaterial * ) this->domain->giveMaterial(this->material);

	IntegrationRule *iRule = integrationRulesArray [ 0 ];
	gp = iRule->getIntegrationPoint(0);
	mat->giveIPValue(answer, gp, type, tStep);

}

void
Tr1Darcy :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                         InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

Interface *
Tr1Darcy :: giveInterface(InterfaceType interface)
{
    if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    }

    return NULL;
}

int
Tr1Darcy :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance ) {
        return 6;
    } else if ( ut == EID_ConservationEquation ) {
        return 3;
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 9;
    } else {
        _error("computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}
}
