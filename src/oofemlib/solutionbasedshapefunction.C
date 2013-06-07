/*
 * solutionbasedshapefunction.C
 *
 *  Created on: May 27, 2013
 *      Author: carl
 */

#include "solutionbasedshapefunction.h"
#include "oofemtxtdatareader.h"
#include "activebc.h"
#include "activedof.h"
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
#include "util.h"
#include "inputrecord.h"
#include "engngm.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"

#include <vector>

namespace oofem {

REGISTER_BoundaryCondition( SolutionbasedShapeFunction );

SolutionbasedShapeFunction::SolutionbasedShapeFunction(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
	// TODO Auto-generated constructor stub

}

SolutionbasedShapeFunction::~SolutionbasedShapeFunction()
{
	// TODO Auto-generated destructor stub
}

IRResultType
SolutionbasedShapeFunction :: initializeFrom(InputRecord *ir)
{

	const char *__proc = "initializeFrom";
	IRResultType result;

	ActiveBoundaryCondition :: initializeFrom(ir);

	// Load and solve problem file
	this->filename="";
	IR_GIVE_FIELD(ir, this->filename, _IFT_SolutionbasedShapeFunction_ShapeFunctionFile);
	useConstantBase = (this->filename == "") ? true : false;
	if (!useConstantBase) {loadProblem();};

	// Set up master dofs
	myNode = new Node(0, this->giveDomain());

	for (int i=1; i<=this->giveDofIDs().giveSize(); i++)
		myNode->appendDof( new MasterDof( 0, myNode, (DofIDItem)this->domain->giveNextFreeDofID() ) );

	init();

	return IRRT_OK;

}

void
SolutionbasedShapeFunction :: init()
{
	Node *n;
	n=this->giveDomain()->giveNode(1);

	maxCoord = *n->giveCoordinates();
	minCoord = *n->giveCoordinates();

	for (int i=1; i<=this->giveDomain()->giveNumberOfDofManagers(); i++) {
		for (int j=1; j<=maxCoord.giveSize(); j++) {
			maxCoord.at(j)=max( this->giveDomain()->giveDofManager(i)->giveCoordinate(j), maxCoord.at(j) );
		}
	}

	maxCoord.printYourself();
	minCoord.printYourself();

}

double
SolutionbasedShapeFunction :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
	return this->giveUnknown(mode, tStep, dof);
}

double
SolutionbasedShapeFunction :: giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
{
	// Return value of pertinent quantity in coordinate given by dof

	FloatArray gamma;

	myNode->giveCompleteUnknownVector( gamma, mode, tStep);
	int index = dofs.findFirstIndexOf(dof->giveDofID());

	double shapeFunctionValue = computeBaseFunctionValueAt(dof->giveDofManager()->giveCoordinates(), dof);

	double out = gamma.at(index) * shapeFunctionValue;

	return out;
}

void
SolutionbasedShapeFunction :: giveLocationArrays(std::vector<IntArray> &rows, std::vector<IntArray> &cols, EquationID eid, CharType type,
		const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
	if (eid != EID_MomentumBalance_ConservationEquation && eid != EID_MomentumBalance)
		return;

	rows.resize(myNode->giveNumberOfDofs());
	cols.resize(myNode->giveNumberOfDofs());

	IntArray dofIDArray, temp;

	myNode->giveDofArray(dofIDArray, temp);

	for (int i=1; i<=myNode->giveNumberOfDofs(); i++) {
		myNode->giveLocationArray(dofIDArray, rows[i], r_s);
		myNode->giveLocationArray(dofIDArray, cols[i], c_s);
	}

}

void
SolutionbasedShapeFunction :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
	double baseFunctionvalue = computeBaseFunctionValueAt(dof->giveDofManager()->giveCoordinates(), dof);
	masterContribs.setValues(1, baseFunctionvalue);
}

Dof *
SolutionbasedShapeFunction :: giveMasterDof(ActiveDof *dof, int mdof)
{
	int index = dofs.findFirstIndexOf(dof->giveDofID());
	return myNode->giveDof(index);
}

void
SolutionbasedShapeFunction :: loadProblem()
{
	OOFEM_LOG_INFO("************************** Instanciating microproblem from file %s\n", filename.c_str());
	OOFEMTXTDataReader drMicro(filename.c_str());

	this->myEngngModel= InstanciateProblem(& drMicro, _processor, 0);
	drMicro.finish();

	this->myEngngModel->checkProblemConsistency();
	this->myEngngModel->initMetaStepAttributes( this->myEngngModel->giveMetaStep( 1 ) );
	thisTimestep = this->myEngngModel->giveNextStep();
	this->myEngngModel->init();
	this->myEngngModel->solveYourselfAt(thisTimestep);
	this->myEngngModel->doStepOutput(thisTimestep);

	OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", this->myEngngModel);
}

double
SolutionbasedShapeFunction :: computeBaseFunctionValueAt(FloatArray *coords, Dof *dof)
{
	if (useConstantBase) {
		return 1.0;
	} else {
		FloatArray closest, lcoords, values;
		std :: vector <FloatArray *> checkcoords;
		std :: vector <int> permuteIndex;
		IntArray coordUnion;
		double outvalue=0.0;
		int n=0;

		this->myEngngModel->giveDomain(1)->giveSpatialLocalizer()->init(false);

		if (this->giveDomain()->giveNumberOfSpatialDimensions()==2) {
			coords->resize(2);
		}

		coordUnion.resize(coords->giveSize());

		for (int i=1; i<=coords->giveSize(); i++) {
			coordUnion.at(i)=0;
			if (fabs(maxCoord.at(i)-coords->at(i))<1e-8) {
				coordUnion.at(i)=1;
				permuteIndex.push_back(i);
				n++;
			}
			if (fabs(minCoord.at(i)-coords->at(i))<1e-8) {
				coordUnion.at(i)=1;
				permuteIndex.push_back(i);
				n++;
			}
		}

		for (int i=0; i<pow(2,n); i++) {

			int mask=n, counter=1;
			FloatArray *newCoord = new (FloatArray)(coords->giveSize());
			*newCoord = *coords;

			while ( mask > 0 ) {
				if ( (i & mask) == 0 ) { // Max
					newCoord->at(permuteIndex.at(counter-1))=minCoord.at( permuteIndex.at(counter-1) )+1e-4;
				} else { // Min
					newCoord->at(permuteIndex.at(counter-1))=maxCoord.at( permuteIndex.at(counter-1) )-1e-4;
				}
				counter ++;
				mask = mask >> 1;
			}
			checkcoords.push_back(newCoord);
		}

		for (size_t i=0; i<checkcoords.size(); i++) {

			Element *elementAtCoords = this->myEngngModel->giveDomain(1)->giveSpatialLocalizer()->giveElementClosestToPoint(lcoords, closest, *checkcoords.at(i), 1);
			printf("Closest element: ");
			elementAtCoords->printYourself();
			EIPrimaryUnknownMapperInterface *em = dynamic_cast<EIPrimaryUnknownMapperInterface*> ( elementAtCoords->giveInterface ( EIPrimaryUnknownMapperInterfaceType ) );

			IntArray dofids;

			em->EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(dofids);
			if (fabs(checkcoords.at(i)->at(1)-1.0)<1e-3) {
				printf("faulty!\n");
			}
			em->EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(VM_Total, thisTimestep, *checkcoords.at(i), values);

			printf("Values at (%f, %f) are ", checkcoords.at(i)->at(1), checkcoords.at(i)->at(2) );
			for (int j = 1; j<values.giveSize(); j++) printf("%f ", values.at(j));
			printf("\n");

			outvalue = outvalue + values.at(dofids.findFirstIndexOf(dof->giveDofID())) / ( (double) n );
		}
		return outvalue;
	}
}

}
