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
#include "dynamicinputrecord.h"
#include "bodyload.h"
#include "boundarycondition.h"

//#include "classfactory.h"

#include <vector>

namespace oofem {

REGISTER_BoundaryCondition( SolutionbasedShapeFunction );

SolutionbasedShapeFunction::SolutionbasedShapeFunction(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    isLoaded=false;
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
    IR_GIVE_OPTIONAL_FIELD(ir, this->filename, _IFT_SolutionbasedShapeFunction_ShapeFunctionFile);
    useConstantBase = (this->filename == "") ? true : false;
    //	if (!useConstantBase) {loadProblem();};

    // Set up master dofs
    myNode = new Node(0, this->giveDomain());

    for (int i=1; i<=this->giveDofIDs().giveSize(); i++) {
        MasterDof *newDof = new MasterDof( 0, myNode, (DofIDItem)this->domain->giveNextFreeDofID() );
        myNode->appendDof( newDof );

        //		DynamicInputRecord ir;
        //		FloatArray gradP;
        //		ir.setRecordKeywordField(_IFT_BoundaryCondition_Name, 1);
        //		ir.setField(1, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
        //		if (i==2) {
        //			ir.setField(1.0, _IFT_BoundaryCondition_PrescribedValue);
        //		} else {
        //			ir.setField(0.0, _IFT_BoundaryCondition_PrescribedValue);
        //		}
        //		int bcID = this->giveDomain()->giveNumberOfBoundaryConditions()+1;
        //		GeneralBoundaryCondition *myBC;
        //		myBC = classFactory.createBoundaryCondition(_IFT_BoundaryCondition_Name, bcID, this->giveDomain());
        //		myBC->initializeFrom(&ir);
        //
        //		this->giveDomain()->setBoundaryCondition(bcID, myBC);
        //
        //		newDof->setBcId(bcID);

    }


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
            minCoord.at(j)=min( this->giveDomain()->giveDofManager(i)->giveCoordinate(j), minCoord.at(j) );
        }
    }

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

    FloatArray shapeFunctionValues;
    computeDofTransformation(dof, shapeFunctionValues);

    FloatArray gamma;
    myNode->giveCompleteUnknownVector( gamma, mode, tStep); // alpha1, alpha2,...

    double out = shapeFunctionValues.dotProduct(gamma);

    return out;
}

void
SolutionbasedShapeFunction :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
    masterContribs.resize(this->giveDomain()->giveNumberOfSpatialDimensions());

    for (int i=1; i<=this->giveDomain()->giveNumberOfSpatialDimensions(); i++) {
        masterContribs.at(i) = computeBaseFunctionValueAt(dof->giveDofManager()->giveCoordinates(), dof, i);
    }

}

Dof *
SolutionbasedShapeFunction :: giveMasterDof(ActiveDof *dof, int mdof)
{
    return myNode->giveDof(mdof);
}

void
SolutionbasedShapeFunction :: loadProblem()
{

    for (int i=0; i<this->domain->giveNumberOfSpatialDimensions(); i++) {

        OOFEM_LOG_INFO("************************** Instanciating microproblem from file %s for dimension %u\n", filename.c_str(), i);
        OOFEMTXTDataReader drMicro(filename.c_str());

        EngngModel *myEngngModel = InstanciateProblem(& drMicro, _processor, 0);

        drMicro.finish();

        myEngngModel->checkProblemConsistency();

        myEngngModel->initMetaStepAttributes( myEngngModel->giveMetaStep( 1 ) );
        thisTimestep = myEngngModel->giveNextStep();
        myEngngModel->init();

        this->setLoads(myEngngModel, i+1);

        myEngngModel->solveYourselfAt(thisTimestep);
        myEngngModel->doStepOutput(thisTimestep);

        myEngngModels.push_back( myEngngModel );

        OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", myEngngModel);
    }

    isLoaded=true;
}

void
SolutionbasedShapeFunction :: setLoads(EngngModel *myEngngModel, int d)
{
    DynamicInputRecord ir;
    FloatArray gradP;

    gradP.resize(this->giveDomain()->giveNumberOfSpatialDimensions());
    gradP.zero();
    gradP.at(d)=1.0;

    ir.setRecordKeywordField("deadweight", 1);
    ir.setField(gradP, _IFT_Load_components);
    ir.setField(1, _IFT_GeneralBoundaryCondition_LoadTimeFunct);

    int bcID = myEngngModel->giveDomain(1)->giveNumberOfBoundaryConditions()+1;
    GeneralBoundaryCondition *myBodyLoad;
    myBodyLoad = classFactory.createBoundaryCondition("deadweight", bcID, myEngngModel->giveDomain(1));
    myBodyLoad->initializeFrom(&ir);
    myEngngModel->giveDomain(1)->setBoundaryCondition(bcID, myBodyLoad);

    for (int i=1; i<=myEngngModel->giveDomain(1)->giveNumberOfElements(); i++) {
        IntArray *blArray;
        blArray = myEngngModel->giveDomain(1)->giveElement(i)->giveBodyLoadArray();
        blArray->resizeWithValues(blArray->giveSize()+1);
        blArray->at(blArray->giveSize())=bcID;
    }
}

double
SolutionbasedShapeFunction :: computeBaseFunctionValueAt(FloatArray *coords, Dof *dof, int d)
{
    // The followind define allows for use of weakly periodic bc to be copied. This does not comply with the theory but is used to check the validity of the code.

    if (useConstantBase) {
        return 1.0;
    } else {
        FloatArray closest, lcoords, values;
        std :: vector <FloatArray *> checkcoords;
        std :: vector <int> permuteIndex;
        double outvalue=0.0;
        int n=0;

        if (!isLoaded) loadProblem();

        this->myEngngModels.at(d-1)->giveDomain(1)->giveSpatialLocalizer()->init(false);

        //printf("******************** Check at (%f, %f, %f)\n", coords->at(1), coords->at(2), coords->at(3) );

        if (this->giveDomain()->giveNumberOfSpatialDimensions()==2) {
            coords->resize(2);
        }

        // Determine if current coordinate is at a max or min point and if so, which type (on a surface, edge or a corner?)
        // n is the number of dimensions at which the point is an extremum, permuteIndex tells in which dimension the coordinate is max/min

        int thisMask=0;

        for (int i=1; i<=coords->giveSize(); i++) {
            if ( ( fabs(maxCoord.at(i)-coords->at(i))<1e-8) || (fabs(minCoord.at(i)-coords->at(i))<1e-8) ){
                permuteIndex.push_back(i);
                n++;
                thisMask = thisMask + pow(2.0,i-1);
            }
        }

        for (int i=0; i<pow(2.0,n); i++) {

            int mask=i, counter=1;
            FloatArray *newCoord = new (FloatArray)(coords->giveSize());
            *newCoord = *coords;

            for (int j=1; j<=n; j++) {
                double d=1e-15;
                if ( (mask & 1) == 0 ) { // Max
                    newCoord->at(permuteIndex.at(counter-1))=minCoord.at( permuteIndex.at(counter-1) )+d;
                } else { // Min
                    newCoord->at(permuteIndex.at(counter-1))=maxCoord.at( permuteIndex.at(counter-1) )-d;
                }
                counter ++;
                mask = mask >> 1;
            }
            checkcoords.push_back(newCoord);
        }

#define USEWPBC 0
#if USEWPBC == 1
        FloatArray *tempCoord = new FloatArray;
        *tempCoord = *coords;
        checkcoords.clear();
        checkcoords.push_back(tempCoord);
#endif

        for (size_t i=0; i<checkcoords.size(); i++) {
            //			checkcoords.at(i)->printYourself();
            for (int j=1; j<=checkcoords.at(i)->giveSize(); j++) {
                if (fabs(checkcoords.at(i)->at(j)-minCoord.at(j))<=0.00001) {
                    checkcoords.at(i)->at(j)=minCoord.at(j)+1e-6;
                }
                if (fabs(checkcoords.at(i)->at(j)-maxCoord.at(j))<=0.00001) {
                    checkcoords.at(i)->at(j)=maxCoord.at(j)-1e-6;
                }
            }

            Element *elementAtCoords = this->myEngngModels.at(d-1)->giveDomain(1)->giveSpatialLocalizer()->giveElementContainingPoint (*checkcoords.at(i));

            if (elementAtCoords==NULL) {
                elementAtCoords = this->myEngngModels.at(d-1)->giveDomain(1)->giveSpatialLocalizer()->giveElementClosestToPoint(lcoords, closest, *checkcoords.at(i), 1);
                if (elementAtCoords==NULL) {
                    printf("Cannot find element closest to point\n");
                    checkcoords.at(i)->pY();
                }
            }

            //			checkcoords.at(i)->pY();

            //			lcoords.printYourself();

            //			printf("Closest element: %u\n", elementAtCoords->giveNumber());

            //			elementAtCoords->giveDofManArray().printYourself();

            //			checkcoords.at(i)->printYourself();
            //			closest.printYourself();

            EIPrimaryUnknownMapperInterface *em = dynamic_cast<EIPrimaryUnknownMapperInterface*> ( elementAtCoords->giveInterface ( EIPrimaryUnknownMapperInterfaceType ) );

            IntArray dofids;

            em->EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(dofids);
            em->EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(VM_Total, thisTimestep, *checkcoords.at(i), values);

            //printf("\t%f, %f, %f : %f\n", checkcoords.at(i)->at(1), checkcoords.at(i)->at(2), checkcoords.at(i)->at(3), values.at(dofids.findFirstIndexOf(dof->giveDofID())));
#if USEWPBC == 1
            outvalue = values.at(dofids.findFirstIndexOf(dof->giveDofID()));
#else
            outvalue = outvalue + values.at(dofids.findFirstIndexOf(dof->giveDofID())) / ( (double) pow(2.0, n));
#endif
        }
        //printf("\tResult = %f\n", outvalue);
        //		printf("Output value %u (%f, %f, %f), dofID %u is %f, d=%u\n", dof->giveDofManNumber(), coords->at(1), coords->at(2), coords->at(3), dof->giveDofID(), outvalue, d);
        //		printf("Output value at (%f, %f), dofID %u is %f, d=%u\n", coords->at(1), coords->at(2), dof->giveDofID(), outvalue, d);
        return outvalue;
    }
}

}
