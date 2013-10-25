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
#include "stokesflow.h"
//#include "classfactory.h"

#include <vector>

namespace oofem {

template <class T>
void logData(T myArray) {
    for (int xyz=1; xyz<=myArray.giveSize(); xyz++) {
        printf("%f ", myArray.at(xyz));
    }
}

template <class T>
void logDataMsg(const char *c, T myArray) {
    printf("%s ", c);
    logData(myArray);
    printf("\n");
}


REGISTER_BoundaryCondition( SolutionbasedShapeFunction );

SolutionbasedShapeFunction::SolutionbasedShapeFunction(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    isLoaded=false;
    TOL=1e-5;
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

    // Load problem file
    this->filename="";
    IR_GIVE_OPTIONAL_FIELD(ir, this->filename, _IFT_SolutionbasedShapeFunction_ShapeFunctionFile);
    useConstantBase = (this->filename == "") ? true : false;

    // Set up master dofs
    myNode = new Node(0, this->giveDomain());

    for (int i=1; i<=this->giveDofIDs().giveSize(); i++) {
        MasterDof *newDof = new MasterDof( 0, myNode, (DofIDItem)this->domain->giveNextFreeDofID() );
        myNode->appendDof( newDof );
        /*if (i==1) {
            setBoundaryConditionOnDof(newDof, 1.0);
        } else {
            setBoundaryConditionOnDof(newDof, 0.0);
        }*/
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

void
SolutionbasedShapeFunction :: computeCorrectionFactors(EngngModel &myEngngModel, IntArray *Dofs, double *Am, double *Ap, double *c0, double *cm, double *cp)
{
    /*
     *Compute c0, cp, cm, Ap and Am
     */

    *Am=0.0;
    *Ap=0.0;
    *c0=0.0;
    *cm=0.0;
    *cp=0.0;

    //modeStruct *mode=this->modes.at(EngngModelID);
    EngngModel *m=&myEngngModel;
    Set *mySet=this->domain->giveSet( this->giveSetNumber() );
    IntArray BoundaryList = mySet->giveBoundaryList();

    for (int i=0; i<BoundaryList.giveSize()/2; i++) {
        int ElementID=BoundaryList(2*i);
        int Boundary=BoundaryList(2*i+1);

        Element *thisElement=m->giveDomain(1)->giveElement(ElementID);
        FEInterpolation *geoInterpolation = thisElement->giveInterpolation();
        IntArray bnodes, bnodesIDs, zNodes, pNodes, mNodes;
        FloatArray nodeValues;

        //printf("*** Element %u (%u), boundary %u (i=%u)\n", ElementID, thisElement->giveNumber(), Boundary, i);

        geoInterpolation->boundaryGiveNodes(bnodes, Boundary);
        bnodesIDs.resize(bnodes.giveSize());

        // Change to global ID for bnodes and identify the intersection of bnodes and the zero boundary
        for (int j=1; j<=bnodes.giveSize(); j++) {

            DofManager *dman=thisElement->giveDofManager(bnodes.at(j));
            bnodesIDs.at(j)= dman->giveNumber();

            //printf("Node %u @ (%f, %f, %f)\n", bnodesIDs.at(j), dman->giveCoordinates()->at(1), dman->giveCoordinates()->at(2), dman->giveCoordinates()->at(3));

            bool isZero = false;
            bool isPlus = false;
            bool isMinus = false;

            for (int k=1; k<=dman->giveCoordinates()->giveSize(); k++) {
                isPlus = isPlus || (fabs(dman->giveCoordinates()->at(k)-maxCoord.at(k))<TOL);
                isMinus = isMinus || (fabs(dman->giveCoordinates()->at(k)-minCoord.at(k))<TOL);
            }

            isZero = isPlus && isMinus;

            if (isZero) {
                zNodes.insertSorted(j);
            } else if (isPlus) {
                pNodes.insertSorted(j);
            } else if (isMinus) {
                mNodes.insertSorted(j);
            }

        }

        GaussIntegrationRule iRule(1, thisElement);

        int n = iRule.getRequiredNumberOfIntegrationPoints(_Triangle, 4);
        iRule.setUpIntegrationPoints(_Triangle, n, _Unknown);

        for (int j=0; j<iRule.giveNumberOfIntegrationPoints(); j++) {

            GaussPoint *gp = iRule.getIntegrationPoint(j);
            FloatArray *lcoords = gp->giveCoordinates();
            FloatArray gcoords, normal, N;

            double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( Boundary, * lcoords, FEIElementGeometryWrapper(thisElement) ) )*gp->giveWeight();

            geoInterpolation->boundaryEvalNormal(normal, Boundary, * lcoords, FEIElementGeometryWrapper(thisElement));
            //logDataMsg("normal: ", normal);
            geoInterpolation->boundaryEvalN(N, Boundary, *lcoords, FEIElementGeometryWrapper(thisElement) );
            geoInterpolation->boundaryLocal2Global( gcoords, Boundary, * lcoords, FEIElementGeometryWrapper(thisElement) );

            //logDataMsg("N:", N);
            //N.times(detJ);
            //logDataMsg("N*detJ:", N);

            FloatArray pPhi, mPhi, zPhi;
            pPhi.resize(Dofs->giveSize()); pPhi.zero();
            mPhi.resize(Dofs->giveSize()); mPhi.zero();
            zPhi.resize(Dofs->giveSize()); zPhi.zero();

            // For each DofID
            for (int k=1; k<=Dofs->giveSize(); k++) {

                // Build array containing node values
                nodeValues.resize(bnodes.giveSize());
                for (int l=1; l<=bnodes.giveSize(); l++) {
                    nodeValues.at(l) = giveValueAtPoint(gcoords, (DofIDItem) Dofs->at(k), myEngngModel);
                }


                // Build zPhi for this DofID
                //printf("zeroNodes=[zeroNodes, ");
                for (int l=1; l<=zNodes.giveSize(); l++) {
                    int nodeID = zNodes.at(l);
                    //printf("%u, ", bnodesIDs.at(nodeID));
                    zPhi.at(k)=zPhi.at(k)+N.at(nodeID)*nodeValues.at(nodeID);
                }
                //printf("];\n");

                // Build pPhi for this DofID
                //printf("plusNodes=[plusNodes, ");
                for (int l=1; l<=pNodes.giveSize(); l++) {
                    int nodeID = pNodes.at(l);
                    //printf("%u, ", bnodesIDs.at(nodeID));
                    pPhi.at(k)=pPhi.at(k)+N.at(nodeID)*nodeValues.at(nodeID);
                }
                //printf("];\n");

                // Build mPhi for this DofID
                //printf("minusNodes=[minusNodes, ");
                for (int l=1; l<=mNodes.giveSize(); l++) {
                    int nodeID = mNodes.at(l);
                    //printf("%u, ", bnodesIDs.at(nodeID));
                    mPhi.at(k)=mPhi.at(k)+N.at(nodeID)*nodeValues.at(nodeID);
                }
                //printf("];\n");

            }

            //logDataMsg("zPhi = ", zPhi);
            //logDataMsg("mPhi = ", mPhi);
            //logDataMsg("pPhi = ", pPhi);

            *c0=*c0+zPhi.dotProduct(normal, 3)*detJ;
            *cp=*cp+pPhi.dotProduct(normal, 3)*detJ;
            *cm=*cm+mPhi.dotProduct(normal, 3)*detJ;
            *Ap=*Ap+pPhi.dotProduct(pPhi, 3)*detJ;
            *Am=*Am+mPhi.dotProduct(mPhi, 3)*detJ;

        }
    }

    printf("c0:%f, cp:%f, cm:%f, Ap:%f, Am:%f\n", *c0, *cp, *cm, *Ap, *Am);

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
    if (!isLoaded) loadProblem();

    masterContribs.resize(this->giveDomain()->giveNumberOfSpatialDimensions());

    for (int i=1; i<=this->giveDomain()->giveNumberOfSpatialDimensions(); i++) {
        masterContribs.at(i) = giveValueAtPoint(*dof->giveDofManager()->giveCoordinates(), dof->giveDofID(), *modes.at(i-1)->myEngngModel);
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

//    for (int i=0; i<this->domain->giveNumberOfSpatialDimensions(); i++) {
    for (int i=2; i>=0; i--) {

        OOFEM_LOG_INFO("************************** Instanciating microproblem from file %s for dimension %u\n", filename.c_str(), i);

        // Set up and solve problem
        OOFEMTXTDataReader drMicro(filename.c_str());
        EngngModel *myEngngModel = InstanciateProblem(& drMicro, _processor, 0, NULL, false);
        drMicro.finish();
        myEngngModel->checkProblemConsistency();
        myEngngModel->initMetaStepAttributes( myEngngModel->giveMetaStep( 1 ) );
        thisTimestep = myEngngModel->giveNextStep();
        myEngngModel->init();
        this->setLoads(myEngngModel, i+1);
        myEngngModel->solveYourselfAt(thisTimestep);

        // Set correct export filename
        std :: string originalFilename;
        originalFilename = myEngngModel->giveOutputBaseFileName();
        if (i==0) originalFilename = originalFilename + "_X";
        if (i==1) originalFilename = originalFilename + "_Y";
        if (i==2) originalFilename = originalFilename + "_Z";
        myEngngModel->letOutputBaseFileNameBe(originalFilename);

        myEngngModel->doStepOutput(thisTimestep);

        modeStruct *mode=new(modeStruct);
        mode->myEngngModel=myEngngModel;

        isLoaded=true;

        // Set unknowns to the mean value of opposite sides of the domain.
        // Loop thru all nodes and compute phi for all degrees of freedom on the boundary. Save phi in a list for later use.
        for (int j=1; j<=myEngngModel->giveDomain(1)->giveNumberOfDofManagers(); j++) {
            DofManager *dman=myEngngModel->giveDomain(1)->giveDofManager(j);
            FloatArray *coordinates=dman->giveCoordinates();

            bool isBoundary=false;

            // Used for indentifying the 0-boundary separating the + and - surfaces
            int maxCount=0, minCount=0;

            //printf("%u, %f, %f, %f\n", dman->giveNumber(), coordinates->at(1), coordinates->at(2), coordinates->at(3));

            for (int k=1; k<=coordinates->giveSize(); k++) {
                if (fabs(coordinates->at(k)-maxCoord.at(k))<TOL) {
                    isBoundary=true;
                    maxCount++;
                }
                if (fabs(coordinates->at(k)-minCoord.at(k))<TOL) {
                    isBoundary=true;
                    minCount++;
                }
            }

            double cSum = (coordinates->at(1)+coordinates->at(2)+coordinates->at(3));

            if ( ( cSum>=1.0 ) && (cSum<=2.0) && ( (maxCount==1) || (maxCount==2) ) && ( (minCount==1) || (minCount==2) ) )  {
                mode->ZeroBoundary.push_back(dman);
            }

            if (isBoundary) {

                for (int k=1; k<=dman->giveNumberOfDofs(); k++) {
                    SurfaceDataStruct *surfaceData=new(SurfaceDataStruct);
                    surfaceData->DofID = dman->giveDof(k)->giveDofID();
                    surfaceData->DofMan = dman;
                    surfaceData->value = computeBaseFunctionValueAt(*coordinates, dman->giveDof(k)->giveDofID(), *myEngngModel);

                    // Check with max/min instead of 1
                    surfaceData->isPlus = false;
                    surfaceData->isMinus = false;
                    for (int l=1; l<=dman->giveCoordinates()->giveSize(); l++) {
                        surfaceData->isPlus = surfaceData->isPlus || (fabs(dman->giveCoordinates()->at(l)-maxCoord.at(l))<TOL);
                        surfaceData->isMinus = surfaceData->isMinus || (fabs(dman->giveCoordinates()->at(l)-minCoord.at(l))<TOL);
                    }
                    surfaceData->isZeroBoundary = (surfaceData->isPlus && surfaceData->isMinus);

                    mode->SurfaceData.push_back(surfaceData);
                }
            }
        }


        modeStruct *m=mode;
        m->myEngngModel->doStepOutput(m->myEngngModel->giveCurrentStep());
        std :: vector <Dof *> nonFreeDofs;

        //Project mode onto descretized surface
        for (size_t j=0; j<m->SurfaceData.size(); j++) {

            DofManager *dman=m->SurfaceData.at(j)->DofMan;
            Dof *d=dman->giveDofWithID(m->SurfaceData.at(j)->DofID);

            //Add BC if dof is a master and has no existing boundary conditions
            MasterDof *dMaster = dynamic_cast<MasterDof *>( d );

            printf("%u @ (%f, %f, %f) DofIDItem %u, ", dman->giveNumber(), dman->giveCoordinates()->at(1), dman->giveCoordinates()->at(2), dman->giveCoordinates()->at(3), m->SurfaceData.at(j)->DofID);

            if (dMaster && dMaster->giveBcId()==0) {
                this->setBoundaryConditionOnDof(dman->giveDofWithID(m->SurfaceData.at(j)->DofID), m->SurfaceData.at(j)->value);
            } else {
                nonFreeDofs.push_back(d);
            }

            printf("\n");

        }

        m->myEngngModel->doStepOutput(m->myEngngModel->giveCurrentStep());

        //Compute correctionfactors
        double Am, Ap, c0, cm, cp;
        //for (int l=0; l<5; l++) {
        computeCorrectionFactors(*m->myEngngModel, &dofs, &Am, &Ap, &c0, &cm, &cp );

        m->ap=-(- Ap*cm*cm + Am*cp*cm + Am*c0*cp)/(Ap*cm*cm + Am*cp*cp);
        m->am=-(- Am*cp*cp + Ap*cm*cp + Ap*c0*cm)/(Ap*cm*cm + Am*cp*cp);

        printf("am=%f, ap=%f, error in incompressibility=%10.10f\n", m->am, m->ap, c0+cm*m->am+cp*m->ap);

        //Update values with correction factors
        for (size_t j=0; j<m->SurfaceData.size(); j++) {
            SurfaceDataStruct *sd = m->SurfaceData.at(j);
            DofManager *dman=sd->DofMan;
            Dof *d=dman->giveDofWithID(sd->DofID);

            double factor=1.0;
            factor = m->SurfaceData.at(j)->isPlus ? m->ap : factor;
            factor = m->SurfaceData.at(j)->isMinus ? m->am : factor;
            factor = m->SurfaceData.at(j)->isZeroBoundary ? 1.0 : factor;

            //factor=1.0;

            //printf("%u @ (%f, %f, %f): isPlus:%u isMinus:%u isZero:%u Factor:%f\n", dman->giveNumber(), dman->giveCoordinates()->at(1),dman->giveCoordinates()->at(2),dman->giveCoordinates()->at(3), m->SurfaceData.at(j)->isPlus,  m->SurfaceData.at(j)->isMinus, m->SurfaceData.at(j)->isZeroBoundary, factor);

            bool isFree=true;
            for (size_t k=0; k<nonFreeDofs.size(); k++) {
                if (nonFreeDofs.at(k) == d) {
                    isFree=false;
                    break;
                }
            }

            if (dynamic_cast<MasterDof *>( d ) && isFree) {
                double value = m->SurfaceData.at(j)->value*factor;
                this->setBoundaryConditionOnDof(dman->giveDofWithID(m->SurfaceData.at(j)->DofID), m->SurfaceData.at(j)->value*factor);
                double value2 = giveValueAtPoint(*d->giveDofManager()->giveCoordinates(), d->giveDofID(), *m->myEngngModel);
                //printf("j=%u, DofManagerID=%u, value1-value2=%15.15f\n", j, d->giveDofManager()->giveNumber(), value-value2);
                //if (d->giveDofID()==V_u) {value=0.01;} else {value=0.0;};
                //this->setBoundaryConditionOnDof(dman->giveDofWithID(m->SurfaceData.at(j)->DofID), value);
            }
        }
        //}


        computeCorrectionFactors(*m->myEngngModel, &dofs, &Am, &Ap, &c0, &cm, &cp );

        m->ap=-(- Ap*cm*cm + Am*cp*cm + Am*c0*cp)/(Ap*cm*cm + Am*cp*cp);
        m->am=-(- Am*cp*cp + Ap*cm*cp + Ap*c0*cm)/(Ap*cm*cm + Am*cp*cp);

        printf("am=%f, ap=%f, error in incompressibility=%10.10f\n", m->am, m->ap, c0+cm*m->am+cp*m->ap);

        m->myEngngModel->doStepOutput(m->myEngngModel->giveCurrentStep());

        //m->myEngngModel->initializeYourself(m->myEngngModel->giveCurrentStep());
        //m->myEngngModel->forceEquationNumbering();
        //m->myEngngModel->solveYourselfAt(m->myEngngModel->giveCurrentStep());
        //m->myEngngModel->doStepOutput(m->myEngngModel->giveCurrentStep());

        modes.push_back(mode);

        OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", myEngngModel);

    }
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
SolutionbasedShapeFunction :: computeBaseFunctionValueAt(FloatArray &coords, DofIDItem dofID, EngngModel &myEngngModel)
{

    if (useConstantBase) {
        return 1.0;
    } else {
        std :: vector <FloatArray *> checkcoords;
        std :: vector <int> permuteIndex;
        double outvalue=0.0;
        int n=0;

        if (!isLoaded) loadProblem();

        myEngngModel.giveDomain(1)->giveSpatialLocalizer()->init(false);

        if (this->giveDomain()->giveNumberOfSpatialDimensions()==2) {
            coords.resize(2);
        }

        // Determine if current coordinate is at a max or min point and if so, which type (on a surface, edge or a corner?)
        // n is the number of dimensions at which the point is an extremum, permuteIndex tells in which dimension the coordinate is max/min

        int thisMask=0;

        for (int i=1; i<=coords.giveSize(); i++) {
            if ( ( fabs(maxCoord.at(i)-coords.at(i))<TOL) || (fabs(minCoord.at(i)-coords.at(i))<TOL) ){
                permuteIndex.push_back(i);
                n++;
                thisMask = thisMask + pow(2.0,i-1);
            }
        }

        for (int i=0; i<pow(2.0,n); i++) {

            int mask=i, counter=1;
            FloatArray *newCoord = new (FloatArray)(coords.giveSize());
            *newCoord = coords;

            for (int j=1; j<=n; j++) {
                double d=0.0;//TOL;
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

        // The followind define allows for use of weakly periodic bc to be copied. This does not comply with the theory but is used to check the validity of the code.
#define USEWPBC 0
#if USEWPBC == 1
        FloatArray *tempCoord = new FloatArray;
        *tempCoord = *coords;
        checkcoords.clear();
        checkcoords.push_back(tempCoord);
#endif

        for (size_t i=0; i<checkcoords.size(); i++) {

            double value = giveValueAtPoint(*checkcoords.at(i), dofID, myEngngModel);
#if USEWPBC == 1
            outvalue = value;
#else
            outvalue = outvalue + value / ( (double) pow(2.0, n));
#endif
        }
        return outvalue;
    }
}

double
SolutionbasedShapeFunction :: giveValueAtPoint(const FloatArray &coords, DofIDItem dofID, EngngModel &myEngngModel)
{

    FloatArray checkCoords;
    checkCoords.resize(coords.giveSize());

    for (int j=1; j<=checkCoords.giveSize(); j++) {
        checkCoords.at(j)=coords.at(j);
        if (fabs(checkCoords.at(j)-minCoord.at(j))<=TOL) {
            //            checkCoords.at(j)=minCoord.at(j);//+TOL;
        }
        if (fabs(checkCoords.at(j)-maxCoord.at(j))<=TOL) {
            //            checkCoords.at(j)=maxCoord.at(j);//-TOL;
        }
    }

    // printf("Find value at (%f, %f, %f) for DofID %u\n", checkCoords.at(1), checkCoords.at(2), checkCoords.at(3), dofID);

    FloatArray closest, lcoords, values;

    Element *elementAtCoords = myEngngModel.giveDomain(1)->giveSpatialLocalizer()->giveElementContainingPoint (checkCoords);

    if (elementAtCoords==NULL) {
        elementAtCoords = myEngngModel.giveDomain(1)->giveSpatialLocalizer()->giveElementClosestToPoint(lcoords, closest, checkCoords, 1);
        if (elementAtCoords==NULL) {
            printf("Cannot find element closest to point\n");
            checkCoords.pY();
        }
    }

    EIPrimaryUnknownMapperInterface *em = dynamic_cast<EIPrimaryUnknownMapperInterface*> ( elementAtCoords->giveInterface ( EIPrimaryUnknownMapperInterfaceType ) );

    IntArray dofids;

    em->EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(dofids);
    em->EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(VM_Total, thisTimestep, checkCoords, values);

    return values.at(dofids.findFirstIndexOf(dofID));
}

void
SolutionbasedShapeFunction :: setBoundaryConditionOnDof(Dof *d, double value)
{
    int bcID = d->giveBcId();

    if (bcID==0) {
        DynamicInputRecord ir;
        ir.setRecordKeywordField("boundarycondition", 1);
        ir.setField(1, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
        ir.setField(value, _IFT_BoundaryCondition_PrescribedValue);

        bcID = d->giveDofManager()->giveDomain()->giveNumberOfBoundaryConditions()+1;

        GeneralBoundaryCondition *myBC;
        myBC = classFactory.createBoundaryCondition("boundarycondition", bcID, d->giveDofManager()->giveDomain());
        myBC->initializeFrom(&ir);
        d->giveDofManager()->giveDomain()->setBoundaryCondition(bcID, myBC);

        d->setBcId(bcID);
    } else {
        BoundaryCondition *bc;
        bc=(BoundaryCondition *) d->giveDofManager()->giveDomain()->giveBc(bcID);
        bc->setPrescribedValue(value);
    }
}

}
