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
        if (dynamic_cast<IntArray*>(&myArray)) {
            printf("%u ", myArray.at(xyz));
        } else {
            printf("%f ", myArray.at(xyz));
        }
    }
}

template <class T>
void logDataMsg(const char *c, T myArray) {
    printf("%s ", c);
    logData(myArray);
    printf("\n");
}

template <class T>
void logDataMsg(const char *c, T myArray, const char *c2) {
    printf("%s ", c);
    logData(myArray);
    printf("%s\n", c2);
}

REGISTER_BoundaryCondition( SolutionbasedShapeFunction );

SolutionbasedShapeFunction::SolutionbasedShapeFunction(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    isLoaded=false;
    order=8;
    TOL=1e-5;
    // TODO Auto-generated constructor stub

}

SolutionbasedShapeFunction::~SolutionbasedShapeFunction()
{

}
void
SolutionbasedShapeFunction :: checkIncompressibility(EngngModel &myEngngModel)
{

    EngngModel *m=&myEngngModel;
    Set *mySet=m->giveDomain(1)->giveSet( externalSet );// this->domain->giveSet( this->giveSetNumber() );
    IntArray BoundaryList = mySet->giveBoundaryList();
    double NetInflow=0.0;

    for (int i=0; i<BoundaryList.giveSize()/2; i++) {
        int ElementID=BoundaryList(2*i);
        int Boundary=BoundaryList(2*i+1);
        Element *e = m->giveDomain(1)->giveElement(ElementID);
        FEInterpolation *interp=e->giveInterpolation();

        GaussIntegrationRule irule(order, e);
        int ngp = irule.getRequiredNumberOfIntegrationPoints(_Triangle, order);
        irule.setUpIntegrationPoints(_Triangle, ngp, _Unknown);


        IntArray bnodes;
        FloatMatrix values;
        interp->boundaryGiveNodes(bnodes, Boundary);
        values.resize(3, bnodes.giveSize());

        for (int j=1; j<=bnodes.giveSize(); j++) {
            for (int k=1; k<=3; k++) {
                values.at(k, j)=e->giveNode(bnodes.at(j))->giveDof(k)->giveUnknown(VM_Total, m->giveCurrentStep());
            }
        }

        for (int j=0; j<irule.giveNumberOfIntegrationPoints(); j++) {
            GaussPoint *gp=irule.getIntegrationPoint(j);
            FloatArray *lcoords = gp->giveCoordinates();
            FloatArray N, normal, v;
            double detJ=interp->boundaryGiveTransformationJacobian(Boundary, *lcoords, FEIElementGeometryWrapper(e))*gp->giveWeight();
            interp->boundaryEvalNormal(normal, Boundary, *lcoords, FEIElementGeometryWrapper(e));
            interp->boundaryEvalN(N, Boundary, *lcoords, FEIElementGeometryWrapper(e));

            v.beProductOf(values, N);
            //v.printYourself();
            NetInflow=+v.dotProduct(normal)*detJ;
            //printf("Net inflow @ (%f, %f):%10.10f\n", gp->giveCoordinates()->at(1), gp->giveCoordinates()->at(2), NetInflow);

        }
    }

    printf("Net inflow:%15.15f\n", NetInflow);

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

    externalSet=-1;
    IR_GIVE_OPTIONAL_FIELD(ir, externalSet, _IFT_SolutionbasedShapeFunction_Externalset);

    // Set up master dofs
    myNode = new Node(0, this->giveDomain());

    for (int i=1; i<=this->giveDofIDs().giveSize(); i++) {
        MasterDof *newDof = new MasterDof( 0, myNode, (DofIDItem)this->domain->giveNextFreeDofID() );
        myNode->appendDof( newDof );
        if (i==1) {
            //setBoundaryConditionOnDof(newDof, 1.0);
        } else {
            //setBoundaryConditionOnDof(newDof, 0.0);
        }
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
SolutionbasedShapeFunction :: computeCorrectionFactors(EngngModel &myEngngModel, IntArray *Dofs, long double *Am, long double *Ap, long double *c0, long double *cm, long double *cp)
{
    /*
     *Compute c0, cp, cm, Ap and Am
     */

    *Am=0.0;
    *Ap=0.0;
    *c0=0.0;
    *cm=0.0;
    *cp=0.0;
    IntArray surfaceNodes;
    long double NetInflow1=0.0, NetInflow2=0.0;
    long double d0=0.0, dm=0.0, dp=0.0;
    IntArray pList, mList, zList;

    //modeStruct *mode=this->modes.at(EngngModelID);
    EngngModel *m=&myEngngModel;
    Set *mySet=m->giveDomain(1)->giveSet( externalSet );// this->domain->giveSet( this->giveSetNumber() );
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
            surfaceNodes.insertOnce(dman->giveNumber());

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
                zList.insertOnce(dman->giveNumber());
            } else if (isPlus) {
                pNodes.insertSorted(j);
                pList.insertOnce(dman->giveNumber());
            } else if (isMinus) {
                mNodes.insertSorted(j);
                mList.insertOnce(dman->giveNumber());
            }
        }

        GaussIntegrationRule iRule(order, thisElement);

        int n = iRule.getRequiredNumberOfIntegrationPoints(_Triangle, order);
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


            *c0=*c0 + (long double) zPhi.dotProduct(normal, 3)*detJ;
            *cp=*cp + (long double) pPhi.dotProduct(normal, 3)*detJ;
            *cm=*cm + (long double) mPhi.dotProduct(normal, 3)*detJ;
            d0=d0 + (long double) zPhi.dotProduct(normal, 3)*detJ;
            dm=dm + (long double) mPhi.dotProduct(normal, 3)*detJ;
            dp=dp + (long double) pPhi.dotProduct(normal, 3)*detJ;

            NetInflow1=+ (zPhi.dotProduct(normal, 3) + pPhi.dotProduct(normal, 3) + mPhi.dotProduct(normal, 3))*detJ;
            NetInflow2=+d0+dp+dm;

            //printf("Net inflow @ (%f, %f):%10.10f *\n", gp->giveCoordinates()->at(1), gp->giveCoordinates()->at(2), NetInflow);

            *Ap=*Ap+pPhi.dotProduct(pPhi, 3)*detJ;
            *Am=*Am+mPhi.dotProduct(mPhi, 3)*detJ;

        }
    }

    printf("c0:%Lf, cp:%Lf, cm:%Lf, Ap:%Lf, Am:%Lf\n", *c0, *cp, *cm, *Ap, *Am);
    printf("Net inflow1=%Lf *\n", NetInflow1);
    printf("Net inflow2=%Lf *\n", NetInflow2);
    printf("c0+cp+cm=%Lf\n", *c0+*cp+*cm);
    logDataMsg("pList=[", pList, "]");
    logDataMsg("mList=[", mList, "]");
    logDataMsg("zList=[", zList, "]");
    /*for (int i=1; i<=surfaceNodes.giveSize(); i++) {
        printf("%u ", surfaceNodes.at(i));
    }
    printf("\n");*/

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

    for (int i=0; i<this->domain->giveNumberOfSpatialDimensions(); i++) {
        //    for (int i=2; i>=0; i--) {

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
        isLoaded=true;

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

        // Set unknowns to the mean value of opposite sides of the domain.
        // Loop thru all nodes and compute phi for all degrees of freedom on the boundary. Save phi in a list for later use.
        initializeSurfaceData(mode);

        mode->myEngngModel->doStepOutput(mode->myEngngModel->giveCurrentStep());

        //Project mode onto descretized surface
        for (size_t j=0; j<mode->SurfaceData.size(); j++) {

            DofManager *dman=mode->SurfaceData.at(j)->DofMan;
            Dof *d=dman->giveDofWithID(mode->SurfaceData.at(j)->DofID);

            //Add BC if dof is a master and has no existing boundary conditions
            MasterDof *dMaster = dynamic_cast<MasterDof *>( d );

            //printf("%u @ (%f, %f, %f) DofIDItem %u, ", dman->giveNumber(), dman->giveCoordinates()->at(1), dman->giveCoordinates()->at(2), dman->giveCoordinates()->at(3), m->SurfaceData.at(j)->DofID);

            mode->SurfaceData.at(j)->isFree=dMaster && dMaster->giveBcId()==0;

            if (mode->SurfaceData.at(j)->isFree) {
                this->setBoundaryConditionOnDof(dman->giveDofWithID(mode->SurfaceData.at(j)->DofID), mode->SurfaceData.at(j)->value);
            }

            //printf("\n");

        }

        mode->myEngngModel->doStepOutput(mode->myEngngModel->giveCurrentStep());

        //Compute correctionfactors
        long double Am, Ap, c0, cm, cp;
        double ap, am;
        checkIncompressibility(*mode->myEngngModel);

        computeCorrectionFactors(*mode->myEngngModel, &dofs, &Am, &Ap, &c0, &cm, &cp );

        ap=-(- Ap*cm*cm + Am*cp*cm + Am*c0*cp)/(Ap*cm*cm + Am*cp*cp);
        am=-(- Am*cp*cp + Ap*cm*cp + Ap*c0*cm)/(Ap*cm*cm + Am*cp*cp);

        mode->ap=ap;
        mode->am=am;

        printf("am=%f, ap=%f, error in incompressibility=%Lf\n", mode->am, mode->ap, c0+cm+cp);

        updateModelWithFactors(mode);

        /*        // Perturbation in ap
        mode->ap=ap+1e-6;
        mode->am=am;
        updateModelWithFactors(mode);
        checkIncompressibility(*mode->myEngngModel);

        // Perturbation in am
        mode->ap=ap;
        mode->am=am-1e-6;
        updateModelWithFactors(mode);
        checkIncompressibility(*mode->myEngngModel);*/

        mode->myEngngModel->doStepOutput(mode->myEngngModel->giveCurrentStep());

        modes.push_back(mode);

        OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", myEngngModel);

    }
}

void
SolutionbasedShapeFunction :: updateModelWithFactors(modeStruct *m)
{
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

        if (dynamic_cast<MasterDof *>( d ) && m->SurfaceData.at(j)->isFree) {
            this->setBoundaryConditionOnDof(dman->giveDofWithID(m->SurfaceData.at(j)->DofID), m->SurfaceData.at(j)->value*factor);
        }
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
        *tempCoord = coords;
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

void
SolutionbasedShapeFunction :: initializeSurfaceData(modeStruct *mode)
{

    EngngModel *m=mode->myEngngModel;
    double TOL2=1e-5;
    IntArray pNodes, mNodes, zNodes;

    Set *mySet=this->domain->giveSet( this->giveSetNumber() );
    /*    IntArray NodeList = mySet->giveNodeList();

    printf("Nodes=[");
    for (int i=1; i<=NodeList.giveSize(); i++) {
        printf("%u, ", NodeList.at(i));
    }
    printf("\n");
*/
    IntArray BoundaryList = mySet->giveBoundaryList();

    // First add all nodes to pNodes or nNodes respectively depending on coordinate and normal.

    for (int i=0; i<BoundaryList.giveSize()/2; i++) {

        int ElementID=BoundaryList(2*i);
        int Boundary=BoundaryList(2*i+1);

        Element *e=m->giveDomain(1)->giveElement(ElementID);
        FEInterpolation *geoInterpolation = e->giveInterpolation();

        // Check all sides of element
        IntArray bnodes;

#define usePoints 1
#if usePoints == 1
        // Check if all nodes are on the boundary
        geoInterpolation->boundaryGiveNodes(bnodes, Boundary);
        for (int k=1; k<=bnodes.giveSize(); k++) {
            DofManager *dman=e->giveDofManager(bnodes.at(k));
            for (int l=1; l<=dman->giveCoordinates()->giveSize(); l++) {
                if (fabs(dman->giveCoordinates()->at(l)-maxCoord.at(l))<TOL2) pNodes.insertOnce(dman->giveNumber());
                if (fabs(dman->giveCoordinates()->at(l)-minCoord.at(l))<TOL2) mNodes.insertOnce(dman->giveNumber());
            }
        }
#else
        // Check normal
        FloatArray lcoords;
        lcoords.resize(2);
        lcoords.at(1)=0.33333;
        lcoords.at(2)=0.33333;

        FloatArray normal;
        geoInterpolation->boundaryEvalNormal(normal, j, lcoords, FEIElementGeometryWrapper(e) );
        geoInterpolation->boundaryGiveNodes(bnodes, j);

        printf("i=%u\tj=%u\t(%f\t%f\t%f)\n", i, j, normal.at(1), normal.at(2), normal.at(3));
        for (int k=1; k<=normal.giveSize(); k++) {

            if (fabs((fabs(normal.at(k))-1))<1e-4) { // Points in x, y or z direction

                addTo=NULL;
                if (normal.at(k)> 0.5) addTo = &pNodes;
                if (normal.at(k)<-0.5) addTo = &mNodes;
                if (addTo!=NULL) {
                    for (int l=1; l<=bnodes.giveSize(); l++) {
                        bool isSurface=false;
                        DofManager *dman=e->giveDofManager(bnodes.at(l));
                        dman->giveCoordinates()->printYourself();
                        for (int m=1; m<=dman->giveCoordinates()->giveSize(); m++) {
                            if ((fabs(dman->giveCoordinates()->at(m)-maxCoord.at(m))<TOL2) || (fabs(dman->giveCoordinates()->at(m)-minCoord.at(m))<TOL2)) {
                                isSurface=true;
                            }
                        }

                        if (isSurface) addTo->insertOnce( e->giveDofManagerNumber(bnodes.at(l)) );
                    }
                }
            }
        }
#endif
        //}
    }

#if 0
    printf("p=[");
    for (int i=1; i<pNodes.giveSize(); i++) {
        printf("%u, ", pNodes.at(i));
    }
    printf("];\n");
    printf("m=[");
    for (int i=1; i<mNodes.giveSize(); i++) {
        printf("%u, ", mNodes.at(i));
    }
    printf("];\n");
#endif
    //The intersection of pNodes and mNodes constitutes zNodes
    {
        int i=1, j=1;
        while (i<=pNodes.giveSize()) {
            j=1;
            while (j<=mNodes.giveSize() && (i<=pNodes.giveSize()) ) {
                //printf("%u == %u?\n", pNodes.at(i), mNodes.at(j));
                if (pNodes.at(i)==mNodes.at(j)) {
                    zNodes.insertOnce(pNodes.at(i));
                    pNodes.erase(i);
                    mNodes.erase(j);
                } else {
                    j++;
                }
            }
            i++;
        }
    }

    // Compute base function values on nodes for dofids
    for (int j=1; j<=dofs.giveSize(); j++) {
        for (int i=1; i<=pNodes.giveSize(); i++) {
            SurfaceDataStruct *surfaceData = new (SurfaceDataStruct);
            copyDofManagerToSurfaceData(surfaceData, pNodes.at(i), dofs.at(j), mode->myEngngModel, true, false, false);
            mode->SurfaceData.push_back(surfaceData);
        }

        for (int i=1; i<=mNodes.giveSize(); i++) {
            SurfaceDataStruct *surfaceData = new (SurfaceDataStruct);
            copyDofManagerToSurfaceData(surfaceData, mNodes.at(i), dofs.at(j), mode->myEngngModel, false, true, false);
            mode->SurfaceData.push_back(surfaceData);
        }

        for (int i=1; i<=zNodes.giveSize(); i++) {
            SurfaceDataStruct *surfaceData = new (SurfaceDataStruct);
            copyDofManagerToSurfaceData(surfaceData, zNodes.at(i), dofs.at(j), mode->myEngngModel, false, false, true);
            mode->SurfaceData.push_back(surfaceData);
        }
    }


#if 0
    printf("p2=[");
    for (int i=1; i<=pNodes.giveSize(); i++) {
        printf("%u, ", pNodes.at(i));
    }
    printf("];\n");
    printf("m2=[");
    for (int i=1; i<=mNodes.giveSize(); i++) {
        printf("%u, ", mNodes.at(i));
    }
    printf("];\n");
    printf("z2=[");
    for (int i=1; i<=zNodes.giveSize(); i++) {
        printf("%u, ", zNodes.at(i));
    }
    printf("];\n");

    printf("pCoords=[");
    for (int i=1; i<=pNodes.giveSize(); i++) {
        FloatArray *coords = m->giveDomain(1)->giveDofManager(pNodes.at(i))->giveCoordinates();
        printf("%f, %f, %f; ", coords->at(1), coords->at(2), coords->at(3));
    }
    printf("]\n");
    printf("mCoords=[");
    for (int i=1; i<=mNodes.giveSize(); i++) {
        FloatArray *coords = m->giveDomain(1)->giveDofManager(mNodes.at(i))->giveCoordinates();
        printf("%f, %f, %f; ", coords->at(1), coords->at(2), coords->at(3));
    }
    printf("]\n");
    printf("zCoords=[");
    for (int i=1; i<=zNodes.giveSize(); i++) {
        FloatArray *coords = m->giveDomain(1)->giveDofManager(zNodes.at(i))->giveCoordinates();
        printf("%f, %f, %f; ", coords->at(1), coords->at(2), coords->at(3));
    }
    printf("];\n");
#endif
}

void
SolutionbasedShapeFunction :: copyDofManagerToSurfaceData(SurfaceDataStruct *surfaceData, int DofManID, int DofType, EngngModel *m, bool isPlus, bool isMinus, bool isZeroBoundary)
{

    surfaceData->DofID = (DofIDItem) DofType;
    surfaceData->DofMan = m->giveDomain(1)->giveDofManager(DofManID);
    surfaceData->isPlus = isPlus;
    surfaceData->isMinus = isMinus;
    surfaceData->isZeroBoundary = isZeroBoundary;
    //    if ((surfaceData->DofMan->giveNumber()==6228) || (surfaceData->DofMan->giveNumber()==5532)) {
    //        printf("DofManager number:%u\n", surfaceData->DofMan->giveNumber());
    //        surfaceData->DofMan->giveCoordinates()->printYourself();
    //    }
    surfaceData->value = computeBaseFunctionValueAt(*surfaceData->DofMan->giveCoordinates(), (DofIDItem) DofType, *m);

}
}
