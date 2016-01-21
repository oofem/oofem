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
#include "dynamicinputrecord.h"
#include "bodyload.h"
#include "boundarycondition.h"
//#include "stokesflow.h"
//#include "classfactory.h"

#include <vector>

namespace oofem {
template< class T >
void logData(T myArray) {
    for ( int xyz = 1; xyz <= myArray.giveSize(); xyz++ ) {
        if ( dynamic_cast< IntArray * >(& myArray) ) {
            printf( "%u ", myArray.at(xyz) );
        } else {
            printf( "%f ", myArray.at(xyz) );
        }
    }
}

template< class T >
void logDataMsg(const char *c, T myArray) {
    printf("%s ", c);
    logData(myArray);
    printf("\n");
}

template< class T >
void logDataMsg(const char *c, T myArray, const char *c2) {
    printf("%s ", c);
    logData(myArray);
    printf("%s\n", c2);
}

REGISTER_BoundaryCondition(SolutionbasedShapeFunction)

SolutionbasedShapeFunction :: SolutionbasedShapeFunction(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    isLoaded = false;
    order = 8;
    TOL = 1e-5;
    bigNorm = 0.0;

    // TODO Auto-generated constructor stub
}

SolutionbasedShapeFunction :: ~SolutionbasedShapeFunction()
{
    // TODO Auto-generated destructor stub
}

IRResultType
SolutionbasedShapeFunction :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    // Load problem file
    this->filename = "";
    IR_GIVE_OPTIONAL_FIELD(ir, this->filename, _IFT_SolutionbasedShapeFunction_ShapeFunctionFile);
    useConstantBase = ( this->filename == "" ) ? true : false;

    externalSet = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, externalSet, _IFT_SolutionbasedShapeFunction_Externalset);

    // use correction factors to ensure incompressibility
    useCorrectionFactors=false;
    IR_GIVE_OPTIONAL_FIELD(ir, useCorrectionFactors, _IFT_SolutionbasedShapeFunction_UseCorrectionFactors);

    dumpSnapshot = false;
    IR_GIVE_OPTIONAL_FIELD(ir, dumpSnapshot, _IFT_SolutionbasedShapeFunction_DumpSnapshots);


    // Set up master dofs
    ///@todo This should be in the constructor:
    myNode = new Node( 1, this->giveDomain() );

    for (int i = 1; i <= this->giveDomain()->giveNumberOfSpatialDimensions(); i++) {
        int DofID = this->domain->giveNextFreeDofID();
        MasterDof *newDof = new MasterDof( myNode, (DofIDItem) DofID );
        myNode->appendDof( newDof );
    }

    init();

    return ActiveBoundaryCondition :: initializeFrom(ir);
}

DofManager *
SolutionbasedShapeFunction :: giveInternalDofManager(int i)
{
    return myNode;
}

bool
SolutionbasedShapeFunction :: isCoeff(ActiveDof *dof)
{
    for ( Dof *myDof: *myNode ) {
        if ( dof == myDof ) {
            return true;
        }
    }
    return false;
}

void
SolutionbasedShapeFunction :: init()
{
    Node *n1 = this->giveDomain()->giveNode(1);

    maxCoord = * n1->giveCoordinates();
    minCoord = * n1->giveCoordinates();

    for ( auto &n :this->giveDomain()->giveDofManagers() ) {
        for ( int j = 1; j <= maxCoord.giveSize(); j++ ) {
            maxCoord.at(j) = max( n->giveCoordinate(j), maxCoord.at(j) );
            minCoord.at(j) = min( n->giveCoordinate(j), minCoord.at(j) );
        }
    }
}

void
SolutionbasedShapeFunction :: computeCorrectionFactors(modeStruct &myMode, IntArray *Dofs, double *am, double *ap)
{
    /*
     * *Compute c0, cp, cm, Bp, Bm, Ap and Am
     */

    double A0p = 0.0, App = 0.0, A0m = 0.0, Amm = 0.0, Bp = 0.0, Bm = 0.0, c0 = 0.0, cp = 0.0, cm = 0.0;

    EngngModel *model = myMode.myEngngModel;
    Set *mySet = model->giveDomain(1)->giveSet(externalSet);

    IntArray BoundaryList = mySet->giveBoundaryList();

    for ( int i = 0; i < BoundaryList.giveSize() / 2; i++ ) {
        int ElementID = BoundaryList(2 * i);
        int Boundary = BoundaryList(2 * i + 1);

        Element *thisElement = model->giveDomain(1)->giveElement(ElementID);
        FEInterpolation *geoInterpolation = thisElement->giveInterpolation();
        IntArray bnodes, zNodes, pNodes, mNodes;
        FloatMatrix nodeValues;

        geoInterpolation->boundaryGiveNodes(bnodes, Boundary);

        nodeValues.resize( this->dofs.giveSize(), bnodes.giveSize() );
        nodeValues.zero();

        // Change to global ID for bnodes and identify the intersection of bnodes and the zero boundary
        splitBoundaryNodeIDs(myMode, * thisElement, bnodes, pNodes, mNodes, zNodes, nodeValues);

        std :: unique_ptr< IntegrationRule >iRule(geoInterpolation->giveBoundaryIntegrationRule(order, Boundary));

        for ( GaussPoint *gp: *iRule ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            FloatArray gcoords, normal, N;
            FloatArray Phi;

            double detJ = fabs( geoInterpolation->boundaryGiveTransformationJacobian( Boundary, lcoords, FEIElementGeometryWrapper(thisElement) ) ) * gp->giveWeight();

            geoInterpolation->boundaryEvalNormal( normal, Boundary, lcoords, FEIElementGeometryWrapper(thisElement) );
            geoInterpolation->boundaryEvalN( N, Boundary, lcoords, FEIElementGeometryWrapper(thisElement) );
            geoInterpolation->boundaryLocal2Global( gcoords, Boundary, lcoords, FEIElementGeometryWrapper(thisElement) );

            FloatArray pPhi, mPhi, zPhi;
            pPhi.resize( Dofs->giveSize() );
            pPhi.zero();
            mPhi.resize( Dofs->giveSize() );
            mPhi.zero();
            zPhi.resize( Dofs->giveSize() );
            zPhi.zero();

            // Build phi (analytical averaging, not projected onto the mesh)
            computeBaseFunctionValueAt(Phi, gcoords, * Dofs, * myMode.myEngngModel);

            // Build zPhi for this DofID
            for ( int l = 1; l <= zNodes.giveSize(); l++ ) {
                int nodeID = zNodes.at(l);
                for ( int m = 1; m <= this->dofs.giveSize(); m++ ) {
                    zPhi.at(m) = zPhi.at(m) + N.at(nodeID) * nodeValues.at(m, nodeID);
                }
            }


            // Build pPhi for this DofID
            for ( int l = 1; l <= pNodes.giveSize(); l++ ) {
                int nodeID = pNodes.at(l);
                for ( int m = 1; m <= this->dofs.giveSize(); m++ ) {
                    pPhi.at(m) = pPhi.at(m) + N.at(nodeID) * nodeValues.at(m, nodeID);
                }
            }

            // Build mPhi for this DofID
            for ( int l = 1; l <= mNodes.giveSize(); l++ ) {
                int nodeID = mNodes.at(l);
                for ( int m = 1; m <= this->dofs.giveSize(); m++ ) {
                    mPhi.at(m) = mPhi.at(m) + N.at(nodeID) * nodeValues.at(m, nodeID);
                }
            }

            c0 = c0 + zPhi.dotProduct(normal, 3) * detJ;
            cp = cp + pPhi.dotProduct(normal, 3) * detJ;
            cm = cm + mPhi.dotProduct(normal, 3) * detJ;

            App = App + pPhi.dotProduct(pPhi, 3) * detJ;
            Amm = Amm + mPhi.dotProduct(mPhi, 3) * detJ;
            A0p = A0p + zPhi.dotProduct(pPhi, 3) * detJ;
            A0m = A0m + zPhi.dotProduct(mPhi, 3) * detJ;

            Bp = Bp + Phi.dotProduct(pPhi, 3) * detJ;
            Bm = Bm + Phi.dotProduct(mPhi, 3) * detJ;
        }
    }

    * am = -( A0m * cp * cp - Bm * cp * cp - A0p * cm * cp + App * c0 * cm + Bp * cm * cp ) / ( App * cm * cm + Amm * cp * cp );
    * ap = -( A0p * cm * cm - Bp * cm * cm - A0m * cm * cp + Amm * c0 * cp + Bm * cm * cp ) / ( App * cm * cm + Amm * cp * cp );
}

void
SolutionbasedShapeFunction :: splitBoundaryNodeIDs(modeStruct &mode, Element &e, IntArray &bnodes, IntArray &pList, IntArray &mList, IntArray &zList, FloatMatrix &nodeValues)
{
    pList.clear();
    mList.clear();
    zList.clear();

    for ( int j = 1; j <= bnodes.giveSize(); j++ ) {
        DofManager *dman = e.giveDofManager( bnodes.at(j) );

        bool isZero = false;
        bool isPlus = false;
        bool isMinus = false;

        whichBoundary(* dman->giveCoordinates(), isPlus, isMinus, isZero);

        if ( isZero ) {
            zList.insertSorted(j);
        } else if ( isPlus ) {
            pList.insertSorted(j);
        } else if ( isMinus ) {
            mList.insertSorted(j);
        }

        // Find global DofManager and fetch nodal values
        for ( size_t k = 0; k < mode.SurfaceData.size(); k++ ) {
            if ( mode.SurfaceData.at(k)->DofMan == dman ) {
                int IndexOfDofIDItem = 0;
                for ( int l = 1; l <= dofs.giveSize(); l++ ) {
                    if ( dofs.at(l) == mode.SurfaceData.at(k)->DofID ) {
                        IndexOfDofIDItem = l;
                        break;
                    }
                }
                nodeValues.at(IndexOfDofIDItem, j) = mode.SurfaceData.at(k)->value;
            }
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
    myNode->giveUnknownVector(gamma, myDofIDs, mode, tStep);  // alpha1, alpha2,...

    double out = shapeFunctionValues.dotProduct(gamma);

    return out;
}

void
SolutionbasedShapeFunction :: computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
{
    if ( !isLoaded ) {
        loadProblem();
    }

    FloatArray values, masterContribs2, d, values2;

    masterContribs.resize( this->giveDomain()->giveNumberOfSpatialDimensions() );
    masterContribs2.resize( this->giveDomain()->giveNumberOfSpatialDimensions() );

    IntArray dofIDs = {dof->giveDofID()};

    bool isPlus, isMinus, isZero, found;
    whichBoundary(* dof->giveDofManager()->giveCoordinates(), isPlus, isMinus, isZero);

    for ( int i = 1; i <= this->giveDomain()->giveNumberOfSpatialDimensions(); i++ ) {
        double factor = 1.0;
        found = false;

        modeStruct *ms = modes.at(i - 1);
        for ( size_t j = 0; j < ms->SurfaceData.size(); j++ ) {
            SurfaceDataStruct *sd = ms->SurfaceData.at(j);
            if ( sd->DofMan->giveNumber() == dof->giveDofManager()->giveNumber() ) {
                if ( sd->DofID == dof->giveDofID() ) {
                    values.resize(1);
                    values.at(1) = sd->value;
                    found = true;
                    break;
                }
            }
        }

        if ( !found ) {
            printf( "%u\n", dof->giveDofManager()->giveNumber() );
            OOFEM_ERROR("Node not found");
        }

        //giveValueAtPoint(values2, *dof->giveDofManager()->giveCoordinates(), dofIDs, *modes.at(i-1)->myEngngModel);
        //printf ("Mode %u, DofManager: %u, DofIDItem %u, value %10.10f\n", i, dof->giveDofManager()->giveNumber(), dof->giveDofID(), values.at(1));

        factor = isPlus  ? modes.at(i - 1)->ap : factor;
        factor = isMinus ? modes.at(i - 1)->am : factor;
        factor = isZero  ? 1.0 : factor;

        masterContribs.at(i) = factor * values.at(1);
    }
}

int
SolutionbasedShapeFunction :: giveNumberOfMasterDofs(ActiveDof *dof)
{
    return this->giveDomain()->giveNumberOfSpatialDimensions();
}

Dof *
SolutionbasedShapeFunction :: giveMasterDof(ActiveDof *dof, int mdof)
{
    return myNode->giveDofWithID(myDofIDs.at(mdof));
}

void
SolutionbasedShapeFunction :: loadProblem()
{
    for ( int i = 0; i < this->domain->giveNumberOfSpatialDimensions(); i++ ) {
        OOFEM_LOG_INFO("************************** Instanciating microproblem from file %s for dimension %u\n", filename.c_str(), i);

        // Set up and solve problem
        OOFEMTXTDataReader drMicro( filename.c_str() );
        EngngModel *myEngngModel = InstanciateProblem(& drMicro, _processor, 0, NULL, false);
        drMicro.finish();
        myEngngModel->checkProblemConsistency();
        myEngngModel->initMetaStepAttributes( myEngngModel->giveMetaStep(1) );
        thisTimestep = myEngngModel->giveNextStep();
        myEngngModel->init();
        this->setLoads(myEngngModel, i + 1);

        // Check
        for ( auto &elem : myEngngModel->giveDomain(1)->giveElements() ) {
            FloatArray centerCoord;
            int vlockCount = 0;
            centerCoord.resize(3);
            centerCoord.zero();

            for ( int k = 1; k <= elem->giveNumberOfDofManagers(); k++ ) {
                DofManager *dman = elem->giveDofManager(k);
                centerCoord.add( * dman->giveCoordinates() );
                for ( Dof *dof: *dman ) {
                    if ( dof->giveBcId() != 0 ) {
                        vlockCount++;
                    }
                }
            }
            if ( vlockCount == 30 ) {
                OOFEM_WARNING("Element over-constrained (%u)! Center coordinate: %f, %f, %f\n", elem->giveNumber(), centerCoord.at(1) / 10, centerCoord.at(2) / 10, centerCoord.at(3) / 10);
            }
        }

        myEngngModel->solveYourselfAt(thisTimestep);
        isLoaded = true;

        // Set correct export filename
        std :: string originalFilename;
        originalFilename = myEngngModel->giveOutputBaseFileName();

        if (i==0) originalFilename = originalFilename + "_X";
        if (i==1) originalFilename = originalFilename + "_Y";
        if (i==2) originalFilename = originalFilename + "_Z";
        if (dumpSnapshot) {
            myEngngModel->letOutputBaseFileNameBe(originalFilename + "_1_Base");
            myEngngModel->doStepOutput(thisTimestep);
        }

        modeStruct *mode = new(modeStruct);
        mode->myEngngModel = myEngngModel;

        // Check elements

        // Set unknowns to the mean value of opposite sides of the domain.
        // Loop thru all nodes and compute phi for all degrees of freedom on the boundary. Save phi in a list for later use.

        initializeSurfaceData(mode);
        // Update with factor
        double am=1.0, ap=1.0;

        if (useCorrectionFactors) {
            computeCorrectionFactors(*mode, &dofs, &am, &ap );
        }

        OOFEM_LOG_INFO("Correction factors: am=%f, ap=%f\n", am, ap);

        mode->ap = ap;
        mode->am = am;

        updateModelWithFactors(mode);

        if (dumpSnapshot) {
            myEngngModel->letOutputBaseFileNameBe(originalFilename + "_2_Updated");
            myEngngModel->doStepOutput(thisTimestep);
        }

        modes.push_back(mode);

        OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", myEngngModel);
    }
}

void
SolutionbasedShapeFunction :: updateModelWithFactors(modeStruct *m)
{
    //Update values with correction factors
    for ( size_t j = 0; j < m->SurfaceData.size(); j++ ) {
        SurfaceDataStruct *sd = m->SurfaceData.at(j);
        DofManager *dman = sd->DofMan;
        Dof *d = dman->giveDofWithID(sd->DofID);

        double factor = 1.0;
        factor = m->SurfaceData.at(j)->isPlus ? m->ap : factor;
        factor = m->SurfaceData.at(j)->isMinus ? m->am : factor;
        factor = m->SurfaceData.at(j)->isZeroBoundary ? 1.0 : factor;

        if ( dynamic_cast< MasterDof * >(d) && m->SurfaceData.at(j)->isFree ) {
            double u = m->SurfaceData.at(j)->value;
            this->setBoundaryConditionOnDof(dman->giveDofWithID(m->SurfaceData.at(j)->DofID), u * factor);
        }
    }
}

void
SolutionbasedShapeFunction :: setLoads(EngngModel *myEngngModel, int d)
{
    DynamicInputRecord ir;
    FloatArray gradP;

    gradP.resize( this->giveDomain()->giveNumberOfSpatialDimensions() );
    gradP.zero();
    gradP.at(d) = 1.0;

    ir.setRecordKeywordField("deadweight", 1);
    ir.setField(gradP, _IFT_Load_components);
    ir.setField(1, _IFT_GeneralBoundaryCondition_timeFunct);

    int bcID = myEngngModel->giveDomain(1)->giveNumberOfBoundaryConditions() + 1;
    GeneralBoundaryCondition *myBodyLoad;
    myBodyLoad = classFactory.createBoundaryCondition( "deadweight", bcID, myEngngModel->giveDomain(1) );
    myBodyLoad->initializeFrom(& ir);
    myEngngModel->giveDomain(1)->setBoundaryCondition(bcID, myBodyLoad);

    for ( auto &elem : myEngngModel->giveDomain(1)->giveElements() ) {
        IntArray *blArray;
        blArray = elem->giveBodyLoadArray();
        blArray->resizeWithValues(blArray->giveSize() + 1);
        blArray->at( blArray->giveSize() ) = bcID;
    }
}

void
SolutionbasedShapeFunction :: computeBaseFunctionValueAt(FloatArray &answer, FloatArray &coords, IntArray &dofIDs, EngngModel &myEngngModel)
{
    answer.resize( dofIDs.giveSize() );
    answer.zero();

    if ( useConstantBase ) {
        for ( int i = 1; i <= answer.giveSize(); i++ ) {
            answer.at(i) = 1;
        }
        ;
    } else {
        std :: vector< FloatArray >checkcoords;
        std :: vector< int >permuteIndex;
        int n = 0;

        if ( !isLoaded ) {
            loadProblem();
        }

        myEngngModel.giveDomain(1)->giveSpatialLocalizer()->init(false);

        if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
            coords.resize(2);
        }

        // Determine if current coordinate is at a max or min point and if so, which type (on a surface, edge or a corner?)
        // n is the number of dimensions at which the point is an extremum, permuteIndex tells in which dimension the coordinate is max/min

        int thisMask = 0;

        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            if ( ( fabs( maxCoord.at(i) - coords.at(i) ) < TOL ) || ( fabs( minCoord.at(i) - coords.at(i) ) < TOL ) ) {
                permuteIndex.push_back(i);
                n++;
                //thisMask = thisMask + pow(2.0, i - 1);   // compiler warning on conversion from double to int
                thisMask = thisMask + ( 0x01 << ( i - 1 ) );
            }
        }
        int _s = 0x01 << n;
        for ( int i = 0; i < _s; i++ ) {
            int mask = i, counter = 1;
            FloatArray newCoord = coords;

            for ( int j = 1; j <= n; j++ ) {
                double d = 0.0; //TOL;
                if ( ( mask & 1 ) == 0 ) { // Max
                    newCoord.at( permuteIndex.at(counter - 1) ) = minCoord.at( permuteIndex.at(counter - 1) ) + d;
                } else { // Min
                    newCoord.at( permuteIndex.at(counter - 1) ) = maxCoord.at( permuteIndex.at(counter - 1) ) - d;
                }
                counter++;
                mask = mask >> 1;
            }
            checkcoords.emplace_back(newCoord);
        }

        // The followind define allows for use of weakly periodic bc to be copied. This does not comply with the theory but is used to check the validity of the code.
#define USEWPBC 0

#if USEWPBC == 1
        FloatArray *tempCoord = new FloatArray;
        * tempCoord = coords;
        checkcoords.clear();
        checkcoords.push_back(tempCoord);
#endif
        FloatArray values;
        for ( size_t i = 0; i < checkcoords.size(); i++ ) {
            giveValueAtPoint(values, checkcoords[i], dofIDs, myEngngModel);
            //printf("Values at (%f, %f, %f) are [%f, %f, %f]\n", checkcoords.at(i)->at(1), checkcoords.at(i)->at(2), checkcoords.at(i)->at(3), values.at(1), values.at(2), values.at(3));
#if USEWPBC == 1
            for ( int j = 1; j <= values.giveSize(); j++ ) {
                answer.at(j) = values.at(j);
            }
#else
            for ( int j = 1; j <= values.giveSize(); j++ ) {
                answer.at(j) += values.at(j) / ( ( double ) pow(2.0, n) );
            }
#endif
        }
    }
}

void
SolutionbasedShapeFunction :: giveValueAtPoint(FloatArray &answer, const FloatArray &coords, IntArray &dofIDs, EngngModel &myEngngModel)
{
    answer.resize( dofIDs.giveSize() );

    FloatArray closest, lcoords, values;

    Element *elementAtCoords = myEngngModel.giveDomain(1)->giveSpatialLocalizer()->giveElementClosestToPoint(lcoords, closest, coords, 1);
    if ( elementAtCoords == NULL ) {
        OOFEM_WARNING("Cannot find element closest to point");
        coords.pY();
        return;
    }

    IntArray eldofids;

    elementAtCoords->giveElementDofIDMask(eldofids);
    elementAtCoords->computeField(VM_Total, thisTimestep, lcoords, values);

    for ( int i = 1; i <= dofIDs.giveSize(); i++ ) {
        for ( int j = 1; j <= eldofids.giveSize(); j++ ) {
            if ( dofIDs.at(i) == eldofids.at(j) ) {
                answer.at(i) = values.at(j);
                break;
            }
        }
    }
}

void
SolutionbasedShapeFunction :: setBoundaryConditionOnDof(Dof *d, double value)
{
    int bcID = d->giveBcId();

    if ( bcID == 0 ) {
        DynamicInputRecord ir;
        ir.setRecordKeywordField("boundarycondition", 1);
        ir.setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
        ir.setField(value, _IFT_BoundaryCondition_PrescribedValue);

        bcID = d->giveDofManager()->giveDomain()->giveNumberOfBoundaryConditions() + 1;

        GeneralBoundaryCondition *myBC;
        myBC = classFactory.createBoundaryCondition( "boundarycondition", bcID, d->giveDofManager()->giveDomain() );
        myBC->initializeFrom(& ir);
        d->giveDofManager()->giveDomain()->setBoundaryCondition(bcID, myBC);

        d->setBcId(bcID);
    } else {
        BoundaryCondition *bc = static_cast< BoundaryCondition * >( d->giveDofManager()->giveDomain()->giveBc(bcID) );
        bc->setPrescribedValue(value);
    }
}

void
SolutionbasedShapeFunction :: initializeSurfaceData(modeStruct *mode)
{

    EngngModel *m=mode->myEngngModel;
    double TOL2=1e-3;

    IntArray pNodes, mNodes, zNodes;

    Set *mySet = this->domain->giveSet( this->giveSetNumber() );
    IntArray BoundaryList = mySet->giveBoundaryList();

    // First add all nodes to pNodes or nNodes respectively depending on coordinate and normal.
    for ( int i = 0; i < BoundaryList.giveSize() / 2; i++ ) {
        int ElementID = BoundaryList(2 * i);
        int Boundary = BoundaryList(2 * i + 1);

        Element *e = m->giveDomain(1)->giveElement(ElementID);
        FEInterpolation *geoInterpolation = e->giveInterpolation();

        // Check all sides of element
        IntArray bnodes;

#define usePoints 1
#if usePoints == 1
        // Check if all nodes are on the boundary
        geoInterpolation->boundaryGiveNodes(bnodes, Boundary);
        for ( int k = 1; k <= bnodes.giveSize(); k++ ) {
            DofManager *dman = e->giveDofManager( bnodes.at(k) );
            for ( int l = 1; l <= dman->giveCoordinates()->giveSize(); l++ ) {
                if ( fabs( dman->giveCoordinates()->at(l) - maxCoord.at(l) ) < TOL2 ) {
                    pNodes.insertOnce( dman->giveNumber() );
                }
                if ( fabs( dman->giveCoordinates()->at(l) - minCoord.at(l) ) < TOL2 ) {
                    mNodes.insertOnce( dman->giveNumber() );
                }
            }
        }
#else
        // Check normal
        FloatArray lcoords;
        lcoords.resize(2);
        lcoords.at(1) = 0.33333;
        lcoords.at(2) = 0.33333;

        FloatArray normal;
        geoInterpolation->boundaryEvalNormal( normal, j, lcoords, FEIElementGeometryWrapper(e) );
        geoInterpolation->boundaryGiveNodes(bnodes, j);

        printf( "i=%u\tj=%u\t(%f\t%f\t%f)\n", i, j, normal.at(1), normal.at(2), normal.at(3) );
        for ( int k = 1; k <= normal.giveSize(); k++ ) {
            if ( fabs( ( fabs( normal.at(k) ) - 1 ) ) < 1e-4 ) { // Points in x, y or z direction
                addTo = NULL;
                if ( normal.at(k) > 0.5 ) {
                    addTo = & pNodes;
                }
                if ( normal.at(k) < -0.5 ) {
                    addTo = & mNodes;
                }
                if ( addTo != NULL ) {
                    for ( int l = 1; l <= bnodes.giveSize(); l++ ) {
                        bool isSurface = false;
                        DofManager *dman = e->giveDofManager( bnodes.at(l) );
                        dman->giveCoordinates()->printYourself();
                        for ( int m = 1; m <= dman->giveCoordinates()->giveSize(); m++ ) {
                            if ( ( fabs( dman->giveCoordinates()->at(m) - maxCoord.at(m) ) < TOL2 ) || ( fabs( dman->giveCoordinates()->at(m) - minCoord.at(m) ) < TOL2 ) ) {
                                isSurface = true;
                            }
                        }

                        if ( isSurface ) {
                            addTo->insertOnce( e->giveDofManagerNumber( bnodes.at(l) ) );
                        }
                    }
                }
            }
        }
#endif
    }

#if 0
    printf("p=[");
    for ( int i = 1; i < pNodes.giveSize(); i++ ) {
        printf( "%u, ", pNodes.at(i) );
    }
    printf("];\n");
    printf("m=[");
    for ( int i = 1; i < mNodes.giveSize(); i++ ) {
        printf( "%u, ", mNodes.at(i) );
    }
    printf("];\n");
#endif
    //The intersection of pNodes and mNodes constitutes zNodes
    {
        int i = 1, j = 1;
        while ( i <= pNodes.giveSize() ) {
            j = 1;
            while ( j <= mNodes.giveSize() && ( i <= pNodes.giveSize() ) ) {
                //printf("%u == %u?\n", pNodes.at(i), mNodes.at(j));
                if ( pNodes.at(i) == mNodes.at(j) ) {
                    zNodes.insertOnce( pNodes.at(i) );
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
    copyDofManagersToSurfaceData(mode, pNodes, true, false, false);
    copyDofManagersToSurfaceData(mode, mNodes, false, true, false);
    copyDofManagersToSurfaceData(mode, zNodes, false, false, true);

#if 0
    printf("p2=[");
    for ( int i = 1; i <= pNodes.giveSize(); i++ ) {
        printf( "%u, ", pNodes.at(i) );
    }
    printf("];\n");
    printf("m2=[");
    for ( int i = 1; i <= mNodes.giveSize(); i++ ) {
        printf( "%u, ", mNodes.at(i) );
    }
    printf("];\n");
    printf("z2=[");
    for ( int i = 1; i <= zNodes.giveSize(); i++ ) {
        printf( "%u, ", zNodes.at(i) );
    }
    printf("];\n");

    printf("pCoords=[");
    for ( int i = 1; i <= pNodes.giveSize(); i++ ) {
        FloatArray *coords = m->giveDomain(1)->giveDofManager( pNodes.at(i) )->giveCoordinates();
        printf( "%f, %f, %f; ", coords->at(1), coords->at(2), coords->at(3) );
    }
    printf("]\n");
    printf("mCoords=[");
    for ( int i = 1; i <= mNodes.giveSize(); i++ ) {
        FloatArray *coords = m->giveDomain(1)->giveDofManager( mNodes.at(i) )->giveCoordinates();
        printf( "%f, %f, %f; ", coords->at(1), coords->at(2), coords->at(3) );
    }
    printf("]\n");
    printf("zCoords=[");
    for ( int i = 1; i <= zNodes.giveSize(); i++ ) {
        FloatArray *coords = m->giveDomain(1)->giveDofManager( zNodes.at(i) )->giveCoordinates();
        printf( "%f, %f, %f; ", coords->at(1), coords->at(2), coords->at(3) );
    }
    printf("];\n");
#endif
}

void
SolutionbasedShapeFunction :: whichBoundary(FloatArray &coord, bool &isPlus, bool &isMinus, bool &isZero)
{
    isPlus = false;
    isMinus = false;
    isZero = false;

    for ( int k = 1; k <= coord.giveSize(); k++ ) {
        isPlus = isPlus || ( fabs( coord.at(k) - maxCoord.at(k) ) < TOL );
        isMinus = isMinus || ( fabs( coord.at(k) - minCoord.at(k) ) < TOL );
    }

    isZero = isPlus && isMinus;
}

void
SolutionbasedShapeFunction :: copyDofManagersToSurfaceData(modeStruct *mode, IntArray nodeList, bool isPlus, bool isMinus, bool isZeroBoundary)
{
    for ( int i = 1; i <= nodeList.giveSize(); i++ ) {
        FloatArray values;

        IntArray DofIDs;
        DofManager *dman = mode->myEngngModel->giveDomain(1)->giveDofManager(nodeList.at(i));

        computeBaseFunctionValueAt(values, *dman->giveCoordinates(), this->dofs, *mode->myEngngModel );

/* <<<<<<< HEAD
=======
        for ( int j = 1; j <= this->dofs.giveSize(); j++ ) {
            SurfaceDataStruct *surfaceData = new(SurfaceDataStruct);
            Dof *d = dman->giveDofWithID( dofs.at(j) );
>>>>>>> 147f565295394adef603dae296a820af5f28d9cd
*/
        // Check that current node contains current DofID
        for (int j=1; j<=this->dofs.giveSize(); j++) {
            for (Dof *d: *dman ){ //int k=1; k<= dman->giveNumberOfDofs(); k++ ) {

                //Dof *d = dman->dofArray.at(k);// giveDof(k);

                if (d->giveDofID() == this->dofs.at(j)) {
                    SurfaceDataStruct *surfaceData = new(SurfaceDataStruct);
                    surfaceData->DofID = (DofIDItem) this->dofs.at(j);
                    surfaceData->DofMan = dman;
                    surfaceData->isPlus = isPlus;
                    surfaceData->isMinus = isMinus;
                    surfaceData->isZeroBoundary = isZeroBoundary;
                    surfaceData->isFree = d->giveBcId() == 0;
                    surfaceData->value = values.at(j);

                    mode->SurfaceData.push_back(surfaceData);
                }
            }
        }
    }
}
}
