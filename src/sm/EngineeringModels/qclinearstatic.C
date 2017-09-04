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

#include "../sm/EngineeringModels/qclinearstatic.h"
#include "../sm/Quasicontinuum/fullsolveddomain.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/structuralelementevaluator.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "sparsemtrx.h"
#include "verbose.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

#include "qcnode.h"
#include "domain.h"
#include "dof.h"
#include "crosssection.h"

#include "../sm/Quasicontinuum/quasicontinuum.h"
#include "../sm/Quasicontinuum/fullsolveddomain.h"
#include "../sm/Quasicontinuum/geometrygenerator.h"
#include "unknownnumberingscheme.h"

#include "../sm/Quasicontinuum/quasicontinuumnumberingscheme.h"

 #include "t3dinterface.h"
#ifdef __T3D
 #include "t3d.h"
#endif

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "communicator.h"
#endif

#include <typeinfo>

namespace oofem {
REGISTER_EngngModel(QClinearStatic);

QClinearStatic :: QClinearStatic(int i, EngngModel *_master) : LinearStatic(i, _master), loadVector(), displacementVector()
{
    qcEquationNumbering = new QuasicontinuumNumberingscheme();
}


QClinearStatic :: ~QClinearStatic()
{
    delete qcEquationNumbering;
}



IRResultType
QClinearStatic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LinearStatic :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, qcApproach, _IFT_QuasiContinuum_approach);

    if ( qcApproach == 0 || qcApproach == 1 || qcApproach == 2 || qcApproach == 3 ) {
        // approach A0 = full solved particle model
        // A1 hanging nodes, interpolation only
        // A2 global homogenization
        // A3 individual homogenization
    } else {
        OOFEM_ERROR("Invalid input field \"qcApproach\". Value: %d is not supported", qcApproach);
    }

    if ( qcApproach > 0 ) { // apply qc
        // Initialize Full-Solved Domain
        this->initializeFullSolvedDomain(ir);
        this->Fullsolveddomain.initializeFrom(ir);

        homogenizationMtrxType = 0;
        if ( qcApproach > 1 ) {
            homogenizationMtrxType = 1; // isotropic is default
            IR_GIVE_OPTIONAL_FIELD(ir, homogenizationMtrxType, _IFT_QuasiContinuum_mtrx_type);
        }


        GeometryGenerator gg;

        // Generate particles
        generateParticles = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, generateParticles, _IFT_QuasiContinuum_generate_Particles);
        if ( generateParticles == 0 ) { // 0 - defined in input file
        } else if ( generateParticles == 1 ) { // 1 - generate + save
            gg.initializeParticleGenerator(ir);
            gg.generateParticles();
        } else if ( generateParticles == 2 ) { // 2 - load
            gg.loadParticles();
        } else {
            OOFEM_ERROR("Invalid input field \"genParticles\". Value: %d is not supported", generateParticles);
        }

        // Generate links
        generateLinks = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, generateLinks, _IFT_QuasiContinuum_generate_Links);
        if ( generateLinks == 0 ) {} else if ( generateLinks == 1 ) { // generate + save
            gg.initializeLinkGenerator(ir);
            gg.generateLinks();
        } else if ( generateLinks == 2 ) { // load
            gg.loadLinks();
        } else {
            OOFEM_ERROR("Invalid input field \"genLinks\". Value: %d is not supported", generateLinks);
        }

        // Generate interpolation elements
        generateInterpolationElements = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, generateInterpolationElements, _IFT_QuasiContinuum_generate_Interpolation_Elements);

        interpolationElementsMaterialNumber = 0;
        if ( generateInterpolationElements == 0 ) {
            IR_GIVE_FIELD(ir, interpolationElementsMaterialNumber, _IFT_QuasiContinuum_interp_Mat_Number);
        } else if ( generateInterpolationElements == 1 ) {
            IR_GIVE_FIELD(ir, defaultT3DMeshSize, _IFT_QuasiContinuum_T3D_Interpolation_Mesh_size);
            IR_GIVE_FIELD(ir, t3dFileName, _IFT_QuasiContinuum_t3d_File_Name);
        } else if ( generateInterpolationElements == 2 ) {
            IR_GIVE_FIELD(ir, t3dFileName, _IFT_QuasiContinuum_t3d_File_Name);
        } else {
            OOFEM_ERROR("Invalid input field \"GenInterpElem\". Value: %d is not supported", generateInterpolationElements);
        }



        // test of combination of input geometry generating
#if DEBUG
        if  ( generateParticles == 1 && generateLinks != 1 ) {
            OOFEM_ERROR("Links cannot be set manualy if particles are generated automaticaly");
        }
        if ( generateParticles == 1 &&  generateInterpolationElements != 1 ) {
            OOFEM_ERROR("Interpolation elements cannot be set manually if particles are generated automaticaly");
        }
#endif
    } else { // if qcApproach == 0
        generateParticles = 0;
        generateLinks = 0;
        generateInterpolationElements = 0;
    }



#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new NodeCommunicator( this, commBuff, this->giveRank(),
                                             this->giveNumberOfProcesses() );
    }

#endif


    return IRRT_OK;
}


void
QClinearStatic :: postInitialize()
{
    // checkout print of dof type before interpolation mesh is generated (1/2 = RN/HN)
    /*
     * for (int i=1; i<=this->giveDomain(1)->giveNumberOfDofManagers(); i++)
     *  {
     *      qcNode *n = dynamic_cast< qcNode * >( this->giveDomain(1)->giveDofManager(i) );
     *      OOFEM_LOG_INFO(" xxx:  node %d type %d \n",n->giveGlobalNumber(),n->giveQcNodeType());
     *  }
     */

    Quasicontinuum qc;
    qc.setNoDimensions( this->giveDomain(1) );

    // interpolation mesh
    numberOfIntepolationElements = 0;
    if ( qcApproach > 0 ) {
        createInterpolationMeshNodes( this->giveDomain(1) );
        setRepNodesInVerticesOfInterpolationMesh( this->giveDomain(1) );
        qc.setupInterpolationMesh(this->giveDomain(1), generateInterpolationElements, interpolationElementsMaterialNumber, & interpolationMeshNodes);
        qc.createInterpolationElements( this->giveDomain(1) );
        qc.addCrosssectionToInterpolationElements( this->giveDomain(1) );
    }

    // checkout print of dof type after interpolation mesh is generated (1/2 = RN/HN)
    /*
     * OOFEM_LOG_INFO(" xxx:  --------------------- \n");
     * for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
     *  //qcNode *n = dynamic_cast< qcNode * >( this->giveDomain(1)->giveDofManager(i) );
     *  Node *n = this->giveDomain(1)->giveNode(i);
     *  OOFEM_LOG_INFO( " xxx:  node %d type %d \n", n->giveGlobalNumber(), n->giveQcNodeType() );
     * }
     */

    EngngModel :: postInitialize();


    //  apply QC

    if ( qcApproach == 0 ) {  //A0 - fully solved model
        // all nodes nad elements are activated
        this->activatedNodeList.resize(this->giveDomain(1)->giveNumberOfDofManagers(), true);
        this->activatedElementList.resize(this->giveDomain(1)->giveNumberOfElements(), true);
    } else if ( qcApproach == 1 ) { // A1 - hanging nodes
        // all nodes nad elements are activated
        this->activatedNodeList.resize(this->giveDomain(1)->giveNumberOfDofManagers(), true);
        this->activatedElementList.resize(this->giveDomain(1)->giveNumberOfElements(), true);

        qc.applyApproach1( this->giveDomain(1) );

        // postinitialize interpolation elements // ??km?? TO DO: giveInterpolationElem...
        // TO DO: in EngngModel::postInitialize() elements are postinitialized again
        for ( auto &el : this->giveDomain(1)->giveElements() ) {
            el->postInitialize();
        }
    } else if ( qcApproach == 2 ) { // A2 - global homogenization
        // elements need to be postintialized before computation of volume.
        // TO DO: in EngngModel::postInitialize() elements are postinitialized again
        // postinitialize interpolation elements // ??km?? TO DO: giveInterpolationElem...
        for ( auto &el : this->giveDomain(1)->giveElements() ) {
            el->postInitialize();
        }

        double volumeOfInterpolationMesh = this->computeTotalVolumeOfInterpolationMesh( this->giveDomain(1) );
        qc.applyApproach2(this->giveDomain(1), homogenizationMtrxType, volumeOfInterpolationMesh);
    } else if ( qcApproach == 3 ) { // A3 - local homogenization
        // elements need to be postintialized before computation of volume.
        // TO DO: in EngngModel::postInitialize() elements are postinitialized again
        // postinitialize interpolation elements // ??km?? TO DO: giveInterpolationElem...
        for ( auto &el : this->giveDomain(1)->giveElements() ) {
            el->postInitialize();
        }
        qc.applyApproach3(this->giveDomain(1), homogenizationMtrxType);
    } else {
        OOFEM_ERROR("unsupported qcApproach number");
    }


    // checkout print of activated elements
    /*
     * OOFEM_LOG_INFO(" xxx:  --------------------- \n");
     * for ( int i = 1; i <= ( int ) this->activatedElementList.size(); i++ ) {
     *  if ( this->activatedElementList.at(i - 1) ) {
     *      OOFEM_LOG_INFO(" active el %d d \n", i);
     *  }
     * }
     */
}



void QClinearStatic :: solveYourself()
{
    LinearStatic :: solveYourself();
}

void QClinearStatic :: solveYourselfAt(TimeStep *tStep)
{
    // initialize node eq numbering (for actived nodes only)
    if ( !qcEquationNumbering->giveIsInitializedFlag() ) {
        qcEquationNumbering->init(this->giveDomain(1), activatedNodeList, tStep);
    }
    LinearStatic :: solveYourselfAt(tStep);
}


IRResultType
QClinearStatic :: initializeFullSolvedDomain(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainNodes, _IFT_FullSolvedDomain_nodes);
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainElements, _IFT_FullSolvedDomain_elements);
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainRadius, _IFT_FullSolvedDomain_radius);
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainBox, _IFT_FullSolvedDomain_box);
    // check input format
#ifdef DEBUG
    if ( FullSolvedDomainRadius.giveSize() != 0 && FullSolvedDomainRadius.giveSize() % 4 != 0 ) {
        OOFEM_ERROR("invalid format of FullSolvedDomainRadius");
    }
    if ( FullSolvedDomainBox.giveSize() != 0 && FullSolvedDomainBox.giveSize() % 6 != 0 ) {
        OOFEM_ERROR("invalid format of FullSolvedDomainBox");
    }
#endif

    return IRRT_OK;
}

bool
QClinearStatic :: nodeInFullSolvedDomainTest(Node *n)
{
    FloatArray *coordinates = n->giveCoordinates();
    // is tested node in FullSolvedDomainNodes
    if ( FullSolvedDomainNodes.giveSize() != 0 ) {
        for ( int i = 1; i <= FullSolvedDomainNodes.giveSize(); i++ ) {
            if ( n->giveNumber() == FullSolvedDomainNodes.at(i) ) {
                return true;
            }
        }
    }

    // is tested node in FullSolvedDomainElements
    if ( FullSolvedDomainElements.giveSize() != 0 ) {
        for ( int i = 1; i <= FullSolvedDomainElements.giveSize(); i++ ) {
            OOFEM_ERROR("Definition of Full Solved Domain by list of interpolation element has not been implemented yet");
            // ??km?? TO DO:
            //if ( "node n is in element with number FullSolvedDomainElements.at(i)" ) {
            // return true;
            //}
        }
    }

    // is tested node in FullSolvedDomainRadius
    if ( FullSolvedDomainRadius.giveSize() != 0 ) {
        for ( int i = 0; i <= FullSolvedDomainRadius.giveSize() / 4 - 1; i++ ) {
            FloatArray vector;
            vector.resize(3);
            vector.at(1) = coordinates->at(1) - FullSolvedDomainRadius.at(4 * i + 1);
            vector.at(2) = coordinates->at(2) - FullSolvedDomainRadius.at(4 * i + 2);
            vector.at(3) = coordinates->at(3) - FullSolvedDomainRadius.at(4 * i + 3);
            if ( vector.computeNorm() <= FullSolvedDomainRadius.at(4 * i + 4) ) {
                return true;
            }
        }
    }

    // is tested node in FullSolvedDomainBox
    if ( FullSolvedDomainBox.giveSize() != 0 ) {
        for ( int i = 0; i <= FullSolvedDomainBox.giveSize() / 6 - 1; i++ ) {
            FloatArray A(3), B(3);
            // left bottom corner coords
            A.at(1) = FullSolvedDomainBox.at(6 * i + 1);
            A.at(2) = FullSolvedDomainBox.at(6 * i + 2);
            A.at(3) = FullSolvedDomainBox.at(6 * i + 3);
            // right top corner coords
            B.at(1) = FullSolvedDomainBox.at(6 * i + 4);
            B.at(2) = FullSolvedDomainBox.at(6 * i + 5);
            B.at(3) = FullSolvedDomainBox.at(6 * i + 4);

            if ( ( A.at(1) <= coordinates->at(1) ) && ( coordinates->at(1) <= B.at(1) ) && \
                 ( A.at(2) <= coordinates->at(2) ) && ( coordinates->at(2) <= B.at(2) ) && \
                 ( A.at(3) <= coordinates->at(3) ) && ( coordinates->at(3) <= B.at(3) ) ) {
                return true;
            }
        }
    }


    return false;
}


void
QClinearStatic :: setRepNodesInVerticesOfInterpolationMesh(Domain *d)
{
    int nol = 0; // number of links
    int noe = d->giveNumberOfElements(); // number of all elements
    for ( int i = 1; i <= noe; i++ ) { // over elements
        int noen = d->giveElement(i)->giveNumberOfDofManagers(); // number of element nodes
        if ( noen > 2 ) { // skip links with 2 nodes only
            for ( int j = 1; j <= noen; j++ ) { // over element nodes
                qcNode *n = dynamic_cast< qcNode * >( d->giveElement(i)->giveDofManager(j) );
                if ( n ) {
                    n->setAsRepnode();
                } else {
                    OOFEM_WARNING( "Node %d of interpolation element %d is not \"qcNode\", quasicontinuum is not applied in this node", d->giveElement(i)->giveGlobalNumber(), d->giveElement(i)->giveInternalDofManager(j)->giveGlobalNumber() );
                }
            }
        } else if ( noen == 2 ) {
            nol++;
        }
    }
    // here we assume that in input file, interpolation elements ale listed continuously in the end of element list (after links )
    numberOfIntepolationElements = noe - nol;
}


void
QClinearStatic :: setQCNodeType(Domain *d)
{
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        qcNode *n = dynamic_cast< qcNode * >( d->giveDofManager(i) );
        if ( n ) {
            if  ( !this->nodeInFullSolvedDomainTest(n) ) {
                n->setAsHanging();
            }
        } else {
            //OOFEM_ERROR("Node %d: Only \"qcNode\" type is supported in quasicontinuum simulation", i);
            OOFEM_WARNING("Node %d is not \"qcNode\", quasicontinuum is not applied in this node", i);
        }
    }
}

void
QClinearStatic :: updateNodeTypes(Domain *d)
{
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        //   for ( Dof *dof: *this ) {
        //for ( qcNode *n: *this ) {
        qcNode *n = dynamic_cast< qcNode * >( d->giveDofManager(i) );
        if ( n ) {
            if  ( this->nodeInFullSolvedDomainTest(n) ) {
                n->setAsRepnode();
            } else {
                n->setAsHanging();
            }
        } else {
            //OOFEM_ERROR("Node %d: Only \"qcNode\" type is supported in quasicontinuum simulation", i);
            OOFEM_WARNING("Node %d is not \"qcNode\", quasicontinuum is not applied in this node", i);
        }
    }
}


void
QClinearStatic :: createInterpolationMeshNodes(Domain *d)
{
    // interpolation mesh
    if ( generateInterpolationElements == 0 ) {
        // interpolation mesh is defined in input file - elements with material number qcMatNumber
    } else if ( generateInterpolationElements == 1 ) {
        this->interpolationMeshNodes = this->generateInterpolationMesh( this->giveDomain(1) );
    } else if ( generateInterpolationElements == 2 ) {
        this->interpolationMeshNodes = this->loadInterpolationMesh( this->giveDomain(1) );
    }
}

std :: vector< IntArray >
QClinearStatic :: generateInterpolationMesh(Domain *d)
{
    std :: vector< IntArray >newMeshNodes;
#ifdef __T3D
    T3DInterface *t3dInterface = new T3DInterface( this->giveDomain(1) );
#endif
    std :: vector< FloatArray >nodeCoords;
    std :: vector< IntArray >meshNodes;
    IntArray cellTypes;

    char options_string [ 128 ] = "";
    const char *t3dInFile; //char *t3dInFile = "interpolationMesh_t3d.in";
    const char *t3dOutFile; //char *t3dOutFile = "interpolationMesh_t3d.out";

    std :: string temp;
    temp = t3dFileName;
    temp.append(".oofem.t3d.out");
    t3dOutFile = temp.c_str();
    t3dInFile = t3dFileName.c_str();
    sprintf(options_string, "-i %s -o %s -S -u %f -Q -W", t3dInFile, t3dOutFile, defaultT3DMeshSize);

#ifdef __T3D
    if ( generateInterpolationElements == 1 ) {
        /* run T3D to generate mesh */
        try {
            t3d_main(NULL, options_string);
        } catch(int exit_code) {
            fprintf(stderr, "T3d was prematuraly terminated with error code %d\n\n", exit_code);
        }
    }
#else
    OOFEM_ERROR("OOFEM is NOT compiled with T3D option but T3D in needed");
#endif

    if ( generateInterpolationElements == 1 || generateInterpolationElements == 2 ) {
        // create interpolation mesh from t3d output file
#ifdef __T3D
        t3dInterface->createQCInterpolationMesh(t3dOutFile, nodeCoords, meshNodes, cellTypes);
#endif

        // map mesh to position of particles
        newMeshNodes = this->transformMeshToParticles(this->giveDomain(1), nodeCoords, meshNodes);

        // degenerated elements have been removed
        // check if mesch is consistent after removal. ??km?? TO DO ?? How to check it?
    }
    return newMeshNodes;
}

std :: vector< IntArray >
QClinearStatic :: loadInterpolationMesh(Domain *d)
{
    std :: vector< IntArray >newMeshNodes;
    //#ifdef __T3D
    T3DInterface *t3dInterface = new T3DInterface( this->giveDomain(1) );
    //#endif
    std :: vector< FloatArray >nodeCoords;
    std :: vector< IntArray >meshNodes;
    IntArray cellTypes;

    const char *t3dOutFile; //char *t3dOutFile = "interpolationMesh_t3d.out";

    t3dOutFile = t3dFileName.c_str();

    // create interpolation mesh from given t3d output file
    //#ifdef __T3D
    t3dInterface->createQCInterpolationMesh(t3dOutFile, nodeCoords, meshNodes, cellTypes);
    //#endif

    // map mesh to position of particles
    newMeshNodes = this->transformMeshToParticles(this->giveDomain(1), nodeCoords, meshNodes);

    // degenerated elements have been removed
    // check if mesch is consistent after removal. ??km?? TO DO ?? How to check it?
    return newMeshNodes;
}


std :: vector< IntArray >
QClinearStatic :: transformMeshToParticles(Domain *d, std :: vector< FloatArray > &nodeCoords, std :: vector< IntArray > &meshNodes)
// all nodes in meshNodes are shifted to the position of neares particle (i.e. node in domain)
// this node is set as repnode
// elements degenerated (by the shift of two vertexes into the same particle) are removed from the list
//
{
    // number of particles (nodes in domain)
    //int nop = this->giveDomain(1)->giveNumberOfDofManagers() ; // TO DO we assume that all DofManagers all nodes (particles)
    // number of (mesh) nodes
    int nomn = nodeCoords.size();
    // number of (mesh) element
    int nome = meshNodes.size();

    //newMeshNodes=meshNodes;
    std :: vector< IntArray >newMeshNodes;
    newMeshNodes.clear();
    //particleOfMeshNode=zeros(1,nomn);
    IntArray newNodeNumbers;
    newNodeNumbers.resize(nomn);

    // loop over all nodes of mesh elments
    for ( int i = 1; i <= nomn; i++ ) {
        // ??km?? TO DO: use octree in "findNearestParticle", now it is not effective
        // find nearest node and set him as repnode (is neares node QCnode?)
        qcNode *nearestParticle = dynamic_cast< qcNode * >( findNearestParticle(d, nodeCoords [ i - 1 ]) );
        if ( nearestParticle ) {
            if  ( nearestParticle->giveQcNodeType() != 1 ) { // if nearest node is not repnode
                nearestParticle->setAsRepnode();
            }
        } else {
            OOFEM_ERROR("Vertex of interpolation mesh is shifted to node which is not \"qcNode\"");
            // ??km?? TO DO: add combination of qcNode and normal node in one input. Then not necessarily qcNode is needed.
            //        Test if node "nearestParticle" has independet DOFs should be enough.
            //        OOFEM_WARNING(" ");
        }

        nodeCoords [ i - 1 ] = * nearestParticle->giveCoordinates();
        newNodeNumbers [ i - 1 ] =  nearestParticle->giveNumber();
    }

    //
    // check and remove degenerated elements
    //


    // for 2D triangles
    if ( meshNodes [ 0 ].giveSize() == 3 ) {
        // TO DO: how to do not consider interpolation elements in user define area
        // ??km??       for example in FSD it is better without interpolation elements
        // check if two or more vertexes were shifted to the same particle -- element degeneration
        for ( int i = 1; i <= nome; i++ ) {
            int n1 = meshNodes [ i - 1 ].at(1); // node numbers
            int n2 = meshNodes [ i - 1 ].at(2);
            int n3 = meshNodes [ i - 1 ].at(3);
            if ( newNodeNumbers.at(n1) == newNodeNumbers.at(n2) ||       \
                 newNodeNumbers.at(n1) == newNodeNumbers.at(n3) ||       \
                 newNodeNumbers.at(n2) == newNodeNumbers.at(n3) ) {
                continue; // skip degenerate element
            } else {                                            // check degeneration to negative volume
                // ??km?? TO DO: how to check negative volume in 2D
                // overlaping elements may cause problems in homogenization
                double detJ = 1.0;
                if ( detJ <= 0 ) {
                    OOFEM_WARNING("%d-th interpolation element degenerates to negative volume", i);
                } else {                                                                                                 // element is ok and will be saved
                    IntArray OkmeshNodes(3);
                    OkmeshNodes.at(1) = newNodeNumbers.at(n1);
                    OkmeshNodes.at(2) = newNodeNumbers.at(n2);
                    OkmeshNodes.at(3) = newNodeNumbers.at(n3);
                    newMeshNodes.push_back(OkmeshNodes);
                }
            }
        } // over elements in 2D

        // for 3D tetras
    } else {
        for ( int i = 1; i <= nome; i++ ) {
            // TO DO: how to do not consider interpolation elements in user define area
            // ??km??       for example in FSD it is better without interpolation elements

            // check if two or more vertexes were shifted to the same particle -- element degeneration
            int n1 = meshNodes [ i - 1 ].at(1); // node numbers
            int n2 = meshNodes [ i - 1 ].at(2);
            int n3 = meshNodes [ i - 1 ].at(3);
            int n4 = meshNodes [ i - 1 ].at(4);
            if ( newNodeNumbers.at(n1) == newNodeNumbers.at(n2) ||       \
                 newNodeNumbers.at(n1) == newNodeNumbers.at(n3) ||       \
                 newNodeNumbers.at(n1) == newNodeNumbers.at(n4) ||       \
                 newNodeNumbers.at(n2) == newNodeNumbers.at(n3) ||       \
                 newNodeNumbers.at(n2) == newNodeNumbers.at(n4) ||       \
                 newNodeNumbers.at(n3) == newNodeNumbers.at(n4) ) {
                continue; // skip degenerate element
            } else {                                            // check degeneration to negative volume
                double x1 = d->giveDofManager( newNodeNumbers.at(n1) )->giveCoordinates()->at(1);
                double y1 = d->giveDofManager( newNodeNumbers.at(n1) )->giveCoordinates()->at(2);
                double z1 = d->giveDofManager( newNodeNumbers.at(n1) )->giveCoordinates()->at(3);
                double x2 = d->giveDofManager( newNodeNumbers.at(n2) )->giveCoordinates()->at(1);
                double y2 = d->giveDofManager( newNodeNumbers.at(n2) )->giveCoordinates()->at(2);
                double z2 = d->giveDofManager( newNodeNumbers.at(n2) )->giveCoordinates()->at(3);
                double x3 = d->giveDofManager( newNodeNumbers.at(n3) )->giveCoordinates()->at(1);
                double y3 = d->giveDofManager( newNodeNumbers.at(n3) )->giveCoordinates()->at(2);
                double z3 = d->giveDofManager( newNodeNumbers.at(n3) )->giveCoordinates()->at(3);
                double x4 = d->giveDofManager( newNodeNumbers.at(n4) )->giveCoordinates()->at(1);
                double y4 = d->giveDofManager( newNodeNumbers.at(n4) )->giveCoordinates()->at(2);
                double z4 = d->giveDofManager( newNodeNumbers.at(n4) )->giveCoordinates()->at(3);
                double detJ = ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) + ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) + ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 );
                if ( detJ <= 0 ) {
                    OOFEM_WARNING("%d-th interpolation element degenerates to negative volume", i);
                    continue;
                } else { // element is ok and will be saved
                    IntArray OkmeshNodes(4);
                    OkmeshNodes.at(1) = newNodeNumbers.at(n1);
                    OkmeshNodes.at(2) = newNodeNumbers.at(n2);
                    OkmeshNodes.at(3) = newNodeNumbers.at(n3);
                    OkmeshNodes.at(4) = newNodeNumbers.at(n4);
                    newMeshNodes.push_back(OkmeshNodes);
                }
            }
        } // over elements
    } // in 3D

    return newMeshNodes;
}

double
QClinearStatic :: computeTotalVolumeOfInterpolationMesh(Domain *d)
// TO DO: computeVolumeAreaOrLength() is not working
{
    int dim =  d->giveNumberOfSpatialDimensions();
    // loop ovel all elements
    double volume = 0.;
    double totalVolume = 0.;

    if ( dim == 2 ) {
        double th;

        for ( int i = 1; i <= d->giveNumberOfElements(); i++ ) {
            // skip links with 2 nodes
            // ??KM?? TO DO: use only list of interpolation elements
            if ( d->giveElement(i)->giveNumberOfDofManagers() > 2 ) {
                volume = d->giveElement(i)->computeVolumeAreaOrLength();
                th = d->giveElement(i)->giveCrossSection()->give(CS_Thickness, NULL);
                totalVolume += volume / th;
            }
        }
    } else if ( dim == 3 ) {
        for ( int i = 1; i <= d->giveNumberOfElements(); i++ ) {
            // skip links with 2 nodes
            // ??KM?? TO DO: use only list of interpolation elements
            if ( d->giveElement(i)->giveNumberOfDofManagers() > 2 ) {
                volume = d->giveElement(i)->computeVolumeAreaOrLength();
                totalVolume += volume;
            }
        }
    } else {
        OOFEM_ERROR("Invalid number of dimensions. Only 2d and 3d domains are supported in QC simulation. \n");
    }

    return totalVolume;
}


DofManager *
QClinearStatic :: findNearestParticle(Domain *d, FloatArray coords)
{
    // TO DO: use octree here
    double minDistance = 1.0e100;
    double distance;
    DofManager *p;
    // loop over all particles (nodes in existing domain)
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        distance = coords.distance( d->giveDofManager(i)->giveCoordinates() );
        if ( distance < minDistance ) {
            minDistance = distance;
            p = d->giveDofManager(i);
        }
    }
    if ( p ) {
        return p;
    } else {
        OOFEM_ERROR( "Neares particle for point [%d, %d] not found", coords.at(1), coords.at(2) );
        return NULL;
    }
}

QCFullsolveddomain *
QClinearStatic :: giveFullSolvedDomain()
{
    QCFullsolveddomain *p;
    p = & Fullsolveddomain;
    return p;
}


void
QClinearStatic :: setActivatedNodeList(IntArray nodeList, Domain *d)
{
    this->activatedNodeList.clear();
    this->activatedNodeList.resize(nodeList.giveSize(), false);
    for ( int i = 1; i <= nodeList.giveSize(); i++ ) {
        // if node is activated or master node
        if ( ( nodeList.at(i) == 1 ) || ( d->giveNode(i)->giveQcNodeType() == 1 ) ) {
            this->activatedNodeList [ i - 1 ] = true;
        }
    }
}

void
QClinearStatic :: setActivatedElementList(IntArray elemList)
{
    this->activatedElementList.clear();
    this->activatedElementList.resize(elemList.giveSize(), false);
    for ( int i = 1; i <= elemList.giveSize(); i++ ) {
        // if element is activated
        if ( elemList.at(i) == 1 ) {
            this->activatedElementList [ i - 1 ] = true;
        }
    }
}
} // end namespace oofem
