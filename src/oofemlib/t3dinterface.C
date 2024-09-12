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

#include "t3dinterface.h"
#include "errorestimator.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "engngm.h"
#include "remeshingcrit.h"
#include "oofemtxtdatareader.h"
#include "dynamicinputrecord.h"

#include "crosssection.h"
#include "classfactory.h"

#include "nonlocalbarrier.h"
#include "initialcondition.h"
#include "classfactory.h"
//#include "loadtimefunction.h"
#include "function.h"
#include "outputmanager.h"

#include <cmath>

namespace oofem {
REGISTER_Mesher(T3DInterface, MPT_T3D);


MesherInterface :: returnCode
T3DInterface :: createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew)
{
    * dNew = NULL;
    if ( this->createInput(this->domain, tStep) ) {
        return MI_NEEDS_EXTERNAL_ACTION;
    } else {
        return MI_FAILED;
    }
}

int
T3DInterface :: createInput(Domain *d, TimeStep *tStep)
{
    int nnodes = d->giveNumberOfDofManagers(), nelem = d->giveNumberOfElements();
    double density;
    FILE *outputStrem;
    Node *inode;
    Element *ielem;
    int edges, trias, quads, tetras, pyrams, wedges, hexas;
    IntArray edgeIdArray, triaIdArray, quadIdArray, tetraIdArray, pyramIdArray, wedgeIdArray, hexaIdArray;
    bool tri_tetra = false;
    char fileName [ 32 ];

    edges = trias = quads = tetras = pyrams = wedges = hexas = 0;
    for ( int i = 1; i <= nelem; i++ ) {
        ielem = d->giveElement(i);
        switch ( ielem->giveGeometryType() ) {
        case EGT_point:
            break;
        case EGT_line_1:
        case EGT_line_2:
            edges++;
            break;
        case EGT_triangle_1:
        case EGT_triangle_2:
            trias++;
            break;
        case EGT_quad_1:
        case EGT_quad_2:
            quads++;
            break;
        case EGT_tetra_1:
        case EGT_tetra_2:
            tetras++;
            break;
        case EGT_wedge_1:
        case EGT_wedge_2:
            wedges++;
            break;
        case EGT_hexa_1:
        case EGT_hexa_2:
            hexas++;
            break;
        default:
            OOFEM_ERROR( "unknown element type (%s)",
                         __Element_Geometry_TypeToString( ielem->giveGeometryType() ) );
        }
    }

    edgeIdArray.resize(edges);
    triaIdArray.resize(trias);
    quadIdArray.resize(quads);
    tetraIdArray.resize(tetras);
    pyramIdArray.resize(pyrams);
    wedgeIdArray.resize(wedges);
    hexaIdArray.resize(hexas);

    edges = trias = quads = tetras = pyrams = wedges = hexas = 0;
    for ( int i = 1; i <= nelem; i++ ) {
        ielem = d->giveElement(i);
        switch ( ielem->giveGeometryType() ) {
        case EGT_point:
            break;
        case EGT_line_1:
        case EGT_line_2:
            edgeIdArray.at(++edges) = i;
            break;
        case EGT_triangle_1:
        case EGT_triangle_2:
            triaIdArray.at(++trias) = i;
            break;
        case EGT_quad_1:
        case EGT_quad_2:
            quadIdArray.at(++quads) = i;
            break;
        case EGT_tetra_1:
        case EGT_tetra_2:
            tetraIdArray.at(++tetras) = i;
            break;
        case EGT_wedge_1:
        case EGT_wedge_2:
            wedgeIdArray.at(++wedges) = i;
            break;
        case EGT_hexa_1:
        case EGT_hexa_2:
            hexaIdArray.at(++hexas) = i;
            break;
        default:
            break;
        }
    }

    if ( quads + hexas + pyrams + wedges == 0 ) {
        tri_tetra = true;
    }

    if ( d->giveEngngModel()->isParallel() ) {
        sprintf( fileName, "%s.%d", BMF_FILENAME, d->giveEngngModel()->giveRank() );
    } else {
        sprintf(fileName, "%s", BMF_FILENAME);
    }

    outputStrem = fopen(fileName, "w");
    if ( tri_tetra == true ) {
        fprintf(outputStrem, "3 1 -1\n");
        fprintf(outputStrem, "%d %d %d %d\n", nnodes, edges, trias, tetras);
    } else {
        fprintf(outputStrem, "7 1 -1\n");
        fprintf(outputStrem, "%d %d %d %d %d %d %d %d\n", nnodes, edges, trias, quads, tetras, pyrams, wedges, hexas);
    }

    // loop over nodes
    for ( int i = 1; i <= nnodes; i++ ) {
        density = d->giveErrorEstimator()->giveRemeshingCrit()->giveRequiredDofManDensity(i, tStep);
        inode = d->giveNode(i);
        fprintf(outputStrem, "%d %e %e %e  %e\n", i, inode->giveCoordinate(1), inode->giveCoordinate(2), inode->giveCoordinate(3), density);
    }

    // loop separately for each type of element
    // since T3ds support only linear elements in bg mesh, only corner nodes are stored
    // alternatively quadratic elements may be splitted to linear (not supported now)

    if ( edges != 0 ) {
        for ( int i = 1; i <= edges; i++ ) {
            ielem = d->giveElement( edgeIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 2; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( trias != 0 ) {
        for ( int i = 1; i <= trias; i++ ) {
            ielem = d->giveElement( triaIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 3; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( quads != 0 ) {
        for ( int i = 1; i <= quads; i++ ) {
            ielem = d->giveElement( quadIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 4; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( tetras != 0 ) {
        for ( int i = 1; i <= tetras; i++ ) {
            ielem = d->giveElement( tetraIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 4; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( pyrams != 0 ) {
        for ( int i = 1; i <= pyrams; i++ ) {
            ielem = d->giveElement( pyramIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 5; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( wedges != 0 ) {
        for ( int i = 1; i <= wedges; i++ ) {
            ielem = d->giveElement( wedgeIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 6; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( hexas != 0 ) {
        for ( int i = 1; i <= hexas; i++ ) {
            ielem = d->giveElement( hexaIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( int j = 1; j <= 8; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    fclose(outputStrem);

    OOFEM_LOG_INFO("t3d.bmf file created\n");
    return 1;
}





int
T3DInterface :: t3d_2_OOFEM(const char *t3dOutFile, Domain **dNew)
{
    std::ifstream inputStream;
    inputStream.open( t3dOutFile );
    if ( !inputStream.is_open() ) {
        OOFEM_ERROR("OOFEMTXTDataReader::OOFEMTXTDataReader: Can't open T3D input stream (%s)", t3dOutFile);
        return 0;
    }
    // create new domain
    ( * dNew ) = new Domain( 2, this->domain->giveSerialNumber() + 1, this->domain->giveEngngModel() );
    ( * dNew )->setDomainType( this->domain->giveDomainType() );

    std::string line;

    //read first line from t3d out file - 4 numbers
    // 2 degree of interpolation
    // other are not important so far 
    std::getline(inputStream, line);
    //convert const char to char in order to use strtok
    char *currentLine = new char[line.size() + 1];
    std::strcpy ( currentLine, line.c_str() );
    // tokenizing line
    char *token = std::strtok(currentLine, " ");
    // set counter to 0
    int i = 0;
    while (token != nullptr) {
        if(i == 0)
            {}//int n1 = atoi(token);
        else if (i == 1)
            {}//int interp = atoi(token);
        else if(i == 2)
            {}//int n3 = atoi(token);
        else if(i == 3)
            {}//int n4 = atoi(token);
        else
            break;
        token = std::strtok(nullptr, " ");
        i++;
    }


    int nnodes = 0, ntriangles = 0, ntetras = 0; // nedges,

    /*read second line from t3d out file
        4 numbers
        1 - number of nodes
        2 - number of edges
        3 - number of triangles
        4 - number of tetras
    */
    std::getline(inputStream, line);
    //convert const char to char in order to use strtok
    currentLine = new char[line.size() + 1];
    std::strcpy ( currentLine, line.c_str() );
    // tokenizing line
    token = std::strtok(currentLine, " ");
    // set counter to 0
    i = 0;
    while (token != nullptr) {
        if(i == 0)
            nnodes = atoi(token);
        else if (i == 1)
            {}//nedges = atoi(token);
        else if(i == 2)
            ntriangles = atoi(token);
        else if(i == 3)
            ntetras = atoi(token);
        else
            break;
        token = std::strtok(nullptr, " ");
        i++;
    }
    // create new domain
    (*dNew)->resizeDofManagers(nnodes);

    //one empty line
    std::getline(inputStream, line);  
    // create nodes
    // read dofs
    const IntArray dofIDArrayPtr = domain->giveDefaultNodeDofIDArry(); ///@todo This is bad idea, dofs can be autogenerated. Don't fix the dofs. This method should be removed.
    int ndofs = dofIDArrayPtr.giveSize();
    // loop over number of nodes, read coordinates and create new nodes
    for ( int inode = 1; inode <= nnodes; inode++ ) {
        FloatArray coords(3);
        std::getline(inputStream, line);  
        //convert const char to char in order to use strtok
        currentLine = new char[line.size() + 1];
        std::strcpy ( currentLine, line.c_str() );
        // tokenizing line
        token = std::strtok(currentLine, " ");
        // set counter to 0
        i = 0;
        while (token != nullptr) {
            if(i == 0)
                {}//int nodeNum = atoi(token);
            else if (i == 1)
                coords.at(1) = atof(token);
            else if(i == 2)
                coords.at(2) = atof(token);
            else if(i == 3)
                coords.at(3) = atof(token);
            else
                break;
            token = std::strtok(NULL, " ");
            i++;
        }
        // newly created node
        auto node = std::make_unique<Node>(inode, * dNew);
        //create new node with default DOFs
        node->setNumberOfDofs(ndofs);
        node->setCoordinates( coords );
        // how to set boundary condition ??
        //      node->setBoundaryFlag( mesh->giveNode(inode)->isBoundary() );
        // ( * dNew )->setDofManager(inode, std::move(node));

    }//end loop over nodes

    //read empty line 
    std::getline(inputStream, line);  


    Element *parentElementPtr = domain->giveElement(1); //??km??
    // loop over triangles, read dofman numbers and create new elements
    for (int itriangle = 1; itriangle <= ntriangles; itriangle++) {
        int elemNumber = 0;
        IntArray dofManagers(3);
        std::getline(inputStream, line); 
        //convert const char to char in order to use strtok
        currentLine = new char[line.size() + 1];
        std::strcpy ( currentLine, line.c_str() );
        // tokenizing line
        token = std::strtok(currentLine, " ");
        // set counter to 0
        i = 0;
        while (token != nullptr) {
            if(i == 0)
                elemNumber = atoi(token);
            else if (i >= 1 && i<=4)
                dofManagers.at(i) = atoi(token);
            else
                break;
            token = std::strtok(nullptr, " ");
            i++;
        } 
        auto elem = classFactory.createElement(parentElementPtr->giveClassName(), elemNumber, * dNew);
        elem->setDofManagers( dofManagers );
        elem->setMaterial( parentElementPtr->giveMaterial()->giveNumber() );
        elem->setCrossSection( parentElementPtr->giveCrossSection()->giveNumber() );
        //
        // ...
        //

    }//end loop over triangles 


    // loop over tetras, read dofman numbers and create new elements
    for (int itetra = 1; itetra <= ntetras; itetra++) {
        int elemNumber = 0;
        IntArray dofManagers(4);
        std::getline(inputStream, line); 
        //convert const char to char in order to use strtok
        currentLine = new char[line.size() + 1];
        std::strcpy ( currentLine, line.c_str() );
        // tokenizing line
        token = std::strtok(currentLine, " ");
        // set counter to 0
        i = 0;
        while (token != nullptr) {
            if(i == 0)
                elemNumber = atoi(token);
            else if (i >= 1 && i<=4)
                dofManagers.at(i) = atoi(token);
            else
                break;
            token = std::strtok(nullptr, " ");
            i++;
        }
        //elem = classFactory.createElement(parentElementPtr->giveClassName(), elemNumber, * dNew);
        auto elem = classFactory.createElement("LTRSpace", elemNumber, * dNew);
        elem->setDofManagers( dofManagers );
        elem->setMaterial( parentElementPtr->giveMaterial()->giveNumber() );
        elem->setCrossSection( parentElementPtr->giveCrossSection()->giveNumber() );
        (*dNew)->setElement(elemNumber, std::move(elem));
    }//end loop over tetras


    std::string name;

    // copy of crossections from old domain
    int ncrosssect = domain->giveNumberOfCrossSectionModels();
    ( * dNew )->resizeCrossSectionModels(ncrosssect);
    for ( int i = 1; i <= ncrosssect; i++ ) {
        DynamicInputRecord ir;
        domain->giveCrossSection(i)->giveInputRecord(ir);
        ir.giveRecordKeywordField(name);

        auto crossSection = classFactory.createCrossSection(name.c_str(), i, * dNew);
        crossSection->initializeFrom(ir);
        ( * dNew )->setCrossSection(i, std::move(crossSection));
    }

    // copy of materials  from old domain
    int nmat = domain->giveNumberOfMaterialModels();
    ( * dNew )->resizeMaterials(nmat);
    for ( int i = 1; i <= nmat; i++ ) {
        DynamicInputRecord ir;
        domain->giveMaterial(i)->giveInputRecord(ir);
        ir.giveRecordKeywordField(name);

        auto mat = classFactory.createMaterial(name.c_str(), i, * dNew);
        mat->initializeFrom(ir);
        ( * dNew )->setMaterial(i, std::move(mat));
    }

    // copy of crossections from old domain
    int nbarriers = domain->giveNumberOfNonlocalBarriers();
    ( * dNew )->resizeNonlocalBarriers(nbarriers);
    for ( int i = 1; i <= nbarriers; i++ ) {
        DynamicInputRecord ir;
        domain->giveNonlocalBarrier(i)->giveInputRecord(ir);
        ir.giveRecordKeywordField(name);

        auto barrier = classFactory.createNonlocalBarrier(name.c_str(), i, * dNew);
        barrier->initializeFrom(ir);
        ( * dNew )->setNonlocalBarrier(i, std::move(barrier));
    }

    // copy of boundary conditions from old domain
    int nbc = domain->giveNumberOfBoundaryConditions();
    ( * dNew )->resizeBoundaryConditions(nbc);
    for ( int i = 1; i <= nbc; i++ ) {
        DynamicInputRecord ir;
        domain->giveBc(i)->giveInputRecord(ir);
        ir.giveRecordKeywordField(name);

        auto bc = classFactory.createBoundaryCondition(name.c_str(), i, * dNew);
        bc->initializeFrom(ir);
        ( * dNew )->setBoundaryCondition(i, std::move(bc));
    }

    // copy of initial conditions from old domain
    int nic = domain->giveNumberOfInitialConditions();
    ( * dNew )->resizeInitialConditions(nic);
    for ( int i = 1; i <= nic; i++ ) {
        DynamicInputRecord ir;
        domain->giveIc(i)->giveInputRecord(ir);
        ir.giveRecordKeywordField(name);

        auto ic = std::make_unique<InitialCondition>(i, *dNew);
        ic->initializeFrom(ir);
        ( * dNew )->setInitialCondition(i, std::move(ic));
    }

    // copy of load time functions from old domain
    //  int nltf = domain->giveNumberOfLoadTimeFunctions();
    int nltf = domain->giveNumberOfFunctions();
    // ( * dNew )->resizeLoadTimeFunctions(nltf);
    ( * dNew )->resizeFunctions(nltf);
    for ( int i = 1; i <= nltf; i++ ) {
        DynamicInputRecord ir;
        //domain->giveLoadTimeFunction(i)->giveInputRecord(ir);
        domain->giveFunction(i)->giveInputRecord(ir);
        ir.giveRecordKeywordField(name);
        
        //ltf = classFactory.createLoadTimeFunction(name.c_str(), i, * dNew);
        auto ltf = classFactory.createFunction(name.c_str(), i, * dNew);
        ltf->initializeFrom(ir);
        //( * dNew )->setLoadTimeFunction(i, ltf);
        ( * dNew )->setFunction(i, std::move(ltf));
    }

    // copy output manager settings from old domain
    ( * dNew )->giveOutputManager()->beCopyOf( domain->giveOutputManager() );  
    return 1;
}

// used by HTS element to create export mesh to vtk

int
T3DInterface :: createInput(Element *e, char *t3dInFile)
{ 
  FILE *inputStrem;
  inputStrem = fopen(t3dInFile, "w");
 
  int nnodes = e->giveNumberOfDofManagers();
  // loop over nodes and write vertexes into t3d input file
  for ( int i = 1; i <= nnodes; i++ ) {
    Node *inode;
    inode = e->giveNode(i);
    fprintf(inputStrem, "vertex %d xyz %e %e %e\n", i, inode->giveCoordinate(1), inode->giveCoordinate(2), inode->giveCoordinate(3));
  }
  for (int i = 1; i<=nnodes; i++) {
    if( i< nnodes)
      fprintf(inputStrem, "curve %d order 2 vertex %d %d\n", i, i, i+1);
    else
      fprintf(inputStrem, "curve %d order 2 vertex %d %d\n", i, i, 1);
  }
  fprintf(inputStrem, "patch 1 normal 0 0 1 boundary curve");
  for (int i = 1; i<=nnodes; i++) {
    fprintf(inputStrem, " %d", i);
  }
  fclose(inputStrem);
  
  return 1;

}







int
T3DInterface :: createVTKExportMesh(const char *t3dOutFile,std::vector<FloatArray> &nodeCoords, std::vector<IntArray> &cellNodes, IntArray &cellTypes )
{

   std::ifstream inputStream;
  inputStream.open( t3dOutFile );
  if ( !inputStream.is_open() ) {
    OOFEM_ERROR("OOFEMTXTDataReader::OOFEMTXTDataReader: Can't open T3D input stream (%s)", t3dOutFile);
    return 0;
  }
  std::string line;
  //read and skip first line from t3d out file - 4 numbers
  std::getline(inputStream, line);
  /*read second line from t3d out file
    4 numbers
    1 - number of nodes
    2 - number of edges
    3 - number of triangles
    4 - number of tetras
  */
  int nnodes = 0, ntriangles = 0; // nedges, ntetras
  std::getline(inputStream, line);
  //convert const char to char in order to use strtok
  char *currentLine = new char[line.size() + 1];
  std::strcpy ( currentLine, line.c_str() );
  // tokenizing line
  char *token = std::strtok(currentLine, " ");
  // set counter to 0
  int i = 0;
  while (token != NULL) {
    if(i == 0)
      nnodes = atoi(token);
    else if (i == 1)
	{}//nedges = atoi(token);
    else if(i == 2)
      ntriangles = atoi(token);
    else if(i == 3)
	{}//ntetras = atoi(token);
    else
      break;
    token = std::strtok(NULL, " ");
    i++;
  }  
  delete [] currentLine;
  //read one empty line
  std::getline(inputStream, line);  


 // loop over number of nodes, read coordinates and put them into nodeCoords Array
  for ( int inode = 1; inode <= nnodes; inode++ ) {
      // int nodeNum;
      FloatArray coords(3);
      std::getline(inputStream, line);  
      //convert const char to char in order to use strtok
      currentLine = new char[line.size() + 1];
      std::strcpy ( currentLine, line.c_str() );
      // tokenizing line
      token = std::strtok(currentLine, " ");
      // set counter to 0
      i = 0;
      while (token != NULL) {
	if(i == 0)
	    {}// nodeNum = atoi(token);
	else if (i == 1)
	  coords.at(1) = atof(token);
	else if(i == 2)
	  coords.at(2) = atof(token);
	else if(i == 3)
	  coords.at(3) = atof(token);
	else
	  break;
	token = std::strtok(NULL, " ");
	i++;
      }  

      nodeCoords.push_back(coords);  
      delete [] currentLine;
  }//end loop over nodes



  //one empty line
  std::getline(inputStream, line);  
  cellTypes.resize(ntriangles);
  // loop over triangles, and fill nodeCoords, cellNodes nad cellTypes arrays
  for (int itriangle = 1; itriangle <= ntriangles; itriangle++) {
    // int elemNumber;
    IntArray dofManagers(3);
    std::getline(inputStream, line); 
    //convert const char to char in order to use strtok
    currentLine = new char[line.size() + 1];
    std::strcpy ( currentLine, line.c_str() );
    // tokenizing line
    token = std::strtok(currentLine, " ");
    // set counter to 0
    i = 0;
    while (token != NULL) {
      if(i == 0)
	  {}//elemNumber = atoi(token);
      else if (i >= 1 && i<=3)
	    dofManagers.at(i) = atoi(token) - 1;
      else
	break;
      token = std::strtok(NULL, " ");
      i++;
    } 
    cellNodes.push_back(dofManagers);
    // in vtk number 5 corresponds to linear triangle
    cellTypes.at(itriangle) = 5;
    delete [] currentLine;
  }//end loop over triangles 
  
  return 1;
}





int
T3DInterface :: createQCInterpolationMesh(const char *t3dOutFile,std::vector<FloatArray> &nodeCoords, std::vector<IntArray> &cellNodes, IntArray &cellTypes )
//
// create interpolation mest from t3dOutFile
// node coordinates and element nodes are saved in "matrix" vector<FloatArray/IntArray> 
//
{

   std::ifstream inputStream;
  inputStream.open( t3dOutFile );
  if ( !inputStream.is_open() ) {
    OOFEM_ERROR("OOFEMTXTDataReader::OOFEMTXTDataReader: Can't open T3D input stream (%s)", t3dOutFile);
    return 0;
  }
  std::string line;
  //read and skip first line from t3d out file - 4 numbers
  std::getline(inputStream, line);
  /*read second line from t3d out file
    4 numbers
    1 - number of nodes
    2 - number of edges
    3 - number of triangles
    4 - number of tetras
  */
  int nnodes = 0, ntriangles = 0, ntetras = 0; //nedges,
  std::getline(inputStream, line);
  //convert const char to char in order to use strtok
  char *currentLine = new char[line.size() + 1];
  std::strcpy ( currentLine, line.c_str() );
  // tokenizing line
  char *token = std::strtok(currentLine, " ");
  // set counter to 0
  int i = 0;
  while (token != NULL) {
    if(i == 0)
      nnodes = atoi(token);
    else if (i == 1)
	{}//nedges = atoi(token);
    else if(i == 2)
      ntriangles = atoi(token);
    else if(i == 3)
      ntetras = atoi(token);
    else
      break;
    token = std::strtok(NULL, " ");
    i++;
  }  
  delete [] currentLine;
  //read one empty line
  std::getline(inputStream, line);  

  // check the number of interpolation element
  if (ntriangles!=0 && ntetras!=0 ) {
    OOFEM_ERROR( "T3DInterface: 2D and 3D interpolation elements are not supported together");  
  } else if (ntriangles!=0) {
    cellTypes.resize(ntriangles);
  } else if (ntetras!=0) {
    cellTypes.resize(ntetras);
  } else { 
	OOFEM_ERROR( "T3DInterface: No interpolation element found in %s", t3dOutFile);  
  } 
  
  // loop over number of nodes, read coordinates and put them into nodeCoords Array
  for ( int inode = 1; inode <= nnodes; inode++ ) {
      FloatArray coords(3);
      // int nodeNum;
      std::getline(inputStream, line);  
      //convert const char to char in order to use strtok
      currentLine = new char[line.size() + 1];
      std::strcpy ( currentLine, line.c_str() );
      // tokenizing line
      token = std::strtok(currentLine, " ");
      // set counter to 0
      i = 0;
      while (token != NULL) {
	if(i == 0)
	    {}// nodeNum = atoi(token);
	else if (i == 1)
	  coords.at(1) = atof(token);
	else if(i == 2)
	  coords.at(2) = atof(token);
	else if(i == 3)
	  coords.at(3) = atof(token);
	else
	  break;
	token = std::strtok(NULL, " ");
	i++;
      }  

      nodeCoords.push_back(coords);  
      delete [] currentLine;
  }//end loop over nodes

  //one empty line
  std::getline(inputStream, line);  

// loop over triangles, and fill nodeCoords, cellNodes nad cellTypes arrays
  for (int itriangle = 1; itriangle <= ntriangles; itriangle++) {
    IntArray dofManagers(3);
    // int elemNumber;
    std::getline(inputStream, line); 
    //convert const char to char in order to use strtok
    currentLine = new char[line.size() + 1];
    std::strcpy ( currentLine, line.c_str() );
    // tokenizing line
    token = std::strtok(currentLine, " ");
    // set counter to 0
    i = 0;
    while (token != NULL) {
      if(i == 0)
	  {}//elemNumber = atoi(token);
      else if (i >= 1 && i<=3)
	    dofManagers.at(i) = atoi(token);
      else
	break;
      token = std::strtok(NULL, " ");
      i++;
    } 
    cellNodes.push_back(dofManagers);
    // in vtk number 5 corresponds to linear triangle
    cellTypes.at(itriangle) = 5;
    delete [] currentLine;
  }//end loop over triangles 

// loop over tetras, read dofman numbers and fill nodeCoords, cellNodes nad cellTypes arrays
  for (int itetra = 1; itetra <= ntetras; itetra++) {
    IntArray dofManagers(4);
    // int elemNumber;
    std::getline(inputStream, line); 
    //convert const char to char in order to use strtok
    currentLine = new char[line.size() + 1];
    std::strcpy ( currentLine, line.c_str() );
    // tokenizing line
    token = std::strtok(currentLine, " ");
       // set counter to 0
    i = 0;
    while (token != NULL) {
      if(i == 0)
	  {}//elemNumber = atoi(token);
      else if (i >= 1 && i<=4)
	    dofManagers.at(i) = atoi(token);
      else
	break;
      token = std::strtok(NULL, " ");
      i++;
    } 
    cellNodes.push_back(dofManagers);
    // in vtk number 5 corresponds to linear triangle  
    cellTypes.at(itetra) = 5;    // ??km?? TODO: number for tetras in vtk = ???
    delete [] currentLine;
  }//end loop over tetras
 
  
  return 1;
}


  
} // end namespace oofem
