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

#include "vtkxmlexportmodule.h"
#include "element.h"
#include "gausspoint.h"
#include "timestep.h"
#include "engngm.h"
#include "node.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "cltypes.h"
#include "material.h"
#include "classfactory.h"
#include "crosssection.h"

//#include "xfemmanager.h"
#include "enrichmentitem.h"
#include "enrichmentdomain.h"

#include <string>
#include <sstream>
#include <fstream>
#include <ctime>

#ifdef __VTK_MODULE
 #include <vtkPoints.h>
 #include <vtkPointData.h>
 #include <vtkDoubleArray.h>
 #include <vtkCellArray.h>
 #include <vtkCellData.h>
 #include <vtkXMLUnstructuredGridWriter.h>
 #include <vtkXMLPUnstructuredGridWriter.h>
 #include <vtkUnstructuredGrid.h>
 #include <vtkSmartPointer.h>
#endif

namespace oofem {
REGISTER_ExportModule(VTKXMLExportModule)

VTKXMLExportModule :: VTKXMLExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{
    primVarSmoother = NULL;
    smoother = NULL;
    redToFull.setValues(6, 1, 5, 9, 6, 3, 2); //position of xx, yy, zz, yz, xz, xy in tensor

#ifdef __VTK_MODULE
    //this->intVarArray = vtkSmartPointer<vtkDoubleArray>::New();
#endif
}


VTKXMLExportModule :: ~VTKXMLExportModule()
{
    if ( this->smoother ) {
        delete this->smoother;
    }

    if ( this->primVarSmoother ) {
        delete this->primVarSmoother;
    }
}


IRResultType
VTKXMLExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    ExportModule :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, _IFT_VTKXMLExportModule_cellvars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_VTKXMLExportModule_vars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_VTKXMLExportModule_primvars); // Macro - see unknowntype.h
    IR_GIVE_OPTIONAL_FIELD(ir, ipInternalVarsToExport, _IFT_VTKXMLExportModule_ipvars); // Macro - see internalstatetype.h
    if(ipInternalVarsToExport.giveSize()){
        this->tstep_substeps_out_flag = true;
    }

    val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_VTKXMLExportModule_stype); // Macro
    stype = ( NodalRecoveryModel :: NodalRecoveryModelType ) val;

    regionsToSkip.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, regionsToSkip, _IFT_VTKXMLExportModule_regionstoskip); // Macro

    this->nvr = 0; // number of virtual regions
    IR_GIVE_OPTIONAL_FIELD(ir, nvr, _IFT_VTKXMLExportModule_nvr); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, vrmap, _IFT_VTKXMLExportModule_vrmap); // Macro

    return IRRT_OK;
}


void
VTKXMLExportModule :: initialize()
{
    if ( this->smoother ) {
        delete this->smoother;
        this->smoother = NULL;
    }
}


void
VTKXMLExportModule :: terminate()
{ }


void
VTKXMLExportModule :: makeFullForm(FloatArray &answer, const FloatArray &reducedForm)
{
    answer.resize(9);
    answer.zero();

    for ( int i = 1; i <= reducedForm.giveSize(); i++ ) {
        answer.at( redToFull.at(i) ) = reducedForm.at(i);
    }

    // Symmetrize
    answer.at(4) = answer.at(2);
    answer.at(7) = answer.at(3);
    answer.at(8) = answer.at(6);
}




std :: string
VTKXMLExportModule :: giveOutputFileName(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".vtu";
}


FILE *
VTKXMLExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;
    std :: string fileName = giveOutputFileName(tStep);
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2( "VTKXMLExportModule::giveOutputStream: failed to open file %s", fileName.c_str() );
    }

    return answer;
}

int
VTKXMLExportModule :: giveCellType(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int vtkCellType = 0;

    if ( elemGT == EGT_point ) {
        vtkCellType = 1;
    } else if ( elemGT == EGT_line_1 ) {
        vtkCellType = 3;
    } else if ( elemGT == EGT_line_2 ) {
        vtkCellType = 21;
    } else if ( elemGT == EGT_triangle_1 ) {
        vtkCellType = 5;
    } else if ( elemGT == EGT_triangle_2 ) {
        vtkCellType = 22;
    } else if ( elemGT == EGT_tetra_1 ) {
        vtkCellType = 10;
    } else if ( elemGT == EGT_tetra_2 ) {
        vtkCellType = 24;
    } else if ( elemGT == EGT_quad_1 || elemGT == EGT_quad_1_interface ) {
        vtkCellType = 9;
    } else if ( elemGT == EGT_quad_21_interface ) {
        vtkCellType = 30;
    } else if ( elemGT == EGT_quad_2 ) {
        vtkCellType = 23;
    } else if ( elemGT == EGT_quad9_2 ) {
        vtkCellType = 23;
    } else if ( elemGT == EGT_hexa_1 ) {
        vtkCellType = 12;
    } else if ( elemGT == EGT_hexa_2 ) {
        vtkCellType = 25;
    } else if ( elemGT == EGT_hexa_27 ) {
        vtkCellType = 29;
    } else if ( elemGT == EGT_wedge_1 ) {
        vtkCellType = 13;
    } else if ( elemGT == EGT_wedge_2 ) {
        vtkCellType = 26;
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported element geometry type on element %d", elem->giveNumber() );
    }

    return vtkCellType;
}

int
VTKXMLExportModule :: giveNumberOfNodesPerCell(int cellType)
{
    switch ( cellType ) {
    case 1:
        return 1;

    case 3:
        return 2;

    case 5:
    case 21:
        return 3;

    case 9:
    case 10:
        return 4;

    case 14:
        return 5;

    case 13:
    case 22:
    case 30:
        return 6;

    case 12:
    case 23:
        return 8;

    case 24:
        return 10;

    case 25:
        return 20;

    case 29:
        return 27;

    default:
        OOFEM_ERROR("VTKXMLExportModule: unsupported cell type ID");
    }

    return 0; // to make compiler happy
}


void
VTKXMLExportModule :: giveElementCell(IntArray &answer, Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int nelemNodes;

    if ( ( elemGT == EGT_point ) ||
         ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
         ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
         ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
         ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
         ( elemGT == EGT_hexa_1 ) ||
         ( elemGT == EGT_wedge_1 ) ) {
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber();
        }
    } else if ( elemGT == EGT_hexa_27 ) {
        int HexaQuadNodeMapping [] = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18, 23, 25, 26, 24, 22, 21, 27
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(HexaQuadNodeMapping [ i - 1 ])->giveNumber();
        }
    } else if ( elemGT == EGT_hexa_2 ) {
        int HexaQuadNodeMapping [] = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(HexaQuadNodeMapping [ i - 1 ])->giveNumber();
        }
    } else if ( elemGT == EGT_wedge_2 ) {
        int WedgeQuadNodeMapping [] = {
            4, 6, 5, 1, 3, 2, 12, 11, 10, 9, 8, 7, 13, 15, 14
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(WedgeQuadNodeMapping [ i - 1 ])->giveNumber();
        }
    } else if ( elemGT == EGT_quad9_2 ) {
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber();
        }
    } else if ( elemGT == EGT_quad_1_interface ) {
        int mapping [] = {
            1, 2, 4, 3
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(mapping [ i - 1 ])->giveNumber();
        }
        //} else if ( elemGT == EGT_quad_21_interface ) {int mapping [] = { 1, 2, 5, 4, 3, 6 };
    } else if ( elemGT == EGT_quad_21_interface ) {
        int mapping [] = {
            1, 3, 2, 5, 6, 4
        };                                                                                /// this is not the same ordering as defined in the VTK reference (typo?)
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(mapping [ i - 1 ])->giveNumber();
        }
    } else {
        OOFEM_ERROR("VTKXMLExportModule: unsupported element geometry type");
    }
}


bool
VTKXMLExportModule :: isElementComposite(Element *elem)
{
    return ( elem->giveGeometryType() == EGT_Composite );
}


int ///@todo intended for composite elements? Currently not in use, remove?
VTKXMLExportModule :: giveNumberOfElementCells(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();

    if ( ( elemGT == EGT_point ) ||
         ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
         ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
         ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
         ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) || ( elemGT == EGT_quad9_2 ) ||
         ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_hexa_2 ) || ( elemGT == EGT_hexa_27 ) ||
         ( elemGT == EGT_wedge_1 ) || ( elemGT == EGT_wedge_2 ) ) {
        return 1;
    } else {
        OOFEM_ERROR("VTKXMLExportModule: unsupported element geometry type");
    }

    return 0;
}


void
VTKXMLExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }



#ifdef __VTK_MODULE
    this->fileStream = vtkSmartPointer< vtkUnstructuredGrid > :: New();
    this->nodes = vtkSmartPointer< vtkPoints > :: New();
    this->elemNodeArray = vtkSmartPointer< vtkIdList > :: New();

#else
    this->fileStream = this->giveOutputStream(tStep);
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);

#endif

    // Write output: VTK header
#ifndef __VTK_MODULE
    fprintf(this->fileStream, "<!-- TimeStep %e Computed %d-%02d-%02d at %02d:%02d:%02d -->\n", tStep->giveIntrinsicTime(), current->tm_year + 1900, current->tm_mon + 1, current->tm_mday, current->tm_hour,  current->tm_min,  current->tm_sec);
    fprintf(this->fileStream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(this->fileStream, "<UnstructuredGrid>\n");
#endif

    this->giveSmoother(); // make sure smoother is created, Necessary? If it doesn't exist it is created /JB

    /* Loop over pieces  ///@todo: this feature has been broken but not checked if it currently works /JB
     * Start default pieces containing all single cell elements. Elements built up from several vtk
     * cells (composite elements) are exported as individual pieces after the default ones.
     */
    int nPiecesToExport = this->smoother->giveNumberOfVirtualRegions(); //old name: region
    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        if ( this->regionsToSkip.contains(pieceNum) ) {
            continue;
        }

        // Fills a data struct (VTKPiece) with all the necessary data.
        this->setupVTKPiece(this->defaultVTKPiece, tStep, pieceNum);


        // Write the VTK piece to file.
        this->writeVTKPiece(this->defaultVTKPiece, tStep);
    }

    /*
     * Output all composite elements - one piece per composite element
     * Each element is responsible of setting up a VTKPiece which can then be exported
     */
    Domain *d  = emodel->giveDomain(1);
    Element *el;
    int nelem = d->giveNumberOfElements();

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        if ( this->regionsToSkip.contains( this->smoother->giveElementVirtualRegionNumber(ielem) ) ) {
            continue;
        }

        el = d->giveElement(ielem);
#ifdef __PARALLEL_MODE
        if ( el->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        if ( this->isElementComposite(el) ) {
#ifndef __VTK_MODULE
            this->exportCompositeElement(this->defaultVTKPiece, el, tStep);
            this->writeVTKPiece(this->defaultVTKPiece, tStep);
#else
            // No support for binary export yet
#endif
        }
    }  // end loop over composite elements

    // Finilize the output:
    std :: string fname = giveOutputFileName(tStep);
#ifdef __VTK_MODULE

 #if 0
    // Code fragment intended for future support of composite elements in binary format
    // Doesn't as well as I would want it to, interface to VTK is to limited to control this.
    // * The PVTU-file is written by every process (seems to be impossible to avoid).
    // * Part files are renamed and time step and everything else is cut off => name collisions
    vtkSmartPointer< vtkXMLPUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLPUnstructuredGridWriter > :: New();
    writer->SetTimeStep(tStep->giveNumber() - 1);
    writer->SetNumberOfPieces( this->emodel->giveNumberOfProcesses() );
    writer->SetStartPiece( this->emodel->giveRank() );
    writer->SetEndPiece( this->emodel->giveRank() );


 #else
    vtkSmartPointer< vtkXMLUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter > :: New();
 #endif

    writer->SetFileName( fname.c_str() );
    writer->SetInput(this->fileStream); // VTK 4
    //writer->SetInputData(this->fileStream); // VTK 6

    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    writer->SetDataModeToAscii();
    writer->Write();
#else
    fprintf(this->fileStream, "</UnstructuredGrid>\n</VTKFile>");
    fclose(this->fileStream);
#endif

    // export raw ip values (if required)
    if ( !this->ipInternalVarsToExport.isEmpty() ) {
        this->exportIntVarsInGpAs(ipInternalVarsToExport, tStep);
    }

    // Write the *.pvd-file. Currently only conatains time step information. It's named "timestep" but is actually the total time.
    // First we check to see that there are more than 1 time steps, otherwise it is redundant;
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() && emodel->giveRank() == 0 ) {
        ///@todo Should use probably use PVTU-files instead. It is starting to get messy.
        // For this to work, all processes must have an identical output file name.
        for ( int i = 0; i < this->emodel->giveNumberOfProcesses(); ++i ) {
            std :: ostringstream pvdEntry;
            char fext [ 100 ];
            if ( this->emodel->giveNumberOfProcesses() > 1 ) {
                sprintf( fext, "_%03d.m%d.%d", i, this->number, tStep->giveNumber() );
            } else {
                sprintf( fext, "m%d.%d", this->number, tStep->giveNumber() );
            }

            pvdEntry << "<DataSet timestep=\"" << tStep->giveIntrinsicTime() << "\" group=\"\" part=\"" << i << "\" file=\""
                     << this->emodel->giveOutputBaseFileName() << fext << ".vtu\"/>";
            this->pvdBuffer.push_back( pvdEntry.str() );
        }

        this->writeVTKCollection();
    } else
#endif
    if ( !emodel->isParallel() && tStep->giveNumber() >= 1 ) { // For non-parallel enabled OOFEM, then we only check for multiple steps.
        std :: ostringstream pvdEntry;
        pvdEntry << "<DataSet timestep=\"" << tStep->giveIntrinsicTime() << "\" group=\"\" part=\"\" file=\"" << fname << "\"/>";
        this->pvdBuffer.push_back( pvdEntry.str() );
        this->writeVTKCollection();
    }
}


void
VTKPiece :: setNumberOfNodes(int numNodes)
{
    this->numNodes = numNodes;
    this->nodeCoords.resize(numNodes);
}

void
VTKPiece :: setNumberOfCells(int numCells)
{
    this->numCells = numCells;
    this->connectivity.resize(numCells);
    this->elCellTypes.resize(numCells);
    this->elOffsets.resize(numCells);
}

void
VTKPiece :: setConnectivity(int cellNum, IntArray &nodes)
{
    this->connectivity [ cellNum - 1 ] = nodes;
}

void
VTKPiece :: setNodeCoords(int nodeNum, FloatArray &coords)
{
    this->nodeCoords [ nodeNum - 1 ] = coords;
}

void
VTKPiece :: setNumberOfPrimaryVarsToExport(int numVars, int numNodes)
{
    this->nodeVars.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeVars [ i - 1 ].resize(numNodes);
    }
}

void
VTKPiece :: setNumberOfInternalVarsToExport(int numVars, int numNodes)
{
    this->nodeVarsFromIS.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeVarsFromIS [ i - 1 ].resize(numNodes);
    }
}

void
VTKPiece :: setNumberOfInternalXFEMVarsToExport(int numVars, int numEnrichmentItems, int numNodes)
{
    this->nodeVarsFromXFEMIS.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeVarsFromXFEMIS [ i - 1 ].resize(numEnrichmentItems);
        for ( int j = 1; j <= numEnrichmentItems; j++ ) {
            this->nodeVarsFromXFEMIS [ i - 1 ] [ j - 1 ].resize(numNodes);
        }
    }
}

void
VTKPiece :: setNumberOfCellVarsToExport(int numVars, int numCells)
{
    this->elVars.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->elVars [ i - 1 ].resize(numCells);
    }
}

void
VTKPiece :: setPrimaryVarInNode(int varNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVars [ varNum - 1 ] [ nodeNum - 1 ] = valueArray;
}

void
VTKPiece :: setInternalVarInNode(int varNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVarsFromIS [ varNum - 1 ] [ nodeNum - 1 ] = valueArray;
}

void
VTKPiece :: setInternalXFEMVarInNode(int varNum, int eiNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVarsFromXFEMIS [ varNum - 1 ] [ eiNum - 1 ] [ nodeNum - 1 ] = valueArray;
}


void
VTKPiece :: setCellVar(int varNum, int cellNum, FloatArray valueArray)
{
    this->elVars [ varNum - 1 ] [ cellNum - 1 ] = valueArray;
}


void
VTKXMLExportModule :: setupVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep, int region)
{
    // Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.

    Domain *d  = emodel->giveDomain(1);
    Element *elem;
    FloatArray *coords;

    this->giveSmoother(); // make sure smoother is created

    // output nodes Region By Region
    int numNodes, numRegionEl;
    IntArray mapG2L, mapL2G;

    // Assemble local->global and global->local region map and get number of
    // single cells to process the composite cells exported individually.
    this->initRegionNodeNumbering(mapG2L, mapL2G, numNodes, numRegionEl, d, region);
#ifndef __PARALLEL_MODE
    if ( numNodes && numRegionEl ) {
#else
    if ( 1 ) {
#endif

        // Export nodes as vtk vertices
        vtkPiece.setNumberOfNodes(numNodes);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            coords = d->giveNode( mapL2G.at(inode) )->giveCoordinates();
            vtkPiece.setNodeCoords(inode, * coords);
        }


        //-------------------------------------------
        // Export all the cell data for the piece
        //-------------------------------------------
        IntArray cellNodes;
        vtkPiece.setNumberOfCells(numRegionEl);

        int offset = 0;
        int cellNum = 0;
        for ( int elNum = 1; elNum <= emodel->giveDomain(1)->giveNumberOfElements(); elNum++ ) {
            elem = d->giveElement(elNum);

            // Skip elements that:
            // are inactivated or of composite type ( these are exported individually later)
            if ( this->isElementComposite(elem) || !elem->isActivated(tStep) ) {
                continue;
            }

#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif

            // skip elements not in active region
            if ( ( region > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(elNum) != region ) ) {
                continue;
            }

            cellNum++;

            // Set the connectivity
            this->giveElementCell(cellNodes, elem);  // node numbering of the cell with according to the VTK format

            // Map from global to local node numbers for the current piece
            //int numElNodes = elem->giveNumberOfNodes();
            int numElNodes = cellNodes.giveSize();
            IntArray connectivity(numElNodes);
            for ( int i = 1; i <= numElNodes; i++ ) {
                connectivity.at(i) = mapG2L.at( cellNodes.at(i) );
            }

            vtkPiece.setConnectivity(cellNum, connectivity);

            vtkPiece.setCellType( cellNum, this->giveCellType(elem) ); // VTK cell type

            //offset += elem->giveNumberOfNodes();
            offset += numElNodes;
            vtkPiece.setOffset(cellNum, offset);
        }


        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(vtkPiece, mapG2L, mapL2G, numNodes, region, tStep);
        this->exportIntVars(vtkPiece, mapG2L, mapL2G, numNodes, region, tStep);

        this->exportCellVars(vtkPiece, numRegionEl, tStep);
    } // end of default piece for simple geometry elements
}


void
VTKXMLExportModule :: writeVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep)
{
    // Write a VTK piece to file. This could be the whole domain (most common case) or it can be a
    // (so-called) composite element consisting of several VTK cells (layered structures, XFEM, etc.).

    if ( !vtkPiece.giveNumberOfCells() ) {
        return;                                  //{ // if there are no elements to output
    }

    // Write output: node coords
    int numNodes = vtkPiece.giveNumberOfNodes();
    int numEl = vtkPiece.giveNumberOfCells();
    FloatArray vtkCoords(3), coords;

#ifdef __VTK_MODULE
    for ( int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPiece.giveNodeCoords(inode);
        vtkCoords.zero();
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            vtkCoords.at(i) = coords.at(i);
        }

        this->nodes->InsertNextPoint( vtkCoords.at(1), vtkCoords.at(2), vtkCoords.at(3) );
        this->fileStream->SetPoints(nodes);
    }

#else
    fprintf(this->fileStream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numNodes, numEl);
    fprintf(this->fileStream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ");

    for ( int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPiece.giveNodeCoords(inode);
        ///@todo move this below into setNodeCoords since it should alwas be 3 components anyway
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            fprintf( this->fileStream, "%e ", coords.at(i) );
        }

        for ( int i = coords.giveSize() + 1; i <= 3; i++ ) {
            fprintf(this->fileStream, "%e ", 0.0);
        }
    }

    fprintf(this->fileStream, "</DataArray>\n</Points>\n");
#endif


    // Write output: connectivity, offsets, cell types

    // output the connectivity data
#ifdef __VTK_MODULE
    this->fileStream->Allocate(numEl);
#else
    fprintf(this->fileStream, "<Cells>\n");
    fprintf(this->fileStream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ");
#endif
    IntArray cellNodes;
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        cellNodes = vtkPiece.giveCellConnectivity(ielem);

#ifdef __VTK_MODULE
        elemNodeArray->Reset();
        elemNodeArray->SetNumberOfIds( cellNodes.giveSize() );
#endif
        for ( int i = 1; i <= cellNodes.giveSize(); i++ ) {
#ifdef __VTK_MODULE
            elemNodeArray->SetId(i - 1, cellNodes.at(i) - 1);
#else
            fprintf(this->fileStream, "%d ", cellNodes.at(i) - 1);
#endif
        }

#ifdef __VTK_MODULE
        this->fileStream->InsertNextCell(vtkPiece.giveCellType(ielem), elemNodeArray);
#else
        fprintf(this->fileStream, " ");
#endif
    }

#ifndef __VTK_MODULE
    fprintf(this->fileStream, "</DataArray>\n");

    // output the offsets (index of individual element data in connectivity array)
    fprintf(this->fileStream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ");

    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        fprintf( this->fileStream, "%d ", vtkPiece.giveCellOffset(ielem) );
    }

    fprintf(this->fileStream, "</DataArray>\n");


    // output cell (element) types
    fprintf(this->fileStream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        fprintf( this->fileStream, "%d ", vtkPiece.giveCellType(ielem) );
    }

    fprintf(this->fileStream, "</DataArray>\n");
    fprintf(this->fileStream, "</Cells>\n");


    ///@todo giveDataHeaders is currently not updated wrt the new structure -> no file names in headers /JB
    std :: string pointHeader, cellHeader;
    this->giveDataHeaders(pointHeader, cellHeader);

    fprintf( this->fileStream, "%s", pointHeader.c_str() );
#endif

    this->writePrimaryVars(vtkPiece);       // Variables availablie in the nodes
    this->writeIntVars(vtkPiece);           // Internal State Type variables smoothed to the nodes

    if ( emodel->giveDomain(1)->hasXfemManager() ) {
        this->writeXFEMVars(vtkPiece);      // XFEM State Type variables associated with XFEM structure
    }

#ifndef __VTK_MODULE
    fprintf(this->fileStream, "</PointData>\n");
    fprintf( this->fileStream, "%s", cellHeader.c_str() );
#endif

    this->writeCellVars(vtkPiece);          // Single cell variables ( if given in the integration points then an average will be exported)

#ifndef __VTK_MODULE
    fprintf(this->fileStream, "</CellData>\n");
    fprintf(this->fileStream, "</Piece>\n");
#endif

    //}


    // Clear object so it can be filled with new data for the next piece
    vtkPiece.clear();
}



#ifndef __VTK_MODULE
void
VTKXMLExportModule :: giveDataHeaders(std :: string &pointHeader, std :: string &cellHeader)
{
    std :: string scalars, vectors, tensors;

    int n = primaryVarsToExport.giveSize();

    UnknownType type;

    for ( int i = 1; i <= n; i++ ) {
        type = ( UnknownType ) primaryVarsToExport.at(i);
        if ( ( type == DisplacementVector ) || ( type == EigenVector ) || ( type == VelocityVector ) || ( type == DirectorField ) ) {
            vectors += __UnknownTypeToString(type);
            vectors.append(" ");
        } else if ( ( type == FluxVector ) || ( type == PressureVector ) || ( type == Temperature ) ) {
            scalars += __UnknownTypeToString(type);
            scalars.append(" ");
        } else {
            OOFEM_ERROR2( "VTKXMLExportModule::exportPrimVarAs: unsupported UnknownType %s", __UnknownTypeToString(type) );
        }
    }

    InternalStateType isttype;
    InternalStateValueType vtype;


    n = internalVarsToExport.giveSize();

    // prepare header
    for ( int i = 1; i <= n; i++ ) {
        isttype = ( InternalStateType ) internalVarsToExport.at(i);
        vtype = giveInternalStateValueType(isttype);

        if ( vtype == ISVT_SCALAR ) {
            scalars += __InternalStateTypeToString(isttype);
            scalars.append(" ");
        } else if ( vtype == ISVT_VECTOR ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else if ( ( vtype == ISVT_TENSOR_S3 ) || ( vtype == ISVT_TENSOR_S3E ) ) {
            tensors += __InternalStateTypeToString(isttype);
            tensors.append(" ");
        } else if ( vtype == ISVT_TENSOR_G ) { //@todo shouldn't this go under tensor?
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else {
            fprintf( stderr, "VTKXMLExportModule::exportIntVars: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
        }
    }

    // print header
    pointHeader = "<PointData Scalars=\"" + scalars + "\" "
                  +  "Vectors=\"" + vectors + "\" "
                  +  "Tensors=\"" + tensors + "\" >\n";


    scalars.clear();
    vectors.clear();
    tensors.clear();
    n = this->cellVarsToExport.giveSize();
    // prepare header
    for ( int i = 1; i <= n; i++ ) {
        isttype = ( InternalStateType ) cellVarsToExport.at(i);
        vtype = giveInternalStateValueType(isttype);

        if ( vtype == ISVT_SCALAR ) {
            scalars += __InternalStateTypeToString(isttype);
            scalars.append(" ");
        } else if ( vtype == ISVT_VECTOR ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else if ( ( vtype == ISVT_TENSOR_S3 ) || ( vtype == ISVT_TENSOR_S3E ) ) {
            tensors += __InternalStateTypeToString(isttype);
            tensors.append(" ");
        } else if ( vtype == ISVT_TENSOR_G ) { //@todo shouldn't this go under tensor?
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else {
            fprintf( stderr, "VTKXMLExportModule::cellVarsToExport: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
        }
    }

    // print header
    cellHeader = "<CellData Scalars=\"" + scalars + "\" "
                 +  "Vectors=\"" + vectors + "\" "
                 +  "Tensors=\"" + tensors + "\" >\n";
}
#endif







//----------------------------------------------------
// Internal variables and XFEM realted fields (keyword "vars" in OOFEM input file)
//----------------------------------------------------
void
VTKXMLExportModule :: exportIntVars(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int numNodes, int region, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    InternalStateType isType;
    FloatArray answer;

    this->giveSmoother()->init(); // Makes sure smoother is up-to-date with potentially new mesh.

    // Export of Internal State Type fields
    vtkPiece.setNumberOfInternalVarsToExport(internalVarsToExport.giveSize(), numNodes);
    for ( int field = 1; field <= internalVarsToExport.giveSize(); field++ ) {
        isType = ( InternalStateType ) internalVarsToExport.at(field);

        for ( int nodeNum = 1; nodeNum <= numNodes; nodeNum++ ) {
            Node *node = d->giveNode( mapL2G.at(nodeNum) );
            this->getNodalVariableFromIS(answer, node, tStep, isType, region);
            vtkPiece.setInternalVarInNode(field, nodeNum, answer);
        }
    }

    // Export of XFEM related quantities
    if ( d->hasXfemManager() ) {
        XfemManager *xFemMan = d->giveXfemManager();
        int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();

        vtkPiece.setNumberOfInternalXFEMVarsToExport(xFemMan->vtkExportFields.giveSize(), nEnrIt, numNodes);
        for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
            XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];

            for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
                for ( int nodeIndx = 1; nodeIndx <= numNodes; nodeIndx++ ) {
                    Node *node = d->giveNode( mapL2G.at(nodeIndx) );
                    getNodalVariableFromXFEMST( answer, node, tStep, xfemstype, region, xFemMan->giveEnrichmentItem(enrItIndex) );
                    vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeIndx, answer);
                }
            }
        }
    }
}


void
VTKXMLExportModule :: getNodalVariableFromIS(FloatArray &answer, Node *node, TimeStep *tStep, InternalStateType type, int ireg)
{
    // Recovers nodal values from Internal States defined in the integration points.
    // Should return an array with proper size supported by VTK (1, 3 or 9)

    this->giveSmoother();
    IntArray redIndx;

    if ( !( type == IST_DisplacementVector || type == IST_MaterialInterfaceVal  ) ) {
        this->smoother->recoverValues(type, tStep);
    }


    const FloatArray *val = NULL;
    FloatArray valueArray;
    InternalStateValueType valType = giveInternalStateValueType(type);

    if ( type == IST_DisplacementVector ) { ///@todo Why does this code exists here? DisplacementVector isn't a internal variable. And if its really desired, then why the special treatment for just displacements?
        ///@todo Is this for elements which does not use displacements as primary variables but e.g. the placement x? Move to exportPrimaryVar? /JB
        valueArray.resize(3);
        val = & valueArray;
        for ( int i = 1; i <= 3; i++ ) {
            valueArray.at(i) = node->giveUpdatedCoordinate(i, tStep, 1.0) - node->giveCoordinate(i);
        }
    } else if ( type == IST_MaterialInterfaceVal ) {
        MaterialInterface *mi = emodel->giveMaterialInterface(1);
        if ( mi ) {
            valueArray.resize(1);
            val = & valueArray;
            valueArray.at(1) = mi->giveNodalScalarRepresentation( node->giveGlobalNumber() );
        }
    } else {
        int found = this->smoother->giveNodalVector(val, node->giveGlobalNumber(), ireg);
        if ( !found ) {
            valueArray.resize( redIndx.giveSize() );
            val = & valueArray;
        }
    }

    int ncomponents = giveInternalStateTypeSize(valType);
    answer.resize(ncomponents);
    int valSize = val->giveSize(); // size of recovered quantity

    // check if valSize corresponds to the expected size otherwise pad with zeros
    if ( valType == ISVT_SCALAR ) {
        answer.at(1) = valSize ? val->at(1) : 0.0;
    } else if ( valType == ISVT_VECTOR ) {
        int isize = min(valSize, 3); // so it will simply truncate larger arrays
        for ( int i = 1; i <= isize; i++ ) {
            answer.at(i) = val->at(i);
        }
    } else if ( valType == ISVT_TENSOR_S3 || valType == ISVT_TENSOR_S3E ) {
        this->makeFullForm(answer, * val);
    } else if ( valType == ISVT_TENSOR_G ) { // export general tensor values as scalars
        int isize = min(val->giveSize(), 9);
        for ( int i = 1; i <= isize; i++ ) {
            answer.at(i) = val->at(i);
        }
    } else {
        OOFEM_ERROR("TKXMLExportModule ::getNodalVariableFromIS - ISVT_UNDEFINED encountered")
    }
}


void
VTKXMLExportModule :: getNodalVariableFromXFEMST(FloatArray &answer, Node *node, TimeStep *tStep, XFEMStateType xfemstype, int ireg, EnrichmentItem *ei)
{
    // Recovers nodal values from XFEM state variables (e.g. levelset function)
    // Should return an array with proper size supported by VTK (1, 3 or 9)
    // This could be moved into EnrichmentItem such that VTKExport doesn't need to know anything about XFEM

    const FloatArray *val = NULL;
    FloatArray valueArray;

    Domain *d = emodel->giveDomain(1);
    XfemManager *xFemMan = d->giveXfemManager();
    InternalStateValueType valType = xFemMan->giveXFEMStateValueType(xfemstype);


    // The xfem level set function is defined in the nodes and recovery of nodal values is trivial.
    if ( xfemstype == XFEMST_LevelSetPhi ) {
        valueArray.resize(1);
        val = & valueArray;
        ei->evalLevelSetNormalInNode( valueArray.at(1), node->giveNumber() );
    } else if ( xfemstype == XFEMST_LevelSetGamma ) {
        valueArray.resize(1);
        val = & valueArray;
        ei->evalLevelSetTangInNode( valueArray.at(1), node->giveNumber() );
    } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
        valueArray.resize(1);
        val = & valueArray;
        ei->evalNodeEnrMarkerInNode( valueArray.at(1), node->giveNumber() );
    } else {
        //OOFEM_WARNING2("VTKXMLExportModule::getNodalVariableFromXFEMST: invalid data in node %d", inode);
    }

    ///@todo duplicated code from getNodalVariableFromIS - uneccessary
    int ncomponents = giveInternalStateTypeSize(valType);
    answer.resize(ncomponents);
    int valSize = val->giveSize(); // size of recovered quantity

    // check if valSize corresponds to the expected size otherwise pad with zeros
    if ( valType == ISVT_SCALAR ) {
        answer.at(1) = valSize ? val->at(1) : 0.0;
    } else if ( valType == ISVT_VECTOR ) {
        int isize = min(valSize, 3); // so it will simply truncate larger arrays
        for ( int i = 1; i <= isize; i++ ) {
            answer.at(i) = val->at(i);
        }
    } else if ( valType == ISVT_TENSOR_S3 || valType == ISVT_TENSOR_S3E ) {
        this->makeFullForm(answer, * val);
    } else if ( valType == ISVT_TENSOR_G ) { // export general tensor values as scalars
        int isize = min(val->giveSize(), 9);
        for ( int i = 1; i <= isize; i++ ) {
            answer.at(i) = val->at(i);
        }
    } else {
        OOFEM_ERROR("TKXMLExportModule ::getNodalVariableFromIS - ISVT_UNDEFINED encountered")
    }
}


void
VTKXMLExportModule :: writeIntVars(VTKPiece &vtkPiece)
{
    int n = internalVarsToExport.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        const char *name = __InternalStateTypeToString(type);
        int numNodes = vtkPiece.giveNumberOfNodes();
        FloatArray valueArray;

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray > :: New();
        varArray->SetName(name);
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.giveInternalVarInNode(i, inode);
            for ( int i = 1; i <= ncomponents; ++i ) {
                varArray->SetComponent( inode - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else
        fprintf(this->fileStream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", name, ncomponents);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.giveInternalVarInNode(i, inode);
            this->writeVTKPointData(valueArray);
        }

#endif





        // Footer
#ifndef __VTK_MODULE
        fprintf(this->fileStream, "</DataArray>\n");
#endif
    }
}

void
VTKXMLExportModule :: writeXFEMVars(VTKPiece &vtkPiece)
{
    Domain *d = emodel->giveDomain(1);
    XfemManager *xFemMan = d->giveXfemManager();
    int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();
    FloatArray valueArray;


    for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
        XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields.at(field);
        const char *namePart = __XFEMStateTypeToString(xfemstype);
        InternalStateValueType valType = xFemMan->giveXFEMStateValueType(xfemstype);
        int ncomponents = giveInternalStateTypeSize(valType);

        int numNodes = vtkPiece.giveNumberOfNodes();
        for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
            // Header
            char name [ 100 ];   // Must I define a fixed size? /JB
            sprintf( name, "%s_%d ", namePart, xFemMan->giveEnrichmentItem(enrItIndex)->giveNumber() );

#ifdef __VTK_MODULE
            vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray > :: New();
            varArray->SetName(name);
            varArray->SetNumberOfComponents(ncomponents);
            varArray->SetNumberOfTuples(numNodes);
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                valueArray = vtkPiece.giveInternalXFEMVarInNode(field, enrItIndex, inode);
                for ( int i = 1; i <= ncomponents; ++i ) {
                    varArray->SetComponent( inode - 1, i - 1, valueArray.at(i) );
                }
            }

            this->writeVTKPointData(name, varArray);
#else
            fprintf(this->fileStream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", name, ncomponents);
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                valueArray = vtkPiece.giveInternalXFEMVarInNode(field, enrItIndex, inode);
                this->writeVTKPointData(valueArray);
            }
            fprintf(this->fileStream, "</DataArray>\n");
#endif
        }
    }
}




//----------------------------------------------------
// Misc. functions
//----------------------------------------------------
#ifdef __VTK_MODULE
void
VTKXMLExportModule :: writeVTKPointData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray)
{
    // Write the data to file
    int ncomponents = varArray->GetNumberOfComponents();
    switch ( ncomponents ) {
    case 1:
        this->fileStream->GetPointData()->SetActiveScalars(name);
        this->fileStream->GetPointData()->SetScalars(varArray);
        break;
    case 3:
        this->fileStream->GetPointData()->SetActiveVectors(name);
        this->fileStream->GetPointData()->SetVectors(varArray);
        break;
    case 9:
        this->fileStream->GetPointData()->SetActiveTensors(name);
        this->fileStream->GetPointData()->SetTensors(varArray);
        break;
    }
}
#else
void
VTKXMLExportModule :: writeVTKPointData(FloatArray &valueArray)
{
    // Write the data to file
    for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
        fprintf( this->fileStream, "%e ", valueArray.at(i) );
    }
}
#endif





#ifdef __VTK_MODULE
void
VTKXMLExportModule :: writeVTKCellData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray)
{
    // Write the data to file
    int ncomponents = varArray->GetNumberOfComponents();
    switch ( ncomponents ) {
    case 1:
        this->fileStream->GetCellData()->SetActiveScalars(name);
        this->fileStream->GetCellData()->SetScalars(varArray);
        break;
    case 3:
        this->fileStream->GetCellData()->SetActiveVectors(name);
        this->fileStream->GetCellData()->SetVectors(varArray);
        break;
    case 9:
        this->fileStream->GetCellData()->SetActiveTensors(name);
        this->fileStream->GetCellData()->SetTensors(varArray);
        break;
    }
}

#else

void
VTKXMLExportModule :: writeVTKCellData(FloatArray &valueArray)
{
    // Write the data to file ///@todo exact copy of writeVTKPointData so remove
    for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
        fprintf( this->fileStream, "%e ", valueArray.at(i) );
    }
}
#endif




int
VTKXMLExportModule :: initRegionNodeNumbering(IntArray &regionG2LNodalNumbers,
                                              IntArray &regionL2GNodalNumbers,
                                              int &regionDofMans, int &regionSingleCells,
                                              Domain *domain, int reg)
{
    // regionG2LNodalNumbers is array with mapping from global numbering to local region numbering.
    // The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).

    // regionL2GNodalNumbers is array with mapping from local to global numbering.
    // The i-th value contains the corresponding global node number.


    int nelem = domain->giveNumberOfElements();
    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    regionG2LNodalNumbers.resize(nnodes);
    regionG2LNodalNumbers.zero();
    regionDofMans = 0;
    regionSingleCells = 0;

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
        if ( ( reg > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != reg ) ) {
            continue;
        }

        if ( this->isElementComposite(element) ) {
            continue;                                    // composite cells exported individually
        }

        if ( !element->isActivated( domain->giveEngngModel()->giveCurrentStep() ) ) {                  //skip inactivated elements
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        regionSingleCells++;
        elemNodes = element->giveNumberOfNodes();
        //  elemSides = element->giveNumberOfSides();

        // determine local region node numbering
        for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            node = element->giveNode(elementNode)->giveNumber();
            if ( regionG2LNodalNumbers.at(node) == 0 ) { // assign new number
                /* mark for assignment. This is done later, as it allows to preserve
                 * natural node numbering.
                 */
                regionG2LNodalNumbers.at(node) = 1;
                regionDofMans++;
            }
        }
    }

    regionL2GNodalNumbers.resize(regionDofMans);

    for ( int i = 1; i <= nnodes; i++ ) {
        if ( regionG2LNodalNumbers.at(i) ) {
            regionG2LNodalNumbers.at(i) = currOffset++;
            regionL2GNodalNumbers.at( regionG2LNodalNumbers.at(i) ) = i;
        }
    }

    return 1;
}








//----------------------------------------------------
// Primary variables - readily available in the nodes
//----------------------------------------------------
void
VTKXMLExportModule :: exportPrimaryVars(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int numNodes, int region, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;
    this->givePrimVarSmoother()->init(); // Makes sure primary smoother is up-to-date with potentially new mesh.

    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numNodes);
    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            DofManager *dman = d->giveNode( mapL2G.at(inode) );

            this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region);
            vtkPiece.setPrimaryVarInNode(i, inode, valueArray);
        }
    }
}


void
VTKXMLExportModule :: getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int ireg)
{
    // This code is not perfect. It should be rewritten to handle all cases more gracefully.
    ///@todo This method needs to be cleaned up - maybe define the common vector types so
    // certain dofid's are associated with them /JB

    IntArray dofIDMask(3);
    int indx, size;
    DofIDItem id;
    const FloatArray *recoveredVal;

    InternalStateType iState = IST_DisplacementVector; // Shouldn't be necessary

    dofIDMask.resize(0);
    if ( ( type == DisplacementVector ) || ( type == EigenVector ) || ( type == VelocityVector ) ) {
        dofIDMask.setValues(3, ( int ) Undef, ( int ) Undef, ( int ) Undef);
        for ( int j = 1; j <= dman->giveNumberOfDofs(); j++ ) {
            id = dman->giveDof(j)->giveDofID();
            if ( ( id == V_u ) || ( id == D_u ) ) {
                dofIDMask.at(1) = id;
            } else if ( ( id == V_v ) || ( id == D_v ) ) {
                dofIDMask.at(2) = id;
            } else if ( ( id == V_w ) || ( id == D_w ) ) {
                dofIDMask.at(3) = id;
            }
        }

        answer.resize(3);
    } else if ( type == FluxVector ) {
        dofIDMask.followedBy(C_1);
        iState = IST_MassConcentration_1;
        answer.resize(1);
    } else if ( type == Temperature ) {
        dofIDMask.followedBy(T_f);
        iState = IST_Temperature;
        answer.resize(1);
    } else if ( type == PressureVector ) {
        dofIDMask.followedBy(P_f);
        iState = IST_Pressure;
        answer.resize(1);
    } else if ( type == DirectorField ) {
        for ( int j = 1; j <= dman->giveNumberOfDofs(); j++ ) {
            id = dman->giveDof(j)->giveDofID();
            if ( ( id == W_u ) || ( id == W_v ) || ( id == W_w ) ) {
                dofIDMask.followedBy(id);
            }

            answer.resize(3);
        }

        iState = IST_DirectorField;
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported unknownType %s", __UnknownTypeToString(type) );
    }

    size = dofIDMask.giveSize();
    answer.zero();

    for ( int j = 1; j <= size; j++ ) {
        id = ( DofIDItem ) dofIDMask.at(j);
        if ( id == Undef ) {
            answer.at(j) = 0.;
        } else if ( iState == IST_DirectorField ) {
            indx = dman->findDofWithDofId( ( DofIDItem ) dofIDMask.at(j) );
            answer.at(j) = dman->giveDof(indx)->giveUnknown(VM_Total, tStep);

            this->givePrimVarSmoother()->recoverValues(iState, tStep); // recover values if not done before
            this->givePrimVarSmoother()->giveNodalVector(recoveredVal, dman->giveNumber(), ireg);
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING2("VTKXMLExportModule :: getDofManPrimaryVariable: recovered variable size mismatch for %d", type);
                answer.at(j) = 0.0;
            }
        } else if ( ( indx = dman->findDofWithDofId(id) ) ) {
            // primary variable available directly in DOF-manager
            answer.at(j) = dman->giveDof(indx)->giveUnknown(VM_Total, tStep);
        } else if ( iState != IST_Undefined ) {
            // primary variable not directly available
            // but equivalent InternalStateType provided
            // in this case use smoother to recover nodal value

            // This can't deal with ValueModeType, and would recover over and over for some vectorial quantities like velocity.
            this->givePrimVarSmoother()->recoverValues(iState, tStep); // recover values if not done before
            this->givePrimVarSmoother()->giveNodalVector(recoveredVal, dman->giveNumber(), ireg);
            // here we have a lack of information about how to convert recovered values to response
            // if the size is compatible we accept it, otherwise give a warning and zero value.
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING2("VTKXMLExportModule :: getDofManPrimaryVariable: recovered variable size mismatch for %d", type);
                answer.at(j) = 0.0;
            }
        }
    }



    InternalStateValueType valType = giveInternalStateValueType(type);
    //rotate back from nodal CS to global CS if applies
    if ( valType == ISVT_VECTOR ) { ///@todo in general, shouldn't this apply for 2nd order tensors as well? /JB
        Node *node = dynamic_cast< Node * >( dman );
        if ( node && node->hasLocalCS() ) {
            answer.rotatedWith(* node->giveLocalCoordinateTriplet(), 't');
        }
    }
}

void
VTKXMLExportModule :: writePrimaryVars(VTKPiece &vtkPiece)
{
    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        int numNodes = vtkPiece.giveNumberOfNodes();
        const char *name = __UnknownTypeToString(type);
        FloatArray valueArray;

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray > :: New();
        varArray->SetName(name);
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.givePrimaryVarInNode(i, inode);
            for ( int i = 1; i <= ncomponents; ++i ) {
                varArray->SetComponent( inode - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else
        fprintf(this->fileStream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", name, ncomponents);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.givePrimaryVarInNode(i, inode);
            this->writeVTKPointData(valueArray);
        }
        fprintf(this->fileStream, "</DataArray>\n");
#endif
    }
}




//----------------------------------------------------
// Cell vars
//----------------------------------------------------

void
VTKXMLExportModule :: exportCellVars(VTKPiece &vtkPiece, int numCells, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;
    InternalStateType type;

    int n = cellVarsToExport.giveSize();
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int field = 1; field <= n; field++ ) {
        type = ( InternalStateType ) cellVarsToExport.at(field);

        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            Element *el = d->giveElement(ielem); ///@todo should be a pointer to an element in the region /JB
#ifdef __PARALLEL_MODE
            if ( el->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            this->getCellVariableFromIS(valueArray, el, type, tStep);
            vtkPiece.setCellVar(field, ielem, valueArray);
        }
    }
}


void
VTKXMLExportModule :: getCellVariableFromIS(FloatArray &answer, Element *el, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix rotMat(3, 3);
    int col = 0;
    FloatArray valueArray, temp;
    IntArray redIndx;

    InternalStateValueType valType = giveInternalStateValueType(type);
    int ncomponents = giveInternalStateTypeSize(valType);

    valueArray.resize(ncomponents);

    switch ( type ) {
    // Special scalars
    case IST_MaterialNumber:
        OOFEM_WARNING1("VTKExportModule - Material numbers are deprecated, outputing cross section number instead...");
    case IST_CrossSectionNumber:
        valueArray.at(1) = ( double ) el->giveCrossSection()->giveNumber();
        break;
    case IST_ElementNumber:
        valueArray.at(1) = ( double ) el->giveNumber();
        break;
    case IST_Pressure: //@todo This case seems redundant, remove? /JB, /// Why this special treatment for pressure? / Mikael
        if ( el->giveNumberOfInternalDofManagers() == 1 ) {
            //IntArray pmask(1); pmask.at(1) = P_f;
            //el->giveInternalDofManager(1)->giveUnknownVector (answer, pmask,EID_ConservationEquation, VM_Total, tStep);
            //valueArray.at(1) = answer.at(1);
        }

        break;

    // Special vectors
    case IST_MaterialOrientation_x:
    case IST_MaterialOrientation_y:
    case IST_MaterialOrientation_z:

        if ( type == IST_MaterialOrientation_x ) {
            col = 1;
        }

        if ( type == IST_MaterialOrientation_y ) {
            col = 2;
        }

        if ( type == IST_MaterialOrientation_z ) {
            col = 3;
        }

        if ( !el->giveLocalCoordinateSystem(rotMat) ) {
            rotMat.zero(); ///@todo shouldn't it be an identity matrix? /JB
        }

        valueArray.beColumnOf(rotMat, col);
        break;

    // Export cell data as average from ip's as default
    default:

        // compute cell average from ip values
        IntegrationRule * iRule = el->giveDefaultIntegrationRulePtr();
        computeIPAverage(temp, iRule, el, type, tStep); // if element has more than one iRule?? /JB

        // Reshape the Voigt vectors to include all components (duplicated if necessary, VTK insists on 9 components for tensors.)
#if 1
        // Is this part necessary now when giveIPValue returns full form? Only need to symmetrize in case of 6 components /JB
        if ( ncomponents == 9 && temp.giveSize() != 9 ) { // If it has 9 components, then it is assumed to be proper already.
            this->makeFullForm(valueArray, temp);
        } else if ( valType == ISVT_VECTOR && temp.giveSize() < 3 ) {
            valueArray.setValues(3,
                                 temp.giveSize() > 1 ? temp.at(1) : 0.0,
                                 temp.giveSize() > 2 ? temp.at(2) : 0.0,
                                 0.0);
        } else if ( ncomponents != temp.giveSize() ) { // Trying to gracefully handle bad cases, just output zeros.
            valueArray.resize(9);
        }

#endif
    }


    answer = valueArray;
}


void
VTKXMLExportModule :: writeCellVars(VTKPiece &vtkPiece)
{
    FloatArray valueArray;
    int numCells = vtkPiece.giveNumberOfCells();
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        const char *name = __InternalStateTypeToString(type);

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >cellVarsArray = vtkSmartPointer< vtkDoubleArray > :: New();
        cellVarsArray->SetName(name);
        cellVarsArray->SetNumberOfComponents(ncomponents);
        cellVarsArray->SetNumberOfTuples(numCells);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(i, ielem);
            for ( int i = 1; i <= ncomponents; ++i ) {
                cellVarsArray->SetComponent( ielem - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKCellData(name, cellVarsArray);

#else
        fprintf(this->fileStream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", name, ncomponents);
        valueArray.resize(ncomponents);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(i, ielem);
            this->writeVTKCellData(valueArray);
        }
        fprintf(this->fileStream, "</DataArray>\n");
#endif
    }
}


void
VTKXMLExportModule :: computeIPAverage(FloatArray &answer, IntegrationRule *iRule, Element *elem, InternalStateType isType, TimeStep *tStep)
{
    // Computes the volume average (over an element) for the quantity defined by isType
    double gptot = 0.0;
    answer.resize(0);
    FloatArray temp;
    if ( iRule ) {
        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); ++i ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(i);
            elem->giveIPValue(temp, ip, isType, tStep);
            gptot += ip->giveWeight();
            answer.add(ip->giveWeight(), temp);
        }

        answer.times(1. / gptot);
    }
}





void
VTKXMLExportModule :: writeVTKCollection()
{
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);
    char buff [ 1024 ];

    std :: string fname = this->emodel->giveOutputBaseFileName() + ".pvd";

    std :: ofstream outfile( fname.c_str() );

    sprintf(buff, "<!-- Computation started %d-%02d-%02d at %02d:%02d:%02d -->\n", current->tm_year + 1900, current->tm_mon + 1, current->tm_mday, current->tm_hour,  current->tm_min,  current->tm_sec);
    //     outfile << buff;

    outfile << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for ( std :: list< std :: string > :: iterator it = this->pvdBuffer.begin(); it != this->pvdBuffer.end(); ++it ) {
        outfile << * it << "\n";
    }

    outfile << "</Collection>\n</VTKFile>";

    outfile.close();
}






// Export of composite elements

void VTKXMLExportModule :: exportCompositeElement(VTKPiece &vtkPiece, Element *el, TimeStep *tStep)
{
    VTKXMLExportModuleElementInterface *interface =
        static_cast< VTKXMLExportModuleElementInterface * >( el->giveInterface(VTKXMLExportModuleElementInterfaceType) );
    if ( interface ) {
        interface->giveCompositeExportData(vtkPiece, this->primaryVarsToExport, this->internalVarsToExport, this->cellVarsToExport, tStep);

        //this->writeVTKPiece(this->defaultVTKPiece, tStep);
    }
}


void
VTKPiece :: clear()
{
    ///@todo Will this give a memory leak? / JB
    numCells = 0;
    numNodes = 0;
    this->connectivity.resize(0);
    this->elCellTypes.resize(0);
    this->elOffsets.resize(0);
    this->elVars.resize(0);
    this->nodeCoords.resize(0);
    this->nodeVars.resize(0);
    this->nodeVarsFromIS.resize(0);
    this->nodeVarsFromXFEMIS.resize(0);
}



NodalRecoveryModel *
VTKXMLExportModule :: giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->smoother == NULL ) {
        this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
        this->smoother->setRecoveryMode(nvr, vrmap);
    }

    return this->smoother;
}


NodalRecoveryModel *
VTKXMLExportModule :: givePrimVarSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->primVarSmoother == NULL ) {
        this->primVarSmoother = classFactory.createNodalRecoveryModel(NodalRecoveryModel :: NRM_NodalAveraging, d);
        this->primVarSmoother->setRecoveryMode(nvr, vrmap);
    }

    return this->primVarSmoother;
}


void
VTKXMLExportModule :: exportIntVarsInGpAs(IntArray valIDs, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nip = 0;
    int j, k, nc = 0;
    int nelem = d->giveNumberOfElements();
    FloatArray *lc, gc, value;
    FILE *stream;
    InternalStateType isttype;
    InternalStateValueType vtype;
    std :: string scalars, vectors, tensors;

    // output nodes Region By Region
    int nregions = this->smoother->giveNumberOfVirtualRegions();
    // open output stream
    std :: string outputFileName = this->giveOutputBaseFileName(tStep) + ".gp.vtu";
    if ( ( stream = fopen(outputFileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2( "VTKXMLExportModule::exportIntVarsInGpAs: failed to open file %s", outputFileName.c_str() );
    }

    fprintf(stream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(stream, "<UnstructuredGrid>\n");

    /* loop over regions */
    for ( int ireg = 1; ireg <= nregions; ireg++ ) {
        if ( this->regionsToSkip.contains(ireg) ) {
            continue;
        }

        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            if ( this->regionsToSkip.contains(d->giveElement(ielem)->giveRegionNumber() ) ) {
                continue;
            }

            nip += d->giveElement(ielem)->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
        }

        //Create one cell per each GP
        fprintf(stream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nip, nip);
        fprintf(stream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ");
        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            if ( this->regionsToSkip.contains(d->giveElement(ielem)->giveRegionNumber() ) ) {
                continue;
            }

            int enip = d->giveElement(ielem)->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
            for ( j = 0; j < enip; j++ ) {
                lc = d->giveElement(ielem)->giveDefaultIntegrationRulePtr()->getIntegrationPoint(j)->giveCoordinates();
                d->giveElement(ielem)->computeGlobalCoordinates(gc, * lc);
                for ( k = 1; k <= gc.giveSize(); k++ ) {
                    fprintf( stream, "%e ", gc.at(k) );
                }

                for ( k = gc.giveSize() + 1; k <= 3; k++ ) {
                    fprintf(stream, "%e ", 0.0);
                }
            }
        }

        fprintf(stream, " </DataArray>\n");
        fprintf(stream, "</Points>\n");
        fprintf(stream, "<Cells>\n");
        fprintf(stream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
        for ( j = 0; j < nip; j++ ) {
            fprintf(stream, "%d ", j);
        }

        fprintf(stream, " </DataArray>\n");
        fprintf(stream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
        for ( j = 1; j <= nip; j++ ) {
            fprintf(stream, "%d ", j);
        }

        fprintf(stream, " </DataArray>\n");
        fprintf(stream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
        for ( j = 1; j <= nip; j++ ) {
            fprintf(stream, "1 ");
        }

        fprintf(stream, " </DataArray>\n");
        fprintf(stream, "</Cells>\n");
        // prepare the data header
        for ( int vi = 1; vi <= valIDs.giveSize(); vi++ ) {
            isttype = ( InternalStateType ) valIDs.at(vi);
            vtype = giveInternalStateValueType(isttype);

            if ( vtype == ISVT_SCALAR ) {
                scalars += __InternalStateTypeToString(isttype);
                scalars.append(" ");
            } else if ( vtype == ISVT_VECTOR ) {
                vectors += __InternalStateTypeToString(isttype);
                vectors.append(" ");
            } else if ( ( vtype == ISVT_TENSOR_S3 ) || ( vtype == ISVT_TENSOR_S3E ) ) {
                tensors += __InternalStateTypeToString(isttype);
                tensors.append(" ");
            } else if ( vtype == ISVT_TENSOR_G ) {
                vectors += __InternalStateTypeToString(isttype);
                vectors.append(" ");
            } else {
                fprintf( stderr, "VTKXMLExportModule::exportIntVarsInGpAs: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }
        }

        // print collected data summary in header
        fprintf( stream, "<PointData Scalars=\"%s\" Vectors=\"%s\" Tensors=\"%s\" >\n",
                 scalars.c_str(), vectors.c_str(), tensors.c_str() );

        // export actual data, loop over individual IDs to export
        for ( int vi = 1; vi <= valIDs.giveSize(); vi++ ) {
            isttype = ( InternalStateType ) valIDs.at(vi);
            vtype = giveInternalStateValueType(isttype);
            if ( vtype == ISVT_SCALAR ) {
                nc = 1;
            } else if ( vtype == ISVT_VECTOR ) {
                nc = 3;
            } else if ( ( vtype == ISVT_TENSOR_S3 ) || ( vtype == ISVT_TENSOR_S3E ) ) {
                nc = 9;
            } else if ( vtype == ISVT_TENSOR_G ) {
                nc = 9;
            } else {
                fprintf( stderr, "VTKXMLExportModule::exportIntVarsInGpAs: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }

            fprintf(stream, "  <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">", __InternalStateTypeToString(isttype), nc);
            for ( int ielem = 1; ielem <= nelem; ielem++ ) {
                if ( this->regionsToSkip.contains(d->giveElement(ielem)->giveRegionNumber() ) ) {
                    continue;
                }

                nip = d->giveElement(ielem)->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
                // loop over default IRule gps
                for ( int ip = 0; ip < nip; ip++ ) {
                    d->giveElement(ielem)->giveIPValue(value, d->giveElement(ielem)->giveDefaultIntegrationRulePtr()->getIntegrationPoint(ip),
                                                       isttype, tStep);

                    if ( ( vtype == ISVT_TENSOR_S3 ) || ( vtype == ISVT_TENSOR_S3E ) ) {
                        FloatArray help = value;
                        this->makeFullForm(value, help);
                    }

                    for ( j = 1; j <= nc; j++ ) {
                        fprintf( stream, "%e ", value.at(j) );
                    }
                } // end loop over IPs
            } // end loop over elements

            fprintf(stream, "  </DataArray>\n");
        } // end loop over values to be exported
    } // end loop over regions

    fprintf(stream, "</PointData>\n</Piece>\n");
    fprintf(stream, "</UnstructuredGrid>\n");
    fprintf(stream, "</VTKFile>\n");
    fclose(stream);
}
} // end namespace oofem
