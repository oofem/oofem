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
#include "dof.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "cltypes.h"
#include "material.h"
#include "classfactory.h"
#include "crosssection.h"
#include "unknownnumberingscheme.h"

#include "xfem/xfemmanager.h"
#include "xfem/enrichmentitem.h"

#ifdef __PFEM_MODULE
 #include "pfem/pfemparticle.h"
#endif

#include <string>
#include <sstream>
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

IntArray VTKXMLExportModule::redToFull = {
    1, 5, 9, 8, 7, 4, 6, 3, 2
};                                                                      //position of xx, yy, zz, yz, xz, xy in tensor


VTKXMLExportModule::VTKXMLExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport() {}


VTKXMLExportModule::~VTKXMLExportModule() { }


void
VTKXMLExportModule::initializeFrom(InputRecord &ir)
{
    ExportModule::initializeFrom(ir);

    int val;

    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, _IFT_VTKXMLExportModule_cellvars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_VTKXMLExportModule_vars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_VTKXMLExportModule_primvars); // Macro - see unknowntype.h
    IR_GIVE_OPTIONAL_FIELD(ir, externalForcesToExport, _IFT_VTKXMLExportModule_externalForces); // Macro - see unknowntype.h
    IR_GIVE_OPTIONAL_FIELD(ir, ipInternalVarsToExport, _IFT_VTKXMLExportModule_ipvars); // Macro - see internalstatetype.h

    val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_VTKXMLExportModule_stype); // Macro
    stype = ( NodalRecoveryModel::NodalRecoveryModelType ) val;

    this->particleExportFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, particleExportFlag, _IFT_VTKXMLExportModule_particleexportflag); // Macro
}


void
VTKXMLExportModule::initialize()
{
    this->smoother = nullptr;
    this->primVarSmoother = nullptr;
    ExportModule::initialize();
}


void
VTKXMLExportModule::terminate()
{ }


void
VTKXMLExportModule::makeFullTensorForm(FloatArray &answer, const FloatArray &reducedForm, InternalStateValueType vtype)
{
    answer.resize(9);
    answer.zero();

    for ( int i = 1; i <= reducedForm.giveSize(); i++ ) {
        answer.at(redToFull.at(i) ) = reducedForm.at(i);
    }

    if ( vtype == ISVT_TENSOR_S3E ) {
        answer.at(4) *= 0.5;
        answer.at(7) *= 0.5;
        answer.at(8) *= 0.5;
    }

    // Symmetrize if needed
    if ( vtype != ISVT_TENSOR_G ) {
        answer.at(2) = answer.at(4);
        answer.at(3) = answer.at(7);
        answer.at(6) = answer.at(8);
    }
}


std::string
VTKXMLExportModule::giveOutputFileName(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".vtu";
}




std::ofstream
VTKXMLExportModule::giveOutputStream(TimeStep *tStep)
{
    std::string fileName = giveOutputFileName(tStep);
    std::ofstream streamF;

    if ( pythonExport ) {
        streamF = std::ofstream(NULL_DEVICE);//do not write anything
    } else {
        streamF = std::ofstream(fileName);
    }

    if ( !streamF.good() ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    streamF.fill('0');//zero padding
    return streamF;
}

int
VTKXMLExportModule::giveCellType(Element *elem)
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
        OOFEM_ERROR("unsupported element geometry type on element %d", elem->giveNumber() );
    }

    return vtkCellType;
}

int
VTKXMLExportModule::giveNumberOfNodesPerCell(int cellType)
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
        OOFEM_ERROR("unsupported cell type ID");
    }

    return 0; // to make compiler happy
}


void
VTKXMLExportModule::giveElementCell(IntArray &answer, Element *elem)
{
    // Gives the node mapping from the order used in OOFEM to that used in VTK

    Element_Geometry_Type elemGT = elem->giveGeometryType();
    IntArray nodeMapping(0);
    if ( ( elemGT == EGT_point ) ||
         ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
         ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
         ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
         ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
         ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_quad9_2 ) ||
         ( elemGT == EGT_wedge_1 ) ) {} else if ( elemGT == EGT_hexa_27 ) {
        nodeMapping = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18, 23, 25, 26, 24, 22, 21, 27
        };
    } else if ( elemGT == EGT_hexa_2 ) {
        nodeMapping = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18
        };
    } else if ( elemGT == EGT_wedge_2 ) {
        nodeMapping = {
            4, 6, 5, 1, 3, 2, 12, 11, 10, 9, 8, 7, 13, 15, 14
        };
    } else if ( elemGT == EGT_quad_1_interface ) {
        nodeMapping = {
            1, 2, 4, 3
        };
    } else if ( elemGT == EGT_quad_21_interface ) {
        nodeMapping = {
            1, 2, 5, 4, 3, 6
        };
    } else {
        OOFEM_ERROR("VTKXMLExportModule: unsupported element geometry type");
    }

    int nelemNodes = elem->giveNumberOfNodes();
    answer.resize(nelemNodes);
    if ( nodeMapping.giveSize() > 0 ) {
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(nodeMapping.at(i) )->giveNumber();
        }
    } else {
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber();
        }
    }
}


bool
VTKXMLExportModule::isElementComposite(Element *elem)
{
    return ( elem->giveGeometryType() == EGT_Composite );
}



void
VTKXMLExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }


#ifdef _PYBIND_BINDINGS
//clear all dictionaries so they can be filled during export
Py_PrimaryVars.clear();
Py_IntVars.clear();
Py_CellVars.clear();
Py_Nodes.clear();
Py_Elements.clear();
#endif
    
    
#ifdef __VTK_MODULE
    this->fileStream = vtkSmartPointer< vtkUnstructuredGrid >::New();
    this->nodes = vtkSmartPointer< vtkPoints >::New();
    this->elemNodeArray = vtkSmartPointer< vtkIdList >::New();

#else
    this->fileStream = this->giveOutputStream(tStep);
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);

#endif

    // Write output: VTK header
#ifndef __VTK_MODULE
    this->fileStream << "<!-- TimeStep " << tStep->giveTargetTime() * timeScale << " Computed " << current->tm_year + 1900 << "-" << setw(2) << current->tm_mon + 1 << "-" << setw(2) << current->tm_mday << " at " << current->tm_hour << ":" << current->tm_min << ":" << setw(2) << current->tm_sec << " -->\n";
    this->fileStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    this->fileStream << "<UnstructuredGrid>\n";
#endif

    this->giveSmoother(); // make sure smoother is created, Necessary? If it doesn't exist it is created /JB


    if ( !this->particleExportFlag ) {
        /* Loop over pieces  ///@todo: this feature has been broken but not checked if it currently works /JB
         * Start default pieces containing all single cell elements. Elements built up from several vtk
         * cells (composite elements) are exported as individual pieces after the default ones.
         */
        int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
        int anyPieceNonEmpty = 0;

        for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
            // Fills a data struct (VTKPiece) with all the necessary data.
            this->setupVTKPiece(this->defaultVTKPiece, tStep, pieceNum);

            // Write the VTK piece to file.
            anyPieceNonEmpty += this->writeVTKPiece(this->defaultVTKPiece, tStep);
        }

        /*
         * Output all composite elements - one piece per composite element
         * Each element is responsible of setting up a VTKPiece which can then be exported
         */
        Domain *d = emodel->giveDomain(1);

        for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
            const IntArray &elements = this->giveRegionSet(pieceNum)->giveElementList();
            for ( int i = 1; i <= elements.giveSize(); i++ ) {
                Element *el = d->giveElement(elements.at(i) );
                if ( this->isElementComposite(el) ) {
                    if ( el->giveParallelMode() != Element_local ) {
                        continue;
                    }

#ifndef __VTK_MODULE
                    //this->exportCompositeElement(this->defaultVTKPiece, el, tStep);
                    this->exportCompositeElement(this->defaultVTKPieces, el, tStep);

                    for ( int j = 0; j < ( int ) this->defaultVTKPieces.size(); j++ ) {
                        anyPieceNonEmpty += this->writeVTKPiece(this->defaultVTKPieces [ j ],  tStep);
                    }
#else
                    // No support for binary export yet
#endif
                }
            }
        } // end loop over composite elements

#ifndef __VTK_MODULE
        if ( anyPieceNonEmpty == 0 ) {
            // write empty piece, Otherwise ParaView complains if the whole vtu file is without <Piece></Piece>
            this->fileStream << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n";
            this->fileStream << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> </DataArray>\n</Cells>\n";
            this->fileStream << "</Piece>\n";
        }
#endif
    } else {     // if (particleExportFlag)
#ifdef __PFEM_MODULE
        // write out the particles (nodes exported as vertices = VTK_VERTEX)
        Domain *d  = emodel->giveDomain(1);
        int nnode = d->giveNumberOfDofManagers();

        int nActiveNode = 0;
        for ( int inode = 1; inode <= nnode; inode++ ) {
            PFEMParticle *particle = dynamic_cast< PFEMParticle * >( d->giveNode(inode) );
            if ( particle ) {
                if ( particle->isActive() ) {
                    nActiveNode++;
                }
            }
        }

        DofManager *node;
        FloatArray *coords;
        this->fileStream << "<Piece NumberOfPoints=\"" << nActiveNode << "\" NumberOfCells=\"" << nActiveNode << "\">\n";
        this->fileStream << "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ";

        for ( int inode = 1; inode <= nnode; inode++ ) {
            node = d->giveNode(inode);
            PFEMParticle *particle = dynamic_cast< PFEMParticle * >( node );
            if ( particle ) {
                if ( particle->isActive() ) {
                    coords = node->giveCoordinates();
                    ///@todo move this below into setNodeCoords since it should alwas be 3 components anyway
                    for ( int i = 1; i <= coords->giveSize(); i++ ) {
                        this->fileStream << scientific << coords->at(i) << " ";
                    }

                    for ( int i = coords->giveSize() + 1; i <= 3; i++ ) {
                        this->fileStream << scientific << 0.0 << " ";
                    }
                }
            }
        }

        this->fileStream << "</DataArray>\n</Points>\n";


        // output the cells connectivity data
        this->fileStream << "<Cells>\n";
        this->fileStream << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ";

        for ( int ielem = 1; ielem <= nActiveNode; ielem++ ) {
            this->fileStream << ielem - 1 << " ";
        }

        this->fileStream << "</DataArray>\n";

        // output the offsets (index of individual element data in connectivity array)
        this->fileStream << " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ";

        for ( int ielem = 1; ielem <= nActiveNode; ielem++ ) {
            this->fileStream << ielem << " ";
        }
        this->fileStream << "</DataArray>\n";


        // output cell (element) types
        this->fileStream << " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ";
        for ( int ielem = 1; ielem <= nActiveNode; ielem++ ) {
            this->fileStream << 1 << " ";
        }

        this->fileStream << "</DataArray>\n";
        this->fileStream << "</Cells>\n";
        this->fileStream << "</Piece>\n";
#endif //__PFEM_MODULE
    }


    // Finalize the output:
    std::string fname = giveOutputFileName(tStep);
#ifdef __VTK_MODULE

 #if 0
    // Code fragment intended for future support of composite elements in binary format
    // Doesn't as well as I would want it to, interface to VTK is to limited to control this.
    // * The PVTU-file is written by every process (seems to be impossible to avoid).
    // * Part files are renamed and time step and everything else is cut off => name collisions
    vtkSmartPointer< vtkXMLPUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLPUnstructuredGridWriter >::New();
    writer->SetTimeStep(tStep->giveNumber() - 1);
    writer->SetNumberOfPieces(this->emodel->giveNumberOfProcesses() );
    writer->SetStartPiece(this->emodel->giveRank() );
    writer->SetEndPiece(this->emodel->giveRank() );


 #else
    vtkSmartPointer< vtkXMLUnstructuredGridWriter >writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter >::New();
 #endif

    writer->SetFileName(fname.c_str() );
    //writer->SetInput(this->fileStream); // VTK 4
    writer->SetInputData(this->fileStream); // VTK 6

    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    writer->SetDataModeToAscii();
    writer->Write();
#else
    this->fileStream << "</UnstructuredGrid>\n</VTKFile>";
    if(this->fileStream){
        this->fileStream.close();
    }
#endif

    // export raw ip values (if required), works only on one domain
    if ( !this->ipInternalVarsToExport.isEmpty() ) {
        this->exportIntVarsInGpAs(ipInternalVarsToExport, tStep);
        if ( !emodel->isParallel() && tStep->giveNumber() >= 1 ) { // For non-parallel enabled OOFEM, then we only check for multiple steps.
            std::ostringstream pvdEntry;
            std::stringstream subStep;
            if ( tstep_substeps_out_flag ) {
                subStep << "." << tStep->giveSubStepNumber();
            }
            pvdEntry << "<DataSet timestep=\"" << tStep->giveTargetTime() * this->timeScale << subStep.str() << "\" group=\"\" part=\"\" file=\"" << this->giveOutputBaseFileName(tStep) + ".gp.vtu" << "\"/>";
            this->gpPvdBuffer.push_back(pvdEntry.str() );
            this->writeGPVTKCollection();
        }
    }

    // Write the *.pvd-file. Currently only contains time step information. It's named "timestep" but is actually the total time.
    // First we check to see that there are more than 1 time steps, otherwise it is redundant;
    if ( emodel->isParallel() && emodel->giveRank() == 0 ) {
        ///@todo Should use probably use PVTU-files instead. It is starting to get messy.
        // For this to work, all processes must have an identical output file name.
        for ( int i = 0; i < this->emodel->giveNumberOfProcesses(); ++i ) {
            std::ostringstream pvdEntry;
            std::stringstream subStep;
            char fext [ 100 ];
            if ( this->emodel->giveNumberOfProcesses() > 1 ) {
                sprintf(fext, "_%03d.m%d.%d", i, this->number, tStep->giveNumber() );
            } else {
                sprintf(fext, "m%d.%d", this->number, tStep->giveNumber() );
            }
            if ( tstep_substeps_out_flag ) {
                subStep << "." << tStep->giveSubStepNumber();
            }
            pvdEntry << "<DataSet timestep=\"" << tStep->giveTargetTime() * this->timeScale << subStep.str() << "\" group=\"\" part=\"" << i << "\" file=\"" << this->emodel->giveOutputBaseFileName() << fext << ".vtu\"/>";
            this->pvdBuffer.push_back(pvdEntry.str() );
        }

        this->writeVTKCollection();
    } else if ( !emodel->isParallel() && tStep->giveNumber() >= 1 ) { // For non-parallel, then we only check for multiple steps.
        std::ostringstream pvdEntry;
        std::stringstream subStep;
        if ( tstep_substeps_out_flag ) {
            subStep << "." << tStep->giveSubStepNumber();
        }
        pvdEntry << "<DataSet timestep=\"" << tStep->giveTargetTime() * this->timeScale << subStep.str() << "\" group=\"\" part=\"\" file=\"" << fname << "\"/>";
        this->pvdBuffer.push_back(pvdEntry.str() );
        this->writeVTKCollection();
    }
}


void
VTKPiece::setNumberOfNodes(int numNodes)
{
    this->numNodes = numNodes;
    this->nodeCoords.resize(numNodes);
}

void
VTKPiece::setNumberOfCells(int numCells)
{
    this->numCells = numCells;
    this->connectivity.resize(numCells);
    this->elCellTypes.resize(numCells);
    this->elOffsets.resize(numCells);
}

void
VTKPiece::setConnectivity(int cellNum, IntArray &nodes)
{
    this->connectivity [ cellNum - 1 ] = nodes;
}

void
VTKPiece::setNodeCoords(int nodeNum, const FloatArray &coords)
{
    this->nodeCoords [ nodeNum - 1 ] = coords;
}

void
VTKPiece::setNumberOfPrimaryVarsToExport(int numVars, int numNodes)
{
    this->nodeVars.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeVars [ i - 1 ].resize(numNodes);
    }
}

void
VTKPiece::setNumberOfLoadsToExport(int numVars, int numNodes)
{
    this->nodeLoads.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeLoads [ i - 1 ].resize(numNodes);
    }
}

void
VTKPiece::setNumberOfInternalVarsToExport(int numVars, int numNodes)
{
    this->nodeVarsFromIS.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeVarsFromIS [ i - 1 ].resize(numNodes);
    }
}

void
VTKPiece::setNumberOfInternalXFEMVarsToExport(int numVars, int numEnrichmentItems, int numNodes)
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
VTKPiece::setNumberOfCellVarsToExport(int numVars, int numCells)
{
    this->elVars.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->elVars [ i - 1 ].resize(numCells);
    }
}

void
VTKPiece::setPrimaryVarInNode(int varNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVars [ varNum - 1 ] [ nodeNum - 1 ] = std::move(valueArray);
}

void
VTKPiece::setLoadInNode(int varNum, int nodeNum, FloatArray valueArray)
{
    this->nodeLoads [ varNum - 1 ] [ nodeNum - 1 ] = std::move(valueArray);
}

void
VTKPiece::setInternalVarInNode(int varNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVarsFromIS [ varNum - 1 ] [ nodeNum - 1 ] = std::move(valueArray);
}

void
VTKPiece::setInternalXFEMVarInNode(int varNum, int eiNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVarsFromXFEMIS [ varNum - 1 ] [ eiNum - 1 ] [ nodeNum - 1 ] = std::move(valueArray);
}


void
VTKPiece::setCellVar(int varNum, int cellNum, FloatArray valueArray)
{
    this->elVars [ varNum - 1 ] [ cellNum - 1 ] = std::move(valueArray);
}


void
VTKXMLExportModule::setupVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep, int region)
{
    // Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.

    Domain *d  = emodel->giveDomain(1);
    Element *elem;

    this->giveSmoother(); // make sure smoother is created

    // output nodes Region By Region
    int numNodes, numRegionEl;
    IntArray mapG2L, mapL2G;

    // Assemble local->global and global->local region map and get number of
    // single cells to process, the composite cells exported individually.
    this->initRegionNodeNumbering(mapG2L, mapL2G, numNodes, numRegionEl, d, tStep, region);
    if ( numNodes > 0 && numRegionEl > 0 ) {
        // Export nodes as vtk vertices
        vtkPiece.setNumberOfNodes(numNodes);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            const auto &coords = d->giveNode(mapL2G.at(inode) )->giveCoordinates();
            vtkPiece.setNodeCoords(inode, coords);
        }


        //-------------------------------------------
        // Export all the cell data for the piece
        //-------------------------------------------
        IntArray cellNodes;
        vtkPiece.setNumberOfCells(numRegionEl);
        IntArray regionElInd;

        int offset = 0;
        int cellNum = 0;
        IntArray elems = this->giveRegionSet(region)->giveElementList();
        for ( int ei = 1; ei <= elems.giveSize(); ei++ ) {
            int elNum = elems.at(ei);
            elem = d->giveElement(elNum);

            // Skip elements that:
            // are inactivated or of composite type ( these are exported individually later)
            if ( this->isElementComposite(elem) || !elem->isActivated(tStep) ) {
                continue;
            }

            //skip materials with casting time > current time
            if ( !elem->isCast(tStep) ) {
                continue;
            }

            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            regionElInd.followedBy(elNum);

            cellNum++;

            // Set the connectivity
            this->giveElementCell(cellNodes, elem);  // node numbering of the cell according to the VTK format

            // Map from global to local node numbers for the current piece
            int numElNodes = cellNodes.giveSize();
            IntArray connectivity(numElNodes);
            for ( int i = 1; i <= numElNodes; i++ ) {
                connectivity.at(i) = mapG2L.at(cellNodes.at(i) );
            }

            vtkPiece.setConnectivity(cellNum, connectivity);

            vtkPiece.setCellType(cellNum, this->giveCellType(elem) );  // VTK cell type

            offset += numElNodes;
            vtkPiece.setOffset(cellNum, offset);
        }


        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(vtkPiece, mapG2L, mapL2G, region, tStep);
        this->exportIntVars(vtkPiece, mapG2L, mapL2G, region, tStep);
        this->exportExternalForces(vtkPiece, mapG2L, mapL2G, region, tStep);

        this->exportCellVars(vtkPiece, regionElInd, tStep);
        
        
#ifdef _PYBIND_BINDINGS    
        if ( pythonExport ) {
        //Export nodes
            py::list vals;
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                py::list node;
                const int numberG = d->giveNode(mapL2G.at(inode))->giveGlobalNumber();
                const int number = d->giveNode(mapL2G.at(inode))->giveNumber();
                const auto &coords = d->giveNode(mapL2G.at(inode))->giveCoordinates();
                node.append(inode);
                node.append(number);
                node.append(numberG);
                node.append(coords);
                vals.append(node);
            }
            std::string s = std::to_string(region);
            this->Py_Nodes[s.c_str()] = vals;//keys as region number (VTKPiece)
       
       //Export elements
            py::list elemVals;
            for ( int ei = 1; ei <= vtkPiece.giveNumberOfCells(); ei++ ) {
                py::list element;
                IntArray &conn = vtkPiece.giveCellConnectivity(ei);
                element.append(ei);
                element.append(conn);
                elemVals.append(element);
            }
            this->Py_Elements[s.c_str()] = elemVals;//keys as region number (VTKPiece)
        }
#endif    

        
    } // end of default piece for simple geometry elements
}


bool
VTKXMLExportModule::writeVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep)
{
    // Write a VTK piece to file. This could be the whole domain (most common case) or it can be a
    // (so-called) composite element consisting of several VTK cells (layered structures, XFEM, etc.).

    /*
     * if ( !vtkPiece.giveNumberOfCells() ) { // handle piece with no elements. Otherwise ParaView complains if the whole vtu file is without <Piece></Piece>
     * //          fprintf(this->fileStream, "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n");
     * //          fprintf(this->fileStream, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> </DataArray>\n</Cells>\n");
     * //          fprintf(this->fileStream, "</Piece>\n");
     *    return;
     * }
     */
    if ( !vtkPiece.giveNumberOfCells() ) {
        return false;
    }


    // Write output: node coords
    int numNodes = vtkPiece.giveNumberOfNodes();
    int numEl = vtkPiece.giveNumberOfCells();
    FloatArray coords;

#ifdef __VTK_MODULE
    FloatArray vtkCoords(3);
    for ( int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPiece.giveNodeCoords(inode);
        vtkCoords.zero();
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            vtkCoords.at(i) = coords.at(i);
        }

        this->nodes->InsertNextPoint(vtkCoords.at(1), vtkCoords.at(2), vtkCoords.at(3) );
        this->fileStream->SetPoints(nodes);
    }

#else
    this->fileStream << "<Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << numEl << "\">\n";
    this->fileStream << "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ";

    for ( int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPiece.giveNodeCoords(inode);
        ///@todo move this below into setNodeCoords since it should alwas be 3 components anyway
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            this->fileStream << scientific << coords.at(i) << " ";
        }

        for ( int i = coords.giveSize() + 1; i <= 3; i++ ) {
            this->fileStream << scientific << 0.0 << " ";
        }
    }

    this->fileStream << "</DataArray>\n</Points>\n";
#endif


    // Write output: connectivity, offsets, cell types

    // output the connectivity data
#ifdef __VTK_MODULE
    this->fileStream->Allocate(numEl);
#else
    this->fileStream << "<Cells>\n";
    this->fileStream << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ";
#endif
    IntArray cellNodes;
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        cellNodes = vtkPiece.giveCellConnectivity(ielem);

#ifdef __VTK_MODULE
        elemNodeArray->Reset();
        elemNodeArray->SetNumberOfIds(cellNodes.giveSize() );
#endif
        for ( int i = 1; i <= cellNodes.giveSize(); i++ ) {
#ifdef __VTK_MODULE
            elemNodeArray->SetId(i - 1, cellNodes.at(i) - 1);
#else
            this->fileStream << cellNodes.at(i) - 1 << " ";
#endif
        }

#ifdef __VTK_MODULE
        this->fileStream->InsertNextCell(vtkPiece.giveCellType(ielem), elemNodeArray);
#else
        this->fileStream << " ";
#endif
    }

#ifndef __VTK_MODULE
    this->fileStream << "</DataArray>\n";

    // output the offsets (index of individual element data in connectivity array)
    this->fileStream << " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ";

    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        this->fileStream << vtkPiece.giveCellOffset(ielem) << " ";
    }

    this->fileStream << "</DataArray>\n";


    // output cell (element) types
    this->fileStream << " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ";
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        this->fileStream << vtkPiece.giveCellType(ielem) << " ";
    }

    this->fileStream << "</DataArray>\n";
    this->fileStream << "</Cells>\n";


    ///@todo giveDataHeaders is currently not updated wrt the new structure -> no file names in headers /JB
    std::string pointHeader, cellHeader;
    this->giveDataHeaders(pointHeader, cellHeader);

    this->fileStream << pointHeader.c_str();
#endif

    this->writePrimaryVars(vtkPiece);       // Primary field
    this->writeIntVars(vtkPiece);           // Internal State Type variables smoothed to the nodes
    this->writeExternalForces(vtkPiece);           // External forces

    if ( emodel->giveDomain(1)->hasXfemManager() ) {
        this->writeXFEMVars(vtkPiece);      // XFEM State Type variables associated with XFEM structure
    }

#ifndef __VTK_MODULE
    this->fileStream << "</PointData>\n";
    this->fileStream << cellHeader.c_str();
#endif

    this->writeCellVars(vtkPiece);          // Single cell variables ( if given in the integration points then an average will be exported)

#ifndef __VTK_MODULE
    this->fileStream << "</CellData>\n";
    this->fileStream << "</Piece>\n";
#endif

    //}
    
    // Clear object so it can be filled with new data for the next piece
    vtkPiece.clear();
    return true;
}



#ifndef __VTK_MODULE
void
VTKXMLExportModule::giveDataHeaders(std::string &pointHeader, std::string &cellHeader)
{
    std::string scalars, vectors, tensors;

    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        if ( type == DisplacementVector || type == EigenVector || type == VelocityVector || type == DirectorField ) {
            vectors += __UnknownTypeToString(type);
            vectors.append(" ");
        } else if ( type == FluxVector || type == PressureVector || type == Temperature || type == Humidity || type == DeplanationFunction ) {
            scalars += __UnknownTypeToString(type);
            scalars.append(" ");
        } else {
            OOFEM_ERROR("unsupported UnknownType %s", __UnknownTypeToString(type) );
        }
    }

    for ( int i = 1; i <= internalVarsToExport.giveSize(); i++ ) {
        InternalStateType isttype = ( InternalStateType ) internalVarsToExport.at(i);
        InternalStateValueType vtype = giveInternalStateValueType(isttype);

        if ( vtype == ISVT_SCALAR ) {
            scalars += __InternalStateTypeToString(isttype);
            scalars.append(" ");
        } else if ( vtype == ISVT_VECTOR ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
            tensors += __InternalStateTypeToString(isttype);
            tensors.append(" ");
        } else {
            OOFEM_ERROR("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
        }
    }

    for ( int i = 1; i <= externalForcesToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) externalForcesToExport.at(i);
        if ( type == DisplacementVector || type == VelocityVector || type == DirectorField ) {
            vectors += std::string("Load") + __UnknownTypeToString(type);
            vectors.append(" ");
        } else if ( type == FluxVector || type == PressureVector || type == Temperature || type == Humidity ) {
            scalars += std::string("Load") + __UnknownTypeToString(type);
            scalars.append(" ");
        } else {
            OOFEM_ERROR("unsupported UnknownType %s", __UnknownTypeToString(type) );
        }
    }

    // print header
    pointHeader = "<PointData Scalars=\"" + scalars + "\" "
                  +  "Vectors=\"" + vectors + "\" "
                  +  "Tensors=\"" + tensors + "\" >\n";


    scalars.clear();
    vectors.clear();
    tensors.clear();
    // prepare header
    for ( int i = 1; i <= this->cellVarsToExport.giveSize(); i++ ) {
        InternalStateType isttype = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType vtype = giveInternalStateValueType(isttype);

        if ( vtype == ISVT_SCALAR ) {
            scalars += __InternalStateTypeToString(isttype);
            scalars.append(" ");
        } else if ( vtype == ISVT_VECTOR ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
            tensors += __InternalStateTypeToString(isttype);
            tensors.append(" ");
        } else {
            OOFEM_WARNING("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
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
VTKXMLExportModule::exportIntVars(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int region, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    InternalStateType isType;
    FloatArray answer;

    this->giveSmoother()->clear(); // Makes sure smoother is up-to-date with potentially new mesh.

    // Export of Internal State Type fields
    vtkPiece.setNumberOfInternalVarsToExport(internalVarsToExport.giveSize(), mapL2G.giveSize() );
    for ( int field = 1; field <= internalVarsToExport.giveSize(); field++ ) {
        isType = ( InternalStateType ) internalVarsToExport.at(field);

        for ( int nodeNum = 1; nodeNum <= mapL2G.giveSize(); nodeNum++ ) {
            Node *node = d->giveNode(mapL2G.at(nodeNum) );
            this->getNodalVariableFromIS(answer, node, tStep, isType, region);
            vtkPiece.setInternalVarInNode(field, nodeNum, answer);
        }
    }

    // Export of XFEM related quantities
    if ( d->hasXfemManager() ) {
        XfemManager *xFemMan = d->giveXfemManager();
        int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();

        vtkPiece.setNumberOfInternalXFEMVarsToExport(xFemMan->vtkExportFields.giveSize(), nEnrIt, mapL2G.giveSize() );
        for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
            XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];

            for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
                for ( int nodeIndx = 1; nodeIndx <= mapL2G.giveSize(); nodeIndx++ ) {
                    Node *node = d->giveNode(mapL2G.at(nodeIndx) );
                    getNodalVariableFromXFEMST(answer, node, tStep, xfemstype, region, xFemMan->giveEnrichmentItem(enrItIndex) );
                    vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeIndx, answer);
                }
            }
        }
    }
}


void
VTKXMLExportModule::getNodalVariableFromIS(FloatArray &answer, Node *node, TimeStep *tStep, InternalStateType type, int ireg)
{
    // Recovers nodal values from Internal States defined in the integration points.
    // Should return an array with proper size supported by VTK (1, 3 or 9)
    // Domain *d = emodel->giveDomain(1);
    this->giveSmoother();
    IntArray redIndx;

    if ( !( type == IST_DisplacementVector || type == IST_MaterialInterfaceVal  ) ) {
        this->smoother->recoverValues(* this->giveRegionSet(ireg), type, tStep);
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
            valueArray.at(1) = mi->giveNodalScalarRepresentation(node->giveNumber() );
        }
    } else {
        int found = this->smoother->giveNodalVector(val, node->giveNumber() );
        if ( !found ) {
            valueArray.resize(redIndx.giveSize() );
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
        answer = * val;
        // bp: hack for BeamForceMomentTensor, which should be splitted into force and momentum vectors
        if ( type == IST_BeamForceMomentTensor ) {
            answer.resizeWithValues(6);
        } else {
            answer.resizeWithValues(3);
        }
    } else if ( valType == ISVT_TENSOR_S3 || valType == ISVT_TENSOR_S3E || valType == ISVT_TENSOR_G ) {
        this->makeFullTensorForm(answer, * val, valType);
    } else {
        OOFEM_ERROR("ISVT_UNDEFINED encountered")
    }
}


void
VTKXMLExportModule::getNodalVariableFromXFEMST(FloatArray &answer, Node *node, TimeStep *tStep, XFEMStateType xfemstype, int ireg, EnrichmentItem *ei)
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
        ei->evalLevelSetNormalInNode(valueArray.at(1), node->giveNumber(), node->giveCoordinates() );
    } else if ( xfemstype == XFEMST_LevelSetGamma ) {
        valueArray.resize(1);
        val = & valueArray;
        ei->evalLevelSetTangInNode(valueArray.at(1), node->giveNumber(), node->giveCoordinates() );
    } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
        valueArray.resize(1);
        val = & valueArray;
        ei->evalNodeEnrMarkerInNode(valueArray.at(1), node->giveNumber() );
    } else {
        //OOFEM_WARNING("invalid data in node %d", inode);
    }

    ///@todo duplicated code from getNodalVariableFromIS - uneccessary
    int ncomponents = giveInternalStateTypeSize(valType);
    answer.resize(ncomponents);
    int valSize = val->giveSize(); // size of recovered quantity

    // check if valSize corresponds to the expected size otherwise pad with zeros
    if ( valType == ISVT_SCALAR ) {
        answer.at(1) = valSize ? val->at(1) : 0.0;
    } else if ( valType == ISVT_VECTOR ) {
        answer = * val;
        answer.resizeWithValues(3);
    } else if ( valType == ISVT_TENSOR_S3 || valType == ISVT_TENSOR_S3E || valType == ISVT_TENSOR_G ) {
        this->makeFullTensorForm(answer, * val, valType);
    } else {
        OOFEM_ERROR("ISVT_UNDEFINED encountered")
    }
}


void
VTKXMLExportModule::writeIntVars(VTKPiece &vtkPiece)
{
    int n = internalVarsToExport.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(i);
        int ncomponents;

        const char *name = __InternalStateTypeToString(type);
        ( void ) name;//silence warning
        int numNodes = vtkPiece.giveNumberOfNodes();
        FloatArray valueArray;
        valueArray = vtkPiece.giveInternalVarInNode(i, 1);
        ncomponents = valueArray.giveSize();
        ( void ) ncomponents;//silence warning

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
        varArray->SetName(name);
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.giveInternalVarInNode(i, inode);
            for ( int i = 1; i <= ncomponents; ++i ) {
                varArray->SetComponent(inode - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else

        this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.giveInternalVarInNode(i, inode);
            this->writeVTKPointData(valueArray);
        }

#endif
        // Footer
#ifndef __VTK_MODULE
        this->fileStream << "</DataArray>\n";
#endif
    
#ifdef _PYBIND_BINDINGS
        if ( pythonExport ) {
            py::list vals;
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                valueArray = vtkPiece.giveInternalVarInNode(i, inode);
                vals.append(valueArray);
            }
            this->Py_IntVars[name] = vals;
        }
#endif
    } //end of for
}


void
VTKXMLExportModule::writeXFEMVars(VTKPiece &vtkPiece)
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
        ( void ) ncomponents; //silence the warning

        int numNodes = vtkPiece.giveNumberOfNodes();
        for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
            // Header
            char name [ 100 ];   // Must I define a fixed size? /JB
            sprintf(name, "%s_%d ", namePart, xFemMan->giveEnrichmentItem(enrItIndex)->giveNumber() );

#ifdef __VTK_MODULE
            vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
            varArray->SetName(name);
            varArray->SetNumberOfComponents(ncomponents);
            varArray->SetNumberOfTuples(numNodes);
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                valueArray = vtkPiece.giveInternalXFEMVarInNode(field, enrItIndex, inode);
                for ( int i = 1; i <= ncomponents; ++i ) {
                    varArray->SetComponent(inode - 1, i - 1, valueArray.at(i) );
                }
            }

            this->writeVTKPointData(name, varArray);
#else
            this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                valueArray = vtkPiece.giveInternalXFEMVarInNode(field, enrItIndex, inode);
                this->writeVTKPointData(valueArray);
            }
            this->fileStream << "</DataArray>\n";
#endif
        }
    }
}




//----------------------------------------------------
// Misc. functions
//----------------------------------------------------
#ifdef __VTK_MODULE
void
VTKXMLExportModule::writeVTKPointData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray)
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
VTKXMLExportModule::writeVTKPointData(FloatArray &valueArray)
{
    // Write the data to file
    for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
        this->fileStream << scientific << valueArray.at(i) << " ";
    }
}
#endif


#ifdef __VTK_MODULE
void
VTKXMLExportModule::writeVTKCellData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray)
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
VTKXMLExportModule::writeVTKCellData(FloatArray &valueArray)
{
    // Write the data to file ///@todo exact copy of writeVTKPointData so remove
    for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
        this->fileStream << valueArray.at(i) << " ";
    }
}
#endif




int
VTKXMLExportModule::initRegionNodeNumbering(IntArray &regionG2LNodalNumbers,
                                            IntArray &regionL2GNodalNumbers,
                                            int &regionDofMans,
                                            int &regionSingleCells,
                                            Domain *domain, TimeStep *tStep, int reg)
{
    // regionG2LNodalNumbers is array with mapping from global numbering to local region numbering.
    // The i-th value contains the corresponding local region number (or zero, if global number is not in region).

    // regionL2GNodalNumbers is array with mapping from local to global numbering.
    // The i-th value contains the corresponding global node number.


    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    regionG2LNodalNumbers.resize(nnodes);
    regionG2LNodalNumbers.zero();
    regionDofMans = 0;
    regionSingleCells = 0;

    IntArray elements = this->giveRegionSet(reg)->giveElementList();
    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        int ielem = elements.at(ie);
        element = domain->giveElement(ielem);

        if ( this->isElementComposite(element) ) {
            continue;                                    // composite cells exported individually
        }

        if ( !element->isActivated(tStep) ) {                    //skip inactivated elements
            continue;
        }

        //skip materials with casting time > current time
        if ( !element->isCast(tStep) ) {
            continue;
        }

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

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
            regionL2GNodalNumbers.at(regionG2LNodalNumbers.at(i) ) = i;
        }
    }

    return 1;
}








//----------------------------------------------------
// Primary variables - readily available in the nodes
//----------------------------------------------------
void
VTKXMLExportModule::exportPrimaryVars(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int region, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;
    this->givePrimVarSmoother()->clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.

    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), mapL2G.giveSize() );
    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

        for ( int inode = 1; inode <= mapL2G.giveSize(); inode++ ) {
            DofManager *dman = d->giveNode(mapL2G.at(inode) );

            this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region);
            vtkPiece.setPrimaryVarInNode(i, inode, std::move(valueArray) );
        }
    }
}


void
VTKXMLExportModule::getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int ireg)
{
    // This code is not perfect. It should be rewritten to handle all cases more gracefully.
    ///@todo This method needs to be cleaned up - maybe define the common vector types so
    /// certain dofid's are associated with them /JB

    IntArray dofIDMask(3);
    int size;
    const FloatArray *recoveredVal;

    InternalStateType iState = IST_DisplacementVector; // Shouldn't be necessary

    dofIDMask.clear();

    if ( type == DisplacementVector ) {
        dofIDMask = {
            ( int ) Undef, ( int ) Undef, ( int ) Undef
        };
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( id == D_u ) {
                dofIDMask.at(1) = id;
            } else if ( id == D_v ) {
                dofIDMask.at(2) = id;
            } else if ( id == D_w ) {
                dofIDMask.at(3) = id;
            }
        }

        answer.resize(3);
    } else if ( type == VelocityVector ) {
        dofIDMask = {
            ( int ) Undef, ( int ) Undef, ( int ) Undef
        };
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( id == V_u ) {
                dofIDMask.at(1) = id;
            } else if ( id == V_v ) {
                dofIDMask.at(2) = id;
            } else if ( id == V_w ) {
                dofIDMask.at(3) = id;
            }
        }

        answer.resize(3);
    } else if ( type == EigenVector ) {
        dofIDMask = {
            ( int ) Undef, ( int ) Undef, ( int ) Undef
        };
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( ( id == V_u ) || ( id == D_u ) ) {
                dofIDMask.at(1) = id;
            } else if ( ( id == V_v ) || ( id == D_v ) ) {
                dofIDMask.at(2) = id;
            } else if ( ( id == V_w ) || ( id == D_w ) ) {
                dofIDMask.at(3) = id;
            }
        }

        answer.resize(3);
    } else if ( type == FluxVector || type == Humidity ) {
        dofIDMask.followedBy(C_1);
        iState = IST_MassConcentration_1;
        answer.resize(1);
    } else if ( type == DeplanationFunction ) {
        dofIDMask.followedBy(Warp_PsiTheta);
        iState = IST_Temperature;
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
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( ( id == W_u ) || ( id == W_v ) || ( id == W_w ) ) {
                dofIDMask.followedBy(id);
            }

            answer.resize(3);
        }

        iState = IST_DirectorField;
    } else {
        OOFEM_ERROR("unsupported unknownType %s", __UnknownTypeToString(type) );
    }

    size = dofIDMask.giveSize();
    answer.zero();

    for ( int j = 1; j <= size; j++ ) {
        DofIDItem id = ( DofIDItem ) dofIDMask.at(j);
        if ( id == Undef ) {
            answer.at(j) = 0.;
        } else if ( iState == IST_DirectorField ) {
            answer.at(j) = dman->giveDofWithID(id)->giveUnknown(VM_Total, tStep);
            // recover values if not done before
            this->givePrimVarSmoother()->recoverValues(* this->giveRegionSet(ireg), iState, tStep);
            this->givePrimVarSmoother()->giveNodalVector(recoveredVal, dman->giveNumber() );
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING("Recovered variable size mismatch for %d for id %d", type, id);
                answer.at(j) = 0.0;
            }
        } else if ( dman->hasDofID(id) ) {
            // primary variable available directly in DOF-manager
            answer.at(j) = dman->giveDofWithID(id)->giveUnknown(VM_Total, tStep);
            //mj - if incremental value needed: answer.at(j) = dman->giveDofWithID(id)->giveUnknown(VM_Incremental, tStep);
        } else if ( iState != IST_Undefined ) {
            // primary variable not directly available
            // but equivalent InternalStateType provided
            // in this case use smoother to recover nodal value

            // This can't deal with ValueModeType, and would recover over and over for some vectorial quantities like velocity
            // recover values if not done before.
            this->givePrimVarSmoother()->recoverValues(* this->giveRegionSet(ireg), iState, tStep);
            this->givePrimVarSmoother()->giveNodalVector(recoveredVal, dman->giveNumber() );
            // here we have a lack of information about how to convert recovered values to response
            // if the size is compatible we accept it, otherwise give a warning and zero value.
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING("Recovered variable size mismatch for \"%s\" for dof id %d. Size is %d, should be %d", __UnknownTypeToString(type), id, recoveredVal->giveSize(), size);
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
VTKXMLExportModule::writePrimaryVars(VTKPiece &vtkPiece)
{
    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        ( void ) ncomponents; //silence the warning
        int numNodes = vtkPiece.giveNumberOfNodes();
        const char *name = __UnknownTypeToString(type);
        ( void ) name; //silence the warning

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
        varArray->SetName(name);
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(i, inode);
            for ( int j = 1; j <= ncomponents; ++j ) {
                varArray->SetComponent(inode - 1, j - 1, valueArray.at(j) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else
        this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(i, inode);
            this->writeVTKPointData(valueArray);
        }
        this->fileStream << "</DataArray>\n";

 #ifdef _PYBIND_BINDINGS
        if ( pythonExport ) {
            py::list vals;
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(i, inode);
                vals.append(valueArray);
            }
            this->Py_PrimaryVars[name] = vals;
        }
 #endif
#endif
    }
}


//----------------------------------------------------
// Load vectors
//----------------------------------------------------
void
VTKXMLExportModule::exportExternalForces(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int region, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    this->givePrimVarSmoother()->clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.

    if ( externalForcesToExport.giveSize() == 0 ) {
        return;
    }

    ///@todo Add a more flexible solution here, ask the Engineering model for the equivalent to this (perhaps as part of the primary field?)
    /// This should be looked into, just as "getNodalVariableFromPrimaryField" is particularly complicated.
    int neq = emodel->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering() );
    int npeq = emodel->giveNumberOfDomainEquations(1, EModelDefaultPrescribedEquationNumbering() );
    FloatArray extForces(neq), extForcesP(npeq);
    emodel->assembleVector(extForces, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), d);
    emodel->assembleVector(extForcesP, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultPrescribedEquationNumbering(), d);

    vtkPiece.setNumberOfLoadsToExport(externalForcesToExport.giveSize(), mapL2G.giveSize() );
    for ( int i = 1; i <= externalForcesToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) externalForcesToExport.at(i);
        ///@todo Have some mapping for UnknownType -> DofID array
        IntArray dofids;
        if ( type == VelocityVector ) {
            dofids = {
                V_u, V_v, V_w
            };
        } else if ( type == DisplacementVector ) {
            dofids = {
                D_u, D_v, D_w
            };
        } else if ( type == PressureVector ) {
            dofids = {
                P_f
            };
        } else {
            OOFEM_WARNING("Unrecognized UnknownType (%d), no external forces exported", type);
        }

        for ( int inode = 1; inode <= mapL2G.giveSize(); inode++ ) {
            DofManager *dman = d->giveNode(mapL2G.at(inode) );

            FloatArray valueArray(dofids.giveSize() );
            for ( int k = 1; k <= dofids.giveSize(); ++k ) {
                Dof *dof = dman->giveDofWithID(dofids.at(k) );
                ///@todo Have to make more assumptions here.. we shouldn't assume EModelDefaultEquationNumbering. Do something nicer than extForces and extForcesP instead.
                int eq;
                if ( ( eq = dof->giveEquationNumber(EModelDefaultEquationNumbering() ) ) > 0 ) {
                    valueArray.at(k) = extForces.at(eq);
                } else if ( ( eq = dof->giveEquationNumber(EModelDefaultPrescribedEquationNumbering() ) ) > 0 ) {
                    valueArray.at(k) = extForcesP.at(eq);
                }
            }
            //this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region);
            vtkPiece.setLoadInNode(i, inode, std::move(valueArray) );
        }
    }
}


void
VTKXMLExportModule::writeExternalForces(VTKPiece &vtkPiece)
{
    for ( int i = 1; i <= externalForcesToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) externalForcesToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        ( void ) ncomponents; //silence the warning
        int numNodes = vtkPiece.giveNumberOfNodes();
        std::string name = std::string("Load") + __UnknownTypeToString(type);

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
        varArray->SetName(name.c_str() );
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.giveLoadInNode(i, inode);
            for ( int j = 1; j <= ncomponents; ++j ) {
                varArray->SetComponent(inode - 1, j - 1, valueArray.at(j) );
            }
        }

        this->writeVTKPointData(name.c_str(), varArray);

#else
        this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.giveLoadInNode(i, inode);
            this->writeVTKPointData(valueArray);
        }
        this->fileStream << "</DataArray>\n";
#endif
    }
}




//----------------------------------------------------
// Cell vars
//----------------------------------------------------

void
VTKXMLExportModule::exportCellVars(VTKPiece &vtkPiece, const IntArray &elems, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;

    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), elems.giveSize() );
    for ( int field = 1; field <= cellVarsToExport.giveSize(); field++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(field);

        for ( int subIndex = 1; subIndex <= elems.giveSize(); ++subIndex ) {
            Element *el = d->giveElement(elems.at(subIndex) );  ///@todo should be a pointer to an element in the region /JB
            if ( el->giveParallelMode() != Element_local ) {
                continue;
            }

            this->getCellVariableFromIS(valueArray, el, type, tStep);
            vtkPiece.setCellVar(field, subIndex, valueArray);
        }
    }
}


void
VTKXMLExportModule::getCellVariableFromIS(FloatArray &answer, Element *el, InternalStateType type, TimeStep *tStep)
{
    InternalStateValueType valType = giveInternalStateValueType(type);
    int ncomponents = giveInternalStateTypeSize(valType);

    answer.resize(ncomponents);

    switch ( type ) {
    // Special scalars
    case IST_MaterialNumber:
        // commented by bp: do what user wants
        //OOFEM_WARNING1("Material numbers are deprecated, outputing cross section number instead...");
        answer.at(1) = ( double ) el->giveMaterial()->giveNumber();
        break;
    case IST_CrossSectionNumber:
        answer.at(1) = ( double ) el->giveCrossSection()->giveNumber();
        break;
    case IST_ElementNumber:
        answer.at(1) = ( double ) el->giveNumber();
        break;
    case IST_Pressure: ///@todo This case seems redundant, remove? /JB, /// Why this special treatment for pressure? / Mikael
        if ( el->giveNumberOfInternalDofManagers() == 1 ) {
            //IntArray pmask(1); pmask.at(1) = P_f;
            //el->giveInternalDofManager(1)->giveUnknownVector (answer, pmask, VM_Total, tStep);
            //answer.at(1) = answer.at(1);
        }

        break;
    case IST_AbaqusStateVector:
    {
        // compute cell average from ip values
        IntegrationRule *iRule = el->giveDefaultIntegrationRulePtr();
        computeIPAverage(answer, iRule, el, type, tStep); // if element has more than one iRule?? /JB
    }
    break;

    // Special vectors
    case IST_MaterialOrientation_x:
    case IST_MaterialOrientation_y:
    case IST_MaterialOrientation_z:
    {
        FloatMatrix rotMat;
        int col = 0;
        if ( type == IST_MaterialOrientation_x ) {
            col = 1;
        } else if ( type == IST_MaterialOrientation_y ) {
            col = 2;
        } else if ( type == IST_MaterialOrientation_z ) {
            col = 3;
        }

        if ( !el->giveLocalCoordinateSystem(rotMat) ) {
            rotMat.resize(3, 3);
            rotMat.beUnitMatrix();
        }

        answer.beColumnOf(rotMat, col);
        break;
    }

    // Export cell data as average from ip's as default
    default:

        // compute cell average from ip values
        IntegrationRule *iRule = el->giveDefaultIntegrationRulePtr();
        computeIPAverage(answer, iRule, el, type, tStep); // if element has more than one iRule?? /JB
        // Reshape the Voigt vectors to include all components (duplicated if necessary, VTK insists on 9 components for tensors.)
        /// @todo Is this part necessary now when giveIPValue returns full form? Only need to symmetrize in case of 6 components /JB
        /// @todo Some material models aren't exporting values correctly (yet) / Mikael
        if ( valType == ISVT_TENSOR_S3 || valType == ISVT_TENSOR_S3E || valType == ISVT_TENSOR_G ) {
            FloatArray temp = answer;
            this->makeFullTensorForm(answer, temp, valType);
        } else if ( valType == ISVT_VECTOR && answer.giveSize() < 3 ) {
            answer.resizeWithValues(3);
        } else if ( ncomponents != answer.giveSize() ) { // Trying to gracefully handle bad cases, just output zeros.
            answer.resizeWithValues(ncomponents);
        }
    }
}


void
VTKXMLExportModule::writeCellVars(VTKPiece &vtkPiece)
{
    FloatArray valueArray;
    int numCells = vtkPiece.giveNumberOfCells();
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        const char *name = __InternalStateTypeToString(type);
        ( void ) name; //silence the warning

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >cellVarsArray = vtkSmartPointer< vtkDoubleArray >::New();
        cellVarsArray->SetName(name);
        cellVarsArray->SetNumberOfComponents(ncomponents);
        cellVarsArray->SetNumberOfTuples(numCells);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(i, ielem);
            for ( int i = 1; i <= ncomponents; ++i ) {
                cellVarsArray->SetComponent(ielem - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKCellData(name, cellVarsArray);

#else
        this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        valueArray.resize(ncomponents);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(i, ielem);
            this->writeVTKCellData(valueArray);
        }
        this->fileStream << "</DataArray>\n";
#endif
    
#ifdef _PYBIND_BINDINGS
        if ( pythonExport ) {
            py::list vals;
            for ( int ielem = 1; ielem <= numCells; ielem++ ) {
                valueArray = vtkPiece.giveCellVar(i, ielem);
                vals.append(valueArray);
            }
            this->Py_CellVars[name] = vals;
        }
 #endif        
    }//end of for
}


void
VTKXMLExportModule::computeIPAverage(FloatArray &answer, IntegrationRule *iRule, Element *elem, InternalStateType isType, TimeStep *tStep)
{
    // Computes the volume average (over an element) for the quantity defined by isType
    double gptot = 0.0;
    answer.clear();
    FloatArray temp;
    if ( iRule ) {
        for ( IntegrationPoint *ip : * iRule ) {
            elem->giveIPValue(temp, ip, isType, tStep);
            gptot += ip->giveWeight();
            answer.add(ip->giveWeight(), temp);
        }

        answer.times(1. / gptot);
    }
}


void
VTKXMLExportModule::writeVTKCollection()
{
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);
    char buff [ 1024 ];
    std::string fname;

    if ( tstep_substeps_out_flag ) {
        fname = this->emodel->giveOutputBaseFileName() + ".m" + std::to_string(this->number) + ".substep.pvd";
    } else {
        fname = this->emodel->giveOutputBaseFileName() + ".m" + std::to_string(this->number) + ".pvd";
    }

    std::ofstream streamP;
    if ( pythonExport ) {
        streamP = std::ofstream(NULL_DEVICE);//do not write anything
    } else {
        streamP = std::ofstream(fname.c_str() );
    }

    if ( !streamP.good() ) {
        OOFEM_ERROR("failed to open file %s", fname.c_str() );
    }

    sprintf(buff, "<!-- Computation started %d-%02d-%02d at %02d:%02d:%02d -->\n", current->tm_year + 1900, current->tm_mon + 1, current->tm_mday, current->tm_hour,  current->tm_min,  current->tm_sec);
    //     outfile << buff;

    streamP << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for ( auto pvd : this->pvdBuffer ) {
        streamP << pvd << "\n";
    }

    streamP << "</Collection>\n</VTKFile>";

    if (streamP){
        streamP.close();
    }
}

void
VTKXMLExportModule::writeGPVTKCollection()
{
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);
    char buff [ 1024 ];
    std::string fname;

    if ( tstep_substeps_out_flag ) {
        fname = this->emodel->giveOutputBaseFileName() + ".m" + std::to_string(this->number) + ".substep.gp.pvd";
    } else {
        fname = this->emodel->giveOutputBaseFileName() + ".m" + std::to_string(this->number) + ".gp.pvd";
    }

    std::ofstream outfile(fname.c_str() );

    sprintf(buff, "<!-- Computation started %d-%02d-%02d at %02d:%02d:%02d -->\n", current->tm_year + 1900, current->tm_mon + 1, current->tm_mday, current->tm_hour,  current->tm_min,  current->tm_sec);
    //     outfile << buff;

    outfile << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for ( auto pvd : this->gpPvdBuffer ) {
        outfile << pvd << "\n";
    }

    outfile << "</Collection>\n</VTKFile>";

    if (outfile){
        outfile.close();
    }
}




// Export of composite elements

void VTKXMLExportModule::exportCompositeElement(VTKPiece &vtkPiece, Element *el, TimeStep *tStep)
{
    VTKXMLExportModuleElementInterface *interface =
        static_cast< VTKXMLExportModuleElementInterface * >( el->giveInterface(VTKXMLExportModuleElementInterfaceType) );
    if ( interface ) {
        interface->giveCompositeExportData(vtkPiece, this->primaryVarsToExport, this->internalVarsToExport, this->cellVarsToExport, tStep);

        //this->writeVTKPiece(this->defaultVTKPiece, tStep);
    }
}

void VTKXMLExportModule::exportCompositeElement(std::vector< VTKPiece > &vtkPieces, Element *el, TimeStep *tStep)
{
    VTKXMLExportModuleElementInterface *interface =
        static_cast< VTKXMLExportModuleElementInterface * >( el->giveInterface(VTKXMLExportModuleElementInterfaceType) );
    if ( interface ) {
        interface->giveCompositeExportData(vtkPieces, this->primaryVarsToExport, this->internalVarsToExport, this->cellVarsToExport, tStep);

        //this->writeVTKPiece(this->defaultVTKPiece, tStep);
    }
}

void
VTKPiece::clear()
{
    ///@todo Will this give a memory leak? / JB
    numCells = 0;
    numNodes = 0;
    this->connectivity.clear();
    this->elCellTypes.clear();
    this->elOffsets.clear();
    this->elVars.clear();
    this->nodeCoords.clear();
    this->nodeVars.clear();
    this->nodeVarsFromIS.clear();
    this->nodeVarsFromXFEMIS.clear();
}


NodalRecoveryModel *
VTKXMLExportModule::giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( !this->smoother ) {
        this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
    }

    return this->smoother.get();
}


NodalRecoveryModel *
VTKXMLExportModule::givePrimVarSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( !this->primVarSmoother ) {
        this->primVarSmoother = classFactory.createNodalRecoveryModel(NodalRecoveryModel::NRM_NodalAveraging, d);
    }

    return this->primVarSmoother.get();
}


void
VTKXMLExportModule::exportIntVarsInGpAs(IntArray valIDs, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nc = 0;
    ( void ) nc; //silence the warning
    FloatArray gc, value;
    std::ofstream stream;
    InternalStateType isttype;
    InternalStateValueType vtype;
    std::string scalars, vectors, tensors;

    // output nodes Region By Region
    int nregions = this->giveNumberOfRegions(); // aka sets
    // open output stream
    std::string outputFileName = this->giveOutputBaseFileName(tStep) + ".gp.vtu";
    std::ofstream streamG;
    if ( pythonExport ) {
        streamG = std::ofstream(NULL_DEVICE);
    } else {
        streamG = std::ofstream(outputFileName);
    }

    if ( !streamG.good() ) {
        OOFEM_ERROR("failed to open file %s", outputFileName.c_str() );
    }

    streamG << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    streamG << "<UnstructuredGrid>\n";

    /* loop over regions */
    for ( int ireg = 1; ireg <= nregions; ireg++ ) {
        const IntArray &elements = this->giveRegionSet(ireg)->giveElementList();
        int nip = 0;
        for ( int i = 1; i <= elements.giveSize(); i++ ) {
            nip += d->giveElement(elements.at(i) )->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
        }

        //Create one cell per each GP
        streamG << "<Piece NumberOfPoints=\"" << nip << "\" NumberOfCells=\"" << nip << "\">\n";
        streamG << "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ";
        for ( int i = 1; i <= elements.giveSize(); i++ ) {
            int ielem = elements.at(i);

            for ( GaussPoint *gp : * d->giveElement(ielem)->giveDefaultIntegrationRulePtr() ) {
                d->giveElement(ielem)->computeGlobalCoordinates(gc, gp->giveNaturalCoordinates() );
                for ( double c : gc ) {
                    ( void ) c; //silence the warning
                    streamG << scientific << c << " ";
                }

                for ( int k = gc.giveSize() + 1; k <= 3; k++ ) {
                    streamG << scientific << 0.0 << " ";
                }
            }
        }

        streamG << " </DataArray>\n";
        streamG << "</Points>\n";
        streamG << "<Cells>\n";
        streamG << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
        for ( int j = 0; j < nip; j++ ) {
            streamG << j << " ";
        }

        streamG << " </DataArray>\n";
        streamG << " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
        for ( int j = 1; j <= nip; j++ ) {
            streamG << j << " ";
        }

        streamG << " </DataArray>\n";
        streamG << " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
        for ( int j = 1; j <= nip; j++ ) {
            streamG << "1 ";
        }

        streamG << " </DataArray>\n";
        streamG << "</Cells>\n";
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
            } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
                tensors += __InternalStateTypeToString(isttype);
                tensors.append(" ");
            } else {
                OOFEM_WARNING("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }
        }

        // print collected data summary in header
        streamG << "<PointData Scalars=\"" << scalars.c_str() << "\" Vectors=\"" << vectors.c_str() << "\" Tensors=\"" << tensors.c_str() << "\" >\n";
        scalars.clear();
        vectors.clear();
        tensors.clear();

        // export actual data, loop over individual IDs to export
        for ( int vi = 1; vi <= valIDs.giveSize(); vi++ ) {
            isttype = ( InternalStateType ) valIDs.at(vi);
            vtype = giveInternalStateValueType(isttype);
            if ( vtype == ISVT_SCALAR ) {
                nc = 1;
            } else if ( vtype == ISVT_VECTOR ) {
                nc = 3;
                if ( isttype == IST_BeamForceMomentTensor ) { //AS: to make the hack work
                    nc = 6;
                }
            } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
                nc = 9;
            } else {
                OOFEM_WARNING("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }

            streamG << "  <DataArray type=\"Float64\" Name=\"" << __InternalStateTypeToString(isttype) << "\" NumberOfComponents=\"" << nc << "\" format=\"ascii\">";
            for ( int i = 1; i <= elements.giveSize(); i++ ) {
                int ielem = elements.at(i);

                // loop over default IRule gps
                for ( GaussPoint *gp : * d->giveElement(ielem)->giveDefaultIntegrationRulePtr() ) {
                    d->giveElement(ielem)->giveIPValue(value, gp, isttype, tStep);

                    if ( vtype == ISVT_VECTOR ) {
                        // bp: hack for BeamForceMomentTensor, which should be splitted into force and momentum vectors
                        if ( isttype == IST_BeamForceMomentTensor ) {
                            value.resizeWithValues(6);
                        } else {
                            value.resizeWithValues(3);
                        }
                    } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
                        FloatArray help = value;
                        this->makeFullTensorForm(value, help, vtype);
                    }

                    for ( double v : value ) {
                        ( void ) v; //silence the warning
                        streamG << scientific << v << " ";
                    }
                } // end loop over IPs
            } // end loop over elements

            streamG << "  </DataArray>\n";
        } // end loop over values to be exported
        streamG << "</PointData>\n</Piece>\n";
    } // end loop over regions

    streamG << "</UnstructuredGrid>\n";
    streamG << "</VTKFile>\n";
    if(streamG){
        streamG.close();
    }
}
} // end namespace oofem
