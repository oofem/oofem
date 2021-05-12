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

#include "nodalaveragingrecoverymodel.h"
#include "zznodalrecoverymodel.h"

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


VTKXMLExportModule::VTKXMLExportModule(int n, EngngModel *e) : VTKBaseExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
#ifdef _PYBIND_BINDINGS
,Py_PrimaryVars(), Py_IntVars(), Py_CellVars(), Py_Nodes(), Py_Elements()
#endif
{}


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
}


void
VTKXMLExportModule::initialize()
{
    this->smoother = nullptr;
    this->primVarSmoother = nullptr;
    VTKBaseExportModule::initialize();
}


void
VTKXMLExportModule::terminate()
{ }

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

    /* Loop over pieces  ///@todo: this feature has been broken but not checked if it currently works /JB
        * Start default pieces containing all single cell elements. Elements built up from several vtk
        * cells (composite elements) are exported as individual pieces after the default ones.
        */
    int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;
    NodalRecoveryModel *smoother = giveSmoother();
    NodalRecoveryModel *primVarSmoother = givePrimVarSmoother();

    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        // Fills a data struct (VTKPiece) with all the necessary data.
        Set* region = this->giveRegionSet(pieceNum);
        this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
        this->writeVTKPieceProlog(this->defaultVTKPiece, tStep); 
        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(this->defaultVTKPiece, *region, primaryVarsToExport, *primVarSmoother, tStep);
        this->exportIntVars(this->defaultVTKPiece, *region, internalVarsToExport, *smoother, tStep);
        this->exportExternalForces(this->defaultVTKPiece, *region, externalForcesToExport, tStep);
        this->exportCellVars(this->defaultVTKPiece, *region, cellVarsToExport, tStep);

        // Write the VTK piece to file.
        anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPiece, tStep);
        this->writeVTKPieceEpilog(this->defaultVTKPiece, tStep);   
        this->defaultVTKPiece.clear();
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
                    this->writeVTKPieceProlog(this->defaultVTKPieces[j], tStep);          
                    anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPieces [ j ],  tStep);
                    this->writeVTKPieceEpilog(this->defaultVTKPieces[j], tStep);  
                    this->defaultVTKPieces [ j ].clear();
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


bool
VTKXMLExportModule::writeVTKPieceProlog(VTKPiece &vtkPiece, TimeStep *tStep)
{
   // Writes a VTK piece header + geometry to file.
   // This could be the whole domain (most common case) or it can be a
   // (so-called) composite element consisting of several VTK cells (layered structures, XFEM, etc.).

    // Write output: node coords
    int numNodes = vtkPiece.giveNumberOfNodes();
    int numEl = vtkPiece.giveNumberOfCells();
    FloatArray coords;

    if ( !vtkPiece.giveNumberOfCells() ) {
      return false;
    }

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
#endif
    return true;

}

bool
VTKXMLExportModule::writeVTKPieceEpilog(VTKPiece &vtkPiece, TimeStep *tStep)
{
    if ( !vtkPiece.giveNumberOfCells() ) {
        return false;
    }

#ifndef __VTK_MODULE
    this->fileStream << "</Piece>\n";
#endif
    return true;
}



bool
VTKXMLExportModule::writeVTKPieceVariables(VTKPiece &vtkPiece, TimeStep *tStep)
{
    // Write a VTK piece variables to file.
    // This could be the whole domain (most common case) or it can be a
    // (so-called) composite element consisting of several VTK cells (layered structures, XFEM, etc.).

    if ( !vtkPiece.giveNumberOfCells() ) {
        return false;
    }


#ifndef __VTK_MODULE 
    ///@todo giveDataHeaders is currently not updated wrt the new structure -> no file names in headers /JB
    std::string pointHeader, cellHeader;
    this->giveDataHeaders(pointHeader, cellHeader);

    this->fileStream << pointHeader.c_str();
#endif

    this->writePrimaryVars(vtkPiece);       // Primary field
    this->writeIntVars(vtkPiece);           // Internal State Type variables smoothed to the nodes
    this->writeExternalForces(vtkPiece);           // External forces

    //if ( emodel->giveDomain(1)->hasXfemManager() ) {
    //    this->writeXFEMVars(vtkPiece);      // XFEM State Type variables associated with XFEM structure
    //}

#ifndef __VTK_MODULE
    this->fileStream << "</PointData>\n";
    this->fileStream << cellHeader.c_str();
#endif
    this->writeCellVars(vtkPiece);          // Single cell variables ( if given in the integration points then an average will be exported)

#ifndef __VTK_MODULE
    this->fileStream << "</CellData>\n";
#endif
    return true;
}

#ifndef __VTK_MODULE
void
VTKXMLExportModule::giveDataHeaders(std::string &pointHeader, std::string &cellHeader)
{
    std::string scalars, vectors, tensors;

    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        if ( type == DisplacementVector || type == EigenVector || type == VelocityVector || type == DirectorField || type == MacroSlipVector ) {
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
        valueArray = vtkPiece.giveInternalVarInNode(type, 1);
        ncomponents = valueArray.giveSize();
        ( void ) ncomponents;//silence warning

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
        varArray->SetName(name);
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.giveInternalVarInNode(type, inode);
            for ( int i = 1; i <= ncomponents; ++i ) {
                varArray->SetComponent(inode - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else

        this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            valueArray = vtkPiece.giveInternalVarInNode(type, inode);
            this->writeVTKPointData(valueArray);
        }

#endif
        // Footer
#ifndef __VTK_MODULE
        this->fileStream << "</DataArray>\n";
#endif
    
#ifdef _PYBIND_BINDINGS
#if 0
        if ( pythonExport ) {
            py::list vals;
            for ( int inode = 1; inode <= numNodes; inode++ ) {
                valueArray = vtkPiece.giveInternalVarInNode(i, inode);
                vals.append(valueArray);
            }
            this->Py_IntVars[name] = vals;
        }
#endif
#endif
    } //end of for
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
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(type, inode);
            for ( int j = 1; j <= ncomponents; ++j ) {
                varArray->SetComponent(inode - 1, j - 1, valueArray.at(j) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else
        this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(type, inode);
            this->writeVTKPointData(valueArray);
        }
        this->fileStream << "</DataArray>\n";

 #ifdef _PYBIND_BINDINGS
   #if 0
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
#endif
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
            valueArray = vtkPiece.giveCellVar(type, ielem);
            this->writeVTKCellData(valueArray);
        }
        this->fileStream << "</DataArray>\n";
#endif
    
#ifdef _PYBIND_BINDINGS
        if ( pythonExport ) {
            py::list vals;
            for ( int ielem = 1; ielem <= numCells; ielem++ ) {
                valueArray = vtkPiece.giveCellVar(type, ielem);
                vals.append(valueArray);
            }
            this->Py_CellVars[name] = vals;
        }
 #endif        
    }//end of for
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
