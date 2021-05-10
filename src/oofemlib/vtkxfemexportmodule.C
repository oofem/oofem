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

#include "vtkxfemexportmodule.h"
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
REGISTER_ExportModule(VTKXMLXFemExportModule)



//
// VTKXMLXFemExport module
//

VTKXMLXFemExportModule::VTKXMLXFemExportModule(int n, EngngModel *e) : VTKXMLExportModule(n, e) {}

VTKXMLXFemExportModule::~VTKXMLXFemExportModule() { }

void
VTKXMLXFemExportModule::initializeFrom(InputRecord &ir)
{
    ExportModule::initializeFrom(ir);
}


std::string
VTKXMLXFemExportModule::giveOutputFileName(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".vtu";
}


std::ofstream
VTKXMLXFemExportModule::giveOutputStream(TimeStep *tStep)
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
VTKXMLXFemExportModule::exportIntVars(VTKPiece &vtkPiece, Set& region, int field, int enrItIndex,  IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray answer;
    //IntArray& mapG2L = vtkPiece.getMapG2L();
    IntArray& mapL2G = vtkPiece.getMapL2G();

    // Export of XFEM related quantities
    if ( d->hasXfemManager() ) {
        XfemManager *xFemMan = d->giveXfemManager();
        vtkPiece.setNumberOfInternalVarsToExport(xFemMan->vtkExportFields, mapL2G.giveSize() );
        XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];

        for ( int nodeIndx = 1; nodeIndx <= mapL2G.giveSize(); nodeIndx++ ) {
          Node *node = d->giveNode(mapL2G.at(nodeIndx) );
          getNodalVariableFromXFEMST(answer, node, tStep, xfemstype, region, xFemMan->giveEnrichmentItem(enrItIndex) );
          vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeIndx, answer);
        }
    }
}

void
VTKXMLXFemExportModule::getNodalVariableFromXFEMST(FloatArray &answer, Node *node, TimeStep *tStep, XFEMStateType xfemstype, Set &region, EnrichmentItem *ei)
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


bool
VTKXMLXFemExportModule::writeXFEMVars(VTKPiece &vtkPiece, int field, int enrItIndex)
{
    Domain *d = emodel->giveDomain(1);
    XfemManager *xFemMan = d->giveXfemManager();
    FloatArray valueArray;
    
    if ( !vtkPiece.giveNumberOfCells() ) {
      return false;
    }

    XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields.at(field);
    const char *namePart = __XFEMStateTypeToString(xfemstype);
    InternalStateValueType valType = xFemMan->giveXFEMStateValueType(xfemstype);
    int ncomponents = giveInternalStateTypeSize(valType);
    ( void ) ncomponents; //silence the warning

    int numNodes = vtkPiece.giveNumberOfNodes();

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
    this->fileStream << "<DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
    for ( int inode = 1; inode <= numNodes; inode++ ) {
      valueArray = vtkPiece.giveInternalXFEMVarInNode(field, enrItIndex, inode);
      this->writeVTKPointData(valueArray);
    }
    this->fileStream << "</DataArray>\n";
#endif
    return true;
}

void
VTKXMLXFemExportModule::giveDataHeaders(std::string &pointHeader, std::string &cellHeader)
{
    std::string scalars, vectors, tensors;

    Domain *d = emodel->giveDomain(1);
    XfemManager *xFemMan = d->giveXfemManager();
    int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();
    for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
        for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
            XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields.at(field);
            const char *namePart = __XFEMStateTypeToString(xfemstype);

            // Header
            char name [ 100 ];   // Must I define a fixed size? /JB
            sprintf(name, "%s_%d ", namePart, xFemMan->giveEnrichmentItem(enrItIndex)->giveNumber() );
            scalars += name;
            scalars.append(" ");
        }
    }
    // print header
    pointHeader = "<PointData Scalars=\"" + scalars + "\" "
                  +  "Vectors=\"" + vectors + "\" "
                  +  "Tensors=\"" + tensors + "\" >\n";

}


void
VTKXMLXFemExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    Domain *d = emodel->giveDomain(1);
    //XfemManager *xFemMan = d->giveXfemManager();
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

    /* Loop over pieces  ///@todo: this feature has been broken but not checked if it currently works /JB
     * Start default pieces containing all single cell elements. Elements built up from several vtk
     * cells (composite elements) are exported as individual pieces after the default ones.
     */
    int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;

    // Export of XFEM related quantities

    if ( d->hasXfemManager() ) {
      XfemManager *xFemMan = d->giveXfemManager();
      int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();

      for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        // Fills a data struct (VTKPiece) with all the necessary data.
        Set* region = this->giveRegionSet(pieceNum);
        this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
        if (this->defaultVTKPiece.giveNumberOfNodes() > 0) {
          this->writeVTKPieceProlog(this->defaultVTKPiece, tStep);   
          this->defaultVTKPiece.setNumberOfInternalVarsToExport(xFemMan->vtkExportFields, this->defaultVTKPiece.getMapL2G().giveSize() );

          std::string pointHeader, cellHeader;
          this->giveDataHeaders(pointHeader, cellHeader);
          this->fileStream << pointHeader.c_str();
          
          for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
            //XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];
            
            for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
              // Fills a data struct (VTKPiece) with all the necessary data.
              //this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
              this->exportIntVars(this->defaultVTKPiece, *region, field, enrItIndex, internalVarsToExport, *smoother, tStep);
              // Write the VTK piece to file.
              anyPieceNonEmpty +=this->writeXFEMVars(this->defaultVTKPiece, field, enrItIndex);
            }
          }
          
          this->fileStream << "</PointData>\n";
          this->writeVTKPieceEpilog(this->defaultVTKPiece, tStep);   
          this->defaultVTKPiece.clear();
        }
      }

      /*
       * Output all composite elements - one piece per composite element
       * Each element is responsible of setting up a VTKPiece which can then be exported
       */
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
              std::string pointHeader, cellHeader;
              this->giveDataHeaders(pointHeader, cellHeader);
              this->fileStream << pointHeader.c_str();
              
              for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
                for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) { 
                  anyPieceNonEmpty += this->writeXFEMVars(this->defaultVTKPieces [ j ],  field, enrItIndex);
                }
              }
              this->fileStream << "</PointData>\n";
              this->writeVTKPieceEpilog(this->defaultVTKPieces[j], tStep); 
            }
#else
            // No support for binary export yet
#endif
          }
        } // end loop over composite elements
      }
    }

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
}

} // end namespace oofem
