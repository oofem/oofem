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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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

REGISTER_ExportModule( VTKXMLExportModule )

VTKXMLExportModule :: VTKXMLExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{
    primVarSmoother = NULL;
    smoother = NULL;
    redToFull.setValues(6, 1, 5, 9, 6, 3, 2);//position of xx, yy, zz, yz, xz, xy in tensor 
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

    val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_VTKXMLExportModule_stype); // Macro
    stype = ( NodalRecoveryModel::NodalRecoveryModelType ) val;

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
VTKXMLExportModule :: makeFullForm(FloatArray &answer, const FloatArray &reducedForm, InternalStateValueType type, const IntArray &redIndx)
{
    answer.resize(9);
    answer.zero();
    if ( type == ISVT_TENSOR_S3 ) {
        for (int i = 1; i <= redIndx.giveSize(); i++) {
            if (redIndx.at(i) > 0) {
                answer.at(redToFull.at(i)) = reducedForm.at(redIndx.at(i));
            }
        }
    } else if ( type == ISVT_TENSOR_S3E ) {
        for (int i = 1; i <= redIndx.giveSize(); i++) {
            if (redIndx.at(i) > 3) {
                answer.at(redToFull.at(i)) = reducedForm.at(i)*0.5;
            } else if (redIndx.at(i) > 0) {
                answer.at(redToFull.at(i)) = reducedForm.at(i);
            }
        }
    }
    // Symmetrize
    answer.at(4) = answer.at(2);
    answer.at(7) = answer.at(3);
    answer.at(8) = answer.at(6);
}

std::string
VTKXMLExportModule :: giveOutputFileName(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".vtu";
}


FILE *
VTKXMLExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;
    std::string fileName = giveOutputFileName(tStep);
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2("VTKXMLExportModule::giveOutputStream: failed to open file %s", fileName.c_str());
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
    } else if ( elemGT == EGT_quad_1 ) {
        vtkCellType = 9;
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
VTKXMLExportModule :: giveElementCell(IntArray &answer, Element *elem, int cell)
//VTKXMLExportModule :: giveElementCell(IntArray &answer, Element_Geometry_Type elemGT, int cell)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int nelemNodes;

    if ( ( elemGT == EGT_point ) ||
        ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
        ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
        ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
        ( elemGT == EGT_hexa_1 ) ||
        (elemGT == EGT_wedge_1) ) {
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber() ;
        }
    } else if ( elemGT == EGT_hexa_27 ) {
        int HexaQuadNodeMapping [] = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18, 23, 25, 26, 24, 22, 21, 27
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(HexaQuadNodeMapping [ i - 1 ])->giveNumber() ;
        }
    } else if ( elemGT == EGT_hexa_2 ) {
        int HexaQuadNodeMapping [] = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(HexaQuadNodeMapping [ i - 1 ])->giveNumber() ;
        }
    } else if ( elemGT == EGT_wedge_2 ) {
        int WedgeQuadNodeMapping [] = { 4, 6, 5, 1, 3, 2, 12, 11, 10, 9, 8, 7, 13, 15,14 };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(WedgeQuadNodeMapping [ i - 1 ])->giveNumber() ;
        }
    } else if ( elemGT == EGT_quad9_2 ) {

        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( int i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber() ;
        }

    } else {
        OOFEM_ERROR("VTKXMLExportModule: unsupported element geometry type");
    }
}


bool
VTKXMLExportModule :: isElementComposite(Element *elem)
{
    //return false;
    return ( elem->giveGeometryType() == EGT_Composite );
}


int
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
VTKXMLExportModule :: doOutput(TimeStep *tStep)
{
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkUnstructuredGrid> stream = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> nodes = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdList> elemNodeArray = vtkSmartPointer<vtkIdList>::New();
#else
    FILE *stream = this->giveOutputStream(tStep);
    struct tm *current;
    time_t now;
    time(&now);
    current = localtime(&now);
    fprintf(stream, "<!-- TimeStep %e Computed %d-%02d-%02d at %02d:%02d:%02d -->\n", tStep->giveIntrinsicTime(), current->tm_year+1900, current->tm_mon+1, current->tm_mday, current->tm_hour,  current->tm_min,  current->tm_sec);
    fprintf(stream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(stream, "<UnstructuredGrid>\n");
#endif


    Domain *d  = emodel->giveDomain(1);
    Element *elem;
    FloatArray *coords;
    int nelem = d->giveNumberOfElements();

    this->giveSmoother(); // make sure smoother is created

    // output nodes Region By Region
    int nregions = this->smoother->giveNumberOfVirtualRegions();
    int regionDofMans, totalcells;
    IntArray mapG2L, mapL2G;

    /* loop over regions */
    for ( int ireg = 1; ireg <= nregions; ireg++ ) {
        if ( ( ireg > 0 ) && ( this->regionsToSkip.contains(ireg) ) ) {
            continue;
        }

        // assemble local->global and global->local region map
        // and get number of single cells to process
        // the composite cells exported individually
        this->initRegionNodeNumbering(mapG2L, mapL2G, regionDofMans, totalcells, d, ireg);

        /* start default piece containing all single cell elements
         * the elements with composite geometry are assumed to be exported in individual pieces
         * after the default one
         */
#ifndef __PARALLEL_MODE
        if ( regionDofMans && totalcells ) {
#else
        if ( 1 ) {
#endif

#ifdef __VTK_MODULE
            for ( int inode = 1; inode <= regionDofMans; inode++ ) {
                coords = d->giveNode(mapL2G.at(inode))->giveCoordinates();
                int dims = coords->giveSize();
                nodes->InsertNextPoint(coords->at(1), dims >= 2 ? coords->at(2) : 0.0, dims >= 3 ? coords->at(3) : 0.0);
            }
            stream->SetPoints(nodes);
#else
            fprintf(stream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", regionDofMans, totalcells);
            // export nodes in region as vtk vertices
            fprintf(stream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ");
            for ( int inode = 1; inode <= regionDofMans; inode++ ) {
                coords = d->giveNode( mapL2G.at(inode) )->giveCoordinates();
                for ( int i = 1; i <= coords->giveSize(); i++ ) {
                    fprintf( stream, "%e ", coords->at(i) );
                }

                for ( int i = coords->giveSize() + 1; i <= 3; i++ ) {
                    fprintf(stream, "%e ", 0.0);
                }
            }
            fprintf(stream, "</DataArray>\n</Points>\n");
#endif


            // output all cells of the piece
            int nelemNodes;
            IntArray cellNodes;
#ifdef __VTK_MODULE
            stream->Allocate(nelem);
#else
            fprintf(stream, "<Cells>\n");
            // output the connectivity data
            fprintf(stream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ");
#endif
            for ( int ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
                if ( ( ireg > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != ireg ) ) {
                    continue;
                }

                if ( this->isElementComposite(elem) ) {
                    continue;                                  // composite cells exported individually
                }

                if ( !elem-> isActivated(tStep) ) {                  //skip inactivated elements
                    continue;
                }
                
#ifdef __PARALLEL_MODE
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }
#endif


                nelemNodes = elem->giveNumberOfNodes();
                this->giveElementCell(cellNodes, elem, 0);
                //Element_Geometry_Type elemGT = elem->giveGeometryType();
                //this->giveElementCell(cellNodes, elemGT, 0);
#ifdef __VTK_MODULE
                elemNodeArray->Reset();
                elemNodeArray->SetNumberOfIds(nelemNodes);
#endif
                for ( int i = 1; i <= nelemNodes; i++ ) {
#ifdef __VTK_MODULE
                    elemNodeArray->SetId(i-1, mapG2L.at( cellNodes.at(i) ) - 1);
#else
                    fprintf(stream, "%d ", mapG2L.at( cellNodes.at(i) ) - 1);
#endif
                }
#ifdef __VTK_MODULE
                stream->InsertNextCell(this->giveCellType(elem), elemNodeArray);
#else
                fprintf(stream, " ");
#endif
            }

#ifndef __VTK_MODULE
            int vtkCellType;
            fprintf(stream, "</DataArray>\n");
            // output the offsets (index of individual element data in connectivity array)
            fprintf(stream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ");
            int offset = 0;
            for ( int ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
                if ( ( ireg > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != ireg ) ) {
                    continue;
                }

#ifdef __PARALLEL_MODE
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif
                offset += elem->giveNumberOfNodes();
                fprintf(stream, "%d ", offset);
            }

            fprintf(stream, "</DataArray>\n");

            // output cell (element) types
            fprintf(stream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
            for ( int ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
                if ( ( ireg > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != ireg ) ) {
                    continue;
                }

                if ( this->isElementComposite(elem) ) {
                    continue;                                  // composite cells exported individually
                }

#ifdef __PARALLEL_MODE
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif
                vtkCellType = this->giveCellType(elem);
                fprintf(stream, "%d ", vtkCellType);
            }

            fprintf(stream, "</DataArray>\n");
            fprintf(stream, "</Cells>\n");
#endif

            // export primary and internal variables
#ifndef __VTK_MODULE
            this->exportPointDataHeader(stream, tStep);
#endif
            this->exportPrimaryVars(stream, mapG2L, mapL2G, regionDofMans, ireg, tStep);
            this->exportIntVars(stream, mapG2L, mapL2G, regionDofMans, ireg, tStep);
#ifndef __VTK_MODULE
            fprintf(stream, "</PointData>\n");
#endif

            //export cell data
            this->exportCellVars(stream, ireg, tStep);

#ifndef __VTK_MODULE
            // end of piece record
            fprintf(stream, "</Piece>\n");
#endif
        } // end of default piece for simple geometry elements

#if 1
        // loop over region elements with multi-cell geometry
        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);

            if ( this->regionsToSkip.contains( this->smoother->giveElementVirtualRegionNumber(ielem) ) ) {
                continue;
            }

 #ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

 #endif

            if ( this->isElementComposite(elem) ) {
#ifndef __VTK_MODULE
                //VTKXMLExportModule :: exportCompositeElement(FILE *stream, VTKXMLExportModule *expModule, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, TimeStep *tStep) 
                this->exportCompositeElement(stream, elem, primaryVarsToExport, internalVarsToExport, tStep);
                ///@todo Not sure how to deal with this.
                // multi cell (composite) elements should support vtkxmlexportmoduleinterface
                // and are exported as individual pieces (see VTKXMLExportModuleElementInterface)
               /*
                VTKXMLExportModuleElementInterface *interface =
                    static_cast< VTKXMLExportModuleElementInterface * >( elem->giveInterface(VTKXMLExportModuleElementInterfaceType) );
                if ( interface ) {
                    // passing this to access general piece related methods like exportPointDataHeader, etc.
                    //interface->_export(stream, this, primaryVarsToExport, internalVarsToExport, tStep);
                    interface->exportCompositeElement(stream, this, primaryVarsToExport, internalVarsToExport, tStep);
                }
                */
#endif
            }
        } // end loop over multi-cell elements
#endif
    } // end loop over regions

    std::string fname = giveOutputFileName(tStep);

#ifdef __VTK_MODULE
#if 0
    // Doesn't as well as I would want it to, interface to VTK is to limited to control this.
    // * The PVTU-file is written by every process (seems to be impossible to avoid).
    // * Part files are renamed and time step and everything else is cut off => name collisions
    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
    writer->SetTimeStep(tStep->giveNumber()-1);
    writer->SetNumberOfPieces( this->emodel->giveNumberOfProcesses() );
    writer->SetStartPiece( this->emodel->giveRank() );
    writer->SetEndPiece( this->emodel->giveRank() );
#else
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
#endif
    writer->SetFileName(fname.c_str());
    writer->SetInput(stream);
    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    //writer->SetDataModeToAscii();
    writer->Write();
#else
    // finish unstructured grid data and vtk file
    fprintf(stream, "</UnstructuredGrid>\n</VTKFile>");
    fclose(stream);
#endif

    // Writing the *.pvd-file. Only time step information for now. It's named "timestep" but is actually the total time.
    // First we check to see that there are more than 1 time steps, otherwise it is redundant;
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() && emodel->giveRank() == 0 ) {
        ///@todo Should use probably use PVTU-files instead. It is starting to get messy.
        // For this to work, all processes must have an identical output file name.
        for (int i = 0; i < this->emodel->giveNumberOfProcesses(); ++i) {
            std::ostringstream pvdEntry;
            char fext[100];
            if (this->emodel->giveNumberOfProcesses() > 1) {
                sprintf( fext, "_%03d.m%d.%d", i, this->number, tStep->giveNumber() );
            } else {
                sprintf( fext, "m%d.%d", this->number, tStep->giveNumber() );
            }
            pvdEntry << "<DataSet timestep=\"" << tStep->giveIntrinsicTime() << "\" group=\"\" part=\"" << i << "\" file=\""
                    << this->emodel->giveOutputBaseFileName() << fext << ".vtu\"/>";
            this->pvdBuffer.push_back(pvdEntry.str());
        }
        this->writeVTKCollection();
    } else
#endif
    if ( !emodel->isParallel() && tStep->giveNumber() >= 1 ) { // For non-parallel enabled OOFEM, then we only check for multiple steps.
        std::ostringstream pvdEntry;
        pvdEntry << "<DataSet timestep=\"" << tStep->giveIntrinsicTime() << "\" group=\"\" part=\"\" file=\"" << fname << "\"/>";
        this->pvdBuffer.push_back(pvdEntry.str());
        this->writeVTKCollection();
    }
}

#ifndef __VTK_MODULE
void
VTKXMLExportModule :: exportPointDataHeader(FILE *stream, TimeStep *tStep)
{
    int n;
    std :: string scalars, vectors, tensors;

    n = primaryVarsToExport.giveSize();

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
        } else if ( vtype == ISVT_TENSOR_G ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else {
            fprintf( stderr, "VTKXMLExportModule::exportIntVars: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
        }
    }

    // print header
    fprintf( stream, "<PointData Scalars=\"%s\" Vectors=\"%s\" Tensors=\"%s\" >\n",
            scalars.c_str(), vectors.c_str(), tensors.c_str() );
}
#endif


//keyword "vars" in OOFEM input file
void
VTKXMLExportModule :: exportIntVars(
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
    FILE *stream,
#endif
    IntArray &mapG2L, IntArray &mapL2G, int regionDofMans, int region, TimeStep *tStep)
{
    int n = internalVarsToExport.giveSize();
    InternalStateType isttype;
    InternalStateValueType vtype;

    this->giveSmoother()->init(); // Makes sure smoother is up-to-date with potentially new mesh.
    // should be performed over regions
    for ( int i = 1; i <= n; i++ ) {
        isttype = ( InternalStateType ) internalVarsToExport.at(i);
        vtype = giveInternalStateValueType(isttype);
        //this->exportIntVarAs(isttype, vtype, mapG2L, mapL2G, regionDofMans, region, stream, tStep);
        this->exportIntVarAs(isttype, mapG2L, mapL2G, regionDofMans, region, stream, tStep);
    }
}


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

        if ( !element-> isActivated(domain->giveEngngModel()->giveCurrentStep()) ) {                  //skip inactivated elements
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




NodalRecoveryModel *
VTKXMLExportModule :: giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->smoother == NULL ) {
      this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
      this->smoother->setRecoveryMode (nvr, vrmap);
    }
    return this->smoother;
}


NodalRecoveryModel *
VTKXMLExportModule :: givePrimVarSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->primVarSmoother == NULL ) {
        this->primVarSmoother = classFactory.createNodalRecoveryModel(NodalRecoveryModel::NRM_NodalAveraging, d);
        this->primVarSmoother->setRecoveryMode (nvr, vrmap);
    }
    return this->primVarSmoother;
}


void
VTKXMLExportModule :: exportPrimaryVars(
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
    FILE *stream,
#endif
    IntArray &mapG2L, IntArray &mapL2G, int regionDofMans, int region, TimeStep *tStep)
{
    // should be performed over regions
    UnknownType type;

    this->givePrimVarSmoother()->init(); // Makes sure primary smoother is up-to-date with potentially new mesh.
    for (int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        type = ( UnknownType ) primaryVarsToExport.at(i);
        this->exportPrimVarAs(type, mapG2L, mapL2G, regionDofMans, region, stream, tStep);
    }
}









void
VTKXMLExportModule :: getPrimaryVariable(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int ireg)
{
    IntArray dofIDMask(3);
    int indx, size;
    DofIDItem id;
    const FloatArray *recoveredVal;

    // This code is not perfect. It should be rewritten to handle all cases more gracefully.

    EquationID eid = EID_MomentumBalance; // Shouldn't be necessary
    InternalStateType iState = IST_DisplacementVector; // Shouldn't be necessary

    dofIDMask.resize(0);
    if ( ( type == DisplacementVector ) || ( type == EigenVector ) || ( type == VelocityVector ) ) {
        dofIDMask.setValues(3, (int)Undef, (int) Undef, (int) Undef);
        for (int j = 1; j <= dman->giveNumberOfDofs(); j++ ) {
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
        eid = EID_ConservationEquation;
        iState = IST_MassConcentration_1;
        answer.resize(1);
    } else if ( type == Temperature ) {
        dofIDMask.followedBy(T_f);
        eid = EID_ConservationEquation;
        iState = IST_Temperature;
        answer.resize(1);
    } else if ( type == PressureVector ) {
        dofIDMask.followedBy(P_f);
        eid = EID_ConservationEquation;
        iState = IST_Pressure;
        answer.resize(1);
    } else if ( type == DirectorField ) {
        
        for (int j = 1; j <= dman->giveNumberOfDofs(); j++ ) {
            id = dman->giveDof(j)->giveDofID();
            if ( ( id == W_u ) || ( id == W_v ) || ( id == W_w )  ) {
                dofIDMask.followedBy(id);
            }
            //eid = EID_ConservationEquation;
            //iState = IST_DirectorField;
            answer.resize(3);
        }
        iState = IST_DirectorField;
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported unknownType %s", __UnknownTypeToString(type) );
    }

    size = dofIDMask.giveSize();
    answer.zero();

    for (int j = 1; j <= size; j++ ) {
        id = (DofIDItem)dofIDMask.at(j);
        if ( id == Undef ) {
            answer.at(j) = 0.;
        } else if ( iState == IST_DirectorField ) {
            indx = dman->findDofWithDofId( (DofIDItem)dofIDMask.at(j) );
            answer.at(j) = dman->giveDof(indx)->giveUnknown(VM_Total, tStep);
            
            this->givePrimVarSmoother()->recoverValues(iState, tStep); // recover values if not done before
            this->givePrimVarSmoother()->giveNodalVector(recoveredVal, dman->giveNumber(), ireg);
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING2("VTKXMLExportModule :: getDofManPrimaryVariable: recovered variable size mismatch for %d", type);
                answer.at(j) = 0.0;
            }
        } else if ( ( indx = dman->findDofWithDofId( id ) ) ) {
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
}




void VTKXMLExportModule :: exportCellVars(
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
    FILE *stream,
#endif
    int region, TimeStep *tStep)
{
    int n = cellVarsToExport.giveSize();
    if ( n == 0 ) return;

#ifndef __VTK_MODULE
    //print header
    fprintf(stream, "<CellData Scalars=\"\" Vectors=\"\" Tensors=\"\">\n");  // should contain a list of InternalStateType
#endif

    InternalStateType type;
    for ( int i = 1; i <= n; i++ ) {
        type = ( InternalStateType ) cellVarsToExport.at(i);
        this->exportCellVarAs(type, region, stream, tStep);
    }

#ifndef __VTK_MODULE
    //print footer
    fprintf(stream, "</CellData>\n");
#endif
}






//Export variables defined by the keyword "vars" in the input file //@todo change to the more verbose internalvars? (no backvar compatability)
void
VTKXMLExportModule :: exportIntVarAs(InternalStateType type, IntArray &mapG2L, IntArray &mapL2G, int regionDofMans, int ireg,
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
    FILE *stream,
#endif
    TimeStep *tStep)
{

    InternalStateValueType valType = giveInternalStateValueType(type);
    int ncomponents = giveInternalStateTypeSize(valType);

    // Header
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkDoubleArray> intVarArray = vtkSmartPointer<vtkDoubleArray>::New();
    intVarArray->SetName(__InternalStateTypeToString(type));
    
    intVarArray->SetNumberOfComponents(ncomponents);
    intVarArray->SetNumberOfTuples(regionDofMans);
    
#else
    fprintf( stream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", __InternalStateTypeToString(type), ncomponents );    
#endif

    this->giveSmoother();
    Domain *d = emodel->giveDomain(1);
    IntArray regionVarMap;

    if ( !( ( type == IST_DisplacementVector ) || ( type == IST_MaterialInterfaceVal ) ) ) {
        this->smoother->recoverValues(type, tStep);
        this->smoother->giveRegionRecordMap(regionVarMap, ireg, type);
    }

    FloatArray valueArray;
    for ( int inode = 1; inode <= regionDofMans; inode++ ) {
        Node *node = d->giveNode( mapL2G.at(inode) );
        this->getNodalVariableFromIS(valueArray, node, regionVarMap, tStep, type, ireg);
        

        // Write the data to file 
#ifdef __VTK_MODULE
        for (int i = 1; i <= ncomponents; ++i) {
            intVarArray->SetComponent(inode-1, i-1, answer.at(i));
        }
        switch ( ncomponents ) {
        case 1:
            stream->GetPointData()->SetActiveScalars(__InternalStateTypeToString(type));
            stream->GetPointData()->SetScalars(intVarArray);
            break;
        case 3
            stream->GetPointData()->SetActiveVectors(__InternalStateTypeToString(type));
            stream->GetPointData()->SetVectors(intVarArray);
        case 9
            stream->GetPointData()->SetActiveTensors(__InternalStateTypeToString(type));
            stream->GetPointData()->SetTensors(intVarArray);
            break;
        }
#else
        for ( int i = 1; i <= ncomponents; i++ ) {
            fprintf( stream, "%e ", valueArray.at(i) );
        }
#endif

    } 


#ifndef __VTK_MODULE
    fprintf(stream, "</DataArray>\n");  // footer
#endif


}


void
VTKXMLExportModule :: getNodalVariableFromIS(FloatArray &answer, Node *node, IntArray &regionVarMap, TimeStep *tStep, InternalStateType type, int ireg) 
{
    // Recovers nodal values from Internal States defined in the integration points. 
    // Should return an array with proper size supported by VTK (1, 3 or 9)
    const FloatArray *val;
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
            valueArray.resize( regionVarMap.giveSize() );
            val = & valueArray;
            //OOFEM_WARNING2("VTKXMLExportModule::exportIntVars: smoothing error: invalid data in node %d", inode);
        }
    }

    int ncomponents = giveInternalStateTypeSize(valType);
    answer.resize(ncomponents);
    int valSize = val->giveSize(); // size of recovered quantity

    // check if valSize corresponds to the expected size otherwise pad with zeros
    if ( valType == ISVT_SCALAR ) {
        answer.at(1) = valSize ? val->at(1) : 0.0;

    } else if ( valType == ISVT_VECTOR ) {
        int isize = min( valSize, 3 ); // so it will simply truncate larger arrays
        for ( int i = 1; i <= isize; i++ ) {
            answer.at(i) = val->at(i);
        }

    } else if ( valType == ISVT_TENSOR_S3 || valType == ISVT_TENSOR_S3E ) {
        this->makeFullForm(answer, *val, valType, regionVarMap);
      
    } else if ( valType == ISVT_TENSOR_G ) { // export general tensor values as scalars
        int isize = min( val->giveSize(), 9 );
        for ( int i = 1; i <= isize; i++ ) {
            answer.at(i) = val->at(i);
        }
    } else {
        OOFEM_ERROR("TKXMLExportModule ::getNodalVariableFromIS - ISVT_UNDEFINED encountered")
    }

}




void
VTKXMLExportModule :: exportPrimVarAs(UnknownType type, IntArray &mapG2L, IntArray &mapL2G,
                                      int regionDofMans, int ireg,
#ifdef __VTK_MODULE
                                      vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
                                      FILE *stream,
#endif
                                      TimeStep *tStep)
{
    

    InternalStateValueType valType = giveInternalStateValueType(type);
    int ncomponents = giveInternalStateTypeSize(valType);

    // Header
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkDoubleArray> primVarArray = vtkSmartPointer<vtkDoubleArray>::New();
    primVarArray->SetName(__UnknownTypeToString(type));
    
    primVarArray->SetNumberOfComponents(ncomponents);
    primVarArray->SetNumberOfTuples(regionDofMans);
    //fprintf( stderr, "VTKXMLExportModule::exportPrimVarAs: unsupported variable type %s\n", __UnknownTypeToString(type) );
#else
    fprintf( stream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", __UnknownTypeToString(type), ncomponents );    
    //fprintf( stderr, "VTKXMLExportModule::exportPrimVarAs: unsupported variable type %s\n", __UnknownTypeToString(type) );
#endif


    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray, valueArrayLCS;
    for ( int inode = 1; inode <= regionDofMans; inode++ ) {
        DofManager *dman = d->giveNode( mapL2G.at(inode) );

        this->getPrimaryVariable(valueArray, dman, tStep, type, ireg);

        if ( valType == ISVT_VECTOR ) {
            //rotate back from nodal CS to global CS if applies
            Node *node = dynamic_cast< Node* >( dman );
            if ( node && node->hasLocalCS() ) {
                valueArrayLCS = valueArray;
                valueArray.beTProductOf(* node->giveLocalCoordinateTriplet(), valueArrayLCS);
            }
        }

        // Write the data to file
#ifdef __VTK_MODULE
        for (int i = 1; i <= ncomponents; ++i) {
            primVarArray->SetComponent(inode-1, i-1, answer.at(i));
        }
        switch ( ncomponents ) {
        case 1:
            stream->GetPointData()->SetActiveScalars(__UnknownTypeToString(type));
            stream->GetPointData()->SetScalars(primVarArray);
            break;
        case 3
            stream->GetPointData()->SetActiveVectors(__UnknownTypeToString(type));
            stream->GetPointData()->SetVectors(primVarArray);
            break;
        }
#else
        for ( int i = 1; i <= ncomponents; i++ ) {
            fprintf( stream, "%e ", valueArray.at(i) );
        }
#endif

    } // end loop over nodes

#ifndef __VTK_MODULE
    fprintf(stream, "</DataArray>\n");  // footer
#endif

}




void
VTKXMLExportModule :: exportCellVarAs(InternalStateType type, int region,
#ifdef __VTK_MODULE
    vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
    FILE *stream,
#endif
    TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();

    Element *elem;
    FloatMatrix mtrx(3, 3);
    IntegrationRule *iRule;
    GaussPoint *gp;
    FloatArray answer, temp;

    FloatArray valueArray;
    IntArray redIndx;

    InternalStateValueType valType = giveInternalStateValueType(type);
    int ncomponents = giveInternalStateTypeSize(valType);

#ifdef __VTK_MODULE
    vtkSmartPointer<vtkDoubleArray> cellVarsArray = vtkSmartPointer<vtkDoubleArray>::New();
    cellVarsArray->SetName(__InternalStateTypeToString(type));
    cellVarsArray->SetNumberOfComponents(ncomponents);
    cellVarsArray->SetNumberOfTuples(nelem);
#else
    fprintf( stream, " <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> ", __InternalStateTypeToString(type), ncomponents );    
#endif

    valueArray.resize(ncomponents);
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        elem = d->giveElement(ielem);
        if ( (( region > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != region ))
                || this->isElementComposite(elem) || !elem-> isActivated(tStep) ) { // composite cells exported individually
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }
#endif

        switch ( type ) {

        // Special scalars
        case IST_MaterialNumber:
            valueArray.at(1) = (double) elem->giveMaterial()->giveNumber();
            break;
        case IST_ElementNumber:
            valueArray.at(1) = (double) elem->giveNumber();
            break;
        case IST_Pressure: //@todo This case  seems redundant, remove? /JB, /// Why this special treatment for pressure? / Mikael
            if (elem->giveNumberOfInternalDofManagers() == 1) {
                //IntArray pmask(1); pmask.at(1) = P_f;
                //elem->giveInternalDofManager(1)->giveUnknownVector (answer, pmask,EID_ConservationEquation, VM_Total, tStep);
                //valueArray.at(1) = answer.at(1); 
            }
            break;

        // Special vectors 
        case IST_MaterialOrientation_x:
        case IST_MaterialOrientation_y:
        case IST_MaterialOrientation_z:
            int pos;
            if ( type == IST_MaterialOrientation_x ) {
                pos = 1;
            }
            if ( type == IST_MaterialOrientation_y ) {
                pos = 2;
            }
            if ( type == IST_MaterialOrientation_z ) {
                pos = 3;
            }
            if ( !d->giveElement(ielem)->giveLocalCoordinateSystem(mtrx) ) { //@todo shouldn't it be an identity matrix?
                mtrx.resize(3, 3);
                mtrx.zero();
            }
            valueArray.setValues( 3, mtrx.at(1, pos), mtrx.at(2, pos), mtrx.at(3, pos) );

            break;
    
            // Export cell data as average from gp's
            default:

            // compute cell average from gp values (volume/area/length average)
            #if 1
            double gptot = 0;
            answer.resize(0);
            iRule = elem->giveDefaultIntegrationRulePtr();
            if (iRule) {
                MaterialMode mmode = _Unknown;
                for (int i = 0; i < iRule->getNumberOfIntegrationPoints(); ++i) {
                    gp = iRule->getIntegrationPoint(i);
                    mmode = gp->giveMaterialMode();
                    elem->giveIPValue(temp, gp, type, tStep);
                    gptot += gp->giveWeight();
                    answer.add(gp->giveWeight(), temp);
                }
                answer.times(1./gptot);
                elem->giveMaterial()->giveIntVarCompFullIndx(redIndx, type, mmode);
            }
            #endif
            
            // Reshape the Voigt vectors to include all components (duplicated if necessary, VTK insists on 9 components for tensors.)
            #if 1
            if ( ncomponents == 9 && answer.giveSize() != 9) { // If it has 9 components, then it is assumed to be proper already.
                FloatArray tmp = answer;
                this->makeFullForm(valueArray, tmp, valType, redIndx);
            } else if ( valType == ISVT_VECTOR && answer.giveSize() < 3) {
                valueArray.setValues(3,
                                 answer.giveSize() > 1 ? answer.at(1) : 0.0,
                                 answer.giveSize() > 2 ? answer.at(2) : 0.0,
                                 0.0);
            } else if ( ncomponents != answer.giveSize() ) { // Trying to gracefully handle bad cases, just output zeros.
                valueArray.resize(9);
            }
            #endif
        }


        // Write the data to file
#ifdef __VTK_MODULE
        switch ( ncomponents ) {
        case 1:
            cellVarsArray->SetTuple1(ielem-1, valueArray.at(1) ); // Should be integer..
            stream->GetCellData()->SetActiveScalars(__InternalStateTypeToString(type));
            stream->GetCellData()->SetScalars(cellVarsArray);
            break;
        case 3
            cellVarsArray->SetTuple3(ielem-1,  valueArray.at(1), valueArray.at(2), valueArray.at(3) );
            stream->GetCellData()->SetActiveVectors(__InternalStateTypeToString(type));
            stream->GetCellData()->SetVectors(cellVarsArray);
            break;
        case 9:
            for (int i = 1; i <= ncomponents; ++i) {
                cellVarsArray->SetComponent(ielem-1, i-1, answer.at(i));
            }
            stream->GetCellData()->SetActiveTensors(__InternalStateTypeToString(type));
            stream->GetCellData()->SetTensors(cellVarsArray);
            break;
        }
#else
        for ( int i = 1; i <= ncomponents; i++ ) {
            fprintf( stream, "%e ", valueArray.at(i) );
        }
#endif

    }

#ifndef __VTK_MODULE
    fprintf(stream, "</DataArray>\n");  // footer
#endif
}



void VTKXMLExportModule::writeVTKCollection()
{
    struct tm *current;
    time_t now;
    time(&now);
    current = localtime(&now);
    char buff[1024];

    std::string fname = this->emodel->giveOutputBaseFileName() + ".pvd";

    std::ofstream outfile(fname.c_str());

    sprintf(buff, "<!-- Computation started %d-%02d-%02d at %02d:%02d:%02d -->\n", current->tm_year+1900, current->tm_mon+1, current->tm_mday, current->tm_hour,  current->tm_min,  current->tm_sec);
//     outfile << buff;

    outfile << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for (std::list< std::string >::iterator it = this->pvdBuffer.begin(); it != this->pvdBuffer.end(); ++it) {
        outfile << *it << "\n";
    }
    outfile << "</Collection>\n</VTKFile>";

    outfile.close();
}






// Export composite elements
void VTKXMLExportModule :: exportCompositeElement(FILE *stream, Element *el, IntArray &primaryVarsToExport,  IntArray &internalVarsToExport, TimeStep *tStep)
{

    VTKXMLExportModuleElementInterface *interface =
    static_cast< VTKXMLExportModuleElementInterface * >( el->giveInterface(VTKXMLExportModuleElementInterfaceType) );
    if ( interface ) {
        
        interface->giveCompositeExportData(this->compositeCell, primaryVarsToExport, internalVarsToExport, tStep);

        int numSubEl = this->compositeCell.numSubEl; 
        int numNodes = this->compositeCell.numTotalNodes; 
    
        fprintf(stream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numNodes, numSubEl);

        // Export nodes in region as vtk vertices
        fprintf(stream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ");

    
        // Export (fictious) node coords
        for ( int cell = 0; cell < numSubEl; cell++ ) {
            Cell &el = this->compositeCell.elements[cell];
            for ( int inode = 1; inode <= el.nodeCoords.size(); inode++ ) {
            
                FloatArray &coords = el.nodeCoords[inode-1]; 
                for ( int i = 1; i <= coords.giveSize(); i++ ) {
                    fprintf( stream, "%e ", coords.at(i) );
                }
                if ( coords.giveSize() < 3 ) { // fix for 1d
                    fprintf( stream, "%e ", 0.0 );
                }
            
            }
        }
        fprintf(stream, "</DataArray>\n</Points>\n");
    
        // output all cells of the piece
        fprintf(stream, "<Cells>\n");
        // output the connectivity data
        fprintf(stream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ");
    
        for ( int cell = 0; cell < numSubEl; cell++ ) {
            Cell &el = this->compositeCell.elements[cell];
            for ( int i = 1; i <= el.connectivity.giveSize(); i++ ) {
                fprintf(stream, "%d ", el.connectivity.at(i) );
            }
            fprintf(stream, " ");
        }


        // Output the offsets (index of individual element data in connectivity array)
        fprintf(stream, "</DataArray>\n");
        fprintf(stream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ");
        for ( int cell = 0; cell < numSubEl; cell++ ) {
            Cell &el = this->compositeCell.elements[cell];
            fprintf(stream, "%d ", el.offset);
        }
        fprintf(stream, "</DataArray>\n");


        // Output cell types
        fprintf(stream, " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ");
        for ( int cell = 0; cell < numSubEl; cell++ ) {
            Cell &el = this->compositeCell.elements[cell];
            fprintf(stream, "%d ", el.cellType);
        }
        fprintf(stream, "</DataArray>\n");
        fprintf(stream, "</Cells>\n");



        // Export primary and internal variables
        int nodeVarNum = 0;
        
        this->exportPointDataHeader(stream, tStep);
        for (int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
            UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
            FloatArray val;
            int varSize = this->compositeCell.elements[0].nodeVars[nodeVarNum][0].giveSize();  // assumes they all have the same size
            val.resize(varSize);
            fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%i\" format=\"ascii\"> ", __UnknownTypeToString(type), varSize );
            for ( int cell = 0; cell < numSubEl; cell++ ) {
                Cell &el = this->compositeCell.elements[cell];
                for ( int j = 1; j <= el.connectivity.giveSize(); j++ ) {
                
                    val = el.nodeVars[nodeVarNum][j-1];
                    for ( int component = 1; component <= val.giveSize(); component++ ) {
                        fprintf( stream, "%e ", val.at(component) );
                    }
                }    
            }
            fprintf(stream, "</DataArray>\n");
            nodeVarNum++;
        }

        for (int i = 1; i <= internalVarsToExport.giveSize(); i++ ) {
            InternalStateType type = ( InternalStateType ) internalVarsToExport.at(i);
            exportNodalVarAs(type, nodeVarNum, stream, tStep);
            nodeVarNum++;
        }
        fprintf(stream, "</PointData>\n");



        // Export cell data
        fprintf(stream, "<CellData Scalars=\"\" Vectors=\"\" Tensors=\"\">\n");  // should contain a list of InternalStateType

        for (int i = 1; i <= internalVarsToExport.giveSize(); i++ ) {
            InternalStateType type = ( InternalStateType ) internalVarsToExport.at(i);
            InternalStateValueType valType = giveInternalStateValueType(type);
            int varSize = giveInternalStateTypeSize(valType);

            fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%i\" format=\"ascii\">", __InternalStateTypeToString(type), varSize );
            for ( int cell = 0; cell < numSubEl; cell++ ) {
                FloatArray &var = this->compositeCell.elements[cell].elVars[i-1];
            
                for ( int component = 1; component <= var.giveSize(); component++ ) {
                    fprintf( stream, "%e ", var.at(component) );
                }
            }

            fprintf(stream, "</DataArray>\n");
        
        }
        fprintf(stream, "</CellData>\n");

    
        // end of piece record
        fprintf(stream, "</Piece>\n");
    


    }
}

void 
VTKXMLExportModule :: exportNodalVarAs(InternalStateType type, int nodeVarNum, FILE *stream, TimeStep *tStep)
{
    FloatArray val;
    int varSize = this->compositeCell.elements[0].nodeVars[nodeVarNum][0].giveSize();  // assumes they all have the same size
    val.resize(varSize);
    fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%i\" format=\"ascii\"> ", __InternalStateTypeToString(type), varSize );
    for ( int cell = 0; cell < this->compositeCell.elements.size(); cell++ ) {
        Cell &el = this->compositeCell.elements[cell];
        for ( int j = 1; j <= el.connectivity.giveSize(); j++ ) {
                
            val = el.nodeVars[nodeVarNum][j-1];
            for ( int component = 1; component <= val.giveSize(); component++ ) {
                fprintf( stream, "%e ", val.at(component) );
            }
        }    
    }
    fprintf(stream, "</DataArray>\n");
    nodeVarNum++;
       
}



void
VTKXMLExportModule :: exportCellVarAs(InternalStateType type, std::vector<FloatArray> &cellVars, FILE *stream, TimeStep *tStep)
{
    InternalStateValueType valType = giveInternalStateValueType(type);
    int ncomponents = giveInternalStateTypeSize(valType);
    
    fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%i\" format=\"ascii\">\n", __InternalStateTypeToString(type), ncomponents );

    for (int i = 1; i <= cellVars.at(1).giveSize(); i++) {
        for ( int cell = 0; cell < cellVars.size(); cell++ ) {
            fprintf( stream, "%e ", cellVars.at(cell).at(i) );
        }
    }
    fprintf( stream, "\n" );
    fprintf(stream, "</DataArray>\n");

}







} // end namespace oofem
