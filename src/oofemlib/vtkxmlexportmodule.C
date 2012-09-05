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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "gausspnt.h"
#include "timestep.h"
#include "engngm.h"
#include "node.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "cltypes.h"
#include "material.h"
#include "usrdefsub.h"

#ifndef __MAKEDEPEND
 #include <string>
 #include <sstream>
 #include <fstream>
#endif

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
VTKXMLExportModule :: VTKXMLExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{
    primVarSmoother = NULL;
    smoother = NULL;
    redToFull.setValues(6, 1, 5, 9, 6, 3, 2);//position of yz, xx, yy, zz, yz, xz, xy in tensor 
}


VTKXMLExportModule :: ~VTKXMLExportModule()
{
    if ( this->smoother ) {
        delete this->smoother;
    }
}


IRResultType
VTKXMLExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    ExportModule :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, IFT_VTKXMLExportModule_cellvars, "cellvars"); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, IFT_VTKXMLExportModule_vars, "vars"); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, IFT_VTKXMLExportModule_primvars, "primvars"); // Macro - see unknowntype.h

    val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_VTKXMLExportModule_stype, "stype"); // Macro
    stype = ( NodalRecoveryModel::NodalRecoveryModelType ) val;

    regionsToSkip.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, regionsToSkip, IFT_VTKXMLExportModule_regionstoskip, "regionstoskip"); // Macro

    this->nvr = 0; // number of virtual regions
    IR_GIVE_OPTIONAL_FIELD(ir, nvr, IFT_VTKXMLExportModule_nvr, "nvr"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, vrmap, IFT_VTKXMLExportModule_vrmap, "vrmap"); // Macro

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
    if (type == ISVT_TENSOR_S3) {
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
    } else if ( elemGT == EGT_hexa_1 ) {
        vtkCellType = 12;
    } else if ( elemGT == EGT_hexa_2 ) {
        vtkCellType = 25;
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

        default:
            OOFEM_ERROR("VTKXMLExportModule: unsupported cell type ID");
    }

    return 0; // to make compiler happy
}


void
VTKXMLExportModule :: giveElementCell(IntArray &answer, Element *elem, int cell)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int i, nelemNodes;

    if ( ( elemGT == EGT_point ) ||
        ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
        ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
        ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
        ( elemGT == EGT_hexa_1 ) ) {
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber() ;
        }
    } else if ( elemGT == EGT_hexa_2 ) {
        int HexaQuadNodeMapping [] = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(HexaQuadNodeMapping [ i - 1 ])->giveNumber() ;
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
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
        ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_hexa_2 ) ) {
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
    fprintf(stream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(stream, "<UnstructuredGrid>\n");
#endif


    Domain *d  = emodel->giveDomain(1);
    Element *elem;
    FloatArray *coords;
    int i, inode;
    int ielem, nelem = d->giveNumberOfElements();

    this->giveSmoother(); // make sure smoother is created

    // output nodes Region By Region
    int ireg, nregions = this->smoother->giveNumberOfVirtualRegions();
    int regionDofMans, totalcells;
    IntArray mapG2L, mapL2G;

    /* loop over regions */
    for ( ireg = 1; ireg <= nregions; ireg++ ) {
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
            for ( inode = 1; inode <= regionDofMans; inode++ ) {
                coords = d->giveNode(mapL2G.at(inode))->giveCoordinates();
                int dims = coords->giveSize();
                nodes->InsertNextPoint(coords->at(1), dims >= 2 ? coords->at(2) : 0.0, dims >= 3 ? coords->at(3) : 0.0);
            }
            stream->SetPoints(nodes);
#else
            fprintf(stream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", regionDofMans, totalcells);
            // export nodes in region as vtk vertices
            fprintf(stream, "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ");
            for ( inode = 1; inode <= regionDofMans; inode++ ) {
                coords = d->giveNode( mapL2G.at(inode) )->giveCoordinates();
                for ( i = 1; i <= coords->giveSize(); i++ ) {
                    fprintf( stream, "%e ", coords->at(i) );
                }

                for ( i = coords->giveSize() + 1; i <= 3; i++ ) {
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
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
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


                nelemNodes = elem->giveNumberOfNodes();
                this->giveElementCell(cellNodes, elem, 0);
#ifdef __VTK_MODULE
                elemNodeArray->Reset();
                elemNodeArray->SetNumberOfIds(nelemNodes);
#endif
                for ( i = 1; i <= nelemNodes; i++ ) {
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
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
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
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
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
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
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
                ///@todo Not sure how to deal with this.
                // multi cell (composite) elements should support vtkxmlexportmoduleinterface
                // and are exported as individual pieces (see VTKXMLExportModuleElementInterface)
                VTKXMLExportModuleElementInterface *interface =
                    ( VTKXMLExportModuleElementInterface * ) elem->giveInterface(VTKXMLExportModuleElementInterfaceType);
                if ( interface ) {
                    // passing this to access general piece related methods like exportPointDataHeader, etc.
                    interface->_export(stream, this, primaryVarsToExport, internalVarsToExport, tStep);
                }
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
            this->pvdBuffer.pushBack(pvdEntry.str());
        }
        this->writeVTKCollection();
    } else
#endif
    if ( !emodel->isParallel() && tStep->giveNumber() > 1 ) { // For non-parallel enabled OOFEM, then we only check for multiple steps.
        std::ostringstream pvdEntry;
        pvdEntry << "<DataSet timestep=\"" << tStep->giveIntrinsicTime() << "\" group=\"\" part=\"\" file=\"" << fname << "\"/>";
        this->pvdBuffer.pushBack(pvdEntry.str());
        this->writeVTKCollection();
    }
}

#ifndef __VTK_MODULE
void
VTKXMLExportModule :: exportPointDataHeader(FILE *stream, TimeStep *tStep)
{
    int i, n;
    std :: string scalars, vectors, tensors;

    n = primaryVarsToExport.giveSize();

    UnknownType type;

    for ( i = 1; i <= n; i++ ) {
        type = ( UnknownType ) primaryVarsToExport.at(i);
        if ( ( type == DisplacementVector ) || ( type == EigenVector ) || ( type == VelocityVector ) ) {
            vectors += __UnknownTypeToString(type);
            vectors.append(" ");
        } else if ( ( type == FluxVector ) || ( type == PressureVector ) || ( type == TemperatureVector ) ) {
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
    for ( i = 1; i <= n; i++ ) {
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
    int i, n = internalVarsToExport.giveSize();
    InternalStateType isttype;
    InternalStateValueType vtype;

    this->giveSmoother()->init(); // Makes sure smoother is up-to-date with potentially new mesh.
    // should be performed over regions
    for ( i = 1; i <= n; i++ ) {
        isttype = ( InternalStateType ) internalVarsToExport.at(i);
        vtype = giveInternalStateValueType(isttype);
        this->exportIntVarAs(isttype, vtype, mapG2L, mapL2G, regionDofMans, region, stream, tStep);
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


    int ielem, nelem = domain->giveNumberOfElements();
    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    regionG2LNodalNumbers.resize(nnodes);
    regionG2LNodalNumbers.zero();
    regionDofMans = 0;
    regionSingleCells = 0;

    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
        if ( ( reg > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != reg ) ) {
            continue;
        }

        if ( this->isElementComposite(element) ) {
            continue;                                    // composite cells exported individually
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

    int i;
    for ( i = 1; i <= nnodes; i++ ) {
        if ( regionG2LNodalNumbers.at(i) ) {
            regionG2LNodalNumbers.at(i) = currOffset++;
            regionL2GNodalNumbers.at( regionG2LNodalNumbers.at(i) ) = i;
        }
    }

    return 1;
}


//keyword "vars" in OOFEM input file
void
VTKXMLExportModule :: exportIntVarAs(InternalStateType valID, InternalStateValueType type,
                                     IntArray &mapG2L, IntArray &mapL2G, int regionDofMans, int ireg,
#ifdef __VTK_MODULE
                                     vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
                                     FILE *stream,
#endif
                                     TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int inode;
    int j, jsize, valSize;
    FloatArray iVal(3), t(9);
    const FloatArray *val;
    IntArray regionVarMap;

    this->giveSmoother();

#ifdef __VTK_MODULE
    vtkSmartPointer<vtkDoubleArray> intVarArray = vtkSmartPointer<vtkDoubleArray>::New();
    intVarArray->SetName(__InternalStateTypeToString(valID));
#endif

    if ( type == ISVT_SCALAR ) {
#ifdef __VTK_MODULE
        intVarArray->SetNumberOfComponents(1);
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"> ", __InternalStateTypeToString(valID) );
#endif
    } else if ( type == ISVT_VECTOR ) {
#ifdef __VTK_MODULE
        intVarArray->SetNumberOfComponents(3);
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> ",
                __InternalStateTypeToString(valID) );
#endif
    } else if ( ( type == ISVT_TENSOR_S3 ) || ( type == ISVT_TENSOR_S3E ) ) {
#ifdef __VTK_MODULE
        intVarArray->SetNumberOfComponents(9);
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"9\" format=\"ascii\"> ",
                __InternalStateTypeToString(valID) );
#endif
    } else if ( type == ISVT_TENSOR_G ) {
#ifdef __VTK_MODULE
        intVarArray->SetNumberOfComponents(9);
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"9\" format=\"ascii\"> ",
                __InternalStateTypeToString(valID) );
#endif
    } else {
        fprintf( stderr, "VTKXMLExportModule::exportIntVarAs: unsupported variable type %s\n", __InternalStateTypeToString(valID) );
    }
#ifdef __VTK_MODULE
    intVarArray->SetNumberOfTuples(regionDofMans);
#endif


    if ( !( ( valID == IST_DisplacementVector ) || ( valID == IST_MaterialInterfaceVal ) ) ) {
        this->smoother->recoverValues(valID, tStep);
        this->smoother->giveRegionRecordMap(regionVarMap, ireg, valID);
    }

    for ( inode = 1; inode <= regionDofMans; inode++ ) {
        if ( valID == IST_DisplacementVector ) {
            iVal.resize(3);
            val = & iVal;
            for ( j = 1; j <= 3; j++ ) {
                iVal.at(j) = d->giveNode( mapL2G.at(inode) )->giveUpdatedCoordinate(j, tStep, EID_MomentumBalance, 1.0) -
                             d->giveNode( mapL2G.at(inode) )->giveCoordinate(j);
            }
        } else if ( valID == IST_MaterialInterfaceVal ) {
            MaterialInterface *mi = emodel->giveMaterialInterface(1);
            if ( mi ) {
                iVal.resize(1);
                val = & iVal;
                iVal.at(1) = mi->giveNodalScalarRepresentation( mapL2G.at(inode) );
            }
        } else {
            int found = this->smoother->giveNodalVector(val, mapL2G.at(inode), ireg);
            if ( !found ) {
                iVal.resize( regionVarMap.giveSize() );
                iVal.zero();
                val = & iVal;
                //OOFEM_WARNING2("VTKXMLExportModule::exportIntVars: smoothing error: invalid data in node %d", inode);
            }
        }

        valSize = val->giveSize();
        if ( type == ISVT_SCALAR ) {
#ifdef __VTK_MODULE
            intVarArray->SetTuple1(inode-1, valSize ? val->at(1) : 0.0);
#else
            fprintf( stream, "%e ", valSize ? val->at(1) : 0.0 );
#endif
        } else if ( type == ISVT_VECTOR ) {
            jsize = min( 3, valSize );
            for ( j = 1; j <= jsize; j++ ) {
#ifdef __VTK_MODULE
                intVarArray->SetComponent(inode-1, j-1, val->at(j));
#else
                fprintf( stream, "%e ", val->at(j) );
#endif
            }

            for ( j = jsize + 1; j <= 3; j++ ) {
#ifdef __VTK_MODULE
                intVarArray->SetComponent(inode-1, j-1, 0.0);
#else
                fprintf(stream, "0.0 ");
#endif
            }
#ifndef __VTK_MODULE
            fprintf(stream, " ");
#endif

        } else if ( type == ISVT_TENSOR_S3 || type == ISVT_TENSOR_S3E ) {
            this->makeFullForm(t, *val, type, regionVarMap);

            for ( j = 1; j <= 9; j++ ) {
#ifdef __VTK_MODULE
                intVarArray->SetComponent(inode-1, j-1, t.at(j));
#else
                fprintf( stream, "%e ", t.at(j) );
#endif
            }
#ifndef __VTK_MODULE
            fprintf(stream, " ");
#endif
        } else if ( type == ISVT_TENSOR_G ) { // export general tensor values as scalars
            jsize = min( 9, val->giveSize() );
            for ( j = 1; j <= jsize; j++ ) {
#ifdef __VTK_MODULE
                intVarArray->SetComponent(inode-1, j-1, val->at(j));
#else
                fprintf( stream, "%e ", val->at(j) );
#endif

            }

            for ( j = jsize + 1; j <= 9; j++ ) {
#ifdef __VTK_MODULE
                intVarArray->SetComponent(inode-1, j-1, 0.0);
#else
                fprintf(stream, "0.0 ");
#endif
            }
#ifndef __VTK_MODULE
            fprintf(stream, " ");
#endif
        }
    } // end loop over dofmans

#ifdef __VTK_MODULE
    if (type == ISVT_SCALAR) {
        stream->GetPointData()->SetActiveScalars(__InternalStateTypeToString(valID));
        stream->GetPointData()->SetScalars(intVarArray);
    } else if (type == ISVT_VECTOR) {
        stream->GetPointData()->SetActiveVectors(__InternalStateTypeToString(valID));
        stream->GetPointData()->SetVectors(intVarArray);
    } else {
        stream->GetPointData()->SetActiveTensors(__InternalStateTypeToString(valID));
        stream->GetPointData()->SetTensors(intVarArray);
    }
#else
    fprintf(stream, "</DataArray>\n");
#endif
}


NodalRecoveryModel *
VTKXMLExportModule :: giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->smoother == NULL ) {
      this->smoother = CreateUsrDefNodalRecoveryModel(this->stype, d);
      this->smoother->setRecoveryMode (nvr, vrmap);
    }
    return this->smoother;
}


NodalRecoveryModel *
VTKXMLExportModule :: givePrimVarSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->primVarSmoother == NULL ) {
        this->primVarSmoother = CreateUsrDefNodalRecoveryModel(NodalRecoveryModel::NRM_NodalAveraging, d);
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
VTKXMLExportModule :: exportPrimVarAs(UnknownType valID, IntArray &mapG2L, IntArray &mapL2G,
                                      int regionDofMans, int ireg,
#ifdef __VTK_MODULE
                                      vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
                                      FILE *stream,
#endif
                                      TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    InternalStateValueType type = ISVT_UNDEFINED;

    if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
        type = ISVT_VECTOR;
    } else if ( ( valID == FluxVector ) || ( valID == PressureVector ) || ( valID == TemperatureVector ) ) {
        type = ISVT_SCALAR;
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule::exportPrimVarAs: unsupported UnknownType %s", __UnknownTypeToString(valID) );
    }

#ifdef __VTK_MODULE
    vtkSmartPointer<vtkDoubleArray> primVarArray = vtkSmartPointer<vtkDoubleArray>::New();
    primVarArray->SetName(__UnknownTypeToString(valID));

    if ( type == ISVT_SCALAR ) {
        primVarArray->SetNumberOfComponents(1);
    } else if ( type == ISVT_VECTOR ) {
        primVarArray->SetNumberOfComponents(3);
    } else {
        fprintf( stderr, "VTKXMLExportModule::exportPrimVarAs: unsupported variable type %s\n", __UnknownTypeToString(valID) );
    }
    primVarArray->SetNumberOfTuples(regionDofMans);
#else
    if ( type == ISVT_SCALAR ) {
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"> ", __UnknownTypeToString(valID) );
    } else if ( type == ISVT_VECTOR ) {
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> ", __UnknownTypeToString(valID) );
    } else {
        fprintf( stderr, "VTKXMLExportModule::exportPrimVarAs: unsupported variable type %s\n", __UnknownTypeToString(valID) );
    }
#endif

    DofManager *dman;
    FloatArray iVal;

    for ( int inode = 1; inode <= regionDofMans; inode++ ) {
        dman = d->giveNode( mapL2G.at(inode) );

        this->getPrimaryVariable(iVal, dman, tStep, valID, ireg);

        if ( type == ISVT_SCALAR ) {
#ifdef __VTK_MODULE
            primVarArray->SetTuple1(inode-1, iVal.at(1));
#else
            fprintf(stream, "%e ", iVal.at(1) );
#endif
        } else if ( type == ISVT_VECTOR ) {
#ifdef __VTK_MODULE
            primVarArray->SetTuple3(inode-1, iVal.at(1), iVal.at(2), iVal.at(3));
#else
            fprintf( stream, "%e %e %e ", iVal.at(1), iVal.at(2), iVal.at(3) );
#endif
        }
    } // end loop over nodes
#ifdef __VTK_MODULE
    if ( type == ISVT_SCALAR ) {
        stream->GetPointData()->SetActiveScalars(__UnknownTypeToString(valID));
        stream->GetPointData()->SetScalars(primVarArray);
    } else if ( type == ISVT_VECTOR ) {
        stream->GetPointData()->SetActiveVectors(__UnknownTypeToString(valID));
        stream->GetPointData()->SetVectors(primVarArray);
    }
#else
    fprintf(stream, "</DataArray>\n");
#endif
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
        for (int j = 1; j <= dman->giveNumberOfDofs(); j++ ) {
            id = dman->giveDof(j)->giveDofID();
            if ( ( id == V_u ) || ( id == D_u ) || ( id == V_v ) || ( id == D_v ) || ( id == V_w ) || ( id == D_w ) ) {
                dofIDMask.followedBy(id);
            }
        }
        answer.resize(3);
    } else if ( type == FluxVector ) {
        dofIDMask.followedBy(C_1);
        eid = EID_ConservationEquation;
        iState = IST_MassConcentration_1;
        answer.resize(1);
    } else if ( type == TemperatureVector ) {
        dofIDMask.followedBy(T_f);
        eid = EID_ConservationEquation;
        iState = IST_Temperature;
        answer.resize(1);
    } else if ( type == PressureVector ) {
        dofIDMask.followedBy(P_f);
        eid = EID_ConservationEquation;
        iState = IST_Pressure;
        answer.resize(1);
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported unknownType %s", __UnknownTypeToString(type) );
    }

    size = dofIDMask.giveSize();
    answer.zero();

    for (int j = 1; j <= size; j++ ) {
        if ( ( indx = dman->findDofWithDofId( (DofIDItem)dofIDMask.at(j) ) ) ) {
            // primary variable available directly in DOF-manager
            answer.at(j) = dman->giveDof(indx)->giveUnknown(eid, VM_Total, tStep);
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
    int i, n = cellVarsToExport.giveSize();
    InternalStateType type;

    if ( n == 0 ) {
        return;
    }

#ifndef __VTK_MODULE
    //print header
    fprintf(stream, "<CellData Scalars=\"\" Vectors=\"\" Tensors=\"\">\n");  // should contain a list of InternalStateType
#endif

    for ( i = 1; i <= n; i++ ) {
        type = ( InternalStateType ) cellVarsToExport.at(i);
        this->exportCellVarAs(type, region, stream, tStep);
    }

#ifndef __VTK_MODULE
    //print footer
    fprintf(stream, "</CellData>\n");
#endif
}


//keyword "cellvars" in OOFEM input file
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
    int ielem, nelem = d->giveNumberOfElements();
    int pos;
    Element *elem;
    FloatMatrix mtrx(3, 3);
    IntegrationRule *iRule;
    GaussPoint *gp;
    FloatArray answer, temp;
    double gptot;
    int ncomponents = 1;

#ifdef __VTK_MODULE
    vtkSmartPointer<vtkDoubleArray> cellVarsArray = vtkSmartPointer<vtkDoubleArray>::New();
    cellVarsArray->SetName(__InternalStateTypeToString(type));
#endif

    switch ( type ) {
    case IST_MaterialNumber:
    case IST_ElementNumber:
    case IST_Pressure:
        // if it wasn't for IST_Pressure,
#ifdef __VTK_MODULE
        cellVarsArray->SetNumberOfComponents(1);
        cellVarsArray->SetNumberOfTuples(nelem);
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", __InternalStateTypeToString(type) );
#endif
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);

            if ( (( region > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != region ))
                    || this->isElementComposite(elem) ) { // composite cells exported individually
                continue;
            }

#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            if ( type == IST_MaterialNumber ) {
#ifdef __VTK_MODULE
                cellVarsArray->SetTuple1(ielem-1, elem->giveMaterial()->giveNumber() ); // Should be integer..
#else
                fprintf( stream, "%d ", elem->giveMaterial()->giveNumber() );
#endif
            } else if ( type == IST_ElementNumber ) {
#ifdef __VTK_MODULE
                cellVarsArray->SetTuple1(ielem-1,  elem->giveNumber() ); // Should be integer..
#else
                fprintf( stream, "%d ", elem->giveNumber() );
#endif
            } else if (type == IST_Pressure) { ///@todo Why this special treatment for pressure? / Mikael
                if (elem->giveNumberOfInternalDofManagers() == 1) {
                    IntArray pmask(1); pmask.at(1) = P_f;
                    elem->giveInternalDofManager(1)->giveUnknownVector (answer, pmask,EID_ConservationEquation, VM_Total, tStep);
#ifdef __VTK_MODULE
                    cellVarsArray->SetTuple1(ielem-1,  answer.at(1) ); // Should be integer..
#else
                    fprintf( stream, "%f ", answer.at(1) );
#endif
                }
            }
        }
#ifdef __VTK_MODULE
        stream->GetCellData()->SetActiveScalars(__InternalStateTypeToString(type));
        stream->GetCellData()->SetScalars(cellVarsArray);
#endif
        break;

    case IST_MaterialOrientation_x:
    case IST_MaterialOrientation_y:
    case IST_MaterialOrientation_z:
#ifdef __VTK_MODULE
        cellVarsArray->SetNumberOfComponents(3);
        cellVarsArray->SetNumberOfTuples(nelem);
        ncomponents = 3;
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", __InternalStateTypeToString(type) );
#endif
        if ( type == IST_MaterialOrientation_x ) {
            pos = 1;
        }

        if ( type == IST_MaterialOrientation_y ) {
            pos = 2;
        }

        if ( type == IST_MaterialOrientation_z ) {
            pos = 3;
        }

        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            ///@todo Should no elements be skipped here? / Mikael
            if ( !d->giveElement(ielem)->giveLocalCoordinateSystem(mtrx) ) {
                mtrx.resize(3, 3);
                mtrx.zero();
            }
#ifdef __VTK_MODULE
            cellVarsArray->SetTuple3(ielem-1,  mtrx.at(1, pos), mtrx.at(2,pos), mtrx.at(3,pos) );
#else
            fprintf( stream, "%f %f %f  ", mtrx.at(1, pos), mtrx.at(2, pos), mtrx.at(3, pos) );
#endif
        }

#ifdef __VTK_MODULE
        stream->GetCellData()->SetActiveVectors(__InternalStateTypeToString(type));
        stream->GetCellData()->SetVectors(cellVarsArray);
#endif
        break;

    default:
        bool reshape = false;
        InternalStateValueType vt = giveInternalStateValueType(type);
        if ( vt == ISVT_SCALAR ) {
            ncomponents = 1;
        } else if ( vt == ISVT_VECTOR ) {
            ncomponents = 3;
        } else {
            ncomponents = 9;
            reshape = true;
        }
#ifdef __VTK_MODULE
        cellVarsArray->SetNumberOfComponents(ncomponents);
        cellVarsArray->SetNumberOfTuples(nelem);
#else
        fprintf( stream, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", __InternalStateTypeToString(type), ncomponents );
#endif

        IntArray redIndx;
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
            if ( (( region > 0 ) && ( this->smoother->giveElementVirtualRegionNumber(ielem) != region ))
                    || this->isElementComposite(elem) ) { // composite cells exported individually
                continue;
            }
#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }
#endif
            gptot = 0;
            answer.resize(0);
            iRule = elem->giveDefaultIntegrationRulePtr();
            if (iRule) {
                for (int i = 0; i < iRule->getNumberOfIntegrationPoints(); ++i) {
                    gp = iRule->getIntegrationPoint(i);
                    elem->giveIPValue(temp, gp, type, tStep);
                    gptot += gp->giveWeight();
                    answer.add(gp->giveWeight(), temp);
                }
                answer.times(1/gptot);
                elem->giveMaterial()->giveIntVarCompFullIndx(redIndx, type, gp->giveMaterialMode());
            }
            // Reshape the Voigt vectors to include all components (duplicated if necessary, VTK insists on 9 components for tensors.)
            if ( reshape && answer.giveSize() != 9) { // If it has 9 components, then it is assumed to be proper already.
                FloatArray tmp = answer;
                this->makeFullForm(answer, tmp, vt, redIndx);
            } else if ( vt == ISVT_VECTOR && answer.giveSize() < 3) {
                answer.setValues(3,
                                 answer.giveSize() > 1 ? answer.at(1) : 0.0,
                                 answer.giveSize() > 2 ? answer.at(2) : 0.0,
                                 0.0);
            } else if ( ncomponents != answer.giveSize() ) { // Trying to gracefully handle bad cases, just output zeros.
                answer.resize(ncomponents);
                answer.zero();
            }
            for (int i = 1; i <= ncomponents; ++i) {
#ifdef __VTK_MODULE
                cellVarsArray->SetComponent(ielem-1, i-1, answer.at(i));
#else
                fprintf( stream, "%e ", answer.at(i) );
#endif
            }
#ifndef __VTK_MODULE
            fprintf( stream, "\n" );
#endif
        }
#ifdef __VTK_MODULE
        if (ncomponents == 1) {
            stream->GetCellData()->SetActiveScalars(__InternalStateTypeToString(type));
            stream->GetCellData()->SetScalars(cellVarsArray);
        } else if (ncomponents == 3) {
            stream->GetCellData()->SetActiveVectors(__InternalStateTypeToString(type));
            stream->GetCellData()->SetVectors(cellVarsArray);
        } else {
            stream->GetCellData()->SetActiveTensors(__InternalStateTypeToString(type));
            stream->GetCellData()->SetTensors(cellVarsArray);
        }
#endif
    }
#ifndef __VTK_MODULE
    fprintf(stream, "</DataArray>\n");
#endif
}

void VTKXMLExportModule::writeVTKCollection()
{
    std::string fname = this->emodel->giveOutputBaseFileName() + ".pvd";

    std::ofstream outfile(fname.c_str());

    outfile << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
    for (dynaList< std::string >::iterator it = this->pvdBuffer.begin(); it != this->pvdBuffer.end(); ++it) {
        outfile << *it << "\n";
    }
    outfile << "</Collection>\n</VTKFile>";

    outfile.close();
}


} // end namespace oofem
