/* $Header: /home/cvs/bp/oofem/sm/src/vtkexportmodule.C,v 1.11.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "timestep.h"
#include "engngm.h"
#include "strreader.h"
#include "node.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "oofem_limits.h"
#include "cltypes.h"
#include "material.h"

#ifndef __MAKEDEPEND
 #include <vector>
 #include <string>
#endif

namespace oofem {
VTKXMLExportModule :: VTKXMLExportModule(EngngModel *e) : ExportModule(e), internalVarsToExport(), primaryVarsToExport()
{
    smoother = NULL;
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

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_VTKXMLExportModule_stype, "stype"); // Macro
    stype = ( VTKEM_SmootherType ) val;

    regionsToSkip.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, regionsToSkip, IFT_VTKXMLExportModule_regionstoskip, "regionstoskip"); // Macro

    return IRRT_OK;
}


void
VTKXMLExportModule :: doOutput(TimeStep *tStep)
{
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    FILE *stream = this->giveOutputStream(tStep);

    fprintf(stream, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(stream, "<UnstructuredGrid>\n");


    Domain *d  = emodel->giveDomain(1);
    Element *elem;
    FloatArray *coords;
    int i, inode;
    int ielem, nelem = d->giveNumberOfElements();

    // output nodes Region By Region
    int ireg, nregions = d->giveNumberOfRegions();
    int regionDofMans, totalcells;
    IntArray mapG2L, mapL2G;

    /* loop over regions */
    for ( ireg = 1; ireg <= nregions; ireg++ ) {
        if ( ( ireg > 0 ) && ( this->regionsToSkip.contains(ireg) ) ) {
            continue;
        }

        // asemble local->global and gloabal->local region map
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
            fprintf(stream, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", regionDofMans, totalcells);

            // export nodes in region as vtk vertices
            fprintf(stream, "<Points>\n <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"> ");
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

            // output all cells of the piece
            int vtkCellType, nelemNodes;
            fprintf(stream, "<Cells>\n");
            // output the connectivity data
            fprintf(stream, " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ");
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
                if ( ( ireg > 0 ) && ( elem->giveRegionNumber() != ireg ) ) {
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
                for ( i = 1; i <= nelemNodes; i++ ) {
                    fprintf(stream, "%d ", mapG2L.at( elem->giveNode(i)->giveNumber() ) - 1);
                }

                fprintf(stream, " ");
            }

            fprintf(stream, "</DataArray>\n");
            // output the offsets (index of individual element data in connectivity array)
            fprintf(stream, " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ");
            int offset = 0;
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
                if ( ( ireg > 0 ) && ( elem->giveRegionNumber() != ireg ) ) {
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
                if ( ( ireg > 0 ) && ( elem->giveRegionNumber() != ireg ) ) {
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


            // export primary and internal varialbles
            this->exportPointDataHeader(stream, tStep);
            this->exportPrimaryVars(stream, mapG2L, mapL2G, regionDofMans, ireg, tStep);
            this->exportIntVars(stream, mapG2L, mapL2G, regionDofMans, ireg, tStep);
            this->exportPointDataFooter(stream, tStep);

            //export cell data
            this->exportCellVars(stream, totalcells, ireg, tStep);

            // end of piece record
            fprintf(stream, "</Piece>\n");
        } // end of default piece for simple geometry elements

#if 1
        // loop over region elements with multi-cell geometry
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);

            if ( this->regionsToSkip.contains( elem->giveRegionNumber() ) ) {
                continue;
            }

 #ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

 #endif
            if ( this->isElementComposite(elem) ) {
                // multi cell (composite) elements should support vtkxmlexportmoduleinterface
                // and are exported as individual pieces (see VTKXMLExportModuleElementInterface)
                VTKXMLExportModuleElementInterface *interface =
                    ( VTKXMLExportModuleElementInterface * ) elem->giveInterface(VTKXMLExportModuleElementInterfaceType);
                if ( interface ) {
                    // passing this to access general piece related methods like exportPointDataHeader, etc.
                    interface->_export(stream, this, primaryVarsToExport, internalVarsToExport, tStep);
                }
            }
        } // end loop over multi-cell elements

#endif
    } // end loop over regions

    // finish unstructured grid data and vtk file
    fprintf(stream, "</UnstructuredGrid>\n</VTKFile>");
    fclose(stream);
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


FILE *
VTKXMLExportModule :: giveOutputStream(TimeStep *tStep)
{
    char baseFileName [ MAX_FILENAME_LENGTH ];
    char fileName [ MAX_FILENAME_LENGTH ];
    FILE *answer;

    emodel->giveOutputBaseFileName(baseFileName, MAX_FILENAME_LENGTH);
#ifdef __PARALLEL_MODE
    sprintf( fileName, "%s.%d.%d.vtu", baseFileName, tStep->giveNumber(), emodel->giveRank() );
#else
    sprintf( fileName, "%s.%d.vtu", baseFileName, tStep->giveNumber() );
#endif
    if ( ( answer = fopen(fileName, "w") ) == NULL ) {
        OOFEM_ERROR2("VTKXMLExportModule::giveOutputStream: failed to open file %s", fileName);
    }

    return answer;
}

int
VTKXMLExportModule :: giveCellType(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int vtkCellType = 0;

    if ( elemGT == EGT_triangle_1 ) {
        vtkCellType = 5;
    } else if ( elemGT == EGT_triangle_2 ) {
        vtkCellType = 22;
    } else if ( elemGT == EGT_tetra_1 ) {
        vtkCellType = 10;
    } else if ( elemGT == EGT_quad_1 ) {
        vtkCellType = 9;
    } else if ( elemGT == EGT_quad_2 ) {
        vtkCellType = 23;
    } else if ( elemGT == EGT_hexa_1 ) {
        vtkCellType = 12;
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported element geometry type on element %d", elem->giveNumber() );
    }

    return vtkCellType;
}

int
VTKXMLExportModule :: giveNumberOfNodesPerCell(int cellType)
{
    if ( cellType == 10 ) {
        return 4;
    } else if ( cellType == 5 ) {
        return 3;
    } else if ( cellType == 9 ) {
        return 4;
    } else if ( cellType == 22 ) {
        return 6;
    } else if  ( ( cellType == 12 ) || ( cellType == 23 ) ) {
        return 8;
    } else {
        OOFEM_ERROR2("VTKXMLExportModule: unknown number of nodes of cellType %d\n", cellType);
    }

    return 0; // to make compiler happy
}


void
VTKXMLExportModule :: giveElementCell(IntArray &answer, Element *elem, int cell)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int i, nelemNodes;

    if ( ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) || ( elemGT == EGT_tetra_1 ) ||
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
        ( elemGT == EGT_hexa_1 ) ) {
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber() - 1;
        }
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported element geometry type on element %d", elem->giveNumber() );
    }

    return;
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

    if ( ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) || ( elemGT == EGT_tetra_1 ) ||
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) || ( elemGT == EGT_hexa_1 ) ) {
        return 1;
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule: unsupported element geometry type on element %d", elem->giveNumber() );
    }

    return 0;
}

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

void
VTKXMLExportModule :: exportPointDataFooter(FILE *stream, TimeStep *tStep)
{
    fprintf(stream, "</PointData>\n");
}

//keyword "vars" in OOFEM input file
void
VTKXMLExportModule :: exportIntVars(FILE *stream, IntArray &mapG2L, IntArray &mapL2G,
                                    int regionDofMans, int region, TimeStep *tStep)
{
    int i, n = internalVarsToExport.giveSize();
    InternalStateType isttype;
    InternalStateValueType vtype;
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
        if ( ( reg > 0 ) && ( element->giveRegionNumber() != reg ) ) {
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
                /* mark for assignement. This is done later, as it allows to preserve
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
                                     IntArray &mapG2L, IntArray &mapL2G,
                                     int regionDofMans, int ireg, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    //int ireg, nregions = d->giveNumberOfRegions();
    int inode;
    int j, jsize;
    FloatArray iVal(3);
    FloatMatrix t(3, 3);
    const FloatArray *val;
    IntArray regionVarMap;

    this->giveSmoother();

    if ( type == ISVT_SCALAR ) {
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\"> ", __InternalStateTypeToString(valID) );
    } else if ( type == ISVT_VECTOR ) {
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> ",
                __InternalStateTypeToString(valID) );
    } else if ( ( type == ISVT_TENSOR_S3 ) || ( type == ISVT_TENSOR_S3E ) ) {
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"9\" format=\"ascii\"> ",
                __InternalStateTypeToString(valID) );
    } else if ( type == ISVT_TENSOR_G ) {
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"9\" format=\"ascii\"> ",
                __InternalStateTypeToString(valID) );
    } else {
        fprintf( stderr, "VTKXMLExportModule::exportIntVarAs: unsupported variable type %s\n", __InternalStateTypeToString(valID) );
    }


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
            this->smoother->giveNodalVector(val, mapL2G.at(inode), ireg);
            if ( val == NULL ) {
	      iVal.resize( regionVarMap.giveSize() );
	      iVal.zero();
	      val = & iVal;
	      //OOFEM_ERROR2("VTKXMLExportModule::exportIntVars: smoothing error: invalid data in node %d", inode);
            }
        }

        if ( type == ISVT_SCALAR ) {
            if ( val->giveSize() ) {
                fprintf( stream, "%e ", val->at(1) );
            } else {
                fprintf(stream, "%e ", 0.0);
            }
        } else if ( type == ISVT_VECTOR ) {
            jsize = min( 3, val->giveSize() );
            for ( j = 1; j <= jsize; j++ ) {
                fprintf( stream, "%e ", val->at(j) );
            }

            for ( j = jsize + 1; j <= 3; j++ ) {
                fprintf(stream, "0.0 ");
            }

            fprintf(stream, " ");
        } else if ( type == ISVT_TENSOR_S3 ) {
            int ii, jj, iii;
            t.zero();
            for ( ii = 1; ii <= regionVarMap.giveSize(); ii++ ) {
                iii = regionVarMap.at(ii);

                if ( ( ii == 1 ) && iii ) {
                    t.at(1, 1) = val->at(iii);
                } else if ( ( ii == 2 ) && iii ) {
                    t.at(2, 2) = val->at( regionVarMap.at(2) );
                } else if ( ( ii == 3 ) && iii ) {
                    t.at(3, 3) = val->at( regionVarMap.at(3) );
                } else if ( ( ii == 4 ) && iii ) {
                    t.at(2, 3) = val->at( regionVarMap.at(4) );
                    t.at(3, 2) = val->at( regionVarMap.at(4) );
                } else if ( ( ii == 5 ) && iii ) {
                    t.at(1, 3) = val->at( regionVarMap.at(5) );
                    t.at(3, 1) = val->at( regionVarMap.at(5) );
                } else if ( ( ii == 6 ) && iii ) {
                    t.at(1, 2) = val->at( regionVarMap.at(6) );
                    t.at(2, 1) = val->at( regionVarMap.at(6) );
                }
            }

            for ( ii = 1; ii <= 3; ii++ ) {
                for ( jj = 1; jj <= 3; jj++ ) {
                    fprintf( stream, "%e ", t.at(ii, jj) );
                }
            }

            fprintf(stream, " ");
        } else if ( type == ISVT_TENSOR_G ) { // export general tensor values as scalars
            jsize = min( 9, val->giveSize() );
            for ( j = 1; j <= jsize; j++ ) {
                fprintf( stream, "%e ", val->at(j) );
            }

            for ( j = jsize + 1; j <= 9; j++ ) {
                fprintf(stream, "0.0 ");
            }

            fprintf(stream, " ");
        }
    } // end loop over dofmans

    fprintf(stream, "</DataArray>\n");
}



NodalRecoveryModel *
VTKXMLExportModule :: giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->smoother == NULL ) {
        if ( this->stype == VTK_Smother_NA ) {
            this->smoother  = new NodalAveragingRecoveryModel(d);
        } else if ( this->stype == VTK_Smoother_ZZ ) {
            this->smoother = new ZZNodalRecoveryModel(d);
        } else if ( this->stype == VTK_Smoother_SPR ) {
            this->smoother = new SPRNodalRecoveryModel(d);
        } else {
            OOFEM_ERROR2("VTKXMLExportModule: unsupported smoother type ID %d", this->stype);
        }
    }

    return this->smoother;
}

//keyword "primvars" in OOFEM input file
void
VTKXMLExportModule :: exportPrimaryVars(FILE *stream, IntArray &mapG2L, IntArray &mapL2G,
                                        int regionDofMans, int region, TimeStep *tStep)
{
    // should be performed over regions

    int i, n = primaryVarsToExport.giveSize();
    //int nnodes;
    //Domain *d = emodel->giveDomain(1);
    UnknownType type;

    if ( n == 0 ) {
        return;
    }

    //nnodes = d->giveNumberOfDofManagers();
    //fprintf(stream,"\n\nPOINT_DATA %d\n", nnodes);

    for ( i = 1; i <= n; i++ ) {
        type = ( UnknownType ) primaryVarsToExport.at(i);
        this->exportPrimVarAs(type, mapG2L, mapL2G, regionDofMans, region, stream, tStep);
    }
}


//keyword "primvars" in OOFEM input file
void
VTKXMLExportModule :: exportPrimVarAs(UnknownType valID, IntArray &mapG2L, IntArray &mapL2G,
                                      int regionDofMans, int ireg, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int inode;
    int j, jsize;
    FloatArray iVal;
    FloatMatrix t(3, 3);
    IntArray regionVarMap;
    InternalStateValueType type = ISVT_UNDEFINED;
    int nScalarComp = 1;

    if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
        type = ISVT_VECTOR;
    } else if ( ( valID == FluxVector ) || ( valID == PressureVector ) || ( valID == TemperatureVector ) ) {
        type = ISVT_SCALAR;
        //nScalarComp = d->giveNumberOfDefaultNodeDofs();
    } else {
        OOFEM_ERROR2( "VTKXMLExportModule::exportPrimVarAs: unsupported UnknownType %s", __UnknownTypeToString(valID) );
    }

    // print header
    if ( type == ISVT_SCALAR ) {
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\"> ", __UnknownTypeToString(valID) );
    } else if ( type == ISVT_VECTOR ) {
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> ", __UnknownTypeToString(valID) );
    } else {
        fprintf( stderr, "VTKXMLExportModule::exportPrimVarAs: unsupported variable type %s\n", __UnknownTypeToString(valID) );
    }

    DofManager *dman;
    DofID id;
    int numberOfDofs;
    IntArray mask(3);
    FloatArray dl, dg;

    for ( inode = 1; inode <= regionDofMans; inode++ ) {
        dman = d->giveNode( mapL2G.at(inode) );
        numberOfDofs = dman->giveNumberOfDofs();


        if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
            /*
             * iVal.resize(3);
             * iVal.zero();
             *
             * for ( j = 1; j <= numberOfDofs; j++ ) {
             *    id = dman->giveDof(j)->giveDofID();
             *    if ( ( id == V_u ) || ( id == D_u ) ) {
             *        iVal.at(1) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
             *    } else if ( ( id == V_v ) || ( id == D_v ) ) {
             *        iVal.at(2) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
             *    } else if ( ( id == V_w ) || ( id == D_w ) ) {
             *        iVal.at(3) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
             *    }
             */
            iVal.resize(3);
            iVal.zero();
            mask.resize(0);

            for ( j = 1; j <= numberOfDofs; j++ ) {
                id = dman->giveDof(j)->giveDofID();
                if ( ( id == V_u ) || ( id == D_u ) || ( id == V_v ) || ( id == D_v ) || ( id == V_w ) || ( id == D_w ) ) {
                    mask.followedBy(id);
                }
            }

            if ( dman->requiresTransformation() || dman->hasAnySlaveDofs() ) {
                // handle local coordinate system, slave dofs, etc in node
                dman->giveUnknownVector(dl, mask, EID_MomentumBalance, VM_Total, tStep);
                dman->computeDofTransformation(t, & mask, _toGlobalCS);
                dg.beProductOf(t, dl);
            } else {
                dman->giveUnknownVector(dg, mask, EID_MomentumBalance, VM_Total, tStep);
            }

            for ( j = 1; j <= mask.giveSize(); j++ ) {
                id = dman->giveDof(j)->giveDofID();
                if ( ( id == V_u ) || ( id == D_u ) ) {
                    iVal.at(1) = dg.at(j);
                } else if ( ( id == V_v ) || ( id == D_v ) ) {
                    iVal.at(2) = dg.at(j);
                } else if ( ( id == V_w ) || ( id == D_w ) ) {
                    iVal.at(3) = dg.at(j);
                }
            }
        } else if ( valID == FluxVector ) {
            iVal.resize(1);

            for ( j = 1; j <= numberOfDofs; j++ ) {
                id = dman->giveDof(j)->giveDofID();
                if ( id == C_1 ) {
                    iVal.at(1) = dman->giveDof(j)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                }
            }
        } else if ( valID == TemperatureVector ) {
            iVal.resize(1);

            for ( j = 1; j <= numberOfDofs; j++ ) {
                id = dman->giveDof(j)->giveDofID();
                if ( id == T_f ) {
                    iVal.at(1) = dman->giveDof(j)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                }
            }
        } else if ( valID == PressureVector ) {
            iVal.resize(1);

            for ( j = 1; j <= numberOfDofs; j++ ) {
                id = dman->giveDof(j)->giveDofID();
                if ( ( id == P_f ) ) {
                    iVal.at(1) = dman->giveDof(j)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                }
            }
        } else {
            OOFEM_ERROR2( "VTKXMLExportModule: unsupported unknownType %s", __UnknownTypeToString(valID) );
            //d->giveDofManager(regionNodalNumbers.at(inode))->giveUnknownVector(iVal, d->giveDefaultNodeDofIDArry(), valID, VM_Total, tStep);
        }

        if ( type == ISVT_SCALAR ) {
            if ( iVal.giveSize() ) {
                for ( j = 1; j <= nScalarComp; j++ ) {
                    fprintf( stream, "%e ", iVal.at(j) );
                }

                fprintf(stream, " ");
            } else {
                fprintf(stream, "%e ", 0.0);
                fprintf(stream, " ");
            }
        } else if ( type == ISVT_VECTOR ) {
            jsize = min( 3, iVal.giveSize() );
            for ( j = 1; j <= jsize; j++ ) {
                fprintf( stream, "%e ", iVal.at(j) );
            }

            for ( j = jsize + 1; j <= 3; j++ ) {
                fprintf(stream, "0.0 ");
            }

            fprintf(stream, " ");
        }
    } // end loop over nodes

    fprintf(stream, "</DataArray>\n");
}


void VTKXMLExportModule :: exportCellVars(FILE *stream, int totalcells, int region, TimeStep *tStep)
{
    int i, n = cellVarsToExport.giveSize();
    InternalStateType type;

    if ( n == 0 ) {
        return;
    }

    //print header
    fprintf(stream, "<CellData Scalars=\"\" Vectors=\"\" Tensors=\"\">\n");  // should contain a list of InternalStateType

    for ( i = 1; i <= n; i++ ) {
        type = ( InternalStateType ) cellVarsToExport.at(i);
        this->exportCellVarAs(type, totalcells, region, stream, tStep);
    }

    //print footer
    fprintf(stream, "</CellData>\n");
}


//keyword "cellvars" in OOFEM input file
void
VTKXMLExportModule :: exportCellVarAs(InternalStateType type, int nelem, int region, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int ielem;
    int pos;
    Element *elem;
    FloatMatrix mtrx(3, 3);
    IntegrationRule *iRule;
    GaussPoint *gp;
    FloatArray answer;

    switch ( type ) {
    case IST_MaterialNumber:
    case IST_ElementNumber:
    case IST_AverageTemperature:
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\"> ", __InternalStateTypeToString(type) );
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            if ( type == IST_MaterialNumber ) {
                fprintf( stream, "%d ", elem->giveMaterial()->giveNumber() );
            } else if ( type == IST_ElementNumber ) {
                fprintf( stream, "%d ", elem->giveNumber() );
            } else if ( type == IST_AverageTemperature ) {    //grab from the first IP
                iRule = elem->giveDefaultIntegrationRulePtr();
                gp  = iRule->getIntegrationPoint(0);
                gp->giveMaterialStatus();
                elem->giveIPValue(answer, gp, IST_AverageTemperature, tStep);
                fprintf( stream, "%f ", answer.at(1) );
            } else {
                OOFEM_ERROR2( "Unsupported Cell variable %s\n", __InternalStateTypeToString(type) );
            }
        }

        break;

    case IST_MaterialOrientation_x:
    case IST_MaterialOrientation_y:
    case IST_MaterialOrientation_z:
        fprintf( stream, "<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\"> ", __InternalStateTypeToString(type) );
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
            if ( !d->giveElement(ielem)->giveLocalCoordinateSystem(mtrx) ) {
                mtrx.resize(3, 3);
                mtrx.zero();
            }

            fprintf( stream, "%f %f %f  ", mtrx.at(1, pos), mtrx.at(2, pos), mtrx.at(3, pos) );
        }

        break;

    default:
        OOFEM_ERROR2( "Unsupported Cell variable %s", __InternalStateTypeToString(type) );
    }


    fprintf(stream, "</DataArray>\n");
}
} // end namespace oofem
