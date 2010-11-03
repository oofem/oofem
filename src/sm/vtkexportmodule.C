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

/* COMPONENTS in TENSORS like stress or strain
 *    x  y  z
 * x  1  6  5
 * y  6  2  4
 * z  5  4  3
 *
 *  PARAVIEW - stresses and strains in global c.s., damage tensor in local c.s.
 *    x  y  z
 * x  0  1  2
 * y  3  4  5
 * z  6  7  8
 *
 */


#include "vtkexportmodule.h"
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
#endif

namespace oofem {
VTKExportModule :: VTKExportModule(EngngModel *e) : ExportModule(e), internalVarsToExport(), primaryVarsToExport()
{
    //this->mode = rbrmode;
    this->mode = wdmode; //preserves node numbering
    this->outMode = rbrmode; //applies only when mode == rbrmode
    smoother = NULL;
}


VTKExportModule :: ~VTKExportModule()
{
    if ( this->smoother ) {
        delete this->smoother;
    }
}


IRResultType
VTKExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    ExportModule :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, IFT_VTKExportModule_cellvars, "cellvars"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, IFT_VTKExportModule_vars, "vars"); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, IFT_VTKExportModule_primvars, "primvars"); // Macro - see dofiditem.h

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_VTKExportModule_stype, "stype"); // Macro
    stype = ( VTKEM_SmootherType ) val;

    regionsToSkip.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, regionsToSkip, IFT_VTKExportModule_regionstoskip, "regionstoskip"); // Macro

    return IRRT_OK;
}


void
VTKExportModule :: doOutput(TimeStep *tStep)
{
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    FILE *stream = this->giveOutputStream(tStep);

    fprintf(stream, "# vtk DataFile Version 2.0\n");
    fprintf(stream, "Output for time %f\n", tStep->giveTime() );
    fprintf(stream, "ASCII\n");

    fprintf(stream, "DATASET UNSTRUCTURED_GRID\n");


    Domain *d  = emodel->giveDomain(1);
    FloatArray *coords;
    int i, inode, nnodes = d->giveNumberOfDofManagers();
    // output points

    if ( ( this->mode == wdmode ) || ( ( this->mode == rbrmode ) && ( this->outMode == wdmode ) ) ) {
        int nnodes = this->giveTotalRBRNumberOfNodes(d);
        fprintf(stream, "POINTS %d float\n", nnodes);
        int ireg = -1;
        int regionDofMans;
        IntArray map( d->giveNumberOfDofManagers() );

        // asemble local->global region map
        this->initRegionNodeNumbering(map, regionDofMans, 0, d, ireg, 1);

        printf("vktexportModule: %d %d\n", nnodes, regionDofMans);
        for ( inode = 1; inode <= regionDofMans; inode++ ) {
            coords = d->giveNode( map.at(inode) )->giveCoordinates();
            for ( i = 1; i <= coords->giveSize(); i++ ) {
                fprintf( stream, "%e ", coords->at(i) );
            }

            for ( i = coords->giveSize() + 1; i <= 3; i++ ) {
                fprintf(stream, "%e ", 0.0);
            }

            fprintf(stream, "\n");
        }

        /*
         * fprintf(stream, "POINTS %d double\n", nnodes);
         * for ( inode = 1; inode <= nnodes; inode++ ) {
         * coords = d->giveNode(inode)->giveCoordinates();
         * for ( i = 1; i <= coords->giveSize(); i++ ) {
         * fprintf( stream, "%e ", coords->at(i) );
         * }
         *
         * for ( i = coords->giveSize() + 1; i <= 3; i++ ) {
         * fprintf(stream, "%e ", 0.0);
         * }
         *
         * fprintf(stream, "\n");
         * }
         */
    } else { // outMode = rbrmode (Region By Region)
        // output nodes Region By Region
        int nnodes = this->giveTotalRBRNumberOfNodes(d);
        fprintf(stream, "POINTS %d float\n", nnodes);
        int ireg, nregions = d->giveNumberOfRegions();
        int regionDofMans;
        IntArray map( d->giveNumberOfDofManagers() );
        for ( ireg = 1; ireg <= nregions; ireg++ ) {
            if ( this->regionsToSkip.contains(ireg) ) {
                continue;
            }

            // asemble local->global region map
            this->initRegionNodeNumbering(map, regionDofMans, 0, d, ireg, 1);

            for ( inode = 1; inode <= regionDofMans; inode++ ) {
                coords = d->giveNode( map.at(inode) )->giveCoordinates();
                for ( i = 1; i <= coords->giveSize(); i++ ) {
                    fprintf( stream, "%e ", coords->at(i) );
                }

                for ( i = coords->giveSize() + 1; i <= 3; i++ ) {
                    fprintf(stream, "%e ", 0.0);
                }

                fprintf(stream, "\n");
            }
        }
    }

    int ielem, nelem = d->giveNumberOfElements(), elemToProcess = 0;
    int ncells, celllistsize = 0;
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        if ( ( this->mode == rbrmode ) && ( this->regionsToSkip.contains( d->giveElement(ielem)->giveRegionNumber() ) ) ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( d->giveElement(ielem)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        elemToProcess++;
        // element composed from same-type cells asumed
        ncells = this->giveNumberOfElementCells( d->giveElement(ielem) );
        celllistsize += ncells + ncells *this-> giveNumberOfNodesPerCell( this->giveCellType ( d->giveElement ( ielem ) ) );
    }

    int nelemNodes;
    Element *elem;
    int vtkCellType;
    IntArray cellNodes;
    // output cells
    fprintf(stream, "\nCELLS %d %d\n", elemToProcess, celllistsize);

    if ( ( this->mode == wdmode ) || ( ( this->mode == rbrmode ) && ( this->outMode == wdmode ) ) ) {
        IntArray regionNodalNumbers(nnodes);
        int regionDofMans = 0, offset = 0;
        int ireg = -1;

        // assemble global->local map
        this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 0);
        offset += regionDofMans;
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            vtkCellType = this->giveCellType(elem);

            nelemNodes = this->giveNumberOfNodesPerCell( vtkCellType ); //elem->giveNumberOfNodes(); // It HAS to be the same size as giveNumberOfNodesPerCell, otherwise the file will be incorrect.
	    this->giveElementCell (cellNodes, elem, 0);
            fprintf(stream, "%d ", nelemNodes);
            for ( i = 1; i <= nelemNodes; i++ ) {
	      fprintf(stream, "%d ", regionNodalNumbers.at( cellNodes.at(i)));
            }

            fprintf(stream, "\n");
        }

        // output cell types
        int vtkCellType;
        fprintf(stream, "\nCELL_TYPES %d\n", elemToProcess);
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = d->giveElement(ielem);
#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            vtkCellType = this->giveCellType(elem);
            fprintf(stream, "%d\n", vtkCellType);
        }

        // output cell data (Material ID ...)
        if ( cellVarsToExport.giveSize() ) {
            exportCellVars(stream, elemToProcess, tStep);
        }

        /*
         * for ( ielem = 1; ielem <= nelem; ielem++ ) {
         * elem = d->giveElement(ielem);
         * vtkCellType = this->giveCellType (elem);
         * nelemNodes = elem->giveNumberOfNodes();
         * fprintf(stream, "%d ", nelemNodes);
         * for ( i = 1; i <= nelemNodes; i++ ) {
         * if(vtkCellType==25){
         * j=HexaQuadNodeMapping[i-1];
         * }
         * else {
         * j=i;
         * }
         * fprintf(stream, "%d ", elem->giveNode(j)->giveNumber() - 1);
         * }
         *
         * fprintf(stream, "\n");
         * }
         *
         * // output cell types
         * fprintf(stream, "\nCELL_TYPES %d\n", elemToProcess);
         * for ( ielem = 1; ielem <= nelem; ielem++ ) {
         * elem = d->giveElement(ielem);
         * vtkCellType = this->giveCellType(elem);
         * fprintf(stream, "%d\n", vtkCellType);
         * }
         */
    } else { // rbr mode
        IntArray regionNodalNumbers(nnodes);
        int regionDofMans = 0, offset = 0;
        int ireg, nregions = d->giveNumberOfRegions();
        for ( ireg = 1; ireg <= nregions; ireg++ ) {
            if ( this->regionsToSkip.contains(ireg) ) {
                continue;
            }

            // assemble global->local map
            this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 0);
            offset += regionDofMans;
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);

                if ( elem->giveRegionNumber() != ireg ) {
                    continue;
                }

                vtkCellType = this->giveCellType(elem);
#ifdef __PARALLEL_MODE
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif

                nelemNodes = elem->giveNumberOfNodes();
		this->giveElementCell (cellNodes, elem, 0);

                fprintf(stream, "%d ", nelemNodes);
                for ( i = 1; i <= nelemNodes; i++ ) {
		  fprintf(stream, "%d ", regionNodalNumbers.at( cellNodes.at(i)));
                }

                fprintf(stream, "\n");
            }
        }

        for ( ireg = 1; ireg <= nregions; ireg++ ) {
            if ( this->regionsToSkip.contains(ireg) ) {
                continue;
            }

            // output cell types
            int vtkCellType;
            fprintf(stream, "\nCELL_TYPES %d\n", elemToProcess);
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
                if ( elem->giveRegionNumber() != ireg ) {
                    continue;
                }

#ifdef __PARALLEL_MODE
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif

                vtkCellType = this->giveCellType(elem);
                fprintf(stream, "%d\n", vtkCellType);
            }

            if ( cellVarsToExport.giveSize() ) {
                exportCellVars(stream, elemToProcess, tStep);
            }
        }
    }

    if ( primaryVarsToExport.giveSize() || internalVarsToExport.giveSize() ) {
        fprintf( stream, "\n\nPOINT_DATA %d\n", this->giveTotalRBRNumberOfNodes(d) );
    }

    this->exportPrimaryVars(stream, tStep);
    this->exportIntVars(stream, tStep);

    fclose(stream);
}

void
VTKExportModule :: initialize()
{
    if ( this->smoother ) {
        delete this->smoother;
        this->smoother = NULL;
    }
}


void
VTKExportModule :: terminate()
{ }


FILE *
VTKExportModule :: giveOutputStream(TimeStep *tStep)
{
    char baseFileName [ MAX_FILENAME_LENGTH ];
    char fileName [ MAX_FILENAME_LENGTH ];
    FILE *answer;

    emodel->giveOutputBaseFileName(baseFileName, MAX_FILENAME_LENGTH);
#ifdef __PARALLEL_MODE
    sprintf( fileName, "%s_%03d.%d.vtk", baseFileName, emodel->giveRank(), tStep->giveNumber() );
#else
    sprintf( fileName, "%s.%d.vtk", baseFileName, tStep->giveNumber() );
#endif
    if ( ( answer = fopen(fileName, "w") ) == NULL ) {
        OOFEM_ERROR2("VTKExportModule::giveOutputStream: failed to open file %s", fileName);
    }

    return answer;
}

int
VTKExportModule :: giveCellType(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int vtkCellType = 0;
    if ( elemGT == EGT_line_1 ) {
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
        OOFEM_ERROR("VTKExportModule: unsupported element gemetry type");
    }

    return vtkCellType;
}

int
VTKExportModule :: giveNumberOfNodesPerCell(int cellType)
{
	switch (cellType) {
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
		OOFEM_ERROR("VTKExportModule: unsupported cell type ID");
	}
    return 0; // to make compiler happy
}


void
VTKExportModule :: giveElementCell(IntArray &answer, Element *elem, int cell)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int i, nelemNodes;

    if ( ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
	 ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
	 ( elemGT == EGT_tetra_1 ) ||
	 ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
	 ( elemGT == EGT_hexa_1 ) ) {
      nelemNodes = elem->giveNumberOfNodes();
      answer.resize(nelemNodes);
      for ( i = 1; i <= nelemNodes; i++ ) {
	answer.at(i) = elem->giveNode(i)->giveNumber() - 1;
      }
    } else if (elemGT == EGT_hexa_2) {
      int HexaQuadNodeMapping [] = {
        5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18
      };
      nelemNodes = elem->giveNumberOfNodes();
      answer.resize(nelemNodes);
      for ( i = 1; i <= nelemNodes; i++ ) {
	answer.at(i) = elem->giveNode(HexaQuadNodeMapping[i-1])->giveNumber() - 1;
      }
    } else {
        OOFEM_ERROR("VTKExportModule: unsupported element geometry type");
    }

    return;
}


int
VTKExportModule :: giveNumberOfElementCells(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();

    if ( ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) || ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) || ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) || ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_hexa_2 ) ) {
        return 1;
    } else {
        OOFEM_ERROR("VTKExportModule: unsupported element geometry type");
    }

    return 0;
}

//keyword "vars" in OOFEM input file
void
VTKExportModule :: exportIntVars(FILE *stream, TimeStep *tStep)
{
    // should be performed over regions

    int i, n = internalVarsToExport.giveSize();
    //int nnodes;
    //Domain *d = emodel->giveDomain(1);
    InternalStateType type;

    if ( n == 0 ) {
        return;
    }

    //nnodes = d->giveNumberOfDofManagers();

    /*
     #ifdef RBR_SUPPORT
     * if (outMode == wdmode)
     * fprintf(stream,"\n\nPOINT_DATA %d\n", nnodes);
     * else
     * fprintf(stream,"\n\nPOINT_DATA %d\n", this->giveTotalRBRNumberOfNodes(d));
     #else
     * fprintf(stream,"\n\nPOINT_DATA %d\n", nnodes);
     #endif
     */

    for ( i = 1; i <= n; i++ ) {
        type = ( InternalStateType ) internalVarsToExport.at(i);
        InternalStateValueType iType = giveInternalStateValueType(type);
        this->exportIntVarAs(type, iType, stream, tStep);
    }
}

void
VTKExportModule :: exportCellVars(FILE *stream, int elemToProcess, TimeStep *tStep) {
    int i, ielem, pos;
    InternalStateType type;
    Element *elem;
    Domain *d  = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    FloatMatrix mtrx(3, 3);

    fprintf(stream, "\nCELL_DATA %d\n", elemToProcess);
    for ( i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        type = ( InternalStateType ) cellVarsToExport.at(i);
        switch ( type ) {
        case IST_MaterialNumber:
        case IST_ElementNumber:
            fprintf( stream, "SCALARS %s int\nLOOKUP_TABLE default\n", __InternalStateTypeToString(type) );
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
                elem = d->giveElement(ielem);
#ifdef __PARALLEL_MODE
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif
                if ( type == IST_MaterialNumber ) {
                    fprintf( stream, "%d\n", elem->giveMaterial()->giveNumber() );
                }

                if ( type == IST_ElementNumber ) {
                    fprintf( stream, "%d\n", elem->giveNumber() );
                }
            }

            break;

        case IST_MaterialOrientation_x:
        case IST_MaterialOrientation_y:
        case IST_MaterialOrientation_z:
            if ( type == IST_MaterialOrientation_x ) {
                pos = 1;
            }

            if ( type == IST_MaterialOrientation_y ) {
                pos = 2;
            }

            if ( type == IST_MaterialOrientation_z ) {
                pos = 3;
            }

            fprintf( stream, "VECTORS %s float\n", __InternalStateTypeToString(type) );
            for ( ielem = 1; ielem <= nelem; ielem++ ) {
                if ( !d->giveElement(ielem)->giveLocalCoordinateSystem(mtrx) ) {
                    mtrx.resize(3, 3);
                    mtrx.zero();
                }

                fprintf( stream, "%f %f %f\n", mtrx.at(1, pos), mtrx.at(2, pos), mtrx.at(3, pos) );
            }

            break;

        default:
            OOFEM_ERROR2( "Quantity %s not defined on cells", __InternalStateTypeToString(type) );
        }

        fprintf(stream, "\n\n");
    }
}


int
VTKExportModule :: giveTotalRBRNumberOfNodes(Domain *d)
//
// Returns total number of nodes with region multiplicity
// (RBR = Region by Region)
// The nodes on region boundaries are for each region.
// We need unique node numbers for each region, to prevent
// vtk from smoothing the nodal values at region boundaries.
//
{
    Element *elem;
    int rbrnodes = 0, nnodes = d->giveNumberOfDofManagers(), nelems = d->giveNumberOfElements();
    std :: vector< char > map(nnodes);
    //char map[nnodes];
    int i, j, nregions = d->giveNumberOfRegions();
    int elemnodes, ielemnode;

    if ( ( this->mode == wdmode ) || ( ( this->mode == rbrmode ) && ( this->outMode == wdmode ) ) ) {
        for ( j = 0; j < nnodes; j++ ) {
            map [ j ] = 0;
        }

        for ( j = 1; j <= nelems; j++ ) {
            elem  = d->giveElement(j);
#ifdef __PARALLEL_MODE
            if ( d->giveElement(j)->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            elemnodes = elem->giveNumberOfNodes();
            for ( ielemnode = 1; ielemnode <= elemnodes; ielemnode++ ) {
                map [ elem->giveNode(ielemnode)->giveNumber() - 1 ] = 1;
            }
        }

        for ( j = 0; j < nnodes; j++ ) {
            rbrnodes += map [ j ];
        }
    } else {
        for ( i = 1; i <= nregions; i++ ) {
            if ( this->regionsToSkip.contains(i) ) {
                continue;
            }

            for ( j = 0; j < nnodes; j++ ) {
                map [ j ] = 0;
            }

            for ( j = 1; j <= nelems; j++ ) {
#ifdef __PARALLEL_MODE
                if ( d->giveElement(j)->giveParallelMode() != Element_local ) {
                    continue;
                }

#endif
                if ( ( elem  = d->giveElement(j) )->giveRegionNumber() == i ) {
                    elemnodes = elem->giveNumberOfNodes();
                    for ( ielemnode = 1; ielemnode <= elemnodes; ielemnode++ ) {
                        map [ elem->giveNode(ielemnode)->giveNumber() - 1 ] = 1;
                    }
                }
            }


            for ( j = 0; j < nnodes; j++ ) {
                rbrnodes += map [ j ];
            }
        }
    }

    return rbrnodes;
}


int
VTKExportModule :: initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans,
                                           int offset, Domain *domain, int reg, int mode)
{
    // if mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
    // The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).

    // if mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
    // The i-th value contains the corresponding global node number.


    int ielem, nelem = domain->giveNumberOfElements();
    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int elementNode, node;
    int currOffset = offset + 1;
    Element *element;

    regionNodalNumbers.resize(nnodes);
    regionNodalNumbers.zero();
    regionDofMans = 0;

    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
        if ( ( reg > 0 ) && ( element->giveRegionNumber() != reg ) ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        elemNodes = element->giveNumberOfNodes();
        //  elemSides = element->giveNumberOfSides();

        // determine local region node numbering
        for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            node = element->giveNode(elementNode)->giveNumber();
            if ( regionNodalNumbers.at(node) == 0 ) { // assign new number
                /* mark for assignement. This is done later, as it allows to preserve
                 * natural node numbering.
                 */
                regionNodalNumbers.at(node) = 1;
                regionDofMans++;
            }
        }
    }

    if ( mode == 1 ) {
        IntArray answer(nnodes);
        int i;
        for ( i = 1; i <= nnodes; i++ ) {
            if ( regionNodalNumbers.at(i) ) {
                regionNodalNumbers.at(i) = currOffset++;
                answer.at( regionNodalNumbers.at(i) ) = i;
            }
        }

        regionNodalNumbers = answer;
    } else {
        int i;
        for ( i = 1; i <= nnodes; i++ ) {
            if ( regionNodalNumbers.at(i) ) {
                regionNodalNumbers.at(i) = currOffset++;
            }
        }
    }

    //  for (elementSide = 1; elementSide<= elemSides; elementSide++) {
    //   node = element->giveSide(elementSide)->giveNumber();
    //   if (regionNodalNumbers.at(node) == 0) {// assign new number
    //    regionNodalNumbers.at(node) = currOffset++;
    //      regionDofMans++;
    //  }
    return 1;
}

//keyword "vars" in OOFEM input file
void
VTKExportModule :: exportIntVarAs(InternalStateType valID, InternalStateValueType type, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int ireg, nregions = d->giveNumberOfRegions();
    int nnodes = d->giveNumberOfDofManagers(), inode;
    int i, j, jsize, mapindx;
    FloatArray iVal(3);
    FloatMatrix t(3, 3);
    const FloatArray *val;
    IntArray regionVarMap;

    this->giveSmoother();

    int nindx = 1;
    if ( type == ISVT_TENSOR_G ) {
        this->smoother->giveRegionRecordMap(regionVarMap, 1, valID);
        nindx = regionVarMap.giveSize();
    }

    for ( int indx = 1; indx <= nindx; indx++ ) {
        // print header
        if ( type == ISVT_SCALAR ) {
            fprintf( stream, "SCALARS %s float 1\n", __InternalStateTypeToString(valID) );
        } else if ( type == ISVT_VECTOR ) {
            fprintf( stream, "VECTORS %s float\n", __InternalStateTypeToString(valID) );
        } else if ( ( type == ISVT_TENSOR_S3 ) || ( type == ISVT_TENSOR_S3E ) ) {
            fprintf( stream, "TENSORS %s float\n", __InternalStateTypeToString(valID) );
        } else if ( type == ISVT_TENSOR_G ) {
            fprintf(stream, "SCALARS %s_%d float 1\n", __InternalStateTypeToString(valID), indx);
        } else {
            fprintf(stderr, "exportIntVarAs: unsupported variable type\n");
        }

        if ( ( type == ISVT_SCALAR ) || ( type == ISVT_TENSOR_G ) ) {
            fprintf(stream, "LOOKUP_TABLE default\n");
        }

        if ( !( ( valID == IST_DisplacementVector ) || ( valID == IST_MaterialInterfaceVal ) ) ) {
            this->smoother->recoverValues(valID, tStep);
        }

        if ( this->mode == wdmode ) { // (Whole Domain)
            IntArray regionNodalNumbers(nnodes);
            int regionDofMans = 0, offset = 0;
            ireg = -1;

            // assemble local->global map
            this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);
            this->smoother->giveRegionRecordMap(regionVarMap, ireg, valID);

            for ( inode = 1; inode <= regionDofMans; inode++ ) {
                if ( valID == IST_DisplacementVector ) {
                    iVal.resize(3);
                    val = & iVal;
                    for ( j = 1; j <= 3; j++ ) {
                        iVal.at(j) = d->giveNode( regionNodalNumbers.at(inode) )->giveUpdatedCoordinate(j, tStep, EID_MomentumBalance, 1.0) -
                                     d->giveNode( regionNodalNumbers.at(inode) )->giveCoordinate(j);
                    }
                } else if ( valID == IST_MaterialInterfaceVal ) {
                    MaterialInterface *mi = emodel->giveMaterialInterface(1);
                    if ( mi ) {
                        iVal.resize(1);
                        val = & iVal;
                        iVal.at(1) = mi->giveNodalScalarRepresentation( regionNodalNumbers.at(inode) );
                    }
                } else {
                    for ( i = 1; i <= nregions; i++ ) {
                        this->smoother->giveNodalVector(val, regionNodalNumbers.at(inode), i);
                        if ( val ) {
                            break;
                        }
                    }

                    if ( val == NULL ) {
                        iVal.resize( regionVarMap.giveSize() );
                        iVal.zero();
                        val = & iVal;
                        //OOFEM_ERROR ("VTKExportModule::exportIntVars: internal error: invalid dofman data");
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

                        fprintf(stream, "\n");
                    }
                } else if ( type == ISVT_TENSOR_G ) { // export general tensor values as scalars
                    if ( ( indx > 0 ) && ( indx <= regionVarMap.giveSize() ) ) {
                        mapindx = regionVarMap.at(indx);
                    } else {
                        mapindx = 0;
                    }

                    if ( ( mapindx > 0 ) && ( mapindx <= val->giveSize() ) ) {
                        fprintf( stream, "%e ", val->at(indx) );
                    } else {
                        fprintf(stream, "0.0 ");
                    }
                }

                fprintf(stream, "\n");
            }
        } else { // RBR mode
            IntArray regionNodalNumbers(nnodes);
            int regionDofMans = 0, offset = 0;
            for ( ireg = 1; ireg <= nregions; ireg++ ) {
                if ( this->regionsToSkip.contains(ireg) ) {
                    continue;
                }

                // assemble local->global map
                this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);
                this->smoother->giveRegionRecordMap(regionVarMap, ireg, valID);

                for ( inode = 1; inode <= regionDofMans; inode++ ) {
                    if ( valID == IST_DisplacementVector ) {
                        iVal.resize(3);
                        val = & iVal;
                        for ( j = 1; j <= 3; j++ ) {
                            iVal.at(j) = d->giveNode( regionNodalNumbers.at(inode) )->giveUpdatedCoordinate(j, tStep, EID_MomentumBalance, 1.0) -
                                         d->giveNode( regionNodalNumbers.at(inode) )->giveCoordinate(j);
                        }
                    } else if ( valID == IST_MaterialInterfaceVal ) {
                        MaterialInterface *mi = emodel->giveMaterialInterface(1);
                        if ( mi ) {
                            iVal.resize(1);
                            val = & iVal;
                            iVal.at(1) = mi->giveNodalScalarRepresentation( regionNodalNumbers.at(inode) );
                        }
                    } else {
                        this->smoother->giveNodalVector(val, regionNodalNumbers.at(inode), ireg);
                        if ( val == NULL ) {
                            OOFEM_ERROR("VTKExportModule::exportIntVars: internal error: invalid dofman data");
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

                            fprintf(stream, "\n");
                        }
                    } else if ( type == ISVT_TENSOR_G ) { // export general tensor values as scalars
                        if ( ( indx > 0 ) && ( indx <= regionVarMap.giveSize() ) ) {
                            mapindx = regionVarMap.at(indx);
                        } else {
                            mapindx = 0;
                        }

                        if ( ( mapindx > 0 ) && ( mapindx <= val->giveSize() ) ) {
                            fprintf( stream, "%e ", val->at(indx) );
                        } else {
                            fprintf(stream, "0.0 ");
                        }
                    }

                    fprintf(stream, "\n");
                }
            }
        }

        fprintf(stream, "\n");
    }
}


NodalRecoveryModel *
VTKExportModule :: giveSmoother()
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
            OOFEM_ERROR("VTKExportModule: unsupported smoother type ID");
        }
    }

    return this->smoother;
}


//keyword "primvars" in OOFEM input file
void
VTKExportModule :: exportPrimaryVars(FILE *stream, TimeStep *tStep)
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
        this->exportPrimVarAs(type, stream, tStep);
    }
}

//keyword "primvars" in OOFEM input file
void
VTKExportModule :: exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int ireg, nregions = d->giveNumberOfRegions();
    int nnodes = d->giveNumberOfDofManagers(), inode;
    int j, jsize;
    FloatArray iVal, iValLCS;
    FloatMatrix t(3, 3);
    IntArray regionVarMap, dofIDMask;
    InternalStateValueType type = ISVT_UNDEFINED;
    int nScalarComp = 1;

    if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
        type = ISVT_VECTOR;
    } else if ( ( valID == FluxVector ) || ( valID == PressureVector ) || ( valID == TemperatureVector ) || ( valID == TemperatureVector ) ) {
        type = ISVT_SCALAR;
        //nScalarComp = d->giveNumberOfDefaultNodeDofs();
    } else {
        OOFEM_ERROR2( "VTKExportModule::exportPrimVarAs: unsupported UnknownType (%s)", __UnknownTypeToString(valID) );
    }

    // print header
    if ( type == ISVT_SCALAR ) {
        fprintf(stream, "SCALARS %s float %d\n",  __UnknownTypeToString(valID), nScalarComp);
    } else if ( type == ISVT_VECTOR ) {
        fprintf( stream, "VECTORS %s float\n", __UnknownTypeToString(valID) );
    } else {
        fprintf(stderr, "exportPrimVarAs: unsupported variable type\n");
    }

    if ( type == ISVT_SCALAR ) {
        fprintf(stream, "LOOKUP_TABLE default\n");
    }

    DofManager *dman;
    DofID id;
    int numberOfDofs;

    if ( this->mode == wdmode ) { // (Whole Domain)
        IntArray regionNodalNumbers(nnodes);
        int regionDofMans = 0, offset = 0;
        ireg = -1;


        // assemble local->global map
        this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);

        // Used for special nodes/elements with mixed number of dofs.
        this->giveSmoother();

        for ( inode = 1; inode <= regionDofMans; inode++ ) {
            dman = d->giveNode( regionNodalNumbers.at(inode) );
            numberOfDofs = dman->giveNumberOfDofs();

            if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
                iVal.resize(3);
                iVal.zero();

                for ( j = 1; j <= numberOfDofs; j++ ) {
                    id = dman->giveDof(j)->giveDofID();
                    if ( ( id == V_u ) || ( id == D_u ) ) {
                        iVal.at(1) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                    } else if ( ( id == V_v ) || ( id == D_v ) ) {
                        iVal.at(2) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                    } else if ( ( id == V_w ) || ( id == D_w ) ) {
                        iVal.at(3) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                    }
                }

                /*
                 * iVal.resize(3);
                 * for (j=1; j<= 3 ; j++) {
                 * iVal.at(j) = d->giveNode(regionNodalNumbers.at(inode))->giveUpdatedCoordinate(j, tStep,EID_MomentumBalance, 1.0) -
                 * d->giveNode(regionNodalNumbers.at(inode))->giveCoordinate(j);
                 * }
                 */
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
	      dofIDMask.resize(1); 
	      dofIDMask.at(1) = P_f;
	      this->getDofManPrimaryVariable (iVal, dman, dofIDMask, EID_ConservationEquation, 
					      VM_Total, tStep, IST_Pressure);

            } else {
                OOFEM_ERROR2( "VTKExportModule: unsupported unknownType (%s)", __UnknownTypeToString(valID) );
                //d->giveDofManager(regionNodalNumbers.at(inode))->giveUnknownVector(iVal, d->giveDefaultNodeDofIDArry(), valID, VM_Total, tStep);
            }

            if ( type == ISVT_SCALAR ) {
                if ( iVal.giveSize() ) {
                    for ( j = 1; j <= nScalarComp; j++ ) {
                        fprintf( stream, "%e ", iVal.at(j) );
                    }
                } else {
                    fprintf(stream, "%e ", 0.0);
                }
            } else if ( type == ISVT_VECTOR ) {
                jsize = min( 3, iVal.giveSize() );
                //rotate back from nodal CS to global CS if applies
                if ( d->giveNode( dman->giveNumber() )->hasLocalCS() ) {
                    iVal.resize(3);
                    iValLCS = iVal;
                    iVal.beTProductOf(* d->giveNode( dman->giveNumber() )->giveLocalCoordinateTriplet(), iValLCS);
                }

                for ( j = 1; j <= jsize; j++ ) {
                    fprintf( stream, "%e ", iVal.at(j) );
                }

                for ( j = jsize + 1; j <= 3; j++ ) {
                    fprintf(stream, "0.0 ");
                }
            }

            fprintf(stream, "\n");
        }
    } else {  // RBR mode
        IntArray regionNodalNumbers(nnodes);
        int regionDofMans = 0, offset = 0;
        for ( ireg = 1; ireg <= nregions; ireg++ ) {
            if ( this->regionsToSkip.contains(ireg) ) {
                continue;
            }

            // assemble local->global map
            this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);

            for ( inode = 1; inode <= regionDofMans; inode++ ) {
                dman = d->giveNode( regionNodalNumbers.at(inode) );
                numberOfDofs = dman->giveNumberOfDofs();

                if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
                    iVal.resize(3);
                    iVal.zero();

                    for ( j = 1; j <= numberOfDofs; j++ ) {
                        id = dman->giveDof(j)->giveDofID();
                        if ( ( id == V_u ) || ( id == D_u ) ) {
                            iVal.at(1) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                        } else if ( ( id == V_v ) || ( id == D_v ) ) {
                            iVal.at(2) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                        } else if ( ( id == V_w ) || ( id == D_w ) ) {
                            iVal.at(3) = dman->giveDof(j)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                        }
                    }

                    /*
                     * iVal.resize(3);
                     * for (j=1; j<= 3 ; j++) {
                     * iVal.at(j) = d->giveNode(regionNodalNumbers.at(inode))->giveUpdatedCoordinate(j, tStep,EID_MomentumBalance, 1.0) -
                     * d->giveNode(regionNodalNumbers.at(inode))->giveCoordinate(j);
                     * }
                     */
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
                    OOFEM_ERROR2( "VTKExportModule: unsupported unknownType (%s)", __UnknownTypeToString(valID) );
                    //d->giveDofManager(regionNodalNumbers.at(inode))->giveUnknownVector(iVal, d->giveDefaultNodeDofIDArry(), valID, VM_Total, tStep);
                }

                if ( type == ISVT_SCALAR ) {
                    if ( iVal.giveSize() ) {
                        for ( j = 1; j <= nScalarComp; j++ ) {
                            fprintf( stream, "%e ", iVal.at(j) );
                        }
                    } else {
                        fprintf(stream, "%e ", 0.0);
                    }
                } else if ( type == ISVT_VECTOR ) {
                    jsize = min( 3, iVal.giveSize() );
                    for ( j = 1; j <= jsize; j++ ) {
                        fprintf( stream, "%e ", iVal.at(j) );
                    }

                    for ( j = jsize + 1; j <= 3; j++ ) {
                        fprintf(stream, "0.0 ");
                    }
                }

                fprintf(stream, "\n");
            }
        }

        fprintf(stream, "\n");
    }

    fprintf(stream, "\n");
}

void
VTKExportModule :: getDofManPrimaryVariable (FloatArray& answer, DofManager* dman, IntArray& dofIDMask, EquationID type, 
					     ValueModeType mode, TimeStep *tStep, InternalStateType iType)
{
  int j, indx, size = dofIDMask.giveSize();
  const FloatArray* recoveredVal;
  answer.resize(size);
  // all values zero by default
  answer.zero();
  
  for ( j = 1; j <= size; j++ ) {
    if ((indx = dman->findDofWithDofId(dofIDMask.at(j)))) {
      // primary variable available directly in dof manager
      answer.at(j) = dman->giveDof(indx)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
    } else if (iType != IST_Undefined) {
      // ok primary variable not directly available
      // but equivalent InternalStateType provided
      // in this case use smoother to recover nodal value
      this->giveSmoother()->recoverValues(iType, tStep); /// recover values if not done before
      this->giveSmoother()->giveNodalVector(recoveredVal, dman->giveNumber(), 1);
      // here we have a lack of information about how to convert recoveredVal to response
      // if the size is compatible we accept it, otherwise an error is thrown
      if (size == recoveredVal->giveSize()) {
	answer.at(j)=recoveredVal->at(j);
      } else {
	OOFEM_WARNING ("VTKExportModule :: getDofManPrimaryVariable: recovered variable size mismatch");
	answer.at(j) = 0.0;
      }
    }
  }
}

} // end namespace oofem
