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
#include "gausspoint.h"
#include "engngm.h"
#include "node.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "cltypes.h"
#include "element.h"
#include "material.h"
#include "classfactory.h"
#include "crosssection.h"
#include "dof.h"

#include <vector>

namespace oofem {
REGISTER_ExportModule(VTKExportModule)

VTKExportModule :: VTKExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, _IFT_VTKExportModule_cellvars);
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_VTKExportModule_vars);
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_VTKExportModule_primvars);

    val = NodalRecoveryModel :: NRM_ZienkiewiczZhu;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_VTKExportModule_stype);
    stype = ( NodalRecoveryModel :: NodalRecoveryModelType ) val;

    return ExportModule :: initializeFrom(ir);
}


void
VTKExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    FILE *stream = this->giveOutputStream(tStep);

    fprintf(stream, "# vtk DataFile Version 2.0\n");
    fprintf( stream, "Output for time %f\n", tStep->giveTargetTime() );
    fprintf(stream, "ASCII\n");

    fprintf(stream, "DATASET UNSTRUCTURED_GRID\n");


    Domain *d  = emodel->giveDomain(1);
    FloatArray *coords;
    int i, inode, nnodes = d->giveNumberOfDofManagers();
    this->giveSmoother(); // make sure smoother is created

    // output points

    fprintf(stream, "POINTS %d double\n", nnodes);
    int ireg = -1;
    int regionDofMans;
    IntArray map( d->giveNumberOfDofManagers() );

    // asemble local->global region map
    this->initRegionNodeNumbering(map, regionDofMans, 0, d, ireg, 1);

    OOFEM_LOG_DEBUG("vktexportModule: %d %d\n", nnodes, regionDofMans);
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

    int elemToProcess = 0;
    int ncells, celllistsize = 0;
    for ( auto &elem : d->giveElements() ) {
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }

        elemToProcess++;
        // element composed from same-type cells asumed
        ncells = this->giveNumberOfElementCells( elem.get() );
        celllistsize += ncells + ncells *this->giveNumberOfNodesPerCell( this->giveCellType ( elem.get() ) );
    }

    int nelemNodes;
    int vtkCellType;
    IntArray cellNodes;
    // output cells
    fprintf(stream, "\nCELLS %d %d\n", elemToProcess, celllistsize);

    IntArray regionNodalNumbers(nnodes);

    // assemble global->local map
    this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, 0, d, ireg, 0);
    for ( auto &elem : d->giveElements() ) {
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }

        vtkCellType = this->giveCellType(elem.get());

        nelemNodes = this->giveNumberOfNodesPerCell(vtkCellType); //elem->giveNumberOfNodes(); // It HAS to be the same size as giveNumberOfNodesPerCell, otherwise the file will be incorrect.
        this->giveElementCell(cellNodes, elem.get(), 0);
        fprintf(stream, "%d ", nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            fprintf(stream, "%d ", regionNodalNumbers.at( cellNodes.at(i) ) - 1);
        }

        fprintf(stream, "\n");
    }

    // output cell types
    fprintf(stream, "\nCELL_TYPES %d\n", elemToProcess);
    for ( auto &elem : d->giveElements() ) {
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }

        vtkCellType = this->giveCellType(elem.get());
        fprintf(stream, "%d\n", vtkCellType);
    }

    // output cell data (Material ID ...)
    if ( cellVarsToExport.giveSize() ) {
        exportCellVars(stream, elemToProcess, tStep);
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
    ExportModule :: initialize();
}


void
VTKExportModule :: terminate()
{ }


FILE *
VTKExportModule :: giveOutputStream(TimeStep *tStep)
{
    std :: string fileName;
    FILE *answer;
    fileName = this->giveOutputBaseFileName(tStep) + ".vtk";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }
    return answer;
}

int
VTKExportModule :: giveCellType(Element *elem)
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
    }  else if ( elemGT == EGT_wedge_1 ) {
        vtkCellType = 13;
    } else if ( elemGT == EGT_wedge_2 ) {
        vtkCellType = 26;
    } else {
        OOFEM_ERROR("unsupported element gemetry type");
    }

    return vtkCellType;
}

int
VTKExportModule :: giveNumberOfNodesPerCell(int cellType)
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
        OOFEM_ERROR("unsupported cell type ID");
    }

    return 0; // to make compiler happy
}


void
VTKExportModule :: giveElementCell(IntArray &answer, Element *elem, int cell)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();
    int i, nelemNodes;

    if ( ( elemGT == EGT_point ) ||
        ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
        ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
        ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
        ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_wedge_1 ) ) {
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(i)->giveNumber();
        }
    } else if ( elemGT == EGT_hexa_2 ) {
        int HexaQuadNodeMapping [] = {
            5, 8, 7, 6, 1, 4, 3, 2, 16, 15, 14, 13, 12, 11, 10, 9, 17, 20, 19, 18
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(HexaQuadNodeMapping [ i - 1 ])->giveNumber();
        }
    } else if ( elemGT == EGT_wedge_2 ) {
        int WedgeQuadNodeMapping [] = {
            4, 6, 5, 1, 3, 2, 12, 11, 10, 9, 8, 7, 13, 15, 14
        };
        nelemNodes = elem->giveNumberOfNodes();
        answer.resize(nelemNodes);
        for ( i = 1; i <= nelemNodes; i++ ) {
            answer.at(i) = elem->giveNode(WedgeQuadNodeMapping [ i - 1 ])->giveNumber();
        }
    } else {
        OOFEM_ERROR("unsupported element geometry type");
    }
}


int
VTKExportModule :: giveNumberOfElementCells(Element *elem)
{
    Element_Geometry_Type elemGT = elem->giveGeometryType();

    if ( ( elemGT == EGT_point ) ||
        ( elemGT == EGT_line_1 ) || ( elemGT == EGT_line_2 ) ||
        ( elemGT == EGT_triangle_1 ) || ( elemGT == EGT_triangle_2 ) ||
        ( elemGT == EGT_tetra_1 ) || ( elemGT == EGT_tetra_2 ) ||
        ( elemGT == EGT_quad_1 ) || ( elemGT == EGT_quad_2 ) ||
        ( elemGT == EGT_hexa_1 ) || ( elemGT == EGT_hexa_2 ) ||
        ( elemGT == EGT_wedge_1 ) || ( elemGT == EGT_wedge_2 ) ) {
        return 1;
    } else {
        OOFEM_ERROR("unsupported element geometry type");
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

    for ( i = 1; i <= n; i++ ) {
        type = ( InternalStateType ) internalVarsToExport.at(i);
        InternalStateValueType iType = giveInternalStateValueType(type);
        this->exportIntVarAs(type, iType, stream, tStep);
    }
}

void
VTKExportModule :: exportCellVars(FILE *stream, int elemToProcess, TimeStep *tStep)
{
    int pos;
    InternalStateType type;
    Domain *d  = emodel->giveDomain(1);
    FloatMatrix mtrx(3, 3);
    FloatArray temp, vec;
    double gptot;

    fprintf(stream, "\nCELL_DATA %d\n", elemToProcess);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        type = ( InternalStateType ) cellVarsToExport.at(i);
        switch ( type ) {
        case IST_MaterialNumber:
        case IST_ElementNumber:
            fprintf( stream, "SCALARS %s int\nLOOKUP_TABLE default\n", __InternalStateTypeToString(type) );
            for ( auto &elem : d->giveElements() ) {
                if ( elem->giveParallelMode() != Element_local ) {
                    continue;
                }

                if ( type == IST_MaterialNumber || type == IST_CrossSectionNumber ) {
                    OOFEM_WARNING("Material numbers are deprecated, outputing cross section number instead...");
                    fprintf( stream, "%d\n", elem->giveCrossSection()->giveNumber() );
                } else if ( type == IST_ElementNumber ) {
                    fprintf( stream, "%d\n", elem->giveNumber() );
                } else {
                    OOFEM_ERROR("Unsupported Cell variable %s\n", __InternalStateTypeToString(type));
                }
            }

            break;

        case IST_MaterialOrientation_x:
        case IST_MaterialOrientation_y:
        case IST_MaterialOrientation_z:
            if ( type == IST_MaterialOrientation_x ) {
                pos = 1;
            } else if ( type == IST_MaterialOrientation_y ) {
                pos = 2;
            } else {
                pos = 3;
            }

            fprintf( stream, "VECTORS %s double\n", __InternalStateTypeToString(type) );
            for ( auto &elem : d->giveElements() ) {
                if ( !elem->giveLocalCoordinateSystem(mtrx) ) {
                    mtrx.resize(3, 3);
                    mtrx.zero();
                }

                fprintf( stream, "%f %f %f\n", mtrx.at(1, pos), mtrx.at(2, pos), mtrx.at(3, pos) );
            }

            break;
        default:
            InternalStateValueType vt = giveInternalStateValueType(type);
            switch ( vt ) {
            case ISVT_TENSOR_S3: // These could be written as tensors as well, if one wants to.
            case ISVT_TENSOR_S3E:
            case ISVT_VECTOR:
            case ISVT_SCALAR:
                if ( vt == ISVT_SCALAR ) {
                    fprintf( stream, "SCALARS %s double\nLOOKUP_TABLE default\n", __InternalStateTypeToString(type) );
                } else {
                    fprintf( stream, "VECTORS %s double\nLOOKUP_TABLE default\n", __InternalStateTypeToString(type) );
                }
                for ( auto &elem : d->giveElements() ) {
                    if ( elem->giveParallelMode() != Element_local ) {
                        continue;
                    }

                    gptot = 0;
                    vec.clear();
                    for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
                        elem->giveIPValue(temp, gp, type, tStep);
                        gptot += gp->giveWeight();
                        vec.add(gp->giveWeight(), temp);
                    }
                    vec.times(1 / gptot);
                    for ( int j = 1; j <= vec.giveSize(); ++j ) {
                        fprintf( stream, "%e ", vec.at(j) );
                    }
                    fprintf(stream, "\n");
                }
                break;
#if 0 // Hardly even worth the effort...
            case ISVT_TENSOR_G:
                for ( int indx = 1; indx < 9; ++indx ) {
                    fprintf(stream, "SCALARS %s_%d double 1\n", __InternalStateTypeToString(valID), indx);

                    for ( ielem = 1; ielem <= nelem; ielem++ ) {
                        elem = d->giveElement(ielem);
                        elem->giveIPValue(vec, elem->giveDefaultIntegrationRulePtr()->getIntegrationPoint(1), type, tStep);
                        fprintf( stream, "%e ", vec.at(indx) );
                    }
                }
#endif
            default:
                OOFEM_ERROR("Quantity %s not handled yet.", __InternalStateTypeToString(type));
            }
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
    int rbrnodes = 0, nnodes = d->giveNumberOfDofManagers();
    std :: vector< char > map(nnodes);
    //char map[nnodes];
    int elemnodes;

    for ( int j = 0; j < nnodes; j++ ) {
        map [ j ] = 0;
    }

    for ( auto &elem : d->giveElements() ) {
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }

        elemnodes = elem->giveNumberOfNodes();
        for ( int ielemnode = 1; ielemnode <= elemnodes; ielemnode++ ) {
             map [ elem->giveNode(ielemnode)->giveNumber() - 1 ] = 1;
        }
    }

    for ( int j = 0; j < nnodes; j++ ) {
        rbrnodes += map [ j ];
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


    int nnodes = domain->giveNumberOfDofManagers();
    int elemNodes;
    int currOffset = offset + 1;

    regionNodalNumbers.resize(nnodes);
    regionNodalNumbers.zero();
    regionDofMans = 0;

    for ( auto &element : domain->giveElements() ) {

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        elemNodes = element->giveNumberOfNodes();
        //  elemSides = element->giveNumberOfSides();

        // determine local region node numbering
        for ( int elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            int node = element->giveNode(elementNode)->giveNumber();
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
        for ( int i = 1; i <= nnodes; i++ ) {
            if ( regionNodalNumbers.at(i) ) {
                regionNodalNumbers.at(i) = currOffset++;
                answer.at( regionNodalNumbers.at(i) ) = i;
            }
        }

        regionNodalNumbers = answer;
    } else {
        for ( int i = 1; i <= nnodes; i++ ) {
            if ( regionNodalNumbers.at(i) ) {
                regionNodalNumbers.at(i) = currOffset++;
            }
        }
    }

    return 1;
}

//keyword "vars" in OOFEM input file
void
VTKExportModule :: exportIntVarAs(InternalStateType valID, InternalStateValueType type, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int ireg;
    int nnodes = d->giveNumberOfDofManagers(), inode;
    int j, jsize;
    FloatArray iVal(3);
    FloatMatrix t(3, 3);
    const FloatArray *val = NULL;


    // create a new set containing all elements
    Set elemSet(0, d);
    elemSet.addAllElements();

    this->giveSmoother();

    int nindx = giveInternalStateTypeSize(type);

    for ( int indx = 1; indx <= nindx; indx++ ) {
        // print header
        if ( type == ISVT_SCALAR ) {
            fprintf( stream, "SCALARS %s double 1\n", __InternalStateTypeToString(valID) );
        } else if ( type == ISVT_VECTOR ) {
            fprintf( stream, "VECTORS %s double\n", __InternalStateTypeToString(valID) );
        } else if ( ( type == ISVT_TENSOR_S3 ) || ( type == ISVT_TENSOR_S3E ) ) {
            fprintf( stream, "TENSORS %s double\n", __InternalStateTypeToString(valID) );
        } else if ( type == ISVT_TENSOR_G ) {
            fprintf(stream, "SCALARS %s_%d double 1\n", __InternalStateTypeToString(valID), indx);
        } else {
            fprintf( stderr, "VTKExportModule :: exportIntVarAs: unsupported variable type %s\n", __InternalStateTypeToString(valID) );
        }

        if ( ( type == ISVT_SCALAR ) || ( type == ISVT_TENSOR_G ) ) {
            fprintf(stream, "LOOKUP_TABLE default\n");
        }

        if ( !( ( valID == IST_DisplacementVector ) || ( valID == IST_MaterialInterfaceVal ) ) ) {
            this->smoother->recoverValues(elemSet, valID, tStep);
        }

        IntArray regionNodalNumbers(nnodes);
        int regionDofMans = 0, offset = 0;
        ireg = -1;
        int defaultSize = 0;

        this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);
        if ( !( ( valID == IST_DisplacementVector ) || ( valID == IST_MaterialInterfaceVal ) ) ) {
            // assemble local->global map
            defaultSize = giveInternalStateTypeSize(type);
        } else {
            regionDofMans = nnodes;
        }

        for ( inode = 1; inode <= regionDofMans; inode++ ) {
            if ( valID == IST_DisplacementVector ) {
                iVal.resize(3);
                val = & iVal;
                for ( j = 1; j <= 3; j++ ) {
                    iVal.at(j) = d->giveNode( regionNodalNumbers.at(inode) )->giveUpdatedCoordinate(j, tStep, 1.0) -
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
                this->smoother->giveNodalVector( val, regionNodalNumbers.at(inode) );
            }

            if ( val == NULL ) {
                iVal.resize(defaultSize);
                iVal.zero();
                val = & iVal;
                //OOFEM_ERROR("internal error: invalid dofman data");
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
            t.zero();
            for ( int ii = 1; ii <= 6; ii++ ) {
                if ( ii == 1 ) {
                    t.at(1, 1) = val->at(ii);
                } else if ( ii == 2 ) {
                    t.at(2, 2) = val->at(ii);
                } else if ( ii == 3 ) {
                    t.at(3, 3) = val->at(ii);
                } else if ( ii == 4 ) {
                    t.at(2, 3) = val->at(ii);
                    t.at(3, 2) = val->at(ii);
                } else if ( ii == 5 ) {
                    t.at(1, 3) = val->at(ii);
                    t.at(3, 1) = val->at(ii);
                } else if ( ii == 6 ) {
                    t.at(1, 2) = val->at(ii);
                    t.at(2, 1) = val->at(ii);
                }
            }

            for ( int ii = 1; ii <= 3; ii++ ) {
                for ( int jj = 1; jj <= 3; jj++ ) {
                    fprintf( stream, "%e ", t.at(ii, jj) );
                }

                fprintf(stream, "\n");
            }
        } else if ( type == ISVT_TENSOR_G ) { // export general tensor values as scalars
            fprintf( stream, "%e ", val->at(indx) );
        }

        fprintf(stream, "\n");
    }
}


NodalRecoveryModel *
VTKExportModule :: giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->smoother == NULL ) {
        this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
    }

    return this->smoother;
}


//keyword "primvars" in OOFEM input file
void
VTKExportModule :: exportPrimaryVars(FILE *stream, TimeStep *tStep)
{
    for ( auto &primvar : primaryVarsToExport ) {
        UnknownType type = ( UnknownType ) primvar;
        this->exportPrimVarAs(type, stream, tStep);
    }
}

//keyword "primvars" in OOFEM input file
void
VTKExportModule :: exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int ireg;
    int nnodes = d->giveNumberOfDofManagers();
    int j, jsize;
    FloatArray iVal, iValLCS;
    IntArray dofIDMask;
    InternalStateValueType type = ISVT_UNDEFINED;
    int nScalarComp = 1;

    if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
        type = ISVT_VECTOR;
    } else if ( ( valID == FluxVector ) || ( valID == PressureVector ) || ( valID == Temperature ) || ( valID == Humidity )) {
        type = ISVT_SCALAR;
        //nScalarComp = d->giveNumberOfDefaultNodeDofs();
    } else {
        OOFEM_ERROR("unsupported UnknownType (%s)", __UnknownTypeToString(valID) );
    }

    // print header
    if ( type == ISVT_SCALAR ) {
        fprintf(stream, "SCALARS %s double %d\n",  __UnknownTypeToString(valID), nScalarComp);
    } else if ( type == ISVT_VECTOR ) {
        fprintf( stream, "VECTORS %s double\n", __UnknownTypeToString(valID) );
    } else {
        fprintf(stderr, "exportPrimVarAs: unsupported variable type\n");
    }

    if ( type == ISVT_SCALAR ) {
        fprintf(stream, "LOOKUP_TABLE default\n");
    }

    DofIDItem id;

    IntArray regionNodalNumbers(nnodes);
    int regionDofMans = 0, offset = 0;
    ireg = -1;


    // assemble local->global map
    this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, offset, d, ireg, 1);

    // Used for special nodes/elements with mixed number of dofs.
    this->giveSmoother();

    for ( int inode = 1; inode <= regionDofMans; inode++ ) {
        DofManager *dman = d->giveNode( regionNodalNumbers.at(inode) );

        if ( ( valID == DisplacementVector ) || ( valID == EigenVector ) || ( valID == VelocityVector ) ) {
            iVal.resize(3);
            iVal.zero();

            for ( Dof *dof: *dman ) {
                id = dof->giveDofID();
                if ( ( id == V_u ) || ( id == D_u ) ) {
                    iVal.at(1) = dof->giveUnknown(VM_Total, tStep);
                } else if ( ( id == V_v ) || ( id == D_v ) ) {
                    iVal.at(2) = dof->giveUnknown(VM_Total, tStep);
                } else if ( ( id == V_w ) || ( id == D_w ) ) {
                    iVal.at(3) = dof->giveUnknown(VM_Total, tStep);
                }
            }
        } else if ( (valID == FluxVector) || (valID == Humidity) ) {
            iVal.resize(1);

            for ( Dof *dof: *dman ) {
                id = dof->giveDofID();
                if ( id == C_1 ) {
                    iVal.at(1) = dof->giveUnknown(VM_Total, tStep);
                }
            }
        } else if ( valID == Temperature ) {
            iVal.resize(1);

            for ( Dof *dof: *dman ) {
                id = dof->giveDofID();
                if ( id == T_f ) {
                    iVal.at(1) = dof->giveUnknown(VM_Total, tStep);
                }
            }
        } else if ( valID == PressureVector ) {
            dofIDMask = {P_f};
            this->getDofManPrimaryVariable(iVal, dman, dofIDMask, VM_Total, tStep, IST_Pressure);
        } else {
            OOFEM_ERROR("unsupported unknownType (%s)", __UnknownTypeToString(valID) );
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

    fprintf(stream, "\n");
}

void
VTKExportModule :: getDofManPrimaryVariable(FloatArray &answer, DofManager *dman, IntArray &dofIDMask,
                                            ValueModeType mode, TimeStep *tStep, InternalStateType iType)
{
    int size = dofIDMask.giveSize();
    const FloatArray *recoveredVal;
    answer.resize(size);
    // all values zero by default
    answer.zero();



    for ( int j = 1; j <= size; j++ ) {
        if ( dman->hasDofID( ( DofIDItem ) dofIDMask.at(j) ) ) {
            // primary variable available directly in dof manager
            answer.at(j) = dman->giveDofWithID( dofIDMask.at(j) )->giveUnknown(VM_Total, tStep);
        } else if ( iType != IST_Undefined ) {
            // ok primary variable not directly available
            // but equivalent InternalStateType provided
            // in this case use smoother to recover nodal value

            // create a new set containing all elements
            Set elemSet( 0, dman->giveDomain() );
            elemSet.addAllElements();

            this->giveSmoother()->recoverValues(elemSet, iType, tStep); /// recover values if not done before
            this->giveSmoother()->giveNodalVector( recoveredVal, dman->giveNumber() );
            // here we have a lack of information about how to convert recoveredVal to response
            // if the size is compatible we accept it, otherwise an error is thrown
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING("recovered variable size mismatch for %d", iType);
                answer.at(j) = 0.0;
            }
        }
    }
}
} // end namespace oofem
