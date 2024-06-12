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

#include "vtkbaseexportmodule.h"
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

#include <string>
#include <sstream>
#include <ctime>

namespace oofem {


// Beginning of VTKBaseExportModule

IntArray VTKBaseExportModule::redToFull = {
    1, 5, 9, 8, 7, 4, 6, 3, 2
};                                                                      //position of xx, yy, zz, yz, xz, xy in tensor

VTKBaseExportModule::VTKBaseExportModule(int n, EngngModel *e) : ExportModule(n, e) 
{}

VTKBaseExportModule::~VTKBaseExportModule() 
{}

void
VTKBaseExportModule::initialize()
{
    ExportModule::initialize();
}

void
VTKBaseExportModule::terminate()
{}

void
VTKBaseExportModule::makeFullTensorForm(FloatArray &answer, const FloatArray &reducedForm, InternalStateValueType vtype)
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

int
VTKBaseExportModule::giveCellType(Element *elem)
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
    } else if ( elemGT == EGT_quad_1) {
        vtkCellType = 9;
    } else if ( elemGT == EGT_quad_21_interface ) {
        vtkCellType = 30;
    } else if ( elemGT == EGT_quad_2 ) {
        vtkCellType = 23;
    } else if ( elemGT == EGT_quad9_2 ) {
        vtkCellType = 23;
    } else if ( elemGT == EGT_hexa_1  || elemGT == EGT_quad_1_interface ) {
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

int VTKBaseExportModule::giveCellType(int num) {
  return giveCellType(this->emodel->giveDomain(1)->giveElement(num));
}


int
VTKBaseExportModule::giveNumberOfNodesPerCell(int cellType)
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
VTKBaseExportModule::giveElementCell(IntArray &answer, Element *elem)
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
            1, 2, 4, 3, 5, 6, 8, 7
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
VTKBaseExportModule::isElementComposite(Element *elem)
{
    return ( elem->giveGeometryType() == EGT_Composite );
}

void
VTKBaseExportModule::setupVTKPiece(ExportRegion &vtkPiece, TimeStep *tStep, Set &region)
{
    Domain *d  = emodel->giveDomain(1);
    Element *elem;

    // output nodes Region By Region
   
    
    // Assemble local->global and global->local region map and get number of
    // single cells to process, the composite cells exported individually.
    this->initRegionNodeNumbering(vtkPiece, d, tStep, region);
    const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();
    const int numNodes = vtkPiece.giveNumberOfNodes();
    const int numRegionEl = vtkPiece.giveNumberOfCells();
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

        int offset = 0;
        int cellNum = 0;
        IntArray elems = region.giveElementList();
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

            vtkPiece.getRegionCells().followedBy(elNum);

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


        
/*        
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
            std::string s = std::to_string(region.giveNumber());
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
*/
        
    } // end of default piece for simple geometry elements
}

//----------------------------------------------------
// Primary variables - readily available in the nodes
//----------------------------------------------------
void
VTKBaseExportModule::exportPrimaryVars(ExportRegion &vtkPiece, Set &region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;
    smoother.clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.

    //const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();

    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport, mapL2G.giveSize() );
    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

        for ( int inode = 1; inode <= mapL2G.giveSize(); inode++ ) {
            DofManager *dman = d->giveNode(mapL2G.at(inode) );

            this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region, smoother);
            vtkPiece.setPrimaryVarInNode(type, inode, std::move(valueArray) );
        }
    }
}


void
VTKBaseExportModule::getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, Set &region, NodalRecoveryModel& smoother)
{
    // This code is not perfect. It should be rewritten to handle all cases more gracefully.
    ///@todo This method needs to be cleaned up - maybe define the common vector types so
    /// certain dofid's are associated with them /JB

    IntArray dofIDMask(3);
    int size;
    const FloatArray *recoveredVal;

    InternalStateType iState = IST_DisplacementVector; // Shouldn't be necessary

    dofIDMask.clear();

    if ( (type == DisplacementVector) || (type == ResidualForce) ) {
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
    } else if ( type == MacroSlipVector ) {
        for ( Dof *dof : * dman ) {
            DofIDItem id = dof->giveDofID();
            if ( ( id == S_u ) || ( id == S_v ) || ( id == S_w ) ) {
                dofIDMask.followedBy(id);
            }
            answer.resize(3);
        }
        iState = IST_MacroSlipVector;
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
            smoother.recoverValues(region, iState, tStep);
            smoother.giveNodalVector(recoveredVal, dman->giveNumber() );
            if ( size == recoveredVal->giveSize() ) {
                answer.at(j) = recoveredVal->at(j);
            } else {
                OOFEM_WARNING("Recovered variable size mismatch for %d for id %d", type, id);
                answer.at(j) = 0.0;
            }
        } else if (type == ResidualForce) {
            answer.at(j) = dman->giveDofWithID(id)->giveUnknown(VM_Residual, tStep);
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
            smoother.recoverValues(region, iState, tStep);
            smoother.giveNodalVector(recoveredVal, dman->giveNumber() );
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


//----------------------------------------------------
// Internal variables and XFEM realted fields (keyword "vars" in OOFEM input file)
//----------------------------------------------------
void
VTKBaseExportModule::exportIntVars(ExportRegion &vtkPiece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    InternalStateType isType;
    FloatArray answer;
    //IntArray& mapG2L = vtkPiece.getMapG2L();
    IntArray& mapL2G = vtkPiece.getMapL2G();

    smoother.clear(); // Makes sure smoother is up-to-date with potentially new mesh.

    // Export of Internal State Type fields
    vtkPiece.setNumberOfInternalVarsToExport(internalVarsToExport, mapL2G.giveSize() );
    for ( int field = 1; field <= internalVarsToExport.giveSize(); field++ ) {
        isType = ( InternalStateType ) internalVarsToExport.at(field);

        for ( int nodeNum = 1; nodeNum <= mapL2G.giveSize(); nodeNum++ ) {
            Node *node = d->giveNode(mapL2G.at(nodeNum) );
            this->getNodalVariableFromIS(answer, node, tStep, isType, region, smoother);
            vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
        }
    }

    /*
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
    */
}


void
VTKBaseExportModule::getNodalVariableFromIS(FloatArray &answer, Node *node, TimeStep *tStep, InternalStateType type, Set& region, NodalRecoveryModel& smoother)
{
    // Recovers nodal values from Internal States defined in the integration points.
    // Should return an array with proper size supported by VTK (1, 3 or 9)
    // Domain *d = emodel->giveDomain(1);
    IntArray redIndx;

    if ( !( type == IST_DisplacementVector || type == IST_MaterialInterfaceVal  ) ) {
        smoother.recoverValues(region, type, tStep);
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
        int found = smoother.giveNodalVector(val, node->giveNumber() );
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



//----------------------------------------------------
// Cell vars
//----------------------------------------------------

void
VTKBaseExportModule::exportCellVars(ExportRegion &vtkPiece, Set& region, IntArray& cellVarsToExport, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray valueArray;
    const IntArray &elems = region.giveElementList();

    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport, elems.giveSize() );
    for ( int field = 1; field <= cellVarsToExport.giveSize(); field++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(field);

        for ( int subIndex = 1; subIndex <= elems.giveSize(); ++subIndex ) {
            Element *el = d->giveElement(elems.at(subIndex) );  ///@todo should be a pointer to an element in the region /JB
            if ( el->giveParallelMode() != Element_local ) {
                continue;
            }

            this->getCellVariableFromIS(valueArray, el, type, tStep);
            vtkPiece.setCellVar(type, subIndex, valueArray);
        }
    }
}


void
VTKBaseExportModule::getCellVariableFromIS(FloatArray &answer, Element *el, InternalStateType type, TimeStep *tStep)
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


int
VTKBaseExportModule::initRegionNodeNumbering(ExportRegion& piece, 
                                            Domain *domain, TimeStep *tStep, Set& region)
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

    IntArray &regionG2LNodalNumbers = piece.getMapG2L();
    IntArray &regionL2GNodalNumbers = piece.getMapL2G();

    regionG2LNodalNumbers.resize(nnodes);
    regionG2LNodalNumbers.zero();
    int regionDofMans = 0;
    int regionSingleCells = 0;

    const IntArray& elements = region.giveElementList();
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

    piece.setNumberOfCells(regionSingleCells);
    piece.setNumberOfNodes(regionDofMans);

    return 1;
}

void
VTKBaseExportModule::computeIPAverage(FloatArray &answer, IntegrationRule *iRule, Element *elem, InternalStateType isType, TimeStep *tStep)
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

//----------------------------------------------------
// Load vectors
//----------------------------------------------------
void
VTKBaseExportModule::exportExternalForces(ExportRegion &vtkPiece, Set& region, IntArray& externalForcesToExport, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    //this->givePrimVarSmoother()->clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.
    //IntArray& mapG2L = vtkPiece.getMapG2L();
    IntArray& mapL2G = vtkPiece.getMapL2G();

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

// End of VTKBaseExportModule implementation



} // end namespace oofem
