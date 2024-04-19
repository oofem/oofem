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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "vtkxmlperiodicexportmodule.h"
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
#include "dof.h"

#ifdef __SM_MODULE
 #include "../sm/Elements/LatticeElements/lattice2dboundary.h"
 #include "../sm/Elements/LatticeElements/lattice3dboundary.h"
 #include "../sm/Elements/LatticeElements/latticelink3dboundary.h"
 #include "../sm/Elements/3D/ltrspaceboundary.h"
 #include "../sm/Elements/Beams/libeam3dboundary.h"
#endif

namespace oofem {
REGISTER_ExportModule(VTKXMLPeriodicExportModule)

VTKXMLPeriodicExportModule :: VTKXMLPeriodicExportModule(int n, EngngModel *e) : VTKXMLExportModule(n, e)
{
    //Find out the element type and elemNodes
}


VTKXMLPeriodicExportModule :: ~VTKXMLPeriodicExportModule()
{}


void
VTKXMLPeriodicExportModule :: initializeFrom(InputRecord &ir)
{
    VTKXMLExportModule :: initializeFrom(ir);
}


void
VTKXMLPeriodicExportModule :: giveSwitches(IntArray &answer, int location) {
    answer.resize(3);
    answer.zero();
    int counter = 1;
    for ( int x = -1; x <  2; x++ ) {
        for ( int y = -1; y <  2; y++ ) {
            for ( int z = -1; z <  2; z++ ) {
                if ( !( z == 0 && y == 0 && x == 0 ) ) {
                    if ( counter == location ) {
                        answer(0) = x;
                        answer(1) = y;
                        answer(2) = z;
                    }
                    counter++;
                }
            }
        }
    }

    return;
}


void
VTKXMLPeriodicExportModule :: setupVTKPiece(ExportRegion &vtkPiece, TimeStep *tStep, Set& region)
{
    // Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.

    Domain *d  = emodel->giveDomain(1);
    Element *elem;

    int nnodes = d->giveNumberOfDofManagers();

    this->giveSmoother(); // make sure smoother is created
    
    // Assemble local->global and global->local region map and get number of
    // single cells to process, the composite cells exported individually.
    this->initRegionNodeNumbering(vtkPiece, d, tStep, region);
    const int numNodes = vtkPiece.giveNumberOfNodes();
    const int numRegionEl = vtkPiece.giveNumberOfCells();
    const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();


    if ( numNodes > 0 && numRegionEl > 0 ) {
        // Export nodes as vtk vertices
        vtkPiece.setNumberOfNodes(numNodes);
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            if ( mapL2G.at(inode) <= nnodes && mapL2G.at(inode) != 0 ) { //DofManagers in domain (input file)
                const auto &coords = d->giveNode(mapL2G.at(inode) )->giveCoordinates();
                vtkPiece.setNodeCoords(inode, coords);
            } else if ( mapL2G.at(inode) > nnodes && mapL2G.at(inode) <= numNodes && mapL2G.at(inode) != 0 ) { //extra image nodes
                FloatArray helpArray(3);
                helpArray.at(1) = uniqueNodeTable.at(mapL2G.at(inode), 1);
                helpArray.at(2) = uniqueNodeTable.at(mapL2G.at(inode), 2);
                helpArray.at(3) = uniqueNodeTable.at(mapL2G.at(inode), 3);
                vtkPiece.setNodeCoords(inode, helpArray);
            } else {
                FloatArray helpArray(3);
                helpArray.zero();
                vtkPiece.setNodeCoords(inode, helpArray);
            }
        }


        //-------------------------------------------
        // Export all the cell data for the piece
        //-------------------------------------------
        IntArray cellNodes;
        vtkPiece.setNumberOfCells(numRegionEl);

        int offset = 0;
        int cellNum = 0;
        IntArray elems = region.giveElementList();
        int helpCounter = 0;
        for ( int ei = 1; ei <= elems.giveSize(); ei++ ) {
            int elNum = elems.at(ei);
            elem = d->giveElement(elNum);

            // Skip elements that:
            // are inactivated or of composite type ( these are exported individually later)
            if ( this->isElementComposite(elem) || !elem->isActivated(tStep) ) {
                continue;
            }

            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            cellNum++;

            // Set the connectivity
#ifdef __SM_MODULE
            if ( dynamic_cast< LTRSpaceBoundary * >( elem ) ) {
                cellNodes.resize(4);
                IntArray loc = elem->giveLocation();
                for ( int ielnode = 1; ielnode <= 4; ielnode++ ) {
                    if ( loc.at(ielnode) != 0 ) {
                        helpCounter++;
                        cellNodes.at(ielnode) = regionToUniqueMap.at(nnodes + helpCounter);
                    } else {
                        cellNodes.at(ielnode) = elem->giveNode(ielnode)->giveNumber();
                    }
                }
            } else if ( dynamic_cast< LIBeam3dBoundary * >( elem ) || dynamic_cast< Lattice3dBoundary * >( elem ) || dynamic_cast< LatticeLink3dBoundary * >( elem ) || dynamic_cast< Lattice2dBoundary * >( elem ) ) {
                cellNodes.resize(2);
                IntArray loc = elem->giveLocation();
                for ( int ielnode = 1; ielnode <= 2; ielnode++ ) {
                    if ( loc.at(ielnode) != 0 ) {
                        helpCounter++;
                        cellNodes.at(ielnode) = regionToUniqueMap.at(nnodes + helpCounter);
                    } else {
                        cellNodes.at(ielnode) = elem->giveNode(ielnode)->giveNumber();
                    }
                }
            } else {  //Standard case
#endif
            this->giveElementCell(cellNodes, elem);
#ifdef __SM_MODULE
        }
#endif
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
        this->exportPrimaryVars(vtkPiece, region, primaryVarsToExport, *primVarSmoother, tStep);
        this->exportIntVars (vtkPiece, region, internalVarsToExport, *smoother, tStep);

        this->exportCellVars(vtkPiece, region, cellVarsToExport, tStep);

    } // end of default piece for simple geometry elements
}

int
VTKXMLPeriodicExportModule :: initRegionNodeNumbering(ExportRegion& vtkPiece,
                                                      Domain *domain, TimeStep *tStep, Set& region)
{
    int nnodes = domain->giveNumberOfDofManagers();
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    int regionDofMans = 0;
    int regionSingleCells = 0;
    IntArray& regionG2LNodalNumbers = vtkPiece.getMapG2L();
    IntArray& regionL2GNodalNumbers = vtkPiece.getMapL2G();


    IntArray elements = region.giveElementList();

    int extraNodes = 0.;
    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        element = domain->giveElement(elements.at(ie) );

#ifdef __SM_MODULE
        if ( dynamic_cast< LTRSpaceBoundary * >( element ) ) {
            IntArray loc = element->giveLocation();
            for ( int ielnode = 1; ielnode <= 4; ielnode++ ) {
                if ( loc.at(ielnode) != 0 ) {
                    extraNodes++;
                }
            }
        } else if ( dynamic_cast< LIBeam3dBoundary * >( element ) || dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) {
            IntArray loc = element->giveLocation();
            for ( int ielnode = 1; ielnode <= 2; ielnode++ ) {
                if ( loc.at(ielnode) != 0 ) {
                    extraNodes++;
                }
            }
        }
#endif
    }

    regionG2LNodalNumbers.resize(nnodes + extraNodes);
    regionG2LNodalNumbers.zero();

    periodicMap.resize(nnodes + extraNodes);
    regionToUniqueMap.resize(nnodes + extraNodes);
    locationMap.resize(nnodes + extraNodes);

    uniqueNodeTable.resize(nnodes + extraNodes, 3);
    uniqueNodeTable.zero();

    int totalNodes = 0;
    int uniqueNodes = 0;

    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        int ielem = elements.at(ie);
        element = domain->giveElement(ielem);

        if ( this->isElementComposite(element) ) {
            continue;                                    // composite cells exported individually
        }

        if ( !element->isActivated(tStep) ) {                    //skip inactivated elements
            continue;
        }

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        regionSingleCells++;
        elemNodes = element->giveNumberOfNodes();

#ifdef __SM_MODULE
        if ( dynamic_cast< LTRSpaceBoundary * >( element ) ) {
            elemNodes = 4;
        } else if ( dynamic_cast< LIBeam3dBoundary * >( element ) || dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) {
            elemNodes = 2;
        }
#endif

        for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
#ifdef __SM_MODULE
            if ( dynamic_cast< LTRSpaceBoundary * >( element ) ) {
                IntArray loc = element->giveLocation();
                FloatArray nodeCoords;
                if ( loc.at(elementNode) != 0 ) { //boundary element with mirrored node
                    totalNodes++;
                    regionDofMans++;
                    regionG2LNodalNumbers.at(nnodes + totalNodes) = 1;
                    node = element->giveNode(elementNode)->giveNumber();
                    if ( regionG2LNodalNumbers.at(node) == 0 ) { // assign new number
                        /* mark for assignment. This is done later, as it allows to preserve
                         * natural node numbering.
                         */

                        regionG2LNodalNumbers.at(node) = 1;
                        regionDofMans++;
                    }
                    if ( regionG2LNodalNumbers.at(nnodes) == 0  ) {
                        regionG2LNodalNumbers.at(nnodes) = 1;
                        regionDofMans++;
                    }

                    periodicMap.at(nnodes + totalNodes) = node;
                    locationMap.at(nnodes + totalNodes) = loc.at(elementNode);

                    element->recalculateCoordinates(elementNode, nodeCoords);

                    //get only the unique nodes
                    int repeatFlag = 0;
                    for ( int j = nnodes + 1; j <= nnodes + uniqueNodes; j++ ) {
                        double dx = fabs(uniqueNodeTable.at(j, 1) - nodeCoords.at(1) );
                        double dy = fabs(uniqueNodeTable.at(j, 2) - nodeCoords.at(2) );
                        double dz = fabs(uniqueNodeTable.at(j, 3) - nodeCoords.at(3) );
                        if ( dx < 1e-9 && dy < 1e-9 && dz < 1e-9 ) {//node already present
                            repeatFlag++;
                            regionToUniqueMap.at(nnodes + totalNodes) = j;
                        }
                    }

                    if ( repeatFlag == 0 ) {//new unique node
                        uniqueNodes++;
                        uniqueNodeTable.at(nnodes + uniqueNodes, 1) = nodeCoords.at(1);
                        uniqueNodeTable.at(nnodes + uniqueNodes, 2) = nodeCoords.at(2);
                        uniqueNodeTable.at(nnodes + uniqueNodes, 3) = nodeCoords.at(3);
                        regionToUniqueMap.at(nnodes + totalNodes) = nnodes + uniqueNodes;
                    }
                } else { //boundary element but node not mirrored
                    node = element->giveNode(elementNode)->giveNumber();
                    if ( regionG2LNodalNumbers.at(node) == 0 ) { // assign new number
                        /* mark for assignment. This is done later, as it allows to preserve
                         * natural node numbering.
                         */
                        regionG2LNodalNumbers.at(node) = 1;
                        regionDofMans++;
                    }
                }
            } else if ( dynamic_cast< LIBeam3dBoundary * >( element ) || dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) { //beam elements - only unique nodes
                IntArray loc = element->giveLocation();
                FloatArray nodeCoords;
                if ( loc.at(elementNode) != 0 ) { //boundary element with mirrored node
                    totalNodes++;
                    regionDofMans++;
                    regionG2LNodalNumbers.at(nnodes + totalNodes) = 1;
                    node = element->giveNode(elementNode)->giveNumber();
                    if ( regionG2LNodalNumbers.at(node) == 0 ) { // assign new number
                        /* mark for assignment. This is done later, as it allows to preserve
                         * natural node numbering.
                         */

                        regionG2LNodalNumbers.at(node) = 1;
                        regionDofMans++;
                    }
                    if ( regionG2LNodalNumbers.at(nnodes) == 0  ) {
                        regionG2LNodalNumbers.at(nnodes) = 1;
                        regionDofMans++;
                    }

                    periodicMap.at(nnodes + totalNodes) = node;
                    locationMap.at(nnodes + totalNodes) = loc.at(elementNode);

                    element->recalculateCoordinates(elementNode, nodeCoords);

                    //get only the unique nodes
                    int repeatFlag = 0;
                    for ( int j = nnodes + 1; j <= nnodes + uniqueNodes; j++ ) {
                        double dx = fabs(uniqueNodeTable.at(j, 1) - nodeCoords.at(1) );
                        double dy = fabs(uniqueNodeTable.at(j, 2) - nodeCoords.at(2) );
                        double dz = fabs(uniqueNodeTable.at(j, 3) - nodeCoords.at(3) );
                        if ( dx < 1e-9 && dy < 1e-9 && dz < 1e-9 ) {//node already present
                            repeatFlag++;
                            regionToUniqueMap.at(nnodes + totalNodes) = j;
                        }
                    }

                    if ( repeatFlag == 0 ) {//new unique node
                        uniqueNodes++;
                        uniqueNodeTable.at(nnodes + uniqueNodes, 1) = nodeCoords.at(1);
                        uniqueNodeTable.at(nnodes + uniqueNodes, 2) = nodeCoords.at(2);
                        uniqueNodeTable.at(nnodes + uniqueNodes, 3) = nodeCoords.at(3);
                        regionToUniqueMap.at(nnodes + totalNodes) = nnodes + uniqueNodes;
                    }
                }
            } else {   //regular element
#endif
            node = element->giveNode(elementNode)->giveNumber();
            if ( regionG2LNodalNumbers.at(node) == 0 ) {     // assign new number
                /* mark for assignment. This is done later, as it allows to preserve
                 * natural node numbering.
                 */
                regionG2LNodalNumbers.at(node) = 1;
                regionDofMans++;
            }
#ifdef __SM_MODULE
        }
#endif
        }
    }

    uniqueNodeTable.resizeWithData(nnodes + uniqueNodes, 3);
    regionDofMans = nnodes + uniqueNodes;
    regionG2LNodalNumbers.resizeWithValues(regionDofMans);
    regionL2GNodalNumbers.resize(regionDofMans);

    vtkPiece.setNumberOfNodes(regionDofMans);
    vtkPiece.setNumberOfCells(regionSingleCells);

    for ( int i = 1; i <= regionDofMans; i++ ) {
        if ( regionG2LNodalNumbers.at(i) ) {
            regionG2LNodalNumbers.at(i) = currOffset++;
            regionL2GNodalNumbers.at(regionG2LNodalNumbers.at(i) ) = i;
        }
    }
    return 1;
}


void VTKXMLPeriodicExportModule :: exportPrimaryVars(ExportRegion &vtkPiece, Set& region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep) 
{
    Domain *d = emodel->giveDomain(1);
    int nnodes = d->giveNumberOfDofManagers();
    FloatArray valueArray;
    smoother.clear(); // Makes sure primary smoother is up-to-date with potentially new mesh.

    //const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport, mapL2G.giveSize() );

    //Get the macroscopic field (deformation gradients, curvatures etc.)
    DofManager *controlNode = d->giveNode(nnodes);   //assuming the control node is last
    IntArray dofIdArray;
    controlNode->giveCompleteMasterDofIDArray(dofIdArray);
    std::vector<int> dofIDVector;
    dofIDVector.assign(dofIdArray.begin(), dofIdArray.end());

    FloatArray macroField(controlNode->giveNumberOfDofs() );
    for ( int j = 1; j <= controlNode->giveNumberOfDofs(); j++ ) {
        macroField.at(j) = controlNode->giveDofWithID(dofIdArray.at(j) )->giveUnknown(VM_Total, tStep);
    }

    std::vector<int> macroTrussIDs = { E_xx };
    std::vector<int> macroMembraneIDs = { E_xx, E_yy, E_xy, E_yx };
    std::vector<int> macroBeamIDs = { E_xx, E_zx, K_xx };
    std::vector<int> macroPlateIDs = { E_xx, E_yy, E_zy, E_zx, E_xy, E_yx, K_xx, K_yy, K_xy, K_yx };
    std::vector<int> macro3DVoigtIDs = { E_xx, E_yy, E_zz, G_yz, G_xz, G_xy };
    std::vector<int> macro3DIDs = { E_xx, E_yy, E_zz, E_yz, E_zy, E_xz, E_zx, E_xy, E_yx };
    std::vector<int> macro2DIDs = { E_xx, E_yy, G_xy  };

    //Get unit cell size
    const auto unitCellSize = controlNode->giveCoordinates();

    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

        for ( int inode = 1; inode <= mapL2G.giveSize(); inode++ ) {
            if ( inode <= nnodes && mapL2G.at(inode) <= nnodes && mapL2G.at(inode) != 0 ) { //no special treatment for master nodes
                DofManager *dman = d->giveNode(mapL2G.at(inode) );

                this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region, smoother);
                vtkPiece.setPrimaryVarInNode(type, inode, std :: move(valueArray) );
            } else { //special treatment for image nodes
                //find the periodic node, enough to find the first occurrence
                int pos = 0;
                if ( mapL2G.at(inode) != 0 ) {
                    pos = regionToUniqueMap.findFirstIndexOf(mapL2G.at(inode) );
                }
                if ( pos ) {
                    DofManager *dman = d->giveNode(periodicMap.at(pos) );
                    IntArray switches;
                    giveSwitches(switches, locationMap.at(pos) );
                    //get the master unknown
                    FloatArray helpArray;
                    this->getNodalVariableFromPrimaryField(helpArray, dman, tStep, type, region, smoother);
                    //recalculate the image unknown
                    if ( type == DisplacementVector ) {
                        if ( dofIDVector == macro3DIDs ) { //Macroscale: 3D SOLID, LTRSpaceBoundary
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(8) + unitCellSize.at(3) * switches.at(3) * macroField.at(6);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(1) * switches.at(1) * macroField.at(9) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(2) + unitCellSize.at(3) * switches.at(3) * macroField.at(4);
                            valueArray.at(3) = helpArray.at(3) + unitCellSize.at(1) * switches.at(1) * macroField.at(7) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(5) + unitCellSize.at(3) * switches.at(3) * macroField.at(3);
                        } else if ( dofIDVector == macroTrussIDs ) { //Macroscale: TRUSS
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1);
                            valueArray.at(2) = helpArray.at(2);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIDVector == macroMembraneIDs ) { //Macroscale: 2D MEMBRANE, LTRSpaceBoundaryMembrane
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(3);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(1) * switches.at(1) * macroField.at(4) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(2);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIDVector == macro2DIDs ) { //Macroscale: 2D plane stress
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(2) * switches.at(2) * macroField.at(2) + unitCellSize.at(1) * switches.at(1) * macroField.at(3);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIDVector == macroBeamIDs )  { //Macroscale: 2D BEAM, LTRSpaceBoundaryBeam
                                valueArray.resize(helpArray.giveSize() );
                                valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) -
                                                   dman->giveCoordinate(3) * unitCellSize.at(1) * switches.at(1) * macroField.at(3);
                                valueArray.at(2) = helpArray.at(2);
                                valueArray.at(3) = helpArray.at(3) + unitCellSize.at(1) * switches.at(1) * macroField.at(2);
                        } else if ( dofIDVector == macroPlateIDs ) { //Macroscale: 2D PLATE, LTRSpaceBoundaryPlate
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(5) -
                                               dman->giveCoordinate(3) * unitCellSize.at(1) * switches.at(1) * macroField.at(7) -
                                               dman->giveCoordinate(3) * unitCellSize.at(2) * switches.at(2) * macroField.at(9);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(1) * switches.at(1) * macroField.at(6) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(2) -
                                               dman->giveCoordinate(3) * unitCellSize.at(2) * switches.at(2) * macroField.at(8) -
                                               dman->giveCoordinate(3) * unitCellSize.at(1) * switches.at(1) * macroField.at(10);
                            valueArray.at(3) = helpArray.at(3) + unitCellSize.at(1) * switches.at(1) * macroField.at(4) +
                                               unitCellSize.at(2) * switches.at(2) * macroField.at(3);
                        } else if ( dofIDVector == macro3DVoigtIDs ) { //Macroscale: 3D SOLID, LTRSpaceBoundaryVoigt, Lattice3dBoundary
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) +
                                               unitCellSize.at(3) * switches.at(3) * macroField.at(5) + unitCellSize.at(2) * switches.at(2) * macroField.at(6);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(2) * switches.at(2) * macroField.at(2) +
                                               unitCellSize.at(3) * switches.at(3) * macroField.at(4);
                            valueArray.at(3) = helpArray.at(3) + unitCellSize.at(3) * switches.at(3) * macroField.at(3);
                        } else {
                            OOFEM_ERROR("Unknown element type\n");
                        }
                    }
                } else {
                    valueArray.resize(3);
                }

                vtkPiece.setPrimaryVarInNode(type, inode, std :: move(valueArray) );
            }
        }
    }
}

void 
VTKXMLPeriodicExportModule :: exportIntVars(ExportRegion &vtkPiece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nnodes = d->giveNumberOfDofManagers();
    InternalStateType isType;
    FloatArray answer;

    smoother.clear(); // Makes sure smoother is up-to-date with potentially new mesh.
    //const IntArray& mapG2L = vtkPiece.getMapG2L();
    const IntArray& mapL2G = vtkPiece.getMapL2G();

    // Export of Internal State Type fields
    vtkPiece.setNumberOfInternalVarsToExport(internalVarsToExport, mapL2G.giveSize() );
    for ( int field = 1; field <= internalVarsToExport.giveSize(); field++ ) {
        isType = ( InternalStateType ) internalVarsToExport.at(field);

        for ( int nodeNum = 1; nodeNum <= mapL2G.giveSize(); nodeNum++ ) {
            if ( nodeNum <= nnodes && mapL2G.at(nodeNum) <= nnodes && mapL2G.at(nodeNum) != 0 ) { //no special treatment for master nodes
                Node *node = d->giveNode(mapL2G.at(nodeNum) );
                this->getNodalVariableFromIS(answer, node, tStep, isType, region, smoother);
                vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
            } else { //special treatment for image nodes
                //find the periodic node, enough to find the first occurrence
                int pos = 0;
                if ( mapL2G.at(nodeNum) != 0 ) {
                    pos = regionToUniqueMap.findFirstIndexOf(mapL2G.at(nodeNum) );
                }

                if ( pos ) {
                    Node *node = d->giveNode(periodicMap.at(pos) );
                    this->getNodalVariableFromIS(answer, node, tStep, isType, region, smoother);
                    vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
                } else { //fill with zeroes
                    InternalStateValueType valType = giveInternalStateValueType(isType);
                    int ncomponents = giveInternalStateTypeSize(valType);
                    answer.resize(ncomponents);

                    if ( isType == IST_BeamForceMomentTensor ) {
                        answer.resize(6);
                    }

                    answer.zero();
                    vtkPiece.setInternalVarInNode(isType, nodeNum, answer);
                }
            }
        }
    }
}
} // end namespace oofem
