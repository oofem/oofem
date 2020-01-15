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

#include "quasicontinuumvtkxmlexportmodule.h"
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


namespace oofem {
REGISTER_ExportModule(QuasicontinuumVTKXMLExportModule)


QuasicontinuumVTKXMLExportModule :: QuasicontinuumVTKXMLExportModule(int n, EngngModel *e) : VTKXMLExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{}


QuasicontinuumVTKXMLExportModule :: ~QuasicontinuumVTKXMLExportModule()
{}


void
QuasicontinuumVTKXMLExportModule :: initializeFrom(InputRecord &ir)
{
    VTKXMLExportModule :: initializeFrom(ir);

    deactivatedElementsExportFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, deactivatedElementsExportFlag, _IFT_QuasicontinuumVTKXMLExportModule_ExportDeactivatedElements); // Macro
}



void
QuasicontinuumVTKXMLExportModule :: setupVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep, int region)
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
            const auto &coords = d->giveNode( mapL2G.at(inode) )->giveCoordinates();
            vtkPiece.setNodeCoords(inode, coords);
        }


        //-------------------------------------------
        // Export all the cell data for the piece
        //-------------------------------------------
        IntArray cellNodes;
        vtkPiece.setNumberOfCells(numRegionEl);

        int offset = 0;
        int cellNum = 0;
        IntArray elems = this->giveRegionSet(region)->giveElementList();
        for ( int ei = 1; ei <= elems.giveSize(); ei++ ) {
            int elNum = elems.at(ei);
            elem = d->giveElement(elNum);


            if ( !deactivatedElementsExportFlag ) {
                // Skip elements that are deactivated by QC
                if ( !emodel->isElementActivated(elem) ) {
                    continue;
                }
            }

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

            cellNum++;

            // Set the connectivity
            this->giveElementCell(cellNodes, elem);  // node numbering of the cell with according to the VTK format

            // Map from global to local node numbers for the current piece
            int numElNodes = cellNodes.giveSize();
            IntArray connectivity(numElNodes);
            for ( int i = 1; i <= numElNodes; i++ ) {
                connectivity.at(i) = mapG2L.at( cellNodes.at(i) );
            }

            vtkPiece.setConnectivity(cellNum, connectivity);

            vtkPiece.setCellType( cellNum, this->giveCellType(elem) ); // VTK cell type

            offset += numElNodes;
            vtkPiece.setOffset(cellNum, offset);
        }


        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(vtkPiece, mapG2L, mapL2G, region, tStep);
        this->exportIntVars(vtkPiece, mapG2L, mapL2G, region, tStep);
        this->exportExternalForces(vtkPiece, mapG2L, mapL2G, region, tStep);

        this->exportCellVars(vtkPiece, elems, tStep);
    } // end of default piece for simple geometry elements
}





int
QuasicontinuumVTKXMLExportModule :: initRegionNodeNumbering(IntArray &regionG2LNodalNumbers,
                                                            IntArray &regionL2GNodalNumbers,
                                                            int &regionDofMans,
                                                            int &regionSingleCells,
                                                            Domain *domain, TimeStep *tStep, int reg)
{
    // regionG2LNodalNumbers is array with mapping from global numbering to local region numbering.
    // The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).

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


        if ( !deactivatedElementsExportFlag ) {
            if ( !domain->giveEngngModel()->isElementActivated(element) ) {
                continue;                                    // skip elements deactovated by QC
            }
        }

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
            regionL2GNodalNumbers.at( regionG2LNodalNumbers.at(i) ) = i;
        }
    }

    return 1;
}
} // end namespace oofem
