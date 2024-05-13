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

#include "vtkxmllatticeexportmodule.h"
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
 #include "../sm/Elements/LatticeElements/latticestructuralelement.h"
 #include "../sm/Elements/LatticeElements/lattice2dboundary.h"
 #include "../sm/Elements/LatticeElements/lattice3dboundary.h"
 #include "../sm/Elements/LatticeElements/latticelink3dboundary.h"
#endif

#ifdef __TM_MODULE
 #include "../tm/Elements/LatticeElements/latticetransportelement.h"
#endif


namespace oofem {
REGISTER_ExportModule(VTKXMLLatticeExportModule)

VTKXMLLatticeExportModule::VTKXMLLatticeExportModule(int n, EngngModel *e) : VTKXMLExportModule(n, e)
{
    //Find out the element type and elemNodes
}


VTKXMLLatticeExportModule::~VTKXMLLatticeExportModule()
{}


void
VTKXMLLatticeExportModule::initializeFrom(InputRecord &ir)
{
    VTKXMLExportModule::initializeFrom(ir);
    this->crossSectionExportFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, this->crossSectionExportFlag, _IFT_VTKXMLLatticeExportModule_cross);
}

std::string
VTKXMLLatticeExportModule::giveOutputFileNameCross(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".cross.vtu";
}

std::ofstream
VTKXMLLatticeExportModule::giveOutputStreamCross(TimeStep *tStep)
{
    std::string fileName = giveOutputFileNameCross(tStep);
    std::ofstream streamC;

    if ( pythonExport ) {
        streamC = std::ofstream(NULL_DEVICE);//do not write anything
    } else {
        streamC = std::ofstream(fileName);
    }

    if ( !streamC.good() ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    streamC.fill('0');//zero padding
    return streamC;
}

void
VTKXMLLatticeExportModule::giveSwitches(IntArray &answer, int location) {
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
VTKXMLLatticeExportModule::setupVTKPiece(ExportRegion &vtkPiece, TimeStep *tStep, Set &region)
{
    // Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.

    Domain *d  = emodel->giveDomain(1);
    Element *elem;

    int nnodes = d->giveNumberOfDofManagers();

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
            if ( dynamic_cast< Lattice3dBoundary * >( elem ) || dynamic_cast< LatticeLink3dBoundary * >( elem ) || dynamic_cast< Lattice2dBoundary * >( elem ) ) {
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
    } // end of default piece for simple geometry elements
}


void
VTKXMLLatticeExportModule::setupVTKPieceCross(ExportRegion &vtkPieceCross, TimeStep *tStep, Set& region)
{
    // Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.

    Domain *domain  = emodel->giveDomain(1);

    IntArray elements = region.giveElementList();
    int numberOfElements = elements.giveSize();

    //Loop over the elements and get crossSectionNodes
    int numberOfCrossSectionNodes = 0;

    IntArray crossSectionTable;
    crossSectionTable.resize(numberOfElements);
    int numberOfNodes = 0;
    for ( int ie = 1; ie <= numberOfElements; ie++ ) {
        if (  dynamic_cast< LatticeStructuralElement * >( domain->giveElement(elements.at(ie) ) ) ) {
            numberOfCrossSectionNodes = ( static_cast< LatticeStructuralElement * >( domain->giveElement(elements.at(ie) ) ) )->giveNumberOfCrossSectionNodes();
        } else if ( dynamic_cast< LatticeTransportElement * >( domain->giveElement(elements.at(ie) ) ) ) {
            numberOfCrossSectionNodes = ( static_cast< LatticeTransportElement * >( domain->giveElement(elements.at(ie) ) ) )->giveNumberOfCrossSectionNodes();
        }
        crossSectionTable.at(ie) = numberOfCrossSectionNodes;
        numberOfNodes += numberOfCrossSectionNodes;
    }

    FloatMatrix nodeTable;
    nodeTable.resize(numberOfNodes, 3);

    FloatArray crossSectionCoordinates;
    FloatArray coords(3);

    //Store node coordinates in table
    int nodeCounter = 0;
    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        if (  dynamic_cast< LatticeStructuralElement * >( domain->giveElement(elements.at(ie) ) ) ) {
            numberOfCrossSectionNodes =  ( static_cast< LatticeStructuralElement * >( domain->giveElement(elements.at(ie) ) ) )->giveNumberOfCrossSectionNodes();
            ( static_cast< LatticeStructuralElement * >( domain->giveElement(elements.at(ie) ) ) )->giveCrossSectionCoordinates(crossSectionCoordinates);
        } else if ( dynamic_cast< LatticeTransportElement * >( domain->giveElement(elements.at(ie) ) ) ) {
            numberOfCrossSectionNodes = ( static_cast< LatticeTransportElement * >( domain->giveElement(elements.at(ie) ) ) )->giveNumberOfCrossSectionNodes();
            ( static_cast< LatticeTransportElement * >( domain->giveElement(elements.at(ie) ) ) )->giveCrossSectionCoordinates(crossSectionCoordinates);
        }

        for ( int is = 0; is < numberOfCrossSectionNodes; is++ ) {
            nodeCounter++;
            nodeTable.at(nodeCounter, 1) = crossSectionCoordinates.at(3 * is + 1);
            nodeTable.at(nodeCounter, 2) = crossSectionCoordinates.at(3 * is + 2);
            nodeTable.at(nodeCounter, 3) = crossSectionCoordinates.at(3 * is + 3);
        }
    }

    if ( numberOfNodes > 0 && numberOfElements > 0 ) {
        // Export nodes as vtk vertices
        vtkPieceCross.setNumberOfNodes(numberOfNodes);
        for ( int inode = 1; inode <= numberOfNodes; inode++ ) {
            coords.at(1) = nodeTable.at(inode, 1);
            coords.at(2) = nodeTable.at(inode, 2);
            coords.at(3) = nodeTable.at(inode, 3);
            vtkPieceCross.setNodeCoords(inode, coords);
        }
    }


    //-------------------------------------------
    // Export all the cell data for the piece
    //-------------------------------------------
    IntArray connectivity;
    vtkPieceCross.setNumberOfCells(numberOfElements);
    int numElNodes;

    int offset = 0;
    for ( int ei = 1; ei <= numberOfElements; ei++ ) {
        numElNodes = crossSectionTable.at(ei);

        connectivity.resize(numElNodes);
        for ( int i = 1; i <= numElNodes; i++ ) {
            connectivity.at(i) = offset + i;
        }
        vtkPieceCross.setConnectivity(ei, connectivity);

        vtkPieceCross.setCellType(ei, 7);
        offset += numElNodes;
        vtkPieceCross.setOffset(ei, offset);
    }

    this->exportCellVars(vtkPieceCross, region, cellVarsToExport, tStep);
}


void
VTKXMLLatticeExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
    this->doOutputNormal(tStep, forcedOutput);
    if ( crossSectionExportFlag ) {
        this->doOutputCross(tStep, forcedOutput);
    }
    this->defaultVTKPiece.clear();

}


void
VTKXMLLatticeExportModule::doOutputCross(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    this->fileStreamCross = this->giveOutputStreamCross(tStep);
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);

    this->fileStreamCross << "<!-- TimeStep " << tStep->giveTargetTime() * timeScale << " Computed " << current->tm_year + 1900 << "-" << setw(2) << current->tm_mon + 1 << "-" << setw(2) << current->tm_mday << " at " << current->tm_hour << ":" << current->tm_min << ":" << setw(2) << current->tm_sec << " -->\n";
    this->fileStreamCross << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    this->fileStreamCross << "<UnstructuredGrid>\n";

    int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;

    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        Set* region = this->giveRegionSet(pieceNum);
        // Fills a data struct (VTKPiece) with all the necessary data.
        this->setupVTKPieceCross(this->defaultVTKPieceCross, tStep, *region);

        // Write the VTK piece to file.
        anyPieceNonEmpty += this->writeVTKPieceCross(this->defaultVTKPieceCross, tStep);
        this->defaultVTKPieceCross.clear();
    }


    if ( anyPieceNonEmpty == 0 ) {
        // write empty piece, Otherwise ParaView complains if the whole vtu file is without <Piece></Piece>
        this->fileStreamCross << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n";
        this->fileStreamCross << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> </DataArray>\n</Cells>\n";
        this->fileStreamCross << "</Piece>\n";
    }

    this->fileStreamCross << "</UnstructuredGrid>\n</VTKFile>";
    this->fileStreamCross.close();
}


void
VTKXMLLatticeExportModule::doOutputNormal(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    this->fileStream = this->giveOutputStream(tStep);
    struct tm *current;
    time_t now;
    time(& now);
    current = localtime(& now);

    this->fileStream << "<!-- TimeStep " << tStep->giveTargetTime() * timeScale << " Computed " << current->tm_year + 1900 << "-" << setw(2) << current->tm_mon + 1 << "-" << setw(2) << current->tm_mday << " at " << current->tm_hour << ":" << current->tm_min << ":" << setw(2) << current->tm_sec << " -->\n";

    this->fileStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    this->fileStream << "<UnstructuredGrid>\n";

    this->giveSmoother(); // make sure smoother is created, Necessary? If it doesn't exist it is created /JB

    int nPiecesToExport = this->giveNumberOfRegions();     //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;

    NodalRecoveryModel *smoother = giveSmoother();
    NodalRecoveryModel *primVarSmoother = givePrimVarSmoother();

    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        // Fills a data struct (VTKPiece) with all the necessary data.
        Set *region = this->giveRegionSet(pieceNum);
        this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
        this->writeVTKPieceProlog(this->defaultVTKPiece, tStep);  
        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(this->defaultVTKPiece, *region, primaryVarsToExport, *primVarSmoother, tStep);
        this->exportIntVars(this->defaultVTKPiece, *region, internalVarsToExport, *smoother, tStep);
        //const IntArray &elements = region.giveElementList();
        this->exportCellVars(this->defaultVTKPiece, *region, cellVarsToExport, tStep);
        // Write the VTK piece to file.
        anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPiece, tStep);
        this->writeVTKPieceEpilog(this->defaultVTKPiece, tStep);  
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

                //this->exportCompositeElement(this->defaultVTKPiece, el, tStep);
                this->exportCompositeElement(this->defaultVTKPieces, el, tStep);

                for ( int j = 0; j < ( int ) this->defaultVTKPieces.size(); j++ ) {
                    this->writeVTKPieceProlog(this->defaultVTKPieces[j], tStep);  
                    anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPieces [ j ], tStep);
                    this->writeVTKPieceEpilog(this->defaultVTKPieces[j], tStep);  
                }
            }
        }
    }     // end loop over composite elements

    if ( anyPieceNonEmpty == 0 ) {
        // write empty piece, Otherwise ParaView complains if the whole vtu file is without <Piece></Piece>
        this->fileStream << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\n";
        this->fileStream << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> </DataArray>\n</Cells>\n";
        this->fileStream << "</Piece>\n";
    }

    this->fileStream << "</UnstructuredGrid>\n</VTKFile>";
    this->fileStream.close();
}




int
VTKXMLLatticeExportModule::initRegionNodeNumbering(ExportRegion& vtkPiece, Domain *domain, TimeStep *tStep, Set& region)
{
    int nnodes = domain->giveNumberOfDofManagers();
    int elementNode, node;
    int currOffset = 1;
    Element *element;

    int  regionDofMans = 0;
    int regionSingleCells = 0;

    IntArray elements = region.giveElementList();

    int extraNodes = 0.;
    for ( int ie = 1; ie <= elements.giveSize(); ie++ ) {
        element = domain->giveElement(elements.at(ie) );

#ifdef __SM_MODULE
        if ( dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) {
            IntArray loc = element->giveLocation();
            for ( int ielnode = 1; ielnode <= 2; ielnode++ ) {
                if ( loc.at(ielnode) != 0 ) {
                    extraNodes++;
                }
            }
        }
#endif
    }
    IntArray& regionG2LNodalNumbers = vtkPiece.getMapG2L();
    IntArray& regionL2GNodalNumbers = vtkPiece.getMapL2G();

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
        if ( dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) {
            elemNodes = 2;
        }
#endif

        for ( elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
#ifdef __SM_MODULE
            if ( dynamic_cast< Lattice3dBoundary * >( element ) || dynamic_cast< LatticeLink3dBoundary * >( element ) || dynamic_cast< Lattice2dBoundary * >( element ) ) { //beam elements - only unique nodes
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
    vtkPiece.setNumberOfNodes(regionDofMans);   
    vtkPiece.setNumberOfCells(regionSingleCells);

    regionG2LNodalNumbers.resizeWithValues(regionDofMans);
    regionL2GNodalNumbers.resize(regionDofMans);

    for ( int i = 1; i <= regionDofMans; i++ ) {
        if ( regionG2LNodalNumbers.at(i) ) {
            regionG2LNodalNumbers.at(i) = currOffset++;
            regionL2GNodalNumbers.at(regionG2LNodalNumbers.at(i) ) = i;
        }
    }
    return 1;
}


void
VTKXMLLatticeExportModule::exportPrimaryVars(ExportRegion &vtkPiece, Set& region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
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

    FloatArray macroField(controlNode->giveNumberOfDofs() );
    for ( int j = 1; j <= controlNode->giveNumberOfDofs(); j++ ) {
        macroField.at(j) = controlNode->giveDofWithID(dofIdArray.at(j) )->giveUnknown(VM_Total, tStep);
    }

    //Get unit cell size
    const auto unitCellSize = controlNode->giveCoordinates();

    for ( int i = 1, n = primaryVarsToExport.giveSize(); i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);

        for ( int inode = 1; inode <= mapL2G.giveSize(); inode++ ) {
            if ( inode <= nnodes && mapL2G.at(inode) <= nnodes && mapL2G.at(inode) != 0 ) { //no special treatment for master nodes
                DofManager *dman = d->giveNode(mapL2G.at(inode) );

                this->getNodalVariableFromPrimaryField(valueArray, dman, tStep, type, region, smoother);
                vtkPiece.setPrimaryVarInNode(type, inode, std::move(valueArray) );
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
                        if ( dofIdArray.giveSize() == 3 && dofIdArray.at(1) == E_xx && dofIdArray.at(2) == E_yy && dofIdArray.at(3) == G_xy ) { //Macroscale: 2D solid using Voigt notation
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(2) * switches.at(2) * macroField.at(2) + unitCellSize.at(1) * switches.at(1) * macroField.at(3);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIdArray.giveSize() == 1 && dofIdArray.at(1) == E_xx ) { //Macroscale: 3D truss. 1 DOF: EXX
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1);
                            valueArray.at(2) = helpArray.at(2);
                            valueArray.at(3) = helpArray.at(3);
                        } else if ( dofIdArray.giveSize() == 6 && dofIdArray.at(1) == E_xx && dofIdArray.at(2) == E_yy && dofIdArray.at(3) == E_zz && dofIdArray.at(4) == G_yz && dofIdArray.at(5) == G_xz && dofIdArray.at(6) == G_xy ) { //Macroscale: 3D solid using Voigt notation
                            valueArray.resize(helpArray.giveSize() );
                            valueArray.at(1) = helpArray.at(1) + unitCellSize.at(1) * switches.at(1) * macroField.at(1) +
                                               unitCellSize.at(3) * switches.at(3) * macroField.at(5) + unitCellSize.at(2) * switches.at(2) * macroField.at(6);
                            valueArray.at(2) = helpArray.at(2) + unitCellSize.at(2) * switches.at(2) * macroField.at(2) +
                                               unitCellSize.at(3) * switches.at(3) * macroField.at(4);
                            valueArray.at(3) = helpArray.at(3) + unitCellSize.at(3) * switches.at(3) * macroField.at(3);
                        } else {
                            OOFEM_ERROR("Unknown periodic element type\n");
                        }
                    }
                } else {
                    valueArray.resize(3);
                }
                vtkPiece.setPrimaryVarInNode(type, inode, std::move(valueArray) );
            }
        }
    }
}


void
VTKXMLLatticeExportModule::exportIntVars(ExportRegion &vtkPiece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep)
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


bool
VTKXMLLatticeExportModule::writeVTKPieceCross(ExportRegion &vtkPieceCross, TimeStep *tStep)
{
    if ( !vtkPieceCross.giveNumberOfCells() ) {
        return false;
    }

    // Write output: node coords
    int numNodes = vtkPieceCross.giveNumberOfNodes();
    int numEl = vtkPieceCross.giveNumberOfCells();
    FloatArray coords;

    this->fileStreamCross << "<Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << numEl << "\">\n";
    this->fileStreamCross << "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ";

    for ( int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPieceCross.giveNodeCoords(inode);
        ///@todo move this below into setNodeCoords since it should alwas be 3 components anyway
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            this->fileStreamCross << scientific << coords.at(i) << " ";
        }

        for ( int i = coords.giveSize() + 1; i <= 3; i++ ) {
            this->fileStreamCross << scientific << 0.0 << " ";
        }
    }

    this->fileStreamCross << "</DataArray>\n</Points>\n";
    this->fileStreamCross << "<Cells>\n";
    this->fileStreamCross << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"> ";

    IntArray cellNodes;
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        cellNodes = vtkPieceCross.giveCellConnectivity(ielem);

        for ( int i = 1; i <= cellNodes.giveSize(); i++ ) {
            this->fileStreamCross << cellNodes.at(i) - 1 << " ";
        }
        this->fileStreamCross << " ";
    }

    this->fileStreamCross << "</DataArray>\n";

    // output the offsets (index of individual element data in connectivity array)
    this->fileStreamCross << " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"> ";

    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        this->fileStreamCross << vtkPieceCross.giveCellOffset(ielem) << " ";
    }

    this->fileStreamCross << "</DataArray>\n";


    // output cell (element) types
    this->fileStreamCross << " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"> ";
    for ( int ielem = 1; ielem <= numEl; ielem++ ) {
        this->fileStreamCross << vtkPieceCross.giveCellType(ielem) << " ";
    }

    this->fileStreamCross << "</DataArray>\n";
    this->fileStreamCross << "</Cells>\n";


    ///@todo giveDataHeaders is currently not updated wrt the new structure -> no file names in headers /JB
    std::string pointHeader, cellHeader;
    this->giveDataHeaders(pointHeader, cellHeader);

    this->fileStreamCross << pointHeader.c_str();
    this->fileStreamCross << "</PointData>\n";
    this->fileStreamCross << cellHeader.c_str();

    this->writeCellVarsCross(vtkPieceCross);

    this->fileStreamCross << "</CellData>\n";
    this->fileStreamCross << "</Piece>\n";

    vtkPieceCross.clear();
    return true;
}

void
VTKXMLLatticeExportModule::writeCellVarsCross(ExportRegion &vtkPiece)
{
    FloatArray valueArray;
    int numCells = vtkPiece.giveNumberOfCells();
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        const char *name = __InternalStateTypeToString(type);
        ( void ) name; //silence the warning

        this->fileStreamCross << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        valueArray.resize(ncomponents);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(type, ielem);
            for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
                this->fileStreamCross << valueArray.at(i) << " ";
            }
        }
        this->fileStreamCross << "</DataArray>\n";

#ifdef _PYBIND_BINDINGS
#if 0
        if ( pythonExport ) {
            py::list vals;
            for ( int ielem = 1; ielem <= numCells; ielem++ ) {
                valueArray = vtkPiece.giveCellVar(i, ielem);
                vals.append(valueArray);
            }
            this->Py_CellVars [ name ] = vals;
        }
#endif
#endif
    }//end of for
}
} // end namespace oofem
