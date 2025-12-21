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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef exportregion_h
#define exportregion_h

#include "exportmodule.h"
#include "intarray.h"
#include "internalstatevaluetype.h"
#include "internalstatetype.h"
#include "unknowntype.h"
#include "element.h"



#include <string>
#include <list>

using namespace std;
namespace oofem {
    ///@todo Rename this to something like "ExportPiece" and move it to a separate file (it doesn't actually contain anything VTK-specific).
/// Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.
class OOFEM_EXPORT ExportRegion
{
public:
    typedef  uint8_t setMembershipGroupType;

    ExportRegion()
    {
        numCells = 0;
        numNodes = 0;
        numSetGroups = 0;
    }

    void clear();

    void setNumberOfNodes(int numNodes);
    int giveNumberOfNodes() { return this->numNodes; }

    void setNumberOfCells(int numCells);
    int giveNumberOfCells() { return this->numCells; }

    int giveNumberOfSetGroups() { return this->numSetGroups; }

    void setConnectivity(int cellNum, IntArray &nodes);
    IntArray &giveCellConnectivity(int cellNum) { return this->connectivity [ cellNum - 1 ]; }

    void setCellType(int cellNum, int type) { this->elCellTypes.at(cellNum) = type; }
    int giveCellType(int cellNum) { return this->elCellTypes.at(cellNum); }

    void setOffset(int cellNum, int offset) { this->elOffsets.at(cellNum) = offset; }
    int giveCellOffset(int cellNum) { return this->elOffsets.at(cellNum); }

    void setNodeCoords(int nodeNum, const FloatArray &coords);
    FloatArray &giveNodeCoords(int nodeNum) { return this->nodeCoords [ nodeNum - 1 ]; }

    void setNumberOfPrimaryVarsToExport(const IntArray& primVars, int numNodes);
    void setNumberOfLoadsToExport(int numVars, int numNodes);
    void setNumberOfInternalVarsToExport(const IntArray& ists, int numNodes);
    void setNumberOfInternalXFEMVarsToExport(int numVars, int numEnrichmentItems, int numNodes);
    void setNumberOfCellVarsToExport(const IntArray& cellVars, int numCells);
    void setNumberOfSetMembershipsToExport(int numSets, int numNodes, int numCells) {
        this->numSetGroups = 1+numSets/8;
        this->vertexSetMembership.resize(numNodes);
        for ( int i = 0; i < numNodes; i++ ) {
            this->vertexSetMembership [ i ].resize(numSetGroups, 0);
        }
        this->cellSetMembership.resize(numCells);
        for ( int i = 0; i < numCells; i++ ) {
            this->cellSetMembership [ i ].resize(numSetGroups, 0);
        }
    }

    void setPrimaryVarInNode(UnknownType  type, int nodeNum, FloatArray valueArray);
    FloatArray &givePrimaryVarInNode(UnknownType type, int nodeNum) { return this->nodeVars [ type ] [ nodeNum - 1 ]; }

    void setLoadInNode(int varNum, int nodeNum, FloatArray valueArray);
    FloatArray &giveLoadInNode(int varNum, int nodeNum) { return this->nodeLoads [ varNum - 1 ] [ nodeNum - 1 ]; }

    void setInternalVarInNode(InternalStateType type, int nodeNum, FloatArray valueArray);
    FloatArray &giveInternalVarInNode (InternalStateType type, int nodeNum) { return this->nodeVarsFromIS [ type ] [ nodeNum - 1 ]; }

    void setInternalXFEMVarInNode(int varNum, int eiNum, int nodeNum, FloatArray valueArray);
    FloatArray &giveInternalXFEMVarInNode(int varNum, int eiNum, int nodeNum) { return this->nodeVarsFromXFEMIS [ varNum - 1 ] [ eiNum - 1 ] [ nodeNum - 1 ]; }

    void setCellVar(InternalStateType type, int cellNum, FloatArray valueArray);
    FloatArray &giveCellVar(InternalStateType type, int cellNum) { return this->cellVars [ type ] [ cellNum - 1 ]; }


    const std::vector<setMembershipGroupType>& getVertexSetMembershipGroup(int inode) { return this->vertexSetMembership [ inode-1 ]; }
    const std::vector<setMembershipGroupType>& getCellSetMembershipGroup(int icell) { return this->cellSetMembership [ icell-1 ]; }
    void setVertexSetMembership(int set, int nodeNum) {
         this->vertexSetMembership [ nodeNum - 1 ][ (set-1) / 8 ]  |= 1 << ((set-1) % 8);   
    }
    void setCellSetMembership(int set, int cellNum) {
         this->cellSetMembership [ cellNum - 1 ][ (set-1) / 8 ]  |= 1 << ((set-1) % 8);   
    }
    IntArray& getMapG2L () {return this->mapG2L;}
    IntArray& getMapL2G () {return this->mapL2G;}
    //void setRegionCells(IntArray& cells) {this->regionElInd = cells;}
    IntArray& getRegionCells () {return this->regionElInd;}

    FloatMatrix getVertices () ;
    IntArray getCellConnectivity ();
    IntArray getCellTypes ();
    FloatMatrix getPrimaryVertexValues (UnknownType u);
    FloatMatrix getInternalVertexValues(InternalStateType u);
    FloatMatrix getCellValues(InternalStateType u);


private:
    int numCells;
    int numNodes;
    int numSetGroups;
    IntArray elCellTypes;
    IntArray elOffsets;
    // dofman local->global and global->local region map 
    IntArray mapG2L, mapL2G;
    // region elements 
    IntArray regionElInd;


    std::vector< FloatArray >nodeCoords;   // all the nodes in the piece [node][coords]
    std::vector< IntArray >connectivity;   // cell connectivity [cell][nodes]
    std::map< UnknownType, std::vector< FloatArray > >nodeVars;     // [field][node][valArray]
    std::vector< std::vector< FloatArray > >nodeLoads;     // [field][node][valArray]
    std::map< InternalStateType, std::vector< FloatArray > >nodeVarsFromIS;     // [field][node][valArray]
    std::vector< std::vector< std::vector< FloatArray > > >nodeVarsFromXFEMIS;       // [field][ei][node][valArray]
    std::map< InternalStateType, std::vector< FloatArray > >cellVars;     // [el][field][valArray]
    // set membership (set membership for each vertex encoded into 8bit uint -? setGroups)
    std::vector< std::vector<setMembershipGroupType> > vertexSetMembership; // [setGroup][vertex][membership]
    std::vector< std::vector<setMembershipGroupType> > cellSetMembership; // [setGroup][cell][membership]

};

} // end namespace oofem
#endif // exportregion_h
