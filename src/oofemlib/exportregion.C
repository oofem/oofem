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

#include "exportregion.h"
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

void
ExportRegion::setNumberOfNodes(int numNodes)
{
    this->numNodes = numNodes;
    this->nodeCoords.resize(numNodes);
}

void
ExportRegion::setNumberOfCells(int numCells)
{
    this->numCells = numCells;
    this->connectivity.resize(numCells);
    this->elCellTypes.resize(numCells);
    this->elOffsets.resize(numCells);
}

void
ExportRegion::setConnectivity(int cellNum, IntArray &nodes)
{
    this->connectivity [ cellNum - 1 ] = nodes;
}

void
ExportRegion::setNodeCoords(int nodeNum, const FloatArray &coords)
{
    this->nodeCoords [ nodeNum - 1 ] = coords;
}

void
ExportRegion::setNumberOfPrimaryVarsToExport(const IntArray& primVars, int numNodes)
{
  for (int it=1; it<= primVars.giveSize(); it++) {
    UnknownType type = (UnknownType) primVars.at(it);
    this->nodeVars [type].resize(numNodes);
  }
}

void
ExportRegion::setNumberOfLoadsToExport(int numVars, int numNodes)
{
    this->nodeLoads.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeLoads [ i - 1 ].resize(numNodes);
    }
}

void
ExportRegion::setNumberOfInternalVarsToExport(const IntArray& ists, int numNodes)
{
    for ( int i = 1; i <= ists.giveSize(); i++ ) {
      InternalStateType itype = (InternalStateType) ists.at(i);
        this->nodeVarsFromIS [ itype ].resize(numNodes);
    }
}

void
ExportRegion::setNumberOfInternalXFEMVarsToExport(int numVars, int numEnrichmentItems, int numNodes)
{
    this->nodeVarsFromXFEMIS.resize(numVars);
    for ( int i = 1; i <= numVars; i++ ) {
        this->nodeVarsFromXFEMIS [ i - 1 ].resize(numEnrichmentItems);
        for ( int j = 1; j <= numEnrichmentItems; j++ ) {
            this->nodeVarsFromXFEMIS [ i - 1 ] [ j - 1 ].resize(numNodes);
        }
    }
}

void
ExportRegion::setNumberOfCellVarsToExport(const IntArray& cellVars, int numCells)
{
     for ( int i = 1; i <= cellVars.giveSize(); i++ ) {
       InternalStateType type = (InternalStateType) cellVars.at(i);
        this->cellVars [ type ].resize(numCells);
    }
}

void
ExportRegion::setPrimaryVarInNode(UnknownType type, int nodeNum, FloatArray valueArray)
{
  this->nodeVars [ type ] [ nodeNum - 1 ] = std::move(valueArray);
}

void
ExportRegion::setLoadInNode(int varNum, int nodeNum, FloatArray valueArray)
{
    this->nodeLoads [ varNum - 1 ] [ nodeNum - 1 ] = std::move(valueArray);
}

void
ExportRegion::setInternalVarInNode(InternalStateType type, int nodeNum, FloatArray valueArray)
{
    this->nodeVarsFromIS [ type ] [ nodeNum - 1 ] = std::move(valueArray);
}

void
ExportRegion::setInternalXFEMVarInNode(int varNum, int eiNum, int nodeNum, FloatArray valueArray)
{
    this->nodeVarsFromXFEMIS [ varNum - 1 ] [ eiNum - 1 ] [ nodeNum - 1 ] = std::move(valueArray);
}


void
ExportRegion::setCellVar(InternalStateType type, int cellNum, FloatArray valueArray)
{
    this->cellVars [ type ] [ cellNum - 1 ] = std::move(valueArray);
}

void
ExportRegion::clear()
{
    ///@todo Will this give a memory leak? / JB
    numCells = 0;
    numNodes = 0;
    this->connectivity.clear();
    this->elCellTypes.clear();
    this->elOffsets.clear();
    this->cellVars.clear();
    this->nodeCoords.clear();
    this->nodeVars.clear();
    this->nodeVarsFromIS.clear();
    this->nodeVarsFromXFEMIS.clear();
}

#ifdef _PYBIND_BINDINGS

py::array_t<double>
ExportRegion::getVertices () {
  double* result = new double[this->numNodes*3];
  for (int i=0; i<this->numNodes; i++) {
    FloatArray &c = this->giveNodeCoords(i+1);
    for (int j=0; j<c.giveSize(); j++) {
      result[i*3+j]=c[j];
    }
    for (int j=c.giveSize(); j<3; j++) {
      result[i*3+j] = 0.0;
    }
  }
  return py::array_t<double>(std::vector<ptrdiff_t>{this->numNodes, 3}, &result[0]);
}

py::array_t<int>
ExportRegion::getCellConnectivity () {
  int reqSize = 0;
  for (int i=0; i<this->numCells; i++) {
    reqSize++;
    reqSize+=this->giveCellConnectivity(i+1).giveSize();
  }
  int* result = new int[reqSize];
  int pos=0;
  for (int i=0; i<this->numCells; i++) {
    IntArray &c = this->giveCellConnectivity(i+1);
    int csize = c.giveSize();
    result[pos++]=csize;
    for (int i = 0; i<csize; i++) {
      result[pos++]=c[i]-1;
    }
  }
  return py::array_t<int>(reqSize, result);
}

py::array_t<int>
ExportRegion::getCellTypes () {
  int* result = new int[this->numCells];
  for (int i=0; i<this->numCells; i++) {
    result[i]=this->giveCellType(i+1); // m.giveCellType(this->regionElInd.at(i+1));
  }
  return py::array_t<int>(this->numCells, result);     
}

py::array_t<double>
ExportRegion::getPrimaryVertexValues (UnknownType u) {

  if (this->nodeVars.find(u) != this->nodeVars.end()) {
    // key exists
    std::vector<FloatArray>& nodalVars = this->nodeVars[u];
    // get size of nodal record
    int recSize = nodalVars[0].giveSize();    
    double* result = new double[this->numNodes*recSize];
    int counter = 0;
    for (int i=0;i<this->numNodes; i++) {
      FloatArray& v = nodalVars[i];
      for (int j=0; j<recSize; j++) {
        result[counter++]=v[j];
      }
    }
    return py::array_t<double>(std::vector<ptrdiff_t>{this->numNodes, recSize}, result);
    //return py::array_t<double>(this->numNodes*recSize, result);     
  } else {
    return py::array_t<double>(0);
  }
}

py::array_t<double>
ExportRegion::getInternalVertexValues(InternalStateType u) {
  if (this->nodeVarsFromIS.find(u)!=this->nodeVarsFromIS.end()) {
    // key exists
    std::vector<FloatArray>& nodalVars = this->nodeVarsFromIS[u];
    // get size of nodal record
    int recSize = nodalVars[0].giveSize();
    double* result = new double[this->numNodes*recSize];
    int counter = 0;
    for (int i=0;i<this->numNodes; i++) {
      FloatArray& v = nodalVars[i];
      for (int j=0; j<recSize; j++) {
        result[counter++]=v[j];
      }
    }
    return py::array_t<double>(std::vector<ptrdiff_t>{this->numNodes, recSize}, result);
  } else {
    return py::array_t<double>(0);   
  }
}
  
py::array_t<double>
ExportRegion::getCellValues(InternalStateType u) {
  if (this->cellVars.find(u)!=this->cellVars.end()) {
    // key exists
    std::vector<FloatArray>& cv = this->cellVars[u];
    // get size of nodal record
    int recSize = cv[0].giveSize();
    double* result = new double[this->numCells*recSize];
    int counter = 0;
    for (int i=0;i<this->numCells; i++) {
      FloatArray& v = cv[i];
      for (int j=0; j<recSize; j++) {
        result[counter++]=v[j];
      }
    }
    return py::array_t<double>(std::vector<ptrdiff_t>{this->numCells, recSize}, result);
  } else {
    return py::array_t<double>(0);   
  }
}

#endif


} // end namespace oofem