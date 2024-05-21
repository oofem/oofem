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
 *               Copyright (C) 1993 - 2024   Borek Patzak
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

#ifndef vtkhdf5reader_h
#define vtkhdf5reader_h

#include "unstructuredgridfield.h"
#include "elementgeometrytype.h"

#include "timestep.h"
#include "intarray.h"
#include "internalstatevaluetype.h"


#include <string>
#include <list>

#ifdef __HDF_MODULE
#include "H5Cpp.h"
#endif


using namespace std;
namespace oofem {

/**
 * Represents VTK (Visualization Toolkit) hdf5 reader. It can read VTK hdf5 file format, Unstructured grid dataset.
 * 
 */
class OOFEM_EXPORT VTKHDF5Reader
{
protected:
    /// File name
    std::string fileName;
#ifdef __HDF_MODULE
    /// File
    H5::H5File *file;
    /// Top group
    H5::Group *topGroup;
     /// Steps group
    H5::Group *stepsGroup;
    /// Number of steps
    unsigned int numSteps;
    /// step values (times)
    FloatArray stepValues;
    /// Number of point data
    unsigned int numPointData;
    /// Number of cell data
    unsigned int numCellData;
#endif
public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKHDF5Reader();
    /// Destructor
    virtual ~VTKHDF5Reader();

    void initialize(std::string &fileName) ;
    void finalize() ;
    void readMesh(UnstructuredGridField&, TimeStep* tStep);
    void readField(UnstructuredGridField&, TimeStep* tStep, const std::string &field_name);
protected:
#ifdef __HDF_MODULE
    void readDataSet (H5::DataSet& dset, int rank, hsize_t* dim, hsize_t* offset, H5::DataType type, void* data);
    void getTimeStepOffsets(TimeStep* tStep, int& nParts, int& pOffset, int& pointOffset, int& cellOffset, int& connIdOffset, int &offsetOfOffset, int& nPoints, int &nCells, int& nconectivities);
    Element_Geometry_Type giveElementGeometryType(int vtkCellType);
#endif

};

} // end namespace oofem
#endif // vtkhdf5exportmodule_h
