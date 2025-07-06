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
 *               Copyright (C) 1993 - 20   Borek Patzak
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

#include "vtkhdf5reader.h"
#include "unstructuredgridfield.h"
#include "elementgeometrytype.h"
#include "floatarray.h"
#include "element.h"
#include "timestep.h"

#ifdef __HDF_MODULE
#include <H5Cpp.h>
#endif

namespace oofem {

#ifndef __HDF_MODULE
VTKHDF5Reader::VTKHDF5Reader()
{
    OOFEM_ERROR("VTKHDF5Reader: HDF5 module not available, cannot read VTKHDF5 files");
}


VTKHDF5Reader::~VTKHDF5Reader() {}
void VTKHDF5Reader::initialize(std::string &fileName) {}
void VTKHDF5Reader::finalize() {}       
void VTKHDF5Reader::readMesh(UnstructuredGridField&, TimeStep* tStep) {}
void VTKHDF5Reader::readField(UnstructuredGridField&, TimeStep* tStep, const std::string &field_name) {}

#else


VTKHDF5Reader::VTKHDF5Reader(): stepValues()
{
    this->fileName = "";
    this->file = NULL;
}

VTKHDF5Reader::~VTKHDF5Reader() {
   this->finalize();
}

void
VTKHDF5Reader::initialize (std::string& filename)
{
    this->fileName = filename;
    // open hdf5 file with fileName
    this->file = new H5::H5File(fileName, H5F_ACC_RDONLY);
    // open VTHHDF group
    // test if VTKHDF group exists
    if ( !file->exists("/VTKHDF") ) {
        OOFEM_ERROR("VTKHDF group not found in file");
    }
    this->topGroup = new H5::Group(file->openGroup("/VTKHDF"));
    // read version attribute in topGroup
    H5::Attribute versionAttr = topGroup->openAttribute("Version");
    int version[2];
    versionAttr.read(versionAttr.getDataType(), &version);

    // Read type attribute
    H5::Attribute typeAttr = topGroup->openAttribute("Type");
    std::string type;
    typeAttr.read(typeAttr.getDataType(), type);

    // check version and type
    if (version[0] != 2 || version[1] != 2 || type != "UnstructuredGrid") {
        OOFEM_ERROR("Unsupported version or type of VTKHDF5 file, version 2.2 and UnstructuredGrid type is required.");
    }

    // read NSteps attribute from Steps group
    this->stepsGroup = new H5::Group(topGroup->openGroup("Steps"));
    H5::Attribute nstepsattr = stepsGroup->openAttribute("NSteps");
    nstepsattr.read(nstepsattr.getDataType(), &numSteps);
    // read step values
    H5::DataSet values = stepsGroup->openDataSet("Values");
    hsize_t dim[1];
    values.getSpace().getSimpleExtentDims(dim);
    this->stepValues.resize(dim[0]);
    values.read(this->stepValues.givePointer(), values.getDataType());
}   

void VTKHDF5Reader::finalize()
{
    if ( this->file ) {
        delete this->topGroup;
        delete this->stepsGroup;
        delete this->file;
        this->file=NULL;
    }
}


void VTKHDF5Reader::readDataSet (H5::DataSet& dset, int rank, hsize_t* dim, hsize_t* offset, H5::DataType type, void* data)
{
    H5::DataSpace dspace = dset.getSpace();
    /* create hyperslab to update point data */
    dspace.selectHyperslab( H5S_SELECT_SET, dim, offset );
    H5::DataSpace mem_space(rank, dim, nullptr);
    dset.read(data, type, mem_space, dspace );
}

void VTKHDF5Reader::getTimeStepOffsets(TimeStep* tStep, int& nParts, int& pOffset, int& pointOffset, int& cellOffset, int& connIdOffset, int& offsetOfOffset, int& nPoints, int &nCells, int& nconectivities) {
    H5::IntType itype(H5::PredType::NATIVE_UINT);
    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);

    // find if tStep time present in array of stepValues, assuming stepValues is sorted
    bool found = false;
    double tt = tStep->giveTargetTime();
    unsigned int indx = 0;
    // simple serach for matching time; could be optimized to take advantage of sorted array
    for (indx = 0; indx < stepValues.giveSize(); indx++) {
        if (fabs(stepValues[indx] - tt) < 1.e-6) {
            found = true;
            break;
        }
    }
    if (found) {
        // tStep time is present in stepValues array
        // read PointOffsets, CellOffsets, and ConnectivityIdOffsets arrays from Steps group
        H5::DataSet npDSet = stepsGroup->openDataSet("NumberOfParts");
        H5::DataSet poDSet = stepsGroup->openDataSet("PartOffsets");
        H5::DataSet pointOffsetsDSet = stepsGroup->openDataSet("PointOffsets");
        H5::DataSet cellOffsetsDSet = stepsGroup->openDataSet("CellOffsets");
        H5::DataSet connIdOffsetsDSet = stepsGroup->openDataSet("ConnectivityIdOffsets");
        H5::DataSet offsetsDSet = topGroup->openDataSet( "Offsets" );
        // read the data with offset corresponding to tStep number
        hsize_t offset[] = {indx};
        hsize_t dims[1] = {1};
        this->readDataSet(npDSet, 1, dims, offset, itype, &nParts);
        this->readDataSet(poDSet, 1, dims, offset, itype, &pOffset);

        this->readDataSet(pointOffsetsDSet, 1, dims, offset, itype, &pointOffset);
        hsize_t     dims11[] = {1,1};
        hsize_t offset11[] = {indx,0};
        this->readDataSet(cellOffsetsDSet, 2, dims11, offset11, itype, &cellOffset);
        this->readDataSet(connIdOffsetsDSet, 2, dims11, offset11, itype, &connIdOffset);
    

        // read NumberOfPoints and NumberOfCells and NumberOfConnectivityIds
        H5::DataSet numberOfPointsDSet = topGroup->openDataSet("NumberOfPoints");
        H5::DataSet numberOfCellsDSet = topGroup->openDataSet("NumberOfCells");
        H5::DataSet numberOfConnectivityIds = topGroup->openDataSet("NumberOfConnectivityIds");

        // read numberOfCellsDSet data
        int stepCells[numSteps];
        numberOfCellsDSet.read(stepCells, numberOfCellsDSet.getDataType());
        // evaluate offset array offset
        offsetOfOffset = 0;
        for (int i = 0; i < indx; i++) {
            offsetOfOffset += stepCells[i]+1;
        }
        nCells = stepCells[indx];

        this->readDataSet(numberOfPointsDSet, 1, dims, offset, itype, &nPoints);
        //this->readDataSet(numberOfCellsDSet, 1, dims, offset, itype, &nCells);
        this->readDataSet(numberOfConnectivityIds, 1, dims, offset, itype, &nconectivities);

    } else {
        // tStep time is not present in stepValues array
        OOFEM_ERROR("Matching time not found, time=%lf, step number %d", tt, tStep->giveNumber());
    }
}

Element_Geometry_Type 
VTKHDF5Reader::giveElementGeometryType(int vtkCellType)
{
    Element_Geometry_Type elemGT = EGT_unknown;


    // implement reverse mapping from vtkCellType to Element_Geometry_Type
    if ( vtkCellType == 1 ) {
        elemGT = EGT_point;
    } else if ( vtkCellType == 3 ) {
        elemGT = EGT_line_1;
    } else if ( vtkCellType == 21 ) {
        elemGT = EGT_line_2;
    } else if ( vtkCellType == 5 ) {
        elemGT = EGT_triangle_1;
    } else if ( vtkCellType == 22 ) {
        elemGT = EGT_triangle_2;
    } else if ( vtkCellType == 10 ) {
        elemGT = EGT_tetra_1;
    } else if ( vtkCellType == 24 ) {
        elemGT = EGT_tetra_2;
    } else if ( vtkCellType == 9 ) {
        elemGT = EGT_quad_1;
    } else if ( vtkCellType == 30 ) {
        elemGT = EGT_quad_21_interface;
    } else if ( vtkCellType == 23 ) {
        elemGT = EGT_quad_2;
    } else if ( vtkCellType == 12 ) {
        elemGT = EGT_hexa_1;
    } else if ( vtkCellType == 25 ) {
        elemGT = EGT_hexa_2;
    } else if ( vtkCellType == 29 ) {
        elemGT = EGT_hexa_27;
    } else if ( vtkCellType == 13 ) {
        elemGT = EGT_wedge_1;
    } else if ( vtkCellType == 26 ) {
        elemGT = EGT_wedge_2;
    } else {
        OOFEM_ERROR("unsupported vtk cell type %d", vtkCellType);
    }
    return elemGT;
}

void VTKHDF5Reader::readMesh(UnstructuredGridField& f, TimeStep* tStep) 
{
    int nParts, pOffset, pointOffset, cellOffset, connIdOffset, offsetOfOffset, nPoints, nCells, nconectivities;
    this->getTimeStepOffsets (tStep, nParts, pOffset, pointOffset, cellOffset, connIdOffset, offsetOfOffset, nPoints, nCells, nconectivities);
    H5::IntType itype(H5::PredType::NATIVE_UINT);
    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
    H5::DataSet pointsDSet = topGroup->openDataSet( "Points" );
    H5::DataSet typesDSet = topGroup->openDataSet( "Types" );
    H5::DataSet connectivityDSet = topGroup->openDataSet( "Connectivity" );
    H5::DataSet offsetsDSet = topGroup->openDataSet( "Offsets" );

    // read points
    hsize_t dim2[]={nPoints,3};
    hsize_t offset2[] = {pointOffset, 0};
    double *points = new double[dim2[0]*3];
    this->readDataSet(pointsDSet, 2, dim2, offset2, dtype, points);
    
    // read cell types
    hsize_t dim1[1]={nCells};
    hsize_t offset1[] = {cellOffset};
    unsigned int *types = new unsigned int[dim1[0]];
    this->readDataSet(typesDSet, 1, dim1, offset1, itype, types);
    
    // read cell connectivity
    dim1[0]=nconectivities;
    offset1[0] = connIdOffset;
    unsigned int *connectivity = new unsigned int[dim1[0]];
    this->readDataSet(connectivityDSet, 1, dim1, offset1, itype, connectivity);
    
    // read cell offsets
    dim1[0]=nCells+1;
    offset1[0] = offsetOfOffset;
    unsigned int *offsets = new unsigned int[dim1[0]];
    this->readDataSet(offsetsDSet, 1, dim1, offset1, itype, offsets);
    
    // set up the field mesh
    f.initialize(nPoints, nCells);
    // define mesh points
    for (int i = 0; i < nPoints; i++) {
        FloatArray coords(3);
        coords.at(1) = points[i*3];
        coords.at(2) = points[i*3+1];
        coords.at(3) = points[i*3+2];
        f.addVertex(i+1, coords); //1-based
    }

    // define mesh elements
    for (int i = 0; i < nCells; i++) {
        int type = types[i];
        int nNodes = offsets[i+1] - offsets[i];
        IntArray nodes(nNodes);
        for (int j = 0; j < nNodes; j++) {
            nodes[j] = connectivity[offsets[i]+j]+1; // 1-based
        }
        f.addCell(i+1, this->giveElementGeometryType(type), nodes) ; //1-based
    }

    delete points;
    delete types;
    delete connectivity;
    delete offsets;

    return;
}

void VTKHDF5Reader::readField(UnstructuredGridField& field, TimeStep* tStep, const std::string &field_name) 
{
    int nParts, pOffset, pointOffset, cellOffset, connIdOffset, offsetOfOffset, nPoints, nCells, nconectivities;
    this->getTimeStepOffsets (tStep, nParts, pOffset, pointOffset, cellOffset, connIdOffset, offsetOfOffset, nPoints, nCells, nconectivities);
    // open PointData group in topGroup
    H5::Group* pointDataGroup = new H5::Group(topGroup->openGroup("PointData"));
    // open field_name dataset in PointData group
    if (pointDataGroup->nameExists(field_name.c_str())) {
        H5::DataSet dSet = pointDataGroup->openDataSet(field_name.c_str());
        int nPoints = field.giveNumberOfVertices();
        // read field data
        H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
        hsize_t dim[2];
        dSet.getSpace().getSimpleExtentDims(dim);
        dim[0] = nPoints;
        double *fieldData = new double[nPoints*dim[1]];
        hsize_t offset2[] = {pointOffset, 0};
        this->readDataSet(dSet, 2, dim, offset2, dtype, fieldData);

        FloatArray valueArray((int)dim[1]); 
        for ( unsigned int inode = 0; inode < nPoints; inode++ ) {
            for (unsigned int i=0; i<dim[1]; i++) {
                valueArray(i) = fieldData[(inode)*dim[1]+i];
            }
            field.setVertexValue(inode+1, valueArray);
        }
        delete fieldData;

    } else {
        OOFEM_ERROR("Field %s not found in VTKHDF5 file", field_name.c_str());
    }
    delete pointDataGroup;
    return;
}

#endif // __HDF_MODULE

} // end namespace oofem
