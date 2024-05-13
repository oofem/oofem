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

#include "vtkhdf5exportmodule.h"
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

#include "xfem/xfemmanager.h"
#include "xfem/enrichmentitem.h"

#include "nodalaveragingrecoverymodel.h"
#include "zznodalrecoverymodel.h"

#ifdef __PFEM_MODULE
 #include "pfem/pfemparticle.h"
#endif

#include <string>
#include <sstream>
#include <ctime>

#ifdef __VTK_MODULE
 #include <vtkPoints.h>
 #include <vtkPointData.h>
 #include <vtkDoubleArray.h>
 #include <vtkCellArray.h>
 #include <vtkCellData.h>
 #include <vtkHDF5UnstructuredGridWriter.h> // @check!
 #include <vtkHDF5PUnstructuredGridWriter.h>
 #include <vtkUnstructuredGrid.h>
 #include <vtkSmartPointer.h>
#endif

#ifdef __HDF_MODULE
#include <H5Cpp.h>
#endif

namespace oofem {
REGISTER_ExportModule(VTKHDF5ExportModule)


VTKHDF5ExportModule::VTKHDF5ExportModule(int n, EngngModel *e) : VTKBaseExportModule(n, e), internalVarsToExport(), primaryVarsToExport()
{}


VTKHDF5ExportModule::~VTKHDF5ExportModule() {}


void
VTKHDF5ExportModule::initializeFrom(InputRecord &ir)
{
    ExportModule::initializeFrom(ir);

    int val;

    IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, _IFT_VTKHDF5ExportModule_cellvars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_VTKHDF5ExportModule_vars); // Macro - see internalstatetype.h
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_VTKHDF5ExportModule_primvars); // Macro - see unknowntype.h
    IR_GIVE_OPTIONAL_FIELD(ir, externalForcesToExport, _IFT_VTKHDF5ExportModule_externalForces); // Macro - see unknowntype.h
    IR_GIVE_OPTIONAL_FIELD(ir, ipInternalVarsToExport, _IFT_VTKHDF5ExportModule_ipvars); // Macro - see internalstatetype.h

    val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_VTKHDF5ExportModule_stype); // Macro
    stype = ( NodalRecoveryModel::NodalRecoveryModelType ) val;
}


void
VTKHDF5ExportModule::initialize()
{
    this->smoother = nullptr;
    this->primVarSmoother = nullptr;
    VTKBaseExportModule::initialize();

#ifdef __HDF_MODULE
    char fext[100];
    sprintf( fext, ".m%d.hdf", this->number);
    this->fileName = this->emodel->giveOutputBaseFileName() + fext;
    /* 
    *  Create the named file, truncating the existing file one if any,
    *  using default create and access property lists.
    */
   try {
    H5::IntType itype(H5::PredType::NATIVE_INT);
    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);

    this->file = new H5::H5File (this->fileName, H5F_ACC_TRUNC);
    /*
     * Create a VTHHDF group 
    */
    this->topGroup = new H5::Group (file->createGroup("/VTKHDF"));
    this->pointDataGroup = new H5::Group(topGroup->createGroup("PointData"));
    this->cellDataGroup = new H5::Group(topGroup->createGroup("CellData"));
    /* Create version attribute*/
    int version[]={2,2};
    hsize_t dims[1]={2};
    H5::DataSpace attrSpace(1, dims);
    H5::Attribute va = topGroup->createAttribute("Version", itype, attrSpace);
    va.write(itype, version);
    
    // Create type attribute
    H5::StrType str_type(H5::PredType::C_S1, 16); /* required by vtkhdf reader Type string with ASCII encoding and fixed length*/
    str_type.setCset(H5T_CSET_ASCII);
    H5::DataSpace attrSpace1(H5S_SCALAR);

    H5::Attribute att = topGroup->createAttribute( "Type", str_type, attrSpace1 );
    att.write( str_type, std::string("UnstructuredGrid"));

    /* Create Unstructured grid top level datasets*/
    hsize_t dim1[]={1};
    hsize_t maxdim1[]={H5S_UNLIMITED};
    H5::DataSpace ds1(1, dim1, maxdim1);
    H5::DSetCreatPropList plist1;
    plist1.setLayout(H5D_CHUNKED);
    hsize_t chunk_dims1[] = {1000};
    plist1.setChunk(1, chunk_dims1);
    H5::DataSet NumberOfConnectivityIds = topGroup->createDataSet( "NumberOfConnectivityIds", itype, ds1, plist1 );
    H5::DataSet numberOfPointsDSet = topGroup->createDataSet( "NumberOfPoints", itype, ds1, plist1 );
    H5::DataSet numberOfCellsDSet = topGroup->createDataSet( "NumberOfCells", itype, ds1, plist1 );
    
    hsize_t dim4[]={1,3};
    hsize_t maxdim4[]={H5S_UNLIMITED,3};
    H5::DataSpace pointsDSpace(2, dim4, maxdim4);
    // Create dataset creation property list to enable chunking
    H5::DSetCreatPropList plist;
    plist.setLayout(H5D_CHUNKED);
    hsize_t chunk_dims2[2] = {1000, 3};
    plist.setChunk(2, chunk_dims2);
    H5::DataSet pointsDSet = topGroup->createDataSet( "Points", dtype, pointsDSpace, plist );
    
    H5::DataType ttype(H5::PredType::STD_U8LE);
    H5::DataSet typesDSet = topGroup->createDataSet( "Types", ttype, ds1, plist1 );
    H5::DataSet connectivityDSet = topGroup->createDataSet( "Connectivity", itype, ds1, plist1 );
    H5::DataSet offsetsDSet = topGroup->createDataSet( "Offsets", itype, ds1, plist1 );



    // Create Steps group (transient data support)
    this->stepsGroup = new H5::Group(topGroup->createGroup("Steps"));
    H5::DataSpace ads(H5S_SCALAR);
    H5::Attribute nsa = this->stepsGroup->createAttribute("NSteps", itype, ads);
    int nsteps = 0;
    nsa.write(itype, &nsteps);
    // create Steps data sets
    hsize_t dim[]={1};
    hsize_t maxdim[]={H5S_UNLIMITED};
    H5::DataSpace stepDSpace(1, dim, maxdim);
    // Create dataset creation property list to enable chunking
    plist1.setLayout(H5D_CHUNKED);
    hsize_t chunk_dims3[1] = {10};
    plist1.setChunk(1, chunk_dims3);
    // create Steps/Values (each entry indicates the time value for the associated time step)
    H5::DataSet valuesDSet = this->stepsGroup->createDataSet( "Values", dtype, stepDSpace, plist1 );
    // create PartOffsets [dims = (NSteps)]: each entry indicates at which part offset to start reading the associated time step
    H5::DataSet poDSet = this->stepsGroup->createDataSet( "PartOffsets", itype, stepDSpace, plist1 );
    H5::DataSet npDSet = this->stepsGroup->createDataSet( "NumberOfParts", itype, stepDSpace, plist1 );
    // PointOffsets [dims = (NSteps)]: each entry indicates where in the VTKHDF/Points data set to start reading point coordinates for the associated time step
    H5::DataSet pointOffsetsDSet = this->stepsGroup->createDataSet( "PointOffsets", itype, stepDSpace, plist1 );
    // CellOffsets [dims = (NSteps, NTopologies)]: each entry indicates by how many cells to offset reading into the connectivity offset structures for the associated time step 
    hsize_t dim2[]={1,1};
    hsize_t maxdim2[]={H5S_UNLIMITED,1};
    H5::DataSpace ds2(2, dim2, maxdim2);
    H5::DSetCreatPropList plist2;
    plist2.setLayout(H5D_CHUNKED);
    hsize_t chunk_dims4[2] = {10,1};
    plist2.setChunk(2, chunk_dims4);
    H5::DataSet cellOffsetsDSet = this->stepsGroup->createDataSet( "CellOffsets", itype, ds2, plist2);
    // ConnectivityIdOffsets [dims = (NSteps, NTopologies)]: each entry indicates by how many values to offset reading into the connectivity indexing structures for the associated time step 
    H5::DataSet connIdOffsetsDSet = this->stepsGroup->createDataSet( "ConnectivityIdOffsets", itype, ds2, plist2);
        

    this->pointCounter = 0;
    this->cellCounter=0;
    this->connCounter=0;
    this->offsetCounter=0;

} catch ( H5::Exception& error) {
    error.getDetailMsg(); // error.printErrorStack();
    return;
}


#endif

}


void
VTKHDF5ExportModule::terminate()
{
#ifdef __HDF_MODULE
    delete topGroup;
    delete pointDataGroup;
    delete cellDataGroup;
    delete stepsGroup;
    delete file;
#endif
}

void
VTKHDF5ExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }
  
#ifdef __VTK_MODULE
    //this->fileStream = vtkSmartPointer< vtkUnstructuredGrid >::New();
    this->nodes = vtkSmartPointer< vtkPoints >::New();
    this->elemNodeArray = vtkSmartPointer< vtkIdList >::New();

#endif
#ifdef __HDF_MODULE
    //int rank = 0; // at present no support fro parallel partitioned output
    H5::IntType itype(H5::PredType::NATIVE_UINT);
    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);

    // update Steps/offset informations
    // increment number of steps
    H5::Attribute nstepsattr = stepsGroup->openAttribute("NSteps");
    unsigned int steps;
    nstepsattr.read(itype, &steps);
    steps ++;
    nstepsattr.write(itype, &steps);
    // resize and update steps datasets
    // update Steps/Values Attribute (step time values)
    hsize_t ssize[] = {steps};
    H5::DataSet values = stepsGroup->openDataSet("Values");
    values.extend( ssize );
    /* create hyperslab to update point data */
    hsize_t     offset[] = {steps-1};
    hsize_t     dims[1] = {1};
    double t = tStep->giveTargetTime();
    this->updateDataSet(values, 1, dims, offset, dtype, &t);

    // Update Steps/PointOffsets
    H5::DataSet pods = stepsGroup->openDataSet("PointOffsets");
    pods.extend( ssize );
    this->updateDataSet(pods, 1, dims, offset, itype, &this->pointCounter);

    // Update Steps/CellOffsets
    hsize_t ssize11[] = {steps,1};
    H5::DataSet codset = stepsGroup->openDataSet("CellOffsets");// CellOffsets
    codset.extend(ssize11);
    /* create hyperslab to update point data */
    hsize_t     offset11[] = {steps-1,0};
    hsize_t     dims11[] = {1,1};
    int tempdata11[1][1]={(int)this->cellCounter};
    this->updateDataSet(codset, 2, dims11, offset11, itype, tempdata11);

    // update Steps/ConnectivityIdOffsets
    H5::DataSet cidset = stepsGroup->openDataSet("ConnectivityIdOffsets");
    cidset.extend(ssize11);
    this->updateDataSet(cidset, 2, dims11, offset11, itype, &this->connCounter);

    // update Steps/NumberOfParts and Steps/PartOffsets
    H5::DataSet npset = stepsGroup->openDataSet("NumberOfParts");
    npset.extend(ssize);
    int tval = 1;
    this->updateDataSet(npset, 1, dims, offset, itype, &tval);
    H5::DataSet poset = stepsGroup->openDataSet("PartOffsets");
    poset.extend(ssize);
    tval =steps-1;
    this->updateDataSet(poset, 1, dims, offset, itype, &tval);



    H5::DataSet pointsDSet = topGroup->openDataSet( "Points" );
    H5::DataSet typesDSet = topGroup->openDataSet( "Types" );
    H5::DataSet connectivityDSet = topGroup->openDataSet( "Connectivity" );
    H5::DataSet offsetsDSet = topGroup->openDataSet( "Offsets" );


    int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
    int anyPieceNonEmpty = 0;
    NodalRecoveryModel *smoother = giveSmoother();
    NodalRecoveryModel *primVarSmoother = givePrimVarSmoother();

    unsigned int stepPointCounter = 0, stepCellCounter = 0, stepConnCounter = 0;
    
    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        // Fills a data struct (VTKPiece) with all the necessary data.
        Set* region = this->giveRegionSet(pieceNum);
        this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
        this->writeVTKPieceProlog(this->defaultVTKPiece, tStep, pointCounter, cellCounter, connCounter, offsetCounter, 
                                    stepPointCounter, stepCellCounter, stepConnCounter, 
                                    pointsDSet, connectivityDSet, offsetsDSet, typesDSet); 

        this->exportPrimaryVars(this->defaultVTKPiece, *region, primaryVarsToExport, *primVarSmoother, tStep);
        this->exportIntVars(this->defaultVTKPiece, *region, internalVarsToExport, *smoother, tStep);
        anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPiece, tStep);

        this->defaultVTKPiece.clear();

    }

    this->pointCounter += stepPointCounter;
    this->cellCounter += stepCellCounter;
    this->connCounter += stepConnCounter;
    this->offsetCounter += (stepCellCounter+1);

    /* Update basic datasets*/
    H5::DataSet numberOfPointsDSet = topGroup->openDataSet("NumberOfPoints");
    H5::DataSet numberOfCellsDSet = topGroup->openDataSet("NumberOfCells");
    H5::DataSet NumberOfConnectivityIds = topGroup->openDataSet("NumberOfConnectivityIds");

    numberOfPointsDSet.extend(ssize);
    numberOfCellsDSet.extend(ssize);
    NumberOfConnectivityIds.extend(ssize);

    this->updateDataSet(numberOfPointsDSet, 1, dims, offset, itype, &stepPointCounter);
    this->updateDataSet(numberOfCellsDSet, 1, dims, offset, itype, &stepCellCounter);
    this->updateDataSet(NumberOfConnectivityIds, 1, dims, offset, itype, &stepConnCounter);

    return;

#endif

#if 0
    
    for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
        // Fills a data struct (VTKPiece) with all the necessary data.
        Set* region = this->giveRegionSet(pieceNum);
        this->setupVTKPiece(this->defaultVTKPiece, tStep, *region);
        /*
        this->writeVTKPieceProlog(this->defaultVTKPiece, tStep); 
        // Export primary, internal and XFEM variables as nodal quantities
        this->exportPrimaryVars(this->defaultVTKPiece, *region, primaryVarsToExport, *primVarSmoother, tStep);
        this->exportIntVars(this->defaultVTKPiece, *region, internalVarsToExport, *smoother, tStep);
        this->exportExternalForces(this->defaultVTKPiece, *region, externalForcesToExport, tStep);
        this->exportCellVars(this->defaultVTKPiece, *region, cellVarsToExport, tStep);

        // Write the VTK piece to file.
        anyPieceNonEmpty += this->writeVTKPieceVariables(this->defaultVTKPiece, tStep);
        this->writeVTKPieceEpilog(this->defaultVTKPiece, tStep);   
        */
        this->defaultVTKPiece.clear();
        }
        
    return;
#endif
}

#ifdef __HDF_MODULE
bool
VTKHDF5ExportModule::writeVTKPieceProlog(ExportRegion &vtkPiece, TimeStep *tStep, unsigned int &pointCounter, unsigned int &cellCounter, unsigned int &connCounter, unsigned int& offsetCounter,
                                        unsigned int& stepPointCounter, unsigned int& stepCellCounter, unsigned int& stepConnCounter,
                                        H5::DataSet &pointsDSet, H5::DataSet &connectivityDSet, H5::DataSet &offsetDSet, H5::DataSet &typeDSet)
{
   // Writes a VTK piece header + geometry to file.
   // This could be the whole domain (most common case) or it can be a
   // (so-called) composite element consisting of several VTK cells (layered structures, XFEM, etc.).

    // Write output: node coords
    unsigned int numNodes = vtkPiece.giveNumberOfNodes();
    unsigned int numEl = vtkPiece.giveNumberOfCells();
    FloatArray coords;

    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
    H5::IntType itype(H5::PredType::NATIVE_INT);
    H5::DataType ttype(H5::PredType::NATIVE_UINT); /* required by vtkhdf reader datset Types of type unsigned int */

    if ( !vtkPiece.giveNumberOfCells() ) {
      return false;
    }

    double *points = new double[numNodes*3];

    for ( unsigned int inode = 1; inode <= numNodes; inode++ ) {
        coords = vtkPiece.giveNodeCoords(inode);
        ///@todo move this below into setNodeCoords since it should alwas be 3 components anyway
        for ( int i = 0; i < coords.giveSize(); i++ ) {
            points[(inode-1)*3+i] = coords[i];
        }
        for ( int i = coords.giveSize() ; i < 3; i++ ) {
            points[(inode-1)*3+i] = 0.0;
        }
    }

    /*
     * Extend the dataset. This call assures that dataset is at least 3 x 3.
     */
      hsize_t psize[] = {pointCounter+numNodes, 3};
      pointsDSet.extend( psize );
    /* create hyperslab to update point data */
      hsize_t     offset[] = {pointCounter, 0};
      hsize_t      dims1[] = { numNodes, 3};            /* data1 dimensions */
      this->updateDataSet(pointsDSet, 2, dims1, offset, dtype, points);

    delete(points);
    

    // Write output: connectivity, offsets, cell types
    // output the connectivity data

    IntArray cellNodes;
    unsigned int connectivitySize=0;
    for ( unsigned int ielem = 1; ielem <= numEl; ielem++ ) {
        connectivitySize+=vtkPiece.giveCellConnectivity(ielem).giveSize();
    }
    IntArray conn(connectivitySize);
    IntArray offsets(numEl+1);
    unsigned int * types = new unsigned int [numEl];
    int cp = 0;

    offsets.at(1)=0;
    for ( unsigned int ielem = 1; ielem <= numEl; ielem++ ) {
        offsets.at(ielem+1) = vtkPiece.giveCellOffset(ielem);
        types[ielem-1] = vtkPiece.giveCellType(ielem);
        cellNodes = vtkPiece.giveCellConnectivity(ielem);
        for ( int i = 0; i < cellNodes.giveSize(); i++ ) {
            conn[cp++]=cellNodes[i]-1;
        }
    }

    H5::DataSpace cspace = connectivityDSet.getSpace();
    H5::DataSpace ospace = offsetDSet.getSpace();
    H5::DataSpace tspace =     typeDSet.getSpace();

    hsize_t csize[] = {connectivitySize+connCounter};
    connectivityDSet.extend( csize );
    hsize_t cellsize[] = {numEl+cellCounter};
    hsize_t cellsize1[] = {numEl+offsetCounter+1};
    offsetDSet.extend( cellsize1 );
    typeDSet.extend(cellsize);

    /* create hyperslab to update point data */
      hsize_t     coffset[] = {connCounter};
      hsize_t     cdims[] = {connectivitySize };            /* data1 dimensions */
      this->updateDataSet(connectivityDSet, 1, cdims, coffset, itype, conn.givePointer());
      


    /* create hyperslab to update cell data */
      hsize_t     celloffset[] = {cellCounter};
      hsize_t     offsetOffset[] = {offsetCounter};
      hsize_t     celldims[] = {numEl };            /* data1 dimensions */
      hsize_t     celldims1[] = {numEl+1};            /* data1 dimensions */
    this->updateDataSet(offsetDSet, 1, celldims1, offsetOffset, itype, offsets.givePointer());
    this->updateDataSet(typeDSet, 1, celldims, celloffset, ttype, types);

    delete(types);
    
    stepPointCounter+=numNodes;
    stepCellCounter = numEl;
    stepConnCounter = connectivitySize;
    return true;


}

bool
VTKHDF5ExportModule::writeVTKPieceEpilog(ExportRegion &vtkPiece, TimeStep *tStep)
{
    if ( !vtkPiece.giveNumberOfCells() ) {
        return false;
    }

#ifndef __VTK_MODULE
    //this->fileStream << "</Piece>\n";
#endif
    return true;
}



bool
VTKHDF5ExportModule::writeVTKPieceVariables(ExportRegion &vtkPiece, TimeStep *tStep)
{
    // Write a VTK piece variables to file.
    // This could be the whole domain (most common case) or it can be a
    // (so-called) composite element consisting of several VTK cells (layered structures, XFEM, etc.).

    if ( !vtkPiece.giveNumberOfCells() ) {
        return false;
    }


#ifndef __VTK_MODULE 
    ///@todo giveDataHeaders is currently not updated wrt the new structure -> no file names in headers /JB
    //std::string pointHeader, cellHeader;
    //this->giveDataHeaders(pointHeader, cellHeader);

    //this->fileStream << pointHeader.c_str();
#endif

    this->writePrimaryVars(vtkPiece);       // Primary field
    this->writeIntVars(vtkPiece);           // Internal State Type variables smoothed to the nodes
    this->writeExternalForces(vtkPiece);           // External forces

    //if ( emodel->giveDomain(1)->hasXfemManager() ) {
    //    this->writeXFEMVars(vtkPiece);      // XFEM State Type variables associated with XFEM structure
    //}

#ifndef __VTK_MODULE
    //this->fileStream << "</PointData>\n";
    //this->fileStream << cellHeader.c_str();
#endif
    this->writeCellVars(vtkPiece);          // Single cell variables ( if given in the integration points then an average will be exported)

#ifndef __VTK_MODULE
    //this->fileStream << "</CellData>\n";
#endif
    return true;
}

#ifdef __VTK_MODULE
void
VTKHDF5ExportModule::giveDataHeaders(std::string &pointHeader, std::string &cellHeader)
{
    std::string scalars, vectors, tensors;

    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        if ( type == DisplacementVector || type == EigenVector || type == VelocityVector || type == DirectorField || type == MacroSlipVector || type == ResidualForce ) {
            vectors += __UnknownTypeToString(type);
            vectors.append(" ");
        } else if ( type == FluxVector || type == PressureVector || type == Temperature || type == Humidity || type == DeplanationFunction ) {
            scalars += __UnknownTypeToString(type);
            scalars.append(" ");
        } else {
            OOFEM_ERROR("unsupported UnknownType %s", __UnknownTypeToString(type) );
        }
    }

    for ( int i = 1; i <= internalVarsToExport.giveSize(); i++ ) {
        InternalStateType isttype = ( InternalStateType ) internalVarsToExport.at(i);
        InternalStateValueType vtype = giveInternalStateValueType(isttype);

        if ( vtype == ISVT_SCALAR ) {
            scalars += __InternalStateTypeToString(isttype);
            scalars.append(" ");
        } else if ( vtype == ISVT_VECTOR ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
            tensors += __InternalStateTypeToString(isttype);
            tensors.append(" ");
        } else {
            OOFEM_ERROR("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
        }
    }

    for ( int i = 1; i <= externalForcesToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) externalForcesToExport.at(i);
        if ( type == DisplacementVector || type == VelocityVector || type == DirectorField ) {
            vectors += std::string("Load") + __UnknownTypeToString(type);
            vectors.append(" ");
        } else if ( type == FluxVector || type == PressureVector || type == Temperature || type == Humidity ) {
            scalars += std::string("Load") + __UnknownTypeToString(type);
            scalars.append(" ");
        } else {
            OOFEM_ERROR("unsupported UnknownType %s", __UnknownTypeToString(type) );
        }
    }

    // print header 
    pointHeader = "<PointData Scalars=\"" + scalars + "\" "
                  +  "Vectors=\"" + vectors + "\" "
                  +  "Tensors=\"" + tensors + "\" >\n";


    scalars.clear();
    vectors.clear();
    tensors.clear();
    // prepare header
    for ( int i = 1; i <= this->cellVarsToExport.giveSize(); i++ ) {
        InternalStateType isttype = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType vtype = giveInternalStateValueType(isttype);

        if ( vtype == ISVT_SCALAR ) {
            scalars += __InternalStateTypeToString(isttype);
            scalars.append(" ");
        } else if ( vtype == ISVT_VECTOR ) {
            vectors += __InternalStateTypeToString(isttype);
            vectors.append(" ");
        } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
            tensors += __InternalStateTypeToString(isttype);
            tensors.append(" ");
        } else {
            OOFEM_WARNING("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
        }
    }

    // print header
    cellHeader = "<CellData Scalars=\"" + scalars + "\" "
                 +  "Vectors=\"" + vectors + "\" "
                 +  "Tensors=\"" + tensors + "\" >\n";
}
#endif


void
VTKHDF5ExportModule::writeIntVars(ExportRegion &vtkPiece)
{

    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
    H5::DataSet dset;
    unsigned int varIndx;
    
    int n = internalVarsToExport.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(i);
        unsigned int ncomponents;

        const char *name = __InternalStateTypeToString(type);
        ( void ) name;//silence warning
        unsigned int numNodes = vtkPiece.giveNumberOfNodes();
        FloatArray valueArray;
        valueArray = vtkPiece.giveInternalVarInNode(type, 1);
        ncomponents = valueArray.giveSize();
        ( void ) ncomponents;//silence warning
     
#ifdef __HDF_MODULE
        if (this->pointDataGroup->nameExists(name)) {
            dset = this->pointDataGroup->openDataSet( name );
            hsize_t dim[2];
            dset.getSpace().getSimpleExtentDims(dim); // get existing dimensions to append the data
            varIndx = dim[0]; 
        } else {
            hsize_t dim4[]={1,ncomponents};
            hsize_t maxdim4[]={H5S_UNLIMITED,ncomponents};
            H5::DataSpace pointsDSpace(2, dim4, maxdim4);
            // Create dataset creation property list to enable chunking
            H5::DSetCreatPropList plist;
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunk_dims[2] = {1000, ncomponents};
            plist.setChunk(2, chunk_dims);
            dset = this->pointDataGroup->createDataSet( name, dtype, pointsDSpace, plist );
            varIndx = 0;    

        }   
        // extend the dataset
        hsize_t psize[] = {varIndx+numNodes, ncomponents};
        dset.extend( psize );

        hsize_t     offset[2];
        hsize_t      dims1[2] = { numNodes, ncomponents};
        H5::DataSpace dspace = dset.getSpace();
        H5::DataSpace mem_space(2, dims1, nullptr);
        double *pdata = new double[numNodes*ncomponents];
        for ( unsigned int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.giveInternalVarInNode(type, inode);
            for (unsigned int i=0; i<ncomponents; i++) {
                pdata[(inode-1)*ncomponents+i]=valueArray(i);
            }
        }
        /* create hyperslab to update point data */
        offset[0] = varIndx;
        offset[1] = 0;
        //dspace.selectNone();
        dspace.selectHyperslab( H5S_SELECT_SET, dims1, offset );
        /*
            * Write the data to the hyperslab.
        */
        //dspace = dset->getSpace();
        dset.write( pdata, dtype, mem_space, dspace );
        delete pdata;
        //this->fileStream << "</DataArray>\n";
#endif
    }
}



//----------------------------------------------------
// Misc. functions
//----------------------------------------------------
#ifdef __VTK_MODULE
void
VTKHDF5ExportModule::writeVTKPointData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray)
{
    // Write the data to file
    int ncomponents = varArray->GetNumberOfComponents();
    switch ( ncomponents ) {
    case 1:
        this->fileStream->GetPointData()->SetActiveScalars(name);
        this->fileStream->GetPointData()->SetScalars(varArray);
        break;
    case 3:
        this->fileStream->GetPointData()->SetActiveVectors(name);
        this->fileStream->GetPointData()->SetVectors(varArray);
        break;
    case 9:
        this->fileStream->GetPointData()->SetActiveTensors(name);
        this->fileStream->GetPointData()->SetTensors(varArray);
        break;
    }
}
#else
void
VTKHDF5ExportModule::writeVTKPointData(FloatArray &valueArray)
{
    // Write the data to file
    for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
        //this->fileStream << scientific << valueArray.at(i) << " ";
    }
}
#endif


#ifdef __VTK_MODULE
void
VTKHDF5ExportModule::writeVTKCellData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray)
{
    // Write the data to file
    int ncomponents = varArray->GetNumberOfComponents();
    switch ( ncomponents ) {
    case 1:
        this->fileStream->GetCellData()->SetActiveScalars(name);
        this->fileStream->GetCellData()->SetScalars(varArray);
        break;
    case 3:
        this->fileStream->GetCellData()->SetActiveVectors(name);
        this->fileStream->GetCellData()->SetVectors(varArray);
        break;
    case 9:
        this->fileStream->GetCellData()->SetActiveTensors(name);
        this->fileStream->GetCellData()->SetTensors(varArray);
        break;
    }
}

#else

void
VTKHDF5ExportModule::writeVTKCellData(FloatArray &valueArray)
{
    // Write the data to file ///@todo exact copy of writeVTKPointData so remove
    for ( int i = 1; i <= valueArray.giveSize(); i++ ) {
        //this->fileStream << valueArray.at(i) << " ";
    }
}
#endif


void
VTKHDF5ExportModule::writePrimaryVars(ExportRegion &vtkPiece)
{
    H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
    H5::DataSet dset;
    int primvarIndx;
    

    for ( int i = 1; i <= primaryVarsToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        unsigned int ncomponents = giveInternalStateTypeSize(valType);
        ( void ) ncomponents; //silence the warning
        unsigned int numNodes = vtkPiece.giveNumberOfNodes();
        const char *name = __UnknownTypeToString(type);
        ( void ) name; //silence the warning
#ifdef __HDF_MODULE
        if (this->pointDataGroup->nameExists(name)) {
            dset = this->pointDataGroup->openDataSet( name );
            hsize_t dim[2];
            dset.getSpace().getSimpleExtentDims(dim); // get existing dimensions to append the data
            primvarIndx = dim[0]; 
        } else {
            hsize_t dim4[]={1,ncomponents};
            hsize_t maxdim4[]={H5S_UNLIMITED,ncomponents};
            H5::DataSpace pointsDSpace(2, dim4, maxdim4);
            // Create dataset creation property list to enable chunking
            H5::DSetCreatPropList plist;
            plist.setLayout(H5D_CHUNKED);
            hsize_t chunk_dims[2] = {1000, ncomponents};
            plist.setChunk(2, chunk_dims);
            dset = this->pointDataGroup->createDataSet( name, dtype, pointsDSpace, plist );
            primvarIndx = 0;    

        }   
        // extend the dataset
        hsize_t psize[] = {primvarIndx+numNodes, ncomponents};
        dset.extend( psize );
#endif
        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
        varArray->SetName(name);
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( unsigned int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(type, inode);
            for ( int j = 1; j <= ncomponents; ++j ) {
                varArray->SetComponent(inode - 1, j - 1, valueArray.at(j) );
            }
        }

        this->writeVTKPointData(name, varArray);

#else
        hsize_t     offset[2];
        hsize_t      dims1[2] = { numNodes, ncomponents};
        H5::DataSpace dspace = dset.getSpace();
        H5::DataSpace mem_space(2, dims1, nullptr);
        double *pdata = new double[numNodes*ncomponents];
        for ( unsigned int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.givePrimaryVarInNode(type, inode);
            for (unsigned int i=0; i<ncomponents; i++) {
                pdata[(inode-1)*ncomponents+i]=valueArray(i);
            }
        }
        /* create hyperslab to update point data */
        offset[0] = primvarIndx;
        offset[1] = 0;
        //dspace.selectNone();
        dspace.selectHyperslab( H5S_SELECT_SET, dims1, offset );
        /*
            * Write the data to the hyperslab.
        */
        //dspace = dset->getSpace();
        dset.write( pdata, dtype, mem_space, dspace );
        delete pdata;
        //this->fileStream << "</DataArray>\n";
#endif
    }
}


void
VTKHDF5ExportModule::writeExternalForces(ExportRegion &vtkPiece)
{
    for ( int i = 1; i <= externalForcesToExport.giveSize(); i++ ) {
        UnknownType type = ( UnknownType ) externalForcesToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        ( void ) ncomponents; //silence the warning
        int numNodes = vtkPiece.giveNumberOfNodes();
        std::string name = std::string("Load") + __UnknownTypeToString(type);

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >varArray = vtkSmartPointer< vtkDoubleArray >::New();
        varArray->SetName(name.c_str() );
        varArray->SetNumberOfComponents(ncomponents);
        varArray->SetNumberOfTuples(numNodes);

        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.giveLoadInNode(i, inode);
            for ( int j = 1; j <= ncomponents; ++j ) {
                varArray->SetComponent(inode - 1, j - 1, valueArray.at(j) );
            }
        }

        this->writeVTKPointData(name.c_str(), varArray);

#else
        //this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        for ( int inode = 1; inode <= numNodes; inode++ ) {
            FloatArray &valueArray = vtkPiece.giveLoadInNode(i, inode);
            this->writeVTKPointData(valueArray);
        }
        //this->fileStream << "</DataArray>\n";
#endif
    }
}

void
VTKHDF5ExportModule::writeCellVars(ExportRegion &vtkPiece)
{
    FloatArray valueArray;
    int numCells = vtkPiece.giveNumberOfCells();
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
        InternalStateValueType valType = giveInternalStateValueType(type);
        int ncomponents = giveInternalStateTypeSize(valType);
        const char *name = __InternalStateTypeToString(type);
        ( void ) name; //silence the warning

        // Header
#ifdef __VTK_MODULE
        vtkSmartPointer< vtkDoubleArray >cellVarsArray = vtkSmartPointer< vtkDoubleArray >::New();
        cellVarsArray->SetName(name);
        cellVarsArray->SetNumberOfComponents(ncomponents);
        cellVarsArray->SetNumberOfTuples(numCells);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(i, ielem);
            for ( int i = 1; i <= ncomponents; ++i ) {
                cellVarsArray->SetComponent(ielem - 1, i - 1, valueArray.at(i) );
            }
        }

        this->writeVTKCellData(name, cellVarsArray);

#else
        //this->fileStream << " <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomponents << "\" format=\"ascii\"> ";
        valueArray.resize(ncomponents);
        for ( int ielem = 1; ielem <= numCells; ielem++ ) {
            valueArray = vtkPiece.giveCellVar(type, ielem);
            this->writeVTKCellData(valueArray);
        }
        //this->fileStream << "</DataArray>\n";
#endif
    
    }//end of for
}

#endif


NodalRecoveryModel *
VTKHDF5ExportModule::giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( !this->smoother ) {
        this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
    }

    return this->smoother.get();
}


NodalRecoveryModel *
VTKHDF5ExportModule::givePrimVarSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( !this->primVarSmoother ) {
        this->primVarSmoother = classFactory.createNodalRecoveryModel(NodalRecoveryModel::NRM_NodalAveraging, d);
    }

    return this->primVarSmoother.get();
}

#ifdef __HDF_MODULE
void
VTKHDF5ExportModule::exportIntVarsInGpAs(IntArray valIDs, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    int nc = 0;
    ( void ) nc; //silence the warning
    FloatArray gc, value;
    std::ofstream stream;
    InternalStateType isttype;
    InternalStateValueType vtype;
    std::string scalars, vectors, tensors;

    // output nodes Region By Region
    int nregions = this->giveNumberOfRegions(); // aka sets
    // open output stream
    std::string outputFileName = this->giveOutputBaseFileName(tStep) + ".gp.vtu";
    std::ofstream streamG;
    if ( pythonExport ) {
        streamG = std::ofstream(NULL_DEVICE);
    } else {
        streamG = std::ofstream(outputFileName);
    }

    if ( !streamG.good() ) {
        OOFEM_ERROR("failed to open file %s", outputFileName.c_str() );
    }

    streamG << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    streamG << "<UnstructuredGrid>\n";

    /* loop over regions */
    for ( int ireg = 1; ireg <= nregions; ireg++ ) {
        const IntArray &elements = this->giveRegionSet(ireg)->giveElementList();
        int nip = 0;
        for ( int i = 1; i <= elements.giveSize(); i++ ) {
            nip += d->giveElement(elements.at(i) )->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
        }

        //Create one cell per each GP
        streamG << "<Piece NumberOfPoints=\"" << nip << "\" NumberOfCells=\"" << nip << "\">\n";
        streamG << "<Points>\n <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> ";
        for ( int i = 1; i <= elements.giveSize(); i++ ) {
            int ielem = elements.at(i);

            for ( GaussPoint *gp : * d->giveElement(ielem)->giveDefaultIntegrationRulePtr() ) {
                d->giveElement(ielem)->computeGlobalCoordinates(gc, gp->giveNaturalCoordinates() );
                for ( double c : gc ) {
                    ( void ) c; //silence the warning
                    streamG << scientific << c << " ";
                }

                for ( int k = gc.giveSize() + 1; k <= 3; k++ ) {
                    streamG << scientific << 0.0 << " ";
                }
            }
        }

        streamG << " </DataArray>\n";
        streamG << "</Points>\n";
        streamG << "<Cells>\n";
        streamG << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
        for ( int j = 0; j < nip; j++ ) {
            streamG << j << " ";
        }

        streamG << " </DataArray>\n";
        streamG << " <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
        for ( int j = 1; j <= nip; j++ ) {
            streamG << j << " ";
        }

        streamG << " </DataArray>\n";
        streamG << " <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
        for ( int j = 1; j <= nip; j++ ) {
            streamG << "1 ";
        }

        streamG << " </DataArray>\n";
        streamG << "</Cells>\n";
        // prepare the data header
        for ( int vi = 1; vi <= valIDs.giveSize(); vi++ ) {
            isttype = ( InternalStateType ) valIDs.at(vi);
            vtype = giveInternalStateValueType(isttype);

            if ( vtype == ISVT_SCALAR ) {
                scalars += __InternalStateTypeToString(isttype);
                scalars.append(" ");
            } else if ( vtype == ISVT_VECTOR ) {
                vectors += __InternalStateTypeToString(isttype);
                vectors.append(" ");
            } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
                tensors += __InternalStateTypeToString(isttype);
                tensors.append(" ");
            } else {
                OOFEM_WARNING("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }
        }

        // print collected data summary in header
        streamG << "<PointData Scalars=\"" << scalars.c_str() << "\" Vectors=\"" << vectors.c_str() << "\" Tensors=\"" << tensors.c_str() << "\" >\n";
        scalars.clear();
        vectors.clear();
        tensors.clear();

        // export actual data, loop over individual IDs to export
        for ( int vi = 1; vi <= valIDs.giveSize(); vi++ ) {
            isttype = ( InternalStateType ) valIDs.at(vi);
            vtype = giveInternalStateValueType(isttype);
            if ( vtype == ISVT_SCALAR ) {
                nc = 1;
            } else if ( vtype == ISVT_VECTOR ) {
                nc = 3;
                if ( isttype == IST_BeamForceMomentTensor ) { //AS: to make the hack work
                    nc = 6;
                }
            } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
                nc = 9;
            } else {
                OOFEM_WARNING("unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
            }

            streamG << "  <DataArray type=\"Float64\" Name=\"" << __InternalStateTypeToString(isttype) << "\" NumberOfComponents=\"" << nc << "\" format=\"ascii\">";
            for ( int i = 1; i <= elements.giveSize(); i++ ) {
                int ielem = elements.at(i);

                // loop over default IRule gps
                for ( GaussPoint *gp : * d->giveElement(ielem)->giveDefaultIntegrationRulePtr() ) {
                    d->giveElement(ielem)->giveIPValue(value, gp, isttype, tStep);

                    if ( vtype == ISVT_VECTOR ) {
                        // bp: hack for BeamForceMomentTensor, which should be splitted into force and momentum vectors
                        if ( isttype == IST_BeamForceMomentTensor ) {
                            value.resizeWithValues(6);
                        } else {
                            value.resizeWithValues(3);
                        }
                    } else if ( vtype == ISVT_TENSOR_S3 || vtype == ISVT_TENSOR_S3E || vtype == ISVT_TENSOR_G ) {
                        FloatArray help = value;
                        this->makeFullTensorForm(value, help, vtype);
                    }

                    for ( double v : value ) {
                        ( void ) v; //silence the warning
                        streamG << scientific << v << " ";
                    }
                } // end loop over IPs
            } // end loop over elements

            streamG << "  </DataArray>\n";
        } // end loop over values to be exported
        streamG << "</PointData>\n</Piece>\n";
    } // end loop over regions

    streamG << "</UnstructuredGrid>\n";
    streamG << "</VTKFile>\n";
    if(streamG){
        streamG.close();
    }
}
#endif

#ifdef __HDF_MODULE
void VTKHDF5ExportModule::updateDataSet (H5::DataSet& dset, int rank, hsize_t* dim, hsize_t* offset, H5::DataType type, const void* data)
{
    H5::DataSpace dspace = dset.getSpace();
    /* create hyperslab to update point data */
    dspace.selectHyperslab( H5S_SELECT_SET, dim, offset );
    H5::DataSpace mem_space(rank, dim, nullptr);
    dset.write(data, type, mem_space, dspace );
}

#endif

} // end namespace oofem
