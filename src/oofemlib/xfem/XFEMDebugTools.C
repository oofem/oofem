/*
 * XFEMDebugTools.C
 *
 *  Created on: Jun 5, 2013
 *      Author: svennine
 */

#include "XFEMDebugTools.h"

namespace oofem {
XFEMDebugTools :: XFEMDebugTools() { }

XFEMDebugTools :: ~XFEMDebugTools() { }

void XFEMDebugTools :: WriteTrianglesToVTK(const std :: string &iName, const std :: vector< Triangle > &iTriangles)
{
    //printf("Entering XFEMDebugTools::WriteTrianglesToVTK().\n");
    size_t numTri = iTriangles.size();


    std :: ofstream file;
    file.open( iName.data() );

    // Write header
    file << "# vtk DataFile Version 2.0\n";
    file << "Geometry of a PolygonLine\n";
    file << "ASCII\n";

    file << "DATASET UNSTRUCTURED_GRID\n";

    int numPoints = numTri * 3;
    // Write points
    file << "POINTS " << numPoints << "double\n";

    for ( size_t i = 0; i < numTri; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            const double &x = iTriangles.at(i).giveVertex(j).at(1);
            const double &y = iTriangles.at(i).giveVertex(j).at(2);
            file << x << " " << y << " 0.0\n";
        }
    }


    // Write segments
    file << "CELLS " << numTri << " " << numTri * 4 << "\n";

    for ( size_t i = 0; i < numTri; i++ ) {
        file << 3 << " " << 3 * i << " " << 3 * i + 1 << " " << 3 * i + 2 << "\n";
    }


    // Write cell types
    file << "CELL_TYPES " << numTri << "\n";
    int vtkCellType = 5;     // triangle
    for ( size_t i = 0; i < numTri; i++ ) {
        file << vtkCellType << "\n";
    }

    file.close();
}

void XFEMDebugTools :: WritePointsToVTK(const std :: string &iName, const std :: vector< FloatArray > &iPoints)
{
    //printf("Entering XFEMDebugTools::WriteTrianglesToVTK().\n");


    std :: ofstream file;
    file.open( iName.data() );

    // Write header
    file << "# vtk DataFile Version 2.0\n";
    file << "Gauss points\n";
    file << "ASCII\n";

    file << "DATASET UNSTRUCTURED_GRID\n";

    int numPoints = iPoints.size();
    // Write points
    file << "POINTS " << numPoints << " double\n";

    for ( int i = 1; i <= numPoints; i++ ) {
        //for(int j = 1; j <= 3; j++) {
        const double &x = iPoints [ i - 1 ].at(1);
        const double &y = iPoints [ i - 1 ].at(2);
        double z = 0.0;
        if ( iPoints [ i - 1 ].giveSize() == 3 ) {
            z = iPoints [ i - 1 ].at(3);
        }

        file << x << " " << y << " " << z << "\n";
        //}
    }


    // Write segments
    file << "CELLS " << numPoints << " " << numPoints * 2 << "\n";

    for ( int i = 0; i < numPoints; i++ ) {
        file << 1 << " " << i << "\n";
    }


    // Write cell types
    file << "CELL_TYPES " << numPoints << "\n";
    int vtkCellType = 1;             // vertex
    for ( int i = 0; i < numPoints; i++ ) {
        file << vtkCellType << "\n";
    }

    file.close();
}

void XFEMDebugTools :: WriteArrayToMatlab(const std :: string &iName, const std :: vector< double > &iX, const std :: vector< double > &iY)
{
    std :: ofstream file;
    file.open( iName.data() );

    file << "x = [";

    for ( size_t i = 0; i < iX.size(); i++ ) {
        file << iX [ i ] << "\n";
    }

    file << "];\n\n";



    file << "y = [";

    for ( size_t i = 0; i < iY.size(); i++ ) {
        file << iY [ i ] << "\n";
    }

    file << "];\n\n";

    file.close();
}

void XFEMDebugTools :: WriteArrayToGnuplot(const std :: string &iName, const std :: vector< double > &iX, const std :: vector< double > &iY)
{
    if ( iX.size() != iY.size() ) {
        OOFEM_ERROR("iX.size() != iY.size().")
    }

    std :: ofstream file;
    file.open( iName.data() );

    // Set some output options
    file << std :: scientific;

    file << "# x y\n";

    for ( size_t i = 0; i < iX.size(); i++ ) {
        file << iX [ i ] << " " << iY [ i ] << "\n";
    }

    file.close();
}
} /* namespace oofem */
