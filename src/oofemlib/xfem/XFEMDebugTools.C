/*
 * XFEMDebugTools.C
 *
 *  Created on: Jun 5, 2013
 *      Author: svennine
 */

#include "XFEMDebugTools.h"

namespace oofem {

XFEMDebugTools::XFEMDebugTools() {

}

XFEMDebugTools::~XFEMDebugTools() {

}

void XFEMDebugTools::WriteTrianglesToVTK( const std::string &iName, const AList< Triangle > &iTriangles )
{
//	printf("Entering XFEMDebugTools::WriteTrianglesToVTK().\n");
	int numTri = iTriangles.giveSize();


	std::ofstream file;
	file.open (iName.data());

	// Write header
	file << "# vtk DataFile Version 2.0\n";
	file << "Geometry of a PolygonLine\n";
	file << "ASCII\n";

	file << "DATASET UNSTRUCTURED_GRID\n";

	int numPoints = numTri*3;
	// Write points
	file << "POINTS " << numPoints << "double\n";

	for(int i = 1; i <= numTri; i++)
	{
		for(int j = 1; j <= 3; j++)
		{
			const double &x = iTriangles.at(i)->giveVertex(j)->at(1);
			const double &y = iTriangles.at(i)->giveVertex(j)->at(2);
			file << x << " " << y << " 0.0\n";
		}
	}


	// Write segments
	file << "CELLS " << numTri << " " << numTri*4 << "\n";

	for(int i = 0; i < numTri; i++)
	{
		file << 3 << " " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
	}


	// Write cell types
	file << "CELL_TYPES " << numTri << "\n";
	int vtkCellType = 5; // triangle
	for(int i = 0; i < numTri; i++)
	{
		file << vtkCellType << "\n";
	}

	file.close();


}


} /* namespace oofem */

