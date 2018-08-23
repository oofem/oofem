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

#ifndef t3dinterface_h
#define t3dinterface_h

#include "mesherinterface.h"
#include "element.h"

#define _IFT_T3DInterface_Name "t3d"

namespace oofem {
#define BMF_FILENAME "t3d.bmf"

class TimeStep;

/**
 * This class represents the interface to t3d mesh generation package.
 * This interface is primarily responsible for two main tasks:
 * - to create input mesher file, containing all information including the mesh density information
 *   based on information from remeshing criteria.
 * - possibly to launch the mesher and transform its output to oofem input (using t3d2oofem)
 */
class OOFEM_EXPORT T3DInterface : public MesherInterface
{
public:
    /// Constructor
    T3DInterface(Domain * d) : MesherInterface(d) { }
    /// Destructor
    virtual ~T3DInterface() { }

    returnCode createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew) override;
    int t3d_2_OOFEM(const char *t3dOutFile, Domain **dNew);
   // Used by HTS elemnt to mesh one element in order to export the HTS element into vtk
    int createInput(Element *e, char *t3dInFile);
    int createVTKExportMesh(const char *t3dOutFile,std::vector<FloatArray> &nodeCoords, std::vector<IntArray> &cellNodes, IntArray &cellTypes );
   // udes in Quasicontinuum    
    int createQCInterpolationMesh(const char *t3dOutFile,std::vector<FloatArray> &nodeCoords, std::vector<IntArray> &cellNodes, IntArray &cellTypes );

protected:
    /// Creates the mesher input, containing the required mesh density information.
    int createInput(Domain *d, TimeStep *tStep);
};
} // end namespace oofem
#endif // t3dinterface_h
