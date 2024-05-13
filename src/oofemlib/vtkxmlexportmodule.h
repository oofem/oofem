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

#ifndef vtkxmlexportmodule_h
#define vtkxmlexportmodule_h

#include "vtkbaseexportmodule.h"
#include "intarray.h"
#include "nodalrecoverymodel.h"
#include "interface.h"
#include "internalstatevaluetype.h"
#include "integrationrule.h"
#include "xfem/xfemmanager.h"
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef __VTK_MODULE
 #include <vtkUnstructuredGrid.h>
 #include <vtkSmartPointer.h>
#endif

#ifdef _PYBIND_BINDINGS
 #include <pybind11/pybind11.h>
 #include <pybind11/stl.h>   //Conversion for lists
 #include "pybind11/numpy.h"
namespace py = pybind11;
#endif

#ifdef _WIN32
 #define NULL_DEVICE "NUL:"
#else
 #define NULL_DEVICE "/dev/null"
#endif


#include <string>
#include <list>

///@name Input fields for VTK XML export module
//@{
#define _IFT_VTKXMLExportModule_Name "vtkxml"
#define _IFT_VTKXMLExportModule_cellvars "cellvars"
#define _IFT_VTKXMLExportModule_vars "vars"
#define _IFT_VTKXMLExportModule_primvars "primvars"
#define _IFT_VTKXMLExportModule_externalForces "externalforces"
#define _IFT_VTKXMLExportModule_ipvars "ipvars"
#define _IFT_VTKXMLExportModule_stype "stype"
//@}

using namespace std;
namespace oofem {
class Node;

/**
 * Represents VTK (Visualization Toolkit) export module. It uses VTK (.vtu) file format, Unstructured grid dataset.
 * The export of data is done on Region By Region basis, possibly taking care about possible nonsmooth character of
 * some internal variables at region boundaries.
 * Each region is usually exported as a single piece. When region contains composite cells, these are assumed to be
 * exported in individual subsequent pieces after the default one for the particular region.
 */
class OOFEM_EXPORT VTKXMLExportModule : public VTKBaseExportModule
{
protected:
    /// List of InternalStateType values, identifying the selected vars for export.
    IntArray internalVarsToExport;
    /// List of primary unknowns to export.
    IntArray primaryVarsToExport;
    /// List of primary unknowns to export.
    IntArray externalForcesToExport;
    /// List of cell data to export.
    IntArray cellVarsToExport;
    /// List of internal variables to export directly in Integration Points (no smoothing to nodes)
    IntArray ipInternalVarsToExport;

    

    /// Smoother type.
    NodalRecoveryModel::NodalRecoveryModelType stype;
    /// Smoother.
    std::unique_ptr< NodalRecoveryModel >smoother;
    /// Smoother for primary variables.
    std::unique_ptr< NodalRecoveryModel >primVarSmoother;

    /// Buffer for earlier time steps exported to *.pvd file.
    std::list< std::string >pvdBuffer;

    /// Buffer for earlier time steps with gauss points exported to *.gp.pvd file.
    std::list< std::string >gpPvdBuffer;


public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKXMLExportModule();

    void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    void initialize() override;
    void terminate() override;
    const char *giveClassName() const override { return "VTKXMLExportModule"; }
    /**
     * Prints point data header.
     */
    void exportPointDataHeader(std::ofstream fileStream, TimeStep *tStep);
    virtual void giveDataHeaders(std::string &pointHeader, std::string &cellHeader);     // returns the headers
    /// Returns the internal smoother.
    NodalRecoveryModel *giveSmoother();
    /// Returns the smoother for primary variables (nodal averaging).
    NodalRecoveryModel *givePrimVarSmoother();
    

#ifdef __VTK_MODULE
    vtkSmartPointer< vtkUnstructuredGrid >fileStream;
    vtkSmartPointer< vtkPoints >nodes;
    vtkSmartPointer< vtkIdList >elemNodeArray;

    vtkSmartPointer< vtkDoubleArray >intVarArray;
    vtkSmartPointer< vtkDoubleArray >primVarArray;
#else
    std::ofstream fileStream;
#endif

    ExportRegion defaultVTKPiece;

    std::vector< ExportRegion >defaultVTKPieces;
    
    ExportRegion& getExportRegions() {return this->defaultVTKPiece;}
protected:

    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    std::ofstream giveOutputStream(TimeStep *tStep);

    void writeIntVars(ExportRegion &vtkPiece);
    void writeXFEMVars(ExportRegion &vtkPiece);
    void writePrimaryVars(ExportRegion &vtkPiece);
    void writeCellVars(ExportRegion &vtkPiece);
    void writeExternalForces(ExportRegion &vtkPiece);
    /**
     * Writes Piece header+geometry
     * @return true if piece is not empty and thus written
     */
    bool writeVTKPieceEpilog(ExportRegion &vtkPiece, TimeStep *tStep);
    /**
     * Writes Piece variables
     * @return true if piece is not empty and thus written
     */
    bool writeVTKPieceVariables(ExportRegion &vtkPiece, TimeStep *tStep);
    /**
     * Writes piece epiloque
     * @return true if piece is not empty and thus written
     */
    bool writeVTKPieceProlog(ExportRegion &vtkPiece, TimeStep *tStep);
    /**
     * Exports given internal variables directly in integration points (raw data, no smoothing)
     * @param valIDs the UnknownType values identifying the internal variables to export
     * @param tStep solution step
     */
    void exportIntVarsInGpAs(IntArray valIDs, TimeStep *tStep);
    /**
     * Writes a VTK collection file where time step data is stored.
     */
    void writeVTKCollection();

    /// Writes a VTK collection file for Gauss points.
    void writeGPVTKCollection();

#ifdef __VTK_MODULE
    void writeVTKPointData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray);
#else
    void writeVTKPointData(FloatArray &valueArray);
#endif

#ifdef __VTK_MODULE
    void writeVTKCellData(const char *name, vtkSmartPointer< vtkDoubleArray >varArray);
#else
    void writeVTKCellData(FloatArray &valueArray);
#endif

    // Export of composite elements (built up from several subcells)
    void exportCompositeElement(ExportRegion &vtkPiece, Element *el, TimeStep *tStep);
    void exportCompositeElement(std::vector< ExportRegion > &vtkPieces, Element *el, TimeStep *tStep);
};

/**
 * Elements with geometry defined as EGT_Composite are exported using individual pieces.
 * The VTKXMLExportModuleElementInterface serves for this purpose, defining abstract
 * export method, responsible for exporting individual element piece in xml vtk syntax.
 * Elements with geometry defined as EGT_Composite should implement this interface.
 */
class OOFEM_EXPORT VTKXMLExportModuleElementInterface : public Interface
{

public:
    VTKXMLExportModuleElementInterface() : Interface() { }
    virtual void giveCompositeExportData(ExportRegion &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep) { }
    virtual void giveCompositeExportData(std::vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep) { }
};
} // end namespace oofem
#endif // vtkxmlexportmodule_h
