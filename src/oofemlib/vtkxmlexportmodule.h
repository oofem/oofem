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

#include "exportmodule.h"
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
//#define _IFT_VTKXMLExportModule_particleexportflag "particleexportflag"
//@}

using namespace std;
namespace oofem {
class Node;

///@todo Rename this to something like "ExportPiece" and move it to a separate file (it doesn't actually contain anything VTK-specific).
/// Stores all neccessary data (of a region) in a VTKPiece so it can be exported later.
class OOFEM_EXPORT VTKPiece
{
public:
    VTKPiece()
    {
        numCells = 0;
        numNodes = 0;
    }

    void clear();

    void setNumberOfNodes(int numNodes);
    int giveNumberOfNodes() { return this->numNodes; }

    void setNumberOfCells(int numCells);
    int giveNumberOfCells() { return this->numCells; }

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

    IntArray& getMapG2L () {return this->mapG2L;}
    IntArray& getMapL2G () {return this->mapL2G;}
    //void setRegionCells(IntArray& cells) {this->regionElInd = cells;}
    IntArray& getRegionCells () {return this->regionElInd;}
private:
    int numCells;
    int numNodes;
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
};

/**
 * Base class for VTK related export modules. Defines commmon methods.
 */
class OOFEM_EXPORT VTKBaseExportModule : public ExportModule
{
protected:
    /// Map from Voigt to full tensor.
    static IntArray redToFull;

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKBaseExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKBaseExportModule();

    void initialize() override;
    void terminate() override;
    const char *giveClassName() const override { return "VTKBaseExportModule"; }
    /**
     * Computes a cell average of an InternalStateType varible based on the weights
     * in the integrationpoints (=> volume/area/length average)
     */
    static void computeIPAverage(FloatArray &answer, IntegrationRule *iRule, Element *elem,  InternalStateType isType, TimeStep *tStep);

    
protected:

    virtual void setupVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep, Set& region);
       
    /**
     * Export primary variables.
     */
    virtual void exportPrimaryVars(VTKPiece &piece, Set& region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep);
    /**
     * Export internal variables by smoothing.
     */
    virtual void exportIntVars(VTKPiece &piece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep);
    /**
     * Export external forces.
     */
    void exportExternalForces(VTKPiece &piece, int region, TimeStep *tStep);
    
    void exportXFEMVarAs(XFEMStateType xfemstype, int regionDofMans, int ireg, TimeStep *tStep, EnrichmentItem *ei);
    ///  Exports cell variables (typically internal variables).
    void exportCellVars(VTKPiece &piece, Set& region, IntArray &cellVarsToExport, TimeStep *tStep);
    /**
     * Export external forces.
     */
    void exportExternalForces(VTKPiece &piece, Set& region, IntArray& externalForcesToExport, TimeStep *tStep);


    //  Tries to find the value of a primary field on the given DofManager.
    //  Some elements have different interpolation of some fields, and requires some additional code to compute node values (if available).
    void getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, Set& region, NodalRecoveryModel& smoother);
    //
    //  Exports single internal variable by smoothing.
    //
    void getNodalVariableFromIS(FloatArray &answer, Node *node, TimeStep *tStep, InternalStateType type, Set& region, NodalRecoveryModel& smoother);
    void getNodalVariableFromXFEMST(FloatArray &answer, Node *node, TimeStep *tStep, XFEMStateType xfemstype,Set& region, EnrichmentItem *ei);
    //
    //  Exports a single cell variable (typically an internal variable).
    //
    void getCellVariableFromIS(FloatArray &answer, Element *el, InternalStateType type, TimeStep *tStep);
 

    /// Gives the full form of given symmetrically stored tensors, missing components are filled with zeros.
    static void makeFullTensorForm(FloatArray &answer, const FloatArray &reducedForm, InternalStateValueType vtype);
    /**
     * Returns corresponding element cell_type.
     * Some common element types are supported, others can be supported via interface concept.
    */
    int giveCellType(Element *element);
    /**
     * Returns number of nodes corresponding to cell type
     */
    int giveNumberOfNodesPerCell(int cellType);
    /**
     * Returns the element cell geometry.
     */
    void giveElementCell(IntArray &answer, Element *elem);
    /**
     * Assembles the region node map. Also computes the total number of nodes in region.
     * The region are numbered starting from offset+1.
     * If mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
     * The i-th value contains the corresponding local region number (or zero, if global number is not in region).
     * If mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
     * The i-th value contains the corresponding global node number.
     */
    virtual int initRegionNodeNumbering(VTKPiece& vtkPiece,
                                        Domain *domain, TimeStep *tStep, Set& region); 

    // Export of composite elements (built up from several subcells)
    bool isElementComposite(Element *elem); /// Returns true if element geometry type is composite (not a single cell).
    void exportCompositeElement(VTKPiece &vtkPiece, Element *el, TimeStep *tStep);
    void exportCompositeElement(std::vector< VTKPiece > &vtkPieces, Element *el, TimeStep *tStep);
};

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

#ifdef _PYBIND_BINDINGS
    ///Dictionaries used for Python export
    py::dict Py_PrimaryVars, Py_IntVars, Py_CellVars, Py_Nodes, Py_Elements;
#endif

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

    VTKPiece defaultVTKPiece;

    std::vector< VTKPiece >defaultVTKPieces;

#ifdef _PYBIND_BINDINGS
    void clearPyList(py::list &list) {py::list tmp; list = tmp; }
    
    py::dict getPrimaryVars(void) { return Py_PrimaryVars; }
    py::dict getInternalVars(void) { return Py_IntVars; }
    py::dict getCellVars(void) { return Py_CellVars; }
    py::dict getNodes(void) { return Py_Nodes; } //from all VTKpieces sorted in dictionary
    py::dict getElementsConnectivity(void) { return Py_Elements; }
    VTKPiece& getVTKPieces() {return this->defaultVTKPiece;}
    

#endif


protected:

    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    std::ofstream giveOutputStream(TimeStep *tStep);

    void writeIntVars(VTKPiece &vtkPiece);
    void writeXFEMVars(VTKPiece &vtkPiece);
    void writePrimaryVars(VTKPiece &vtkPiece);
    void writeCellVars(VTKPiece &vtkPiece);
    void writeExternalForces(VTKPiece &vtkPiece);
    /**
     * Writes Piece header+geometry
     * @return true if piece is not empty and thus written
     */
    bool writeVTKPieceEpilog(VTKPiece &vtkPiece, TimeStep *tStep);
    /**
     * Writes Piece variables
     * @return true if piece is not empty and thus written
     */
    bool writeVTKPieceVariables(VTKPiece &vtkPiece, TimeStep *tStep);
    /**
     * Writes piece epiloque
     * @return true if piece is not empty and thus written
     */
    bool writeVTKPieceProlog(VTKPiece &vtkPiece, TimeStep *tStep);
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
    void exportCompositeElement(VTKPiece &vtkPiece, Element *el, TimeStep *tStep);
    void exportCompositeElement(std::vector< VTKPiece > &vtkPieces, Element *el, TimeStep *tStep);
};

class OOFEM_EXPORT VTKPFEMXMLExportModule : public VTKBaseExportModule
{
protected:
    std::ofstream fileStream;
public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKPFEMXMLExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKPFEMXMLExportModule();

    void initializeFrom(InputRecord &ir) override;

    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
protected:
    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    std::ofstream giveOutputStream(TimeStep *tStep);

};

/** 
 * VTK-like export class intended mainly to provide easy to use Pythonic interface by providing 
 * acccess to VTKPieces.
 */  
class OOFEM_EXPORT VTKMemoryExportModule : public VTKBaseExportModule
{
protected:
    std::vector< VTKPiece > vtkPieces;
public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKMemoryExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKMemoryExportModule();
    void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    std::vector< VTKPiece>& getVTKPieces(); 
protected:
};

///@name Input fields for VTK XML export module
//@{
#define _IFT_VTKXMLXFemExportModule_Name "vtkxmlxfem"
//@}

/**
 * Represents VTK (Visualization Toolkit) export module for Xfem. 
 * It uses VTK (.vtu) file format, Unstructured grid dataset.
 * The export of data is done on Region By Region basis, possibly taking care about possible
 * nonsmooth character of some internal variables at region boundaries.
 * The exported variables are dermined by XFemManager (FemMan->vtkExportFields).
 * Each region and each enrichment item is exported as a single piece. 
 * When region contains composite cells, these are assumed to be
 * exported in individual subsequent pieces after the default one for the particular region.
 */
class OOFEM_EXPORT VTKXMLXFemExportModule : public VTKXMLExportModule
{
protected:

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLXFemExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKXMLXFemExportModule();

    void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    void terminate() override {}
    const char *giveClassName() const override { return "VTKXMLXFemExportModule"; }

protected:

    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    std::ofstream giveOutputStream(TimeStep *tStep);

    bool writeXFEMVars(VTKPiece &vtkPiece, int field, int enrItIndex);
    void getNodalVariableFromXFEMST(FloatArray &answer, Node *node, TimeStep *tStep, XFEMStateType xfemstype, Set &region, EnrichmentItem *ei);
    void exportIntVars(VTKPiece &vtkPiece, Set& region, int field, int enrItIndex,  IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep);
    void giveDataHeaders(std::string &pointHeader, std::string &cellHeader) override;     // returns the headers
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
    virtual void giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep) { }
    virtual void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep) { }
};
} // end namespace oofem
#endif // vtkxmlexportmodule_h
