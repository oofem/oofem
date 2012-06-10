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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef vtkxmlexportmodule_h
#define vtkxmlexportmodule_h

#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"
#include "intarray.h"
#include "nodalrecoverymodel.h"
#include "interface.h"
#include "internalstatevaluetype.h"

#ifdef __VTK_MODULE
 #include <vtkUnstructuredGrid.h>
 #include <vtkSmartPointer.h>
#endif


namespace oofem {
/**
 * Represents VTK (Visualization Toolkit) export module. It uses VTK (.vtu) file format, Unstructured grid dataset.
 * The export of data is done on Region By Region basis, possibly taking care about possible nonsmooth character of
 * some internal variables at region boundaries.
 * Each region is usually exported as a single piece. When region contains composite cells, these are assumed to be
 * exported in individual subsequent pieces after the default one for the particular region.
 */
class VTKXMLExportModule : public ExportModule
{
protected:
    /// List of InternalStateType values, identifying the selected vars for export.
    IntArray internalVarsToExport;
    /// List of primary unknowns to export.
    IntArray primaryVarsToExport;
    /// List of cell data to export.
    IntArray cellVarsToExport;

    /// Smoother type.
    NodalRecoveryModel::NodalRecoveryModelType stype;
    /// Smoother.
    NodalRecoveryModel *smoother;
    /// Smoother for primary variables.
    NodalRecoveryModel *primVarSmoother;
    /// List of regions to skip.
    IntArray regionsToSkip;
    /// Number of virtual regions.
    int nvr;
    /// Real->virtual region map.
    IntArray vrmap;

    /// Buffer for earlier time steps exported to *.pvd file.
    dynaList< std::string > pvdBuffer;

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLExportModule(int n, EngngModel *e);
    /// Destructor
    ~VTKXMLExportModule();

    virtual IRResultType initializeFrom(InputRecord *ir);
    void doOutput(TimeStep *tStep);
    void initialize();
    void terminate();
    virtual const char *giveClassName() const { return "VTKXMLExportModule"; }

protected:
    /// Returns the internal smoother.
    NodalRecoveryModel *giveSmoother();
    /// Returns the smoother for primary variables (nodal averaging).
    NodalRecoveryModel *givePrimVarSmoother();

    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    FILE *giveOutputStream(TimeStep *tStep);
    /**
     * Returns corresponding element cell_type.
     * Some common element types are supported, others can be supported via interface concept.
     */
    int giveCellType(Element *element);
    /**
     * Returns the number of elements vtk cells.
     */
    int giveNumberOfElementCells(Element *element);
    /**
     * Returns number of nodes corresponding to cell type
     */
    int giveNumberOfNodesPerCell(int cellType);
    /**
     * Returns the element cell geometry.
     */
    void giveElementCell(IntArray &answer, Element *elem, int cell);
#ifndef __VTK_MODULE
    /**
     * Prints point data header.
     */
    void exportPointDataHeader(FILE *stream, TimeStep *tStep);
#endif
    /**
     * Export internal variables by smoothing.
     */
    void exportIntVars(
#ifdef __VTK_MODULE
        vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
        FILE *stream,
#endif
        IntArray &mapG2L, IntArray &mapL2G, int regionDofMans, int ireg, TimeStep *tStep);
    /**
     * Export primary variables.
     */
    void exportPrimaryVars(
#ifdef __VTK_MODULE
        vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
        FILE *stream,
#endif
        IntArray &mapG2L, IntArray &mapL2G, int regionDofMans, int region, TimeStep *tStep);
    /**
     * Tries to find the value of a primary field on the given DofManager.
     * Some elements have different interpolation of some fields, and requires some additional code to compute node values (if available).
     */
    void getPrimaryVariable(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int ireg);
    /**
     * Exports single internal variable by smoothing.
     */
    void exportIntVarAs(InternalStateType valID, InternalStateValueType type, IntArray &mapG2L, IntArray &mapL2G,
                        int regionDofMans, int ireg,
#ifdef __VTK_MODULE
                        vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
                        FILE *stream,
#endif
                        TimeStep *tStep);
    /**
     * Exports single primary variable.
     */
    void exportPrimVarAs(UnknownType valID, IntArray &mapG2L, IntArray &mapL2G,
                         int regionDofMans, int region,
#ifdef __VTK_MODULE
                         vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
                         FILE *stream,
#endif
                         TimeStep *tStep);

    /**
     * Exports cell variables (typically internal variables).
     */
    void exportCellVars(
#ifdef __VTK_MODULE
        vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
        FILE *stream,
#endif
        int region, TimeStep *tStep);
    /**
     * Exports a single cell variable (typically an internal variable).
     */
    void exportCellVarAs(InternalStateType type, int region,
#ifdef __VTK_MODULE
                         vtkSmartPointer<vtkUnstructuredGrid> &stream,
#else
                         FILE *stream,
#endif
                         TimeStep *tStep);

    /**
     * Assembles the region node map. Also computes the total number of nodes in region.
     * The region are numbered starting from offset+1.
     * If mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
     * The i-th value contains the corresponding local region number (or zero, if global number is not in region).
     * If mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
     * The i-th value contains the corresponding global node number.
     */
    int initRegionNodeNumbering(IntArray &mapG2L, IntArray &mapL2G,
                                int &regionDofMans, int &totalcells,
                                Domain *domain, int reg);

    /// Returns true if element geometry type is composite (not a single cell).
    bool isElementComposite(Element *elem);

    /**
     * Writes a VTK collection file where time step data is stored.
     */
    void writeVTKCollection();
};

/**
 * Elements with geometry defined as EGT_Composite are exported using individual pieces.
 * The VTKXMLExportModuleElementInterface serves for this purpose, defining abstract
 * export method, responsible for exporting individual element piece in xml vtk syntax.
 * Elements with geometry defined as EGT_Composite should implement this interface.
 */
class VTKXMLExportModuleElementInterface : public Interface
{
public:
    VTKXMLExportModuleElementInterface() : Interface() {}
    const char *giveClassName() const { return "VTKXMLExportModuleElementInterface"; }
    virtual void _export(FILE *stream, VTKXMLExportModule *m, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, TimeStep *tStep) = 0;
};
} // end namespace oofem
#endif // vtkxmlexportmodule_h
