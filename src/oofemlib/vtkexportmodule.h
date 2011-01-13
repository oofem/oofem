/* $Header: /home/cvs/bp/oofem/sm/src/vtkexportmodule.h,v 1.9 2003/04/06 14:08:32 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//
// class vtkExportModule
//

#ifndef vtkexportmodule_h
#define vtkexportmodule_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif
#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"
#include "intarray.h"
#include "nodalrecoverymodel.h"

namespace oofem {
/**
 * Represents VTK (Visualization Toolkit) export module. It uses vtk file format, Unstructured grid dataset.
 * There is built in support for Region By Region otput, taking care about possible nonsmooth character of
 * some internal variables at region boundaries. This, hovewer, is rathrer complication and since application
 * of vtk is naturally in 3d, the corresponding sections are commented out.
 *
 */
class VTKExportModule : public ExportModule
{
protected:

    /// list of InternalStateType values, identifying the selected vars for export
    IntArray internalVarsToExport;
    /// list of primary unknowns to export
    IntArray primaryVarsToExport;
    /// list of cell data to export
    IntArray cellVarsToExport;

    enum modeType { wdmode, rbrmode }; // WholeDomain or RegionByRegion output
    modeType outMode;
    modeType mode;

    /// smoother type
    NodalRecoveryModel::NodalRecoveryModelType stype;
    /// smoother
    NodalRecoveryModel *smoother;
    /// list of regions to skip
    IntArray regionsToSkip;


public:

    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKExportModule(EngngModel *e);
    /// Destructor
    ~VTKExportModule();
    /// Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Writes the output. Abstract service.
     * @param tStep time step.
     */
    void              doOutput(TimeStep *tStep);
    /**
     * Initializes receiver.
     * The init file messages should be printed.
     */
    void              initialize();
    /**
     * Terminates the receiver.
     * The terminating messages should be printed.
     * All the streams should be closed.
     */
    void              terminate();
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "VTKExportModule"; }




protected:
    /// returns the internal smoother
    NodalRecoveryModel *giveSmoother();

    /// returns the output stream for given solution step
    FILE *giveOutputStream(TimeStep *);
    /**
     * Returns corresponding element cell_type.
     * Some common element types are supported, others can be supported via interface concept.
     */
    int   giveCellType(Element *);
    /**
     * Returns the number of elements vtk cells
     */
    int giveNumberOfElementCells(Element *);
    /**
     * Returns number of nodes correpsonding to cell type
     */
    int giveNumberOfNodesPerCell(int cellType);
    /**
     * Returns the element cell geometry.
     */
    void giveElementCell(IntArray &answer, Element *elem, int cell);
    /**
     * export internal variables
     */
    void exportIntVars(FILE *stream, TimeStep *tStep);
    /**
     * export primary variables
     */
    void exportPrimaryVars(FILE *stream, TimeStep *tStep);
    /** exports single variable */
    void exportIntVarAs(InternalStateType valID, InternalStateValueType type, FILE *stream, TimeStep *tStep);
    /** exports single variable */
    void exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep);
    /** export variables defined on cells */
    void exportCellVars(FILE *stream, int elemToProcess, TimeStep *tStep);

    /**
     * Assembles the region node map. Also computes the total number of nodes in region.
     * The region are numbered starting from offset+1.
     * if mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
     * The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).
     * if mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
     * The i-th value contains the corresponding global node number.
     */
    int initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans,
                                int offset, Domain *domain, int reg, int mode);
    /**
     * Computes total number of nodes (summed Region by Region, nodes on region boundaries are
     * added multiple times.
     */
    int giveTotalRBRNumberOfNodes(Domain *d);
    /**
     * Returns the value of Primary variable at given dof manager.
     * If such variable not directly available in dofman dofs, can use
     * smoother to recover this nodal value.
     */
    void getDofManPrimaryVariable(FloatArray &answer, DofManager *dman, IntArray &dofIDMask, EquationID type,
                                  ValueModeType mode, TimeStep *tStep, InternalStateType iType);
};
} // end namespace oofem
#endif // vtkexportmodule_h
