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
// class vtkXMLExportModule
//

#ifndef vtkxmlexportmodule_h
#define vtkxmlexportmodule_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"
#include "intarray.h"
#include "nodalrecoverymodel.h"
#include "interface.h"


namespace oofem {

/**
 * Represents VTK (Visualization Toolkit) export module. It uses vtk file format, Unstructured grid dataset.
 * The export of data is done on Region By Region basis, possibly taking care about possible nonsmooth character of
 * some internal variables at region boundaries. 
 * Each region is usually exported as a single piece. When region contains composite cells, these are assumed to be 
 * exported in individual subsequent pices after the default one for the particular region.
 *
 */
class VTKXMLExportModule : public ExportModule
{
protected:

    /// list of InternalStateType values, identifying the selected vars for export
    IntArray internalVarsToExport;
    /// list of primary unknowns to export
    IntArray primaryVarsToExport;
    /// list of cell data to export
    IntArray cellVarsToExport;

    /// smoother type
    enum VTKEM_SmootherType { VTK_Smother_NA, VTK_Smoother_ZZ, VTK_Smoother_SPR } stype;
    /// smotther
    NodalRecoveryModel *smoother;
    /// list of regions to skip
    IntArray regionsToSkip;


public:

    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLExportModule(EngngModel *e);
    /// Destructor
    ~VTKXMLExportModule();
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
    virtual const char *giveClassName() const { return "VTKXMLExportModule"; }




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
    /// Prints point data header 
    void exportPointDataHeader (FILE* stream, TimeStep* tStep);
    /// Prints point data footer 
    void exportPointDataFooter (FILE* stream, TimeStep* tStep);
    /**
     * export internal variables
     */
    void exportIntVars(FILE *stream, IntArray& mapG2L, IntArray& mapL2G,
               int regionDofMans, int ireg, TimeStep *tStep);
    /**
     * export primary variables
     */
    void exportPrimaryVars(FILE *stream, IntArray& mapG2L, IntArray& mapL2G,
               int regionDofMans, int region, TimeStep *tStep);
    /** exports single variable */
    void exportIntVarAs(InternalStateType valID, InternalStateValueType type, IntArray& mapG2L, IntArray& mapL2G,
            int regionDofMans, int ireg, FILE *stream, TimeStep *tStep);
    /** exports single variable */
    void exportPrimVarAs(UnknownType valID, IntArray& mapG2L, IntArray& mapL2G,
             int regionDofMans, int region, FILE *stream, TimeStep *tStep);
    
    /** exports cell variables */
    void exportCellVars(FILE *stream, int totalcells, int region, TimeStep *tStep);
    void exportCellVarAs(InternalStateType type, int nelem, int region, FILE *stream, TimeStep *tStep);

    /**
     * Assembles the region node map. Also computes the total number of nodes in region.
     * The region are numbered starting from offset+1.
     * if mode == 0 then regionNodalNumbers is array with mapping from global numbering to local region numbering.
     * The i-th value contains the corresponding local region number (or zero, if global numbar is not in region).
     * if mode == 1 then regionNodalNumbers is array with mapping from local to global numbering.
     * The i-th value contains the corresponding global node number.
     */
    int initRegionNodeNumbering(IntArray &mapG2L, IntArray& mapL2G, 
                int &regionDofMans, int &totalcells,
                                Domain *domain, int reg);

    /// Returns true if element geometry type is composite (not a single cell)
    bool isElementComposite(Element* elem);
};

/**
   Elements with geometry defined as EGT_Composite are exported using individual pieces.
   The VTKXMLExportModuleElementInterface serves for this purpose, defining abstract
   export method, responsible for exporting individual element piece in xml vtk syntax.
   Elements with geometry defined as EGT_Composite should implement this interface.
 */
class VTKXMLExportModuleElementInterface : public Interface {
 public:
  VTKXMLExportModuleElementInterface () : Interface() {}
  const char *giveClassName() const { return "VTKXMLExportModuleElementInterface"; }
  virtual void _export (FILE * stream, VTKXMLExportModule* m, IntArray& primaryVarsToExport, IntArray& internalVarsToExport, TimeStep* tStep ) = 0;
};

} // end namespace oofem
#endif // vtkxmlexportmodule_h
