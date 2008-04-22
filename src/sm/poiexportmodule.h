/* $Header: /home/cvs/bp/oofem/sm/src/poiexportmodule.h,v 1.1 2003/05/19 13:04:00 bp Exp $ */
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
// class poiExportModule
//

#ifndef poiexportmodule_h
#define poiexportmodule_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"
#include "intarray.h"
#include "materialmappingalgorithm.h"
#include "dynalist.h"

/**
 * Represents POI (Point Of Interest) export module.
 * It is able to perform output on given points, which are inside domain but can have arbitrary position.
 */
class POIExportModule : public ExportModule
{
protected:
    /// POIs data structure
    struct POI_dataType {
        int id;
        double x, y, z;
        int region;
    };

    /// list of InternalStateType values, identifying the selected vars for export
    IntArray internalVarsToExport;
    /// list of primary unknowns to export
    IntArray primaryVarsToExport;
    /// list of POIs
    dynaList< POI_dataType >POIList;

    /// smoother type
    enum POIEM_MapperType { POI_CPT, POI_SFT, POI_LST } mtype;
    /// Mapper
    MaterialMappingAlgorithm *mapper;

public:

    /// Constructor. Creates empty Output Manager. By default all components are selected.
    POIExportModule(EngngModel *e);
    /// Destructor
    ~POIExportModule();
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
    virtual const char *giveClassName() const { return "POIExportModule"; }


protected:
    void readPOIFile(char *poiFileName);
    /// returns the output stream for given solution step
    FILE *giveOutputStream(TimeStep *);
    /**
     * export internal variables
     */
    void exportIntVars(FILE *stream, TimeStep *tStep);
    /**
     * export primary variables
     */
    void exportPrimaryVars(FILE *stream, TimeStep *tStep);
    /** exports single variable */
    void exportIntVarAs(InternalStateType valID, FILE *stream, TimeStep *tStep);
    /** exports single variable */
    void exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep);
    MaterialMappingAlgorithm *giveMapper();
};

#endif // poiexportmodule_h



