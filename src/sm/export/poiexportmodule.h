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

#ifndef poiexportmodule_h
#define poiexportmodule_h

#include "exportmodule.h"
#include "intarray.h"
#include "internalstatetype.h"
#include "unknowntype.h"

#include <list>
#include <memory>

///@name Input fields for Point-of-Interest export module
//@{
#define _IFT_POIExportModule_Name "poi"
#define _IFT_POIExportModule_vars "vars"
#define _IFT_POIExportModule_primvars "primvars"
#define _IFT_POIExportModule_mtype "mtype"
#define _IFT_POIExportModule_poifilename "poifilename"
//@}

namespace oofem {
class MaterialMappingAlgorithm;

/**
 * Represents POI (Point Of Interest) export module.
 * It is able to perform output on given points, which are inside domain but can have arbitrary position.
 */
class POIExportModule : public ExportModule
{
protected:
    /// POIs data structure.
    struct POI_dataType {
        int id;
        double x, y, z;
        // bp: region is now the set number (which is more flexible)
        int region;
    };

    /// List of InternalStateType values, identifying the selected vars for export.
    IntArray internalVarsToExport;
    /// List of primary unknowns to export.
    IntArray primaryVarsToExport;
    /// List of POIs.
    std :: list< POI_dataType >POIList;

    /// Smoother type.
    enum POIEM_MapperType { POI_CPT, POI_SFT, POI_LST } mtype;
    /// Mapper.
    std :: unique_ptr< MaterialMappingAlgorithm > mapper;
    //

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    POIExportModule(int n, EngngModel * e);
    /// Destructor
    virtual ~POIExportModule();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveInputRecordName() const { return _IFT_POIExportModule_Name; }
    virtual const char *giveClassName() const { return "POIExportModule"; }

protected:
    void readPOIFile(const std :: string &poiFileName);
    /// Returns the output stream for given solution step
    FILE *giveOutputStream(TimeStep *tStep);
    /**
     * Export internal variables.
     */
    void exportIntVars(FILE *stream, TimeStep *tStep);
    /**
     * Export primary variables.
     */
    void exportPrimaryVars(FILE *stream, TimeStep *tStep);
    /** Exports single variable */
    void exportIntVarAs(InternalStateType valID, FILE *stream, TimeStep *tStep);
    /** Exports single variable */
    void exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep);
    MaterialMappingAlgorithm *giveMapper();
};
} // end namespace oofem
#endif // poiexportmodule_h
