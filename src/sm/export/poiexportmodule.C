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

#include "poiexportmodule.h"
#include "timestep.h"
#include "domain.h"
#include "engngm.h"
#include "materialmappingalgorithm.h"
#include "mmaclosestiptransfer.h"
#include "mmaleastsquareprojection.h"
#include "mmashapefunctprojection.h"
#include "spatiallocalizer.h"
#include "internalstatevaluetype.h"
#include "element.h"
#include "classfactory.h"

#include <string>
#include <fstream>
#include <ios>

namespace oofem {
REGISTER_ExportModule(POIExportModule)

POIExportModule :: POIExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), primaryVarsToExport(), POIList()
{
}


POIExportModule :: ~POIExportModule()
{
}


IRResultType
POIExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_POIExportModule_vars);
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_POIExportModule_primvars);

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_POIExportModule_mtype);
    mtype = ( POIEM_MapperType ) val;

    std :: string poiFileName;
    IR_GIVE_OPTIONAL_FIELD(ir, poiFileName, _IFT_POIExportModule_poifilename);
    this->readPOIFile(poiFileName); // parse poi file

    return ExportModule :: initializeFrom(ir);
}

void
POIExportModule :: readPOIFile(const std :: string &poiFileName)
{
    POI_dataType poi;
    int nPOI;
    // Open the file;
    std :: ifstream file(poiFileName.c_str(), std :: ios :: in);
    if ( !file.is_open() ) {
        OOFEM_ERROR("Failed to open time data file: %s\n", poiFileName.c_str() );
    }

    file >> nPOI; // Not actually needed.

    while ( file >>  poi.id >> poi.x >> poi.y >> poi.z >> poi.region ) {
        POIList.push_back(poi);
    }
}


void
POIExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    FILE *stream = this->giveOutputStream(tStep);

    fprintf(stream, "# POI DataFile\n");
    fprintf( stream, "Output for time %f\n", tStep->giveTargetTime() );

    this->exportPrimaryVars(stream, tStep);
    this->exportIntVars(stream, tStep);

    fclose(stream);
}

void
POIExportModule :: initialize()
{ 
    ExportModule :: initialize();
}


void
POIExportModule :: terminate()
{ }


FILE *
POIExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;
    std :: string fileName = this->giveOutputBaseFileName(tStep) + ".poi";

    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    return answer;
}

void
POIExportModule :: exportIntVars(FILE *stream, TimeStep *tStep)
{
    int i, n = internalVarsToExport.giveSize();
    InternalStateType type;
    FloatArray poiCoords(3);

    if ( n == 0 ) {
        return;
    }

    // loop over POIs
    POI_dataType &poi = *POIList.begin();
    poiCoords.at(1) = poi.x;
    poiCoords.at(2) = poi.y;
    poiCoords.at(3) = poi.z;
    //int region = poi.region;

    for ( i = 1; i <= n; i++ ) {
        type = ( InternalStateType ) internalVarsToExport.at(i);
        fprintf(stream, "\n\nPOI_INTVAR_DATA %d\n", type);
        this->exportIntVarAs(type, stream, tStep);
    }

    this->giveMapper()->finish(tStep);
}


void
POIExportModule :: exportIntVarAs(InternalStateType valID, FILE *stream, TimeStep *tStep)
{
    int region;
    IntArray toMap(1);
    Domain *d = emodel->giveDomain(1);
    FloatArray poiCoords(3);
    FloatArray val;

    toMap.at(1) = ( int ) valID;

    // loop over POIs
    for ( auto &poi: POIList ) {
        poiCoords.at(1) = poi.x;
        poiCoords.at(2) = poi.y;
        poiCoords.at(3) = poi.z;
        region = poi.region;

        this->giveMapper()->__init(d, toMap, poiCoords, * d->giveSet(region), tStep);
        if ( !this->giveMapper()->__mapVariable(val, poiCoords, valID, tStep) ) {
            OOFEM_WARNING("Failed to map variable");
            val.clear();
        }
        fprintf(stream, "%10d ", poi.id);
        for ( auto &x : val ) {
            fprintf( stream, " %15e", x );
        }

        fprintf(stream, "\n");
    }
}


MaterialMappingAlgorithm *
POIExportModule :: giveMapper()
{
    if ( !this->mapper ) {
        if ( this->mtype == POI_CPT ) {
            this->mapper.reset( new MMAClosestIPTransfer() );
        } else if ( this->mtype == POI_SFT ) {
            this->mapper.reset( new MMAShapeFunctProjection() );
        } else if ( this->mtype == POI_LST ) {
            this->mapper.reset( new MMALeastSquareProjection() );
        } else {
            OOFEM_ERROR("unsupported smoother type ID");
        }
    }

    return this->mapper.get();
}


void
POIExportModule :: exportPrimaryVars(FILE *stream, TimeStep *tStep)
{
    // should be performed over regions

    int n = primaryVarsToExport.giveSize();
    Domain *d = emodel->giveDomain(1);

    if ( n == 0 ) {
        return;
    }

    fprintf(stream, "\n\nPOINT_DATA %d\n", d->giveNumberOfDofManagers());

    for ( int i = 1; i <= n; i++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(i);
        this->exportPrimVarAs(type, stream, tStep);
    }
}


void
POIExportModule :: exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray pv, coords(3), lcoords, closest;
    InternalStateValueType type = ISVT_UNDEFINED;

    if ( valID == DisplacementVector ) {
        type = ISVT_VECTOR;
    } else if ( valID == FluxVector || valID == Humidity ) {
        type = ISVT_SCALAR;
    } else {
        OOFEM_ERROR("unsupported UnknownType");
    }

    // print header
    if ( type == ISVT_SCALAR ) {
        fprintf(stream, "SCALARS prim_scalar_%d\n", ( int ) valID);
    } else if ( type == ISVT_VECTOR ) {
        fprintf(stream, "VECTORS vector_%d float\n", ( int ) valID);
    } else {
        OOFEM_ERROR("unsupported variable type");
    }


    SpatialLocalizer *sl = d->giveSpatialLocalizer();
    // loop over POIs
    for ( auto &poi: POIList ) {
        coords.at(1) = poi.x;
        coords.at(3) = poi.z;
        //region = poi.region;

        Element *source = sl->giveElementClosestToPoint(lcoords, closest, coords);
        if ( source ) {
            // ask interface
            source->computeField(VM_Total, tStep, lcoords, pv);

            fprintf(stream, "%10d ", poi.id);
            for ( auto &p : pv ) {
                fprintf( stream, " %15e ", p );
            }

            fprintf(stream, "\n");
        } else {
            OOFEM_ERROR("no element containing POI(%e,%e,%e) found",
                        coords.at(1), coords.at(2), coords.at(3) );
        }
    }
}
} // end namespace oofem
