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
#include "engngm.h"
#include "mmaclosestiptransfer.h"
#include "mmaleastsquareprojection.h"
#include "mmashapefunctprojection.h"
#include "spatiallocalizer.h"
#include "eleminterpmapperinterface.h"
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
    mapper = NULL;
}


POIExportModule :: ~POIExportModule()
{
    if ( this->mapper ) {
        delete this->mapper;
    }
}


IRResultType
POIExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    ExportModule :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_POIExportModule_vars);
    IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_POIExportModule_primvars);

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_POIExportModule_mtype);
    mtype = ( POIEM_MapperType ) val;

    std :: string poiFileName;
    IR_GIVE_OPTIONAL_FIELD(ir, poiFileName, _IFT_POIExportModule_poifilename);
    this->readPOIFile(poiFileName); // parse poi file

    return IRRT_OK;
}

void
POIExportModule :: readPOIFile(const std :: string &poiFileName)
{
    POI_dataType poi;
    int nPOI;
    // Open the file;
    std :: ifstream file(poiFileName.c_str(), std :: ios :: in);
    if ( !file.is_open() ) {
        OOFEM_ERROR2( "POIExportModule :: readPOIFile - Failed to open time data file: %s\n", poiFileName.c_str() );
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
{ }


void
POIExportModule :: terminate()
{ }


FILE *
POIExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;
    std :: string fileName = this->giveOutputBaseFileName(tStep) + ".poi";

    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2( "POIExportModule::giveOutputStream: failed to open file %s", fileName.c_str() );
    }

    return answer;
}

void
POIExportModule :: exportIntVars(FILE *stream, TimeStep *tStep)
{
    int i, region, n = internalVarsToExport.giveSize();
    Domain *d = emodel->giveDomain(1);
    InternalStateType type;
    FloatArray poiCoords(3);

    if ( n == 0 ) {
        return;
    }

    // loop over POIs
    std :: list< POI_dataType > :: iterator PoiIter = POIList.begin();
    poiCoords.at(1) = ( * PoiIter ).x;
    poiCoords.at(2) = ( * PoiIter ).y;
    poiCoords.at(3) = ( * PoiIter ).z;
    region = ( * PoiIter ).region;

    this->giveMapper()->__init(d, internalVarsToExport, poiCoords, region, tStep);

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
    int i, region;
    IntArray toMap(1);
    Domain *d = emodel->giveDomain(1);
    FloatArray poiCoords(3);
    FloatArray val;

    toMap.at(1) = ( int ) valID;

    // loop over POIs
    std :: list< POI_dataType > :: iterator PoiIter;
    for ( PoiIter = POIList.begin(); PoiIter != POIList.end(); ++PoiIter ) {
        poiCoords.at(1) = ( * PoiIter ).x;
        poiCoords.at(2) = ( * PoiIter ).y;
        poiCoords.at(3) = ( * PoiIter ).z;
        region = ( * PoiIter ).region;

        this->giveMapper()->__init(d, toMap, poiCoords, region, tStep);
        if ( !this->giveMapper()->__mapVariable(val, poiCoords, valID, tStep) ) {
            OOFEM_WARNING("POIExportModule :: exportIntVarAs - Failed to map variable");
            val.clear();
        }
        fprintf(stream, "%10d ", ( * PoiIter ).id);
        for ( i = 1; i <= val.giveSize(); i++ ) {
            fprintf( stream, " %15e", val.at(i) );
        }

        fprintf(stream, "\n");
    }
}


MaterialMappingAlgorithm *
POIExportModule :: giveMapper()
{
    if ( this->mapper == NULL ) {
        if ( this->mtype == POI_CPT ) {
            this->mapper  = new MMAClosestIPTransfer();
        } else if ( this->mtype == POI_SFT ) {
            this->mapper = new MMAShapeFunctProjection();
        } else if ( this->mtype == POI_LST ) {
            this->mapper = new MMALeastSquareProjection();
        } else {
            OOFEM_ERROR("POIExportModule: unsupported smoother type ID");
        }
    }

    return this->mapper;
}


void
POIExportModule :: exportPrimaryVars(FILE *stream, TimeStep *tStep)
{
    // should be performed over regions

    int n = primaryVarsToExport.giveSize();
    int nnodes;
    Domain *d = emodel->giveDomain(1);
    UnknownType type;

    if ( n == 0 ) {
        return;
    }

    nnodes = d->giveNumberOfDofManagers();

    fprintf(stream, "\n\nPOINT_DATA %d\n", nnodes);

    for ( int i = 1; i <= n; i++ ) {
        type = ( UnknownType ) primaryVarsToExport.at(i);
        this->exportPrimVarAs(type, stream, tStep);
    }
}


void
POIExportModule :: exportPrimVarAs(UnknownType valID, FILE *stream, TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(1);
    FloatArray pv, coords(3);
    InternalStateValueType type = ISVT_UNDEFINED;

    if ( valID == DisplacementVector ) {
        type = ISVT_VECTOR;
    } else if ( valID == FluxVector ) {
        type = ISVT_SCALAR;
    } else {
        OOFEM_ERROR("POIExportModule::exportPrimVarAs: unsupported UnknownType");
    }

    // print header
    if ( type == ISVT_SCALAR ) {
        fprintf(stream, "SCALARS prim_scalar_%d\n", ( int ) valID);
    } else if ( type == ISVT_VECTOR ) {
        fprintf(stream, "VECTORS vector_%d float\n", ( int ) valID);
    } else {
        fprintf(stderr, "POIExportModule::exportPrimVarAs: unsupported variable type\n");
    }


    SpatialLocalizer *sl = d->giveSpatialLocalizer();
    // loop over POIs
    std :: list< POI_dataType > :: iterator PoiIter;
    for ( PoiIter = POIList.begin(); PoiIter != POIList.end(); ++PoiIter ) {
        coords.at(1) = ( * PoiIter ).x;
        coords.at(2) = ( * PoiIter ).y;
        coords.at(3) = ( * PoiIter ).z;
        //region = (*PoiIter).region;

        ElementGeometry *source = sl->giveElementContainingPoint(coords, NULL);
        if ( source ) {
            // ask interface
            EIPrimaryUnknownMapperInterface *interface =
                static_cast< EIPrimaryUnknownMapperInterface * >( source->giveInterface(EIPrimaryUnknownMapperInterfaceType) );
            if ( interface ) {
                interface->EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(VM_Total, tStep, coords, pv);
            } else {
                pv.clear();
                OOFEM_WARNING2( "POIExportModule::exportPrimVarAs: element %d with no EIPrimaryUnknownMapperInterface support",
                               source->giveNumber() );
            }

            fprintf(stream, "%10d ", ( * PoiIter ).id);
            if ( pv.giveSize() ) {
                for ( int j = 1; j <= pv.giveSize(); j++ ) {
                    fprintf( stream, " %15e ", pv.at(j) );
                }
            }

            fprintf(stream, "\n");
        } else {
            OOFEM_ERROR4( "POIExportModule::exportPrimVarAs: no element containing POI(%e,%e,%e) found",
                         coords.at(1), coords.at(2), coords.at(3) );
        }
    }
}
} // end namespace oofem
