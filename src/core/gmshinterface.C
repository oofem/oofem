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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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


#include "gmshinterface.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "node.h"
#include "errorestimator.h"
#include "remeshingcrit.h"
#include "timestep.h"
#include "classfactory.h"
#include "error.h"
#include "logger.h"

#include <cstdio>
#include <queue>
#include <vector>

#ifndef OOFEM_GMSH2OOFEM_SCRIPT
 #define OOFEM_GMSH2OOFEM_SCRIPT "gmsh2oofem.py"
#endif

namespace oofem {
REGISTER_Mesher(GmshInterface, MPT_GMSH);

GmshInterface :: GmshInterface(Domain *d) : MesherInterface(d),
    gmshCmd("gmsh"), pythonCmd("python3"), script(OOFEM_GMSH2OOFEM_SCRIPT), grading(0.)
{ }


void
GmshInterface :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    MesherInterface :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, geoFile, _IFT_GmshInterface_geo);
    IR_GIVE_OPTIONAL_FIELD(ir, ctrlFile, _IFT_GmshInterface_ctrl);
    IR_GIVE_OPTIONAL_FIELD(ir, gmshCmd, _IFT_GmshInterface_gmshcmd);
    IR_GIVE_OPTIONAL_FIELD(ir, pythonCmd, _IFT_GmshInterface_pythoncmd);
    IR_GIVE_OPTIONAL_FIELD(ir, script, _IFT_GmshInterface_script);
    IR_GIVE_OPTIONAL_FIELD(ir, scriptOpts, _IFT_GmshInterface_opts);
    IR_GIVE_OPTIONAL_FIELD(ir, grading, _IFT_GmshInterface_grading);
}


MesherInterface :: returnCode
GmshInterface :: createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew)
{
    * dNew = nullptr;
    std :: string bgmFile = this->domain->giveEngngModel()->giveOutputBaseFileName() + ".bgm." + std :: to_string(domainSerNum) + ".pos";
    if ( !this->createBackgroundMesh(bgmFile, tStep) ) {
        return MI_FAILED;
    }

    if ( remeshCmd.empty() && !geoFile.empty() && !ctrlFile.empty() ) {
        // default remesh command: gmsh + unv2oofem conversion
        int dim = this->domain->giveNumberOfSpatialDimensions();
        remeshCmd = "\"" + pythonCmd + "\" \"" + script + "\" --gmsh \"" + gmshCmd + "\" --dim " + std :: to_string(dim) +
                    " " + scriptOpts + " --bgm \"%d\" --din \"" + geoFile + "\" \"" + ctrlFile + "\" \"%o\"";
    }

    return this->remeshExternally(bgmFile, tStep, domainNumber, domainSerNum, dNew);
}


/// Returns gmsh view cell type and number of corner nodes for given element (false if not supported).
static bool giveGmshCell(Element *elem, const char * &cell, int &ncorners)
{
    // corner nodes only (gmsh view cell types: ST triangle, SQ quad, SS tetra, SH hexa, SI prism)
    switch ( elem->giveGeometryType() ) {
    case EGT_triangle_1:
    case EGT_triangle_2:
        cell = "ST";
        ncorners = 3;
        return true;
    case EGT_quad_1:
    case EGT_quad_2:
    case EGT_quad9_2:
        cell = "SQ";
        ncorners = 4;
        return true;
    case EGT_tetra_1:
    case EGT_tetra_2:
        cell = "SS";
        ncorners = 4;
        return true;
    case EGT_hexa_1:
    case EGT_hexa_2:
    case EGT_hexa_27:
        cell = "SH";
        ncorners = 8;
        return true;
    case EGT_wedge_1:
    case EGT_wedge_2:
        cell = "SI";
        ncorners = 6;
        return true;
    default:
        return false; // not a volume/area element (e.g. line, interface)
    }
}


void
GmshInterface :: limitSizeGradation(FloatArray &h)
{
    // Limits the gradation of the size field: h(i) <= h(j) + grading * |x_i - x_j| for all node pairs
    // connected by elements of the current mesh. Solved exactly by Dijkstra-like propagation from
    // the smallest sizes (sizes can only decrease).
    int nnodes = this->domain->giveNumberOfDofManagers();
    std :: vector< std :: vector< std :: pair< int, double > > >adj(nnodes + 1);
    const char *cell;
    int ncorners;
    for ( auto &elem : this->domain->giveElements() ) {
        if ( !giveGmshCell(elem.get(), cell, ncorners) ) {
            continue;
        }
        for ( int i = 1; i <= ncorners; i++ ) {
            for ( int j = i + 1; j <= ncorners; j++ ) {
                int a = elem->giveNode(i)->giveNumber(), b = elem->giveNode(j)->giveNumber();
                double l = distance( elem->giveNode(i)->giveCoordinates(), elem->giveNode(j)->giveCoordinates() );
                adj [ a ].emplace_back(b, l);
                adj [ b ].emplace_back(a, l);
            }
        }
    }

    typedef std :: pair< double, int >Entry;
    std :: priority_queue< Entry, std :: vector< Entry >, std :: greater< Entry > >queue;
    for ( int i = 1; i <= nnodes; i++ ) {
        if ( h.at(i) > 0. ) {
            queue.emplace(h.at(i), i);
        }
    }
    int nchanged = 0;
    while ( !queue.empty() ) {
        Entry e = queue.top();
        queue.pop();
        if ( e.first > h.at(e.second) ) {
            continue; // outdated entry
        }
        for ( auto &nb : adj [ e.second ] ) {
            double limit = e.first + this->grading * nb.second;
            if ( h.at(nb.first) > limit ) {
                h.at(nb.first) = limit;
                queue.emplace(limit, nb.first);
                nchanged++;
            }
        }
    }
    OOFEM_LOG_INFO("GmshInterface: size gradation limited to %g, %d nodal sizes reduced\n", this->grading, nchanged);
}


int
GmshInterface :: createBackgroundMesh(const std :: string &fileName, TimeStep *tStep)
{
    RemeshingCriteria *rc = this->domain->giveErrorEstimator()->giveRemeshingCrit();
    int nnodes = this->domain->giveNumberOfDofManagers();

    // required nodal sizes
    FloatArray h(nnodes);
    for ( int i = 1; i <= nnodes; i++ ) {
        h.at(i) = rc->giveRequiredDofManDensity(i, tStep);
    }
    if ( this->grading > 0. ) {
        this->limitSizeGradation(h);
    }

    FILE *file = fopen(fileName.c_str(), "w");
    if ( !file ) {
        OOFEM_WARNING("Can not open background mesh file %s", fileName.c_str() );
        return 0;
    }

    fprintf(file, "View \"oofem required element size\" {\n");
    int nwritten = 0;
    const char *cell;
    int ncorners;
    for ( auto &elem : this->domain->giveElements() ) {
        if ( !giveGmshCell(elem.get(), cell, ncorners) ) {
            continue;
        }

        fprintf(file, "%s(", cell);
        for ( int i = 1; i <= ncorners; i++ ) {
            const auto &c = elem->giveNode(i)->giveCoordinates();
            for ( int j = 0; j < 3; j++ ) {
                fprintf(file, "%s%.12e", ( i == 1 && j == 0 ) ? "" : ",", c [ j ]);
            }
        }
        fprintf(file, "){");
        for ( int i = 1; i <= ncorners; i++ ) {
            fprintf(file, "%s%.12e", i == 1 ? "" : ",", h.at( elem->giveNode(i)->giveNumber() ) );
        }
        fprintf(file, "};\n");
        nwritten++;
    }
    fprintf(file, "};\n");
    fclose(file);

    if ( nwritten == 0 ) {
        OOFEM_WARNING("No supported element written into background mesh %s", fileName.c_str() );
        return 0;
    }
    OOFEM_LOG_INFO("Gmsh background mesh file %s created\n", fileName.c_str() );
    return 1;
}
} // end namespace oofem
