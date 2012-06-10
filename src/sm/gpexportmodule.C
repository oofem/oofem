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

#include "gpexportmodule.h"
#include "gausspnt.h"
#include "material.h"
#include "element.h"
#include "integrationrule.h"
#include "timestep.h"
#include "engngm.h"

namespace oofem {
GPExportModule :: GPExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
    ncoords = -1; // means: export as many coordinates as available
}


GPExportModule :: ~GPExportModule()
{}


IRResultType
GPExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    ExportModule :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, vartypes, IFT_GPExportModule_vartypes, "vars");
    IR_GIVE_OPTIONAL_FIELD(ir, ncoords, IFT_GPExportModule_ncoords, "ncoords");
    return IRRT_OK;
}


void
GPExportModule :: doOutput(TimeStep *tStep)
{
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    int ielem, j, ic, nc, iv, nv, nvars;
    double weight;
    Element *elem;
    GaussPoint *gp;
    FloatArray gcoords, intvar;
    InternalStateType vartype;
    IntegrationRule *iRule;

    Domain *d  = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    FILE *stream = this->giveOutputStream(tStep);

    // print the header
    fprintf(stream, "# gauss point data file\n");
    fprintf( stream, "# output for time %g\n", tStep->giveTargetTime() );
    fprintf(stream, "# variables: ");
    nvars = vartypes.giveSize();
    fprintf(stream, "%d  ", nvars);
    for ( iv = 1; iv <= nvars; iv++ ) {
        fprintf( stream, "%d ", vartypes.at(iv) );
    }

    fprintf(stream, "\n# for interpretation see internalstatetype.h\n");

    // loop over elements
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        elem = d->giveElement(ielem);
        iRule = elem->giveDefaultIntegrationRulePtr();

        // loop over Gauss points
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            gp = iRule->getIntegrationPoint(j);
            // export:
            // 1) element number
            // 2) material number
            // 3) Gauss point number
            // 4) contributing volume around Gauss point
            weight = gp->giveElement()->computeVolumeAround(gp);
            fprintf(stream, "%d %d %d %.6e ", elem->giveNumber(), elem->giveMaterial()->giveNumber(), j + 1, weight);

            // export Gauss point coordinates
            if ( ncoords ) { // no coordinates exported if ncoords==0
                elem->computeGlobalCoordinates( gcoords, * ( gp->giveCoordinates() ) );
                nc = gcoords.giveSize();
                if ( ncoords >= 0 ) {
                    fprintf(stream, "%d ", ncoords);
                } else {
                    fprintf(stream, "%d ", nc);
                }

                if ( ncoords > 0 && ncoords < nc ) {
                    nc = ncoords;
                }

                for ( ic = 1; ic <= nc; ic++ ) {
                    fprintf( stream, "%.6e ", gcoords.at(ic) );
                }

                for ( ic = nc + 1; ic <= ncoords; ic++ ) {
                    fprintf(stream, "%g ", 0.0);
                }
            }

            // export internal variables
            for ( iv = 1; iv <= nvars; iv++ ) {
                vartype = ( InternalStateType ) vartypes.at(iv);
                elem->giveIPValue(intvar, gp, vartype, tStep);
                nv = intvar.giveSize();
                fprintf(stream, "%d ", nv);
                for ( ic = 1; ic <= nv; ic++ ) {
                    fprintf( stream, "%.6e ", intvar.at(ic) );
                }
            }

            fprintf(stream, "\n");
        }

#if 0
        // for CST elements write also nodal coordinates
        // (non-standard part, used only exceptionally)
        int nnode = elem->giveNumberOfNodes();
        if (nnode==3){
            for (int inod=1; inod<=3; inod++)
                fprintf (stream, "%f %f ",elem->giveNode(inod)->giveCoordinate(1),elem->giveNode(inod)->giveCoordinate(2));
        }
#endif
    }

    fclose(stream);
}

void
GPExportModule :: initialize()
{}


void
GPExportModule :: terminate()
{}


FILE *
GPExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    std::string fileName = this->giveOutputBaseFileName(tStep) + ".gp";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2("GPExportModule::giveOutputStream: failed to open file %s", fileName.c_str());
    }

    return answer;
}
} // namespace oofem
