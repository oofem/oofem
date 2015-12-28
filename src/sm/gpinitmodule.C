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

// Initialization module reading data related to Gauss points from a specified file

#include "gpinitmodule.h"
#include "element.h"
#include "integrationrule.h"
#include "material.h"
#include "floatarray.h"
#include "domain.h"
#include "engngm.h"
#include "classfactory.h"
#include "gausspoint.h"
#include "crosssection.h"
#include <cassert>

namespace oofem {
REGISTER_InitModule(GPInitModule)

GPInitModule :: GPInitModule(int n, EngngModel *e) : InitModule(n, e)
{ }


GPInitModule :: ~GPInitModule()
{ }


IRResultType
GPInitModule :: initializeFrom(InputRecord *ir)
{
    return InitModule :: initializeFrom(ir);
}


void
GPInitModule :: doInit()
{
    int ielem, nelem, ie, ig, nv, iv, nc, ic, varsize, vt;
    InternalStateType vartype;
    FloatArray value;

    double coords [ 3 ];

    Domain *d = emodel->giveDomain(1);
    nelem = d->giveNumberOfElements();

    // loop over elements
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        Element *elem = d->giveElement(ielem);
        // loop over Gauss points
        for ( auto &gp: *elem->giveDefaultIntegrationRulePtr() ) {
            Material *mat = elem->giveCrossSection()->giveMaterial(gp);
            MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
            if ( fscanf(initStream, "%d %d", & ie, & ig) != 2 ) {
                OOFEM_ERROR("initStream reading error");
            }

            // check whether the element and GP number agree
            assert(ielem == ie);
            assert( gp->giveNumber() == ig );
            // read coordinates
            if ( fscanf(initStream, "%d", & nc) != 1 ) {
                OOFEM_ERROR("initStream reading error");
            }

            assert(nc >= 0 && nc <= 3);
            for ( ic = 0; ic < nc; ic++ ) {
                if ( fscanf(initStream, "%lg", & coords [ ic ]) != 1 ) {
                    OOFEM_ERROR("initStream reading error");
                }
            }

            if ( fscanf(initStream, "%d", & nv) != 1 ) {
                OOFEM_ERROR("initStream reading error");
            }

            assert(nv >= 0);
            for ( iv = 1; iv <= nv; iv++ ) {
                if ( fscanf(initStream, "%d %d", & vt, & varsize) != 2 ) {
                    OOFEM_ERROR("initStream reading error");
                }

                vartype = ( InternalStateType ) vt;
                value.resize(varsize);
                for ( ic = 1; ic <= varsize; ic++ ) {
                    if ( fscanf( initStream, "%lg", & value.at(ic) ) != 1 ) {
                        OOFEM_ERROR("initStream reading error");
                    }
                }

                mat->setIPValue(value, gp, vartype);
            }
            // restore consistency (compute dependent internal variables)
            status->restoreConsistency();

        }
    }

    fclose(initStream);
}
} // namespace oofem
