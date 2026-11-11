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


#include "mesherinterface.h"
#include "domain.h"
#include "engngm.h"
#include "timestep.h"
#include "oofemtxtdatareader.h"
#include "error.h"
#include "logger.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>

namespace oofem {
void
MesherInterface :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    IR_GIVE_OPTIONAL_FIELD(ir, remeshCmd, _IFT_MesherInterface_remeshcmd);
}


MesherInterface :: returnCode
MesherInterface :: remeshExternally(const std :: string &densityFile, TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew)
{
    * dNew = nullptr;
    if ( remeshCmd.empty() ) {
        return MI_NEEDS_EXTERNAL_ACTION;
    }

    EngngModel *emodel = this->domain->giveEngngModel();
    std :: string domainFile = emodel->giveDomainFileName(1, domainSerNum);

    // substitute placeholders
    std :: string cmd;
    for ( size_t i = 0; i < remeshCmd.size(); i++ ) {
        if ( remeshCmd [ i ] == '%' && i + 1 < remeshCmd.size() ) {
            char c = remeshCmd [ ++i ];
            if ( c == 'd' ) {
                cmd += densityFile;
            } else if ( c == 'o' ) {
                cmd += domainFile;
            } else if ( c == 's' ) {
                cmd += std :: to_string(domainSerNum);
            } else if ( c == 'n' ) {
                cmd += std :: to_string( tStep->giveNumber() );
            } else if ( c == '%' ) {
                cmd += '%';
            } else {
                cmd += '%';
                cmd += c;
            }
        } else {
            cmd += remeshCmd [ i ];
        }
    }

    // remove old domain file, so that the stale file is never used
    std :: remove( domainFile.c_str() );

    OOFEM_LOG_INFO("Running remesh command: %s\n", cmd.c_str() );
    // flush buffered output, so that it is not interleaved with output of the command
    fflush(nullptr);
    int ret = std :: system( cmd.c_str() );
    if ( ret != 0 ) {
        OOFEM_WARNING("Remesh command failed (return code %d): %s", ret, cmd.c_str() );
        return MI_FAILED;
    }

    if ( !std :: ifstream(domainFile).good() ) {
        OOFEM_WARNING("Remesh command did not create domain file %s", domainFile.c_str() );
        return MI_FAILED;
    }

    OOFEM_LOG_INFO("Reading new domain from %s\n", domainFile.c_str() );
    auto d = std :: make_unique< Domain >( domainNumber, domainSerNum, emodel );
    try {
        OOFEMTXTDataReader dr(domainFile, true);
        if ( !d->instanciateYourself(dr) ) {
            OOFEM_WARNING("Instanciation of new domain from %s failed", domainFile.c_str() );
            return MI_FAILED;
        }
        d->postInitialize();
    } catch ( const std :: exception &e ) {
        OOFEM_WARNING("Instanciation of new domain from %s failed: %s", domainFile.c_str(), e.what() );
        return MI_FAILED;
    }
    * dNew = d.release();
    return MI_OK;
}
} // end namespace oofem
