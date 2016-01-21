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

#include "homexportmodule.h"
#include "timestep.h"
#include "element.h"
#include "gausspoint.h"
#include "engngm.h"
#include "material.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ExportModule(HOMExportModule)

HOMExportModule :: HOMExportModule(int n, EngngModel *e) : ExportModule(n, e) { }

HOMExportModule :: ~HOMExportModule() { }

IRResultType
HOMExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    IR_GIVE_FIELD(ir, this->ists, _IFT_HOMExportModule_ISTs);
    this->scale = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->scale, _IFT_HOMExportModule_scale);
    //this->matnum.clear();
    //IR_GIVE_OPTIONAL_FIELD(ir, this->matnum, _IFT_HOMExportModule_matnum);
    return ExportModule :: initializeFrom(ir);
}

void
HOMExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    bool volExported = false;
    fprintf(this->stream, "%.3e  ", tStep->giveIntrinsicTime()*this->timeScale);
    IntArray elements;
    //assemble list of eligible elements. Elements can be present more times in a list but averaging goes just once over each element.
    elements.resize(0);
    for ( int ireg = 1; ireg <= this->giveNumberOfRegions(); ireg++ ) {
        elements.followedBy(this->giveRegionSet(ireg)->giveElementList());
    }
//     elements.printYourself();

    for ( int ist: ists ) {
        FloatArray ipState, avgState;
        double VolTot = 0.;
        Domain *d = emodel->giveDomain(1);
        for ( auto &elem : d->giveElements() ) {
            if ( elements.contains(elem -> giveGlobalNumber()) ){
                for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
                    double dV = elem->computeVolumeAround(gp);
                    VolTot += dV;
                    elem->giveGlobalIPValue(ipState, gp, (InternalStateType)ist, tStep);
                    avgState.add(dV, ipState);
                }
            }
        }

        if ( !volExported ) {
            fprintf(this->stream, "%.3e    ", VolTot);
            volExported = true;
        }
        avgState.times( 1. / VolTot * this->scale );
        fprintf(this->stream, "%d ", avgState.giveSize());
        for ( auto s: avgState ) {
            fprintf(this->stream, "%e ", s);
        }
        fprintf(this->stream, "    ");
    }
    fprintf(this->stream, "\n" );
    fflush(this->stream);
}

void
HOMExportModule :: initialize()
{
    std :: string fileName = emodel->giveOutputBaseFileName() + ".hom";
    if ( ( this->stream = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR( "failed to open file %s", fileName.c_str() );
    }

    fprintf(this->stream, "#Time      Volume       ");
    for ( int var: this->ists ) {
        fprintf(this->stream, "%s    ", __InternalStateTypeToString( ( InternalStateType ) var) );
    }
    fprintf(this->stream, "\n" );
    fflush(this->stream);


    initializeElementSet();

    ExportModule :: initialize();
}

void
HOMExportModule :: terminate()
{
    fclose(this->stream);
}
} // end namespace oofem
