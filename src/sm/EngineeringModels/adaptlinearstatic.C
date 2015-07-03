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

#include "../sm/EngineeringModels/adaptlinearstatic.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"
#include "errorestimator.h"
#include "domain.h"
#include "classfactory.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_EngngModel(AdaptiveLinearStatic);

void
AdaptiveLinearStatic :: updateYourself(TimeStep *tStep)
{
    LinearStatic :: updateYourself(tStep);
    // perform error evaluation
    // evaluate error of the reached solution
    this->defaultErrEstimator->estimateError(temporaryEM, tStep);
    // this->defaultErrEstimator->estimateError (equilibratedEM, this->giveCurrentStep());
    RemeshingStrategy strategy = this->defaultErrEstimator->giveRemeshingCrit()->giveRemeshingStrategy(tStep);

    if ( strategy == NoRemeshing_RS ) {
        return;
    } else {
        // do remeshing
        std :: unique_ptr< MesherInterface >mesher( classFactory.createMesherInterface( meshPackage, this->giveDomain(1) ) );
        Domain *newDomain;

        MesherInterface :: returnCode result =
            mesher->createMesh(tStep, 1, this->giveDomain(1)->giveSerialNumber() + 1, & newDomain);

        if ( result == MesherInterface :: MI_OK ) { } else if ( result == MesherInterface :: MI_NEEDS_EXTERNAL_ACTION ) {
            // terminate step
            //this->terminate( tStep );
            //this->terminateAnalysis();
            //exit(1);
        } else {
            OOFEM_ERROR("createMesh failed");
        }
    }
}

void
AdaptiveLinearStatic :: terminate(TimeStep *tStep)
{
    LinearStatic :: terminate(tStep);
    //
    // print estimated error
    //
    fprintf(outputStream, "\nRelative error estimate: %5.2f%%\n", this->defaultErrEstimator->giveValue(relativeErrorEstimateEEV, tStep) * 100.0);
}




int
AdaptiveLinearStatic :: initializeAdaptive(int tStepNumber)
{
    /*
     * Due to linear character of the problem,
     * the whole analysis is restarted from beginning.
     * The solution steps represent the adaptive steps and for each adaptive step
     * new domain with corresponding domainSerNum is generated.
     */
    int result = 1;
    /*
     * this -> initStepIncrements();
     *
     * int sernum = tStepNumber + 1;
     * printf ("\nrestoring domain %d.%d\n", 1, sernum);
     * Domain* dNew = new Domain (1, sernum, this);
     * FILE* domainInputFile;
     * this->giveDomainFile (&domainInputFile, 1, sernum, contextMode_read);
     * if (!dNew -> instanciateYourself(domainInputFile)) OOFEM_ERROR("domain Instanciation failed");
     * fclose (domainInputFile);
     *
     * printf ("\ndeleting old domain\n");
     * delete domainList->at(1);
     * domainList->put(1, dNew);
     *
     * // init equation numbering
     * this->forceEquationNumbering();
     *
     * // set time step
     * this->giveCurrentStep()->setTime(tStepNumber+1);
     *
     * // init equation numbering
     * // this->forceEquationNumbering();
     * this->giveNumericalMethod(giveCurrentStep())->setDomain (dNew);
     * this->ee->setDomain (dNew);
     */
    return result;
}


contextIOResultType AdaptiveLinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = LinearStatic :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

IRResultType
AdaptiveLinearStatic :: initializeFrom(InputRecord *ir)
// input from inputString
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int meshPackageId = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, meshPackageId, _IFT_AdaptiveLinearStatic_meshpackage);

    if ( meshPackageId == 1 ) {
        meshPackage = MPT_TARGE2;
    } else if ( meshPackageId == 2 ) {
        meshPackage = MPT_FREEM;
    } else {
        meshPackage = MPT_T3D;
    }

    return LinearStatic :: initializeFrom(ir);
}


void
AdaptiveLinearStatic :: updateDomainLinks()
{
    LinearStatic :: updateDomainLinks();
    // associate ee to possibly newly restored mesh
    this->defaultErrEstimator->setDomain( this->giveDomain(1) );
}
} // end namespace oofem
