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

#include "adaptlinearstatic.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"
#include "errorestimator.h"
#include "usrdefsub.h"
#include "contextioerr.h"

namespace oofem {
void
AdaptiveLinearStatic::updateYourself(TimeStep *stepN)
{

   LinearStatic :: updateYourself(stepN);
   // perform error evaluation
   // evaluate error of the reached solution
   this->defaultErrEstimator->estimateError( temporaryEM, stepN );
   // this->defaultErrEstimator->estimateError (equilibratedEM, this->giveCurrentStep());
   RemeshingStrategy strategy = this->defaultErrEstimator->giveRemeshingCrit()->giveRemeshingStrategy( stepN );

    if ( strategy == NoRemeshing_RS ) {
        return;
    } else {
        // do remeshing
        MesherInterface *mesher = CreateUsrDefMesherInterface( meshPackage, this->giveDomain(1) );
        Domain *newDomain;

        MesherInterface :: returnCode result =
        mesher->createMesh(stepN, 1, this->giveDomain(1)->giveSerialNumber() + 1, & newDomain);

        if ( result == MesherInterface :: MI_OK ) {} else if ( result == MesherInterface :: MI_NEEDS_EXTERNAL_ACTION ) {
            // terminate step
            this->terminate( stepN );
            this->terminateAnalysis();
            exit(1);
        } else {
            _error("solveYourselfAt: MesherInterface::createMesh failed");
        }
    }
}


int
AdaptiveLinearStatic :: initializeAdaptive(int stepNumber)
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
     * int sernum = stepNumber + 1;
     * printf ("\nrestoring domain %d.%d\n", 1, sernum);
     * Domain* dNew = new Domain (1, sernum, this);
     * FILE* domainInputFile;
     * this->giveDomainFile (&domainInputFile, 1, sernum, contextMode_read);
     * if (!dNew -> instanciateYourself(domainInputFile)) _error ("initializeAdaptive: domain Instanciation failed");
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
     * this->giveCurrentStep()->setTime(stepNumber+1);
     *
     * // init equation numbering
     * // this->forceEquationNumbering();
     * this->giveNumericalMethod(giveCurrentStep())->setDomain (dNew);
     * this->ee->setDomain (dNew);
     */
    return result;
}


contextIOResultType AdaptiveLinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj) {
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LinearStatic :: initializeFrom(ir);

    int meshPackageId = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, meshPackageId, IFT_AdaptiveLinearStatic_meshpackage, "meshpackage"); // Macro

    if ( meshPackageId == 1 ) {
        meshPackage = MPT_TARGE2;
    } else if ( meshPackageId == 2 ) {
        meshPackage = MPT_FREEM;
    } else {
        meshPackage = MPT_T3D;
    }

    return IRRT_OK;
}


void
AdaptiveLinearStatic :: updateDomainLinks()
{
    LinearStatic :: updateDomainLinks();
    // associate ee to possibly newly restored mesh
    this->defaultErrEstimator->setDomain( this->giveDomain(1) );
}
} // end namespace oofem
