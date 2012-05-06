/* $Header: /home/cvs/bp/oofem/main/src/maindebug.C,v 1.2.4.1 2004/04/05 15:19:41 bp Exp $ */
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

//  MAIN for debugging
//  Solves finite element problems.
//
//  The DEBUG option MUST be used (check in file 'debug.def').
//  See also file 'main2.c'.

#ifdef __OOFEM_MAINDEBUG

#include "engngm.h"
#include "freestor.h"
#include "compiler.h"

#include "plaintextdatareader.h"
#include "util.h"
#include "oofemdef.h"
#include "usrdefsub.h"
#include "error.h"
#include "logger.h"
#ifdef __PETSC_MODULE
 #include "petsc.h"
#endif
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string.h>
 #include <new>
#endif

#include "domain.h"
#include "material.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "timestep.h"
#include "structuralmaterial.h"
#include "mathfem.h"

using namespace oofem;

/* Defaul oofem loggers */
Logger oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);


void debugProcedure(EngngModel *);


int
main(int argc, char *argv[])
{
#ifdef BORLAND
    set_new_handler(& freeStoreError);    // prevents memory overflow
#endif
    int contextFlag = 0;
    EngngModel *problem;


    PlainTextDataReader dr("oofem.debug.in");
    problem = :: InstanciateProblem(& dr, _processor, contextFlag);
    dr.finish();

    debugProcedure(problem);

    delete problem;

    return 1;
}



void
debugProcedure(EngngModel *problem)
//
// procedure for debuging purposes
//
{
    int i;
    FloatArray *coord1, strainIncr(2), strain(2), stress(2);
    Domain *mesh = problem->giveDomain(1);
    TimeStep t(1, problem, 1, 1.0, 1.0, 1);
    coord1             = new FloatArray(1);
    coord1->at(1)    = 0.;
    GaussPoint *gp     = new GaussPoint(mesh->giveElement(1), 1, coord1, 0.5, _2dInterface);

    /*
     * // pure tension
     * StructuralMaterial *mat = (StructuralMaterial*) (mesh -> giveMaterial (1)) ;
     * // loading
     * strainIncr.at(1) = 0.0055 ;
     * for (i = 0; i<= 60; i++) {
     * strain.add(strainIncr);
     * mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     * gp->updateYourself (&t);
     * //printf ("\n");
     * //mat->giveStatus(gp)->printOutputAt (stdout, &t);
     * printf ("\n%e %e %e %e",strain.at(1), strain.at(2), stress.at(1), stress.at(2));
     * }
     */

    /*
     * // tension with small constant shear
     * StructuralMaterial *mat = (StructuralMaterial*) (mesh -> giveMaterial (1)) ;
     * // loading
     * strainIncr.at(1) = 0.0055 ;
     * strain.at(2) = 0.0010;
     * for (i = 0; i<= 60; i++) {
     * strain.add(strainIncr);
     * mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     * gp->updateYourself (&t);
     * //printf ("\n");
     * //mat->giveStatus(gp)->printOutputAt (stdout, &t);
     * printf ("\n%e %e %e %e",strain.at(1), strain.at(2), stress.at(1), stress.at(2));
     * }
     */

    /*
     * // tension with constant constant shear (1+2 yield functs activated)
     * StructuralMaterial *mat = (StructuralMaterial*) (mesh -> giveMaterial (1)) ;
     * // loading
     * strainIncr.at(1) = 0.0055 ;
     * strainIncr.at(2) = 0.0002 ;
     * strain.at(2) = 0.002;
     * for (i = 0; i<= 60; i++) {
     * strain.add(strainIncr);
     * mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     * gp->updateYourself (&t);
     * //printf ("\n");
     * //mat->giveStatus(gp)->printOutputAt (stdout, &t);
     * printf ("\n%e %e %e %e",strain.at(1), strain.at(2), stress.at(1), stress.at(2));
     * }
     */


    /*
     * // pure shear
     * StructuralMaterial *mat = (StructuralMaterial*) (mesh -> giveMaterial (1)) ;
     * // loading
     * strainIncr.at(2) = 0.0002 ;
     * for (i = 0; i<= 400; i++) {
     * strain.add(strainIncr);
     * mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     * gp->updateYourself (&t);
     * printf ("\n%e %e %e %e",strain.at(1), strain.at(2), stress.at(1), stress.at(2));
     * }
     */

    /*
     * // pure compression
     * StructuralMaterial *mat = (StructuralMaterial*) (mesh -> giveMaterial (1)) ;
     * // loading
     * strainIncr.at(1) = -0.004 ;
     * for (i = 0; i<= 300; i++) {
     * strain.add(strainIncr);
     * mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     * gp->updateYourself (&t);
     * printf ("\n%e %e %e %e",strain.at(1), strain.at(2), stress.at(1), stress.at(2));
     * }
     */

    // pure compression
    StructuralMaterial *mat = ( StructuralMaterial * ) ( mesh->giveMaterial(1) );
    // loading
    strainIncr.at(1) = -3.787e-2;
    strainIncr.at(2) =  8.111e-3;
    for ( i = 0; i <= 1; i++ ) {
        strain.add(strainIncr);
        mat->giveRealStressVector(stress, ReducedForm, gp, strain, & t);
        gp->updateYourself(& t);
        printf( "\n%e %e %e %e", strain.at(1), strain.at(2), stress.at(1), stress.at(2) );
    }




    /*
     * // shear under prescribed confiment
     * StructuralMaterial *mat = (StructuralMaterial*) (mesh -> giveMaterial (1)) ;
     * FloatMatrix d, di;
     * mat->giveCharacteristicMatrix(d, ReducedForm, ElasticStiffness, gp, &t);
     * di.beInverseOf(d);
     * double res, conf = -0.35;
     * // loading
     * for (i = 0; i<= 400; i++) {
     * strainIncr.at(1) = 0.0;
     * strainIncr.at(2) = 0.002 ;
     * strain.add(strainIncr);
     * mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     * res = conf - stress.at(1);
     *
     * while (fabs(res) > 1.e-4) {
     *  //printf ("%e\n", res);
     *  strainIncr.at(1) = di.at(1,1)*res;
     *  strainIncr.at(2) = 0.0;
     *  strain.add(strainIncr);
     *  mat->giveRealStressVector(stress, ReducedForm, gp, strain, &t);
     *  res = conf - stress.at(1);
     * };
     *
     * //printf ("%e\n", res);
     * gp->updateYourself (&t);
     *
     * //
     * //printf ("\n");
     * //mat->giveStatus(gp)->printOutputAt (stdout, &t);
     * //
     * printf ("\n%e %e %e %e",strain.at(1), strain.at(2), stress.at(1), stress.at(2));
     *
     * }
     */
}


#endif
