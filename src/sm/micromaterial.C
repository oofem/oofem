/* $Header: /home/cvs/bp/oofem/oofemlib/src/isolinearelasticmaterial.C,v 1.8 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "micromaterial.h"
#include "material.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "structuralmaterial.h"
#include "domain.h"
#include "flotmtrx.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#endif

MicroMaterial :: MicroMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{
 /// Constructor
}


IRResultType MicroMaterial :: initializeFrom(InputRecord *ir)
{
    int i;
    double value;
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

  
  IR_GIVE_FIELD2(ir, inputFileNameMicro, IFT_MicroMaterial_tmp, "file", MAX_FILENAME_LENGTH);
  
  OOFEM_LOG_INFO( "** Instanciating microproblem from file %s\n", inputFileNameMicro);
  OOFEMTXTDataReader drMicro(inputFileNameMicro);
  problemMicro = :: InstanciateProblem(& drMicro, _processor, 0);
  drMicro.finish();
  problemMicro->checkProblemConsistency();
  OOFEM_LOG_INFO( "** Microproblem at address %p instanciated\n", problemMicro);
}




//pure virtual function has to be declared here
void MicroMaterial :: giveRealStressVector (FloatArray& answer,  MatResponseForm, GaussPoint*, const FloatArray&, TimeStep*){
//
}
