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
#include "structuralms.h"
#include "gausspnt.h"
#include "domain.h"
#include "flotmtrx.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#endif

// constructor
//strainVector, tempStrainVector, stressVector, tempStressVector are defined on StructuralMaterialStatus

MicroMaterialStatus :: MicroMaterialStatus(int n, Domain *d, GaussPoint *gp) : StructuralMaterialStatus(n, d, gp){

}

MicroMaterialStatus :: ~MicroMaterialStatus(){

}

void MicroMaterialStatus :: initTempStatus()
{
  StructuralMaterialStatus :: initTempStatus();
}

void MicroMaterialStatus :: updateYourself(TimeStep *atTime)
{
  StructuralMaterialStatus :: updateYourself(atTime);
}

void MicroMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{

}

contextIOResultType MicroMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{

}

contextIOResultType MicroMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{

}


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




//original pure virtual function has to be redeclared here
void MicroMaterial :: giveRealStressVector (FloatArray& answer, MatResponseForm form, GaussPoint* gp, const FloatArray &totalStrain, TimeStep *atTime){
//perform average over microproblem
int i, j, index, ielem;
Element *elem;
double dV, VolTot = 0.;
double scale = 1.;
FloatArray VecStrain, VecStress, SumStrain(6), SumStress(6);
IntArray Mask;
GaussPoint *gpL;
IntegrationRule *iRule;
Domain *microDomain = problemMicro->giveDomain(1);//from engngm.h

int nelem = microDomain->giveNumberOfElements();
int nnodes = microDomain->giveNumberOfDofManagers();

for ( ielem = 1; ielem <= nelem; ielem++ ) {
        elem = microDomain->giveElement(ielem);
        iRule = elem->giveDefaultIntegrationRulePtr();
        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gpL  = iRule->getIntegrationPoint(i);
            gpL->giveCoordinate(1);
            dV  = elem->computeVolumeAround(gpL);
            VolTot += dV;
            //OOFEM_LOG_INFO("Element %d GP %d Vol %f\n", elem->giveNumber(), gp->giveNumber(), dV);
            //fprintf(this->stream, "Element %d GP %d stress %f\n", elem->giveNumber(), gp->giveNumber(), 0.0);
            //((StructuralCrossSection*) gp->giveCrossSection())->giveFullCharacteristicVector(helpVec, gp, strainVector);
            elem->giveIPValue(VecStrain, gpL, IST_StrainTensor, atTime);
            elem->giveIPValue(VecStress, gpL, IST_StressTensor, atTime);
            elem->giveIntVarCompFullIndx(Mask, IST_StrainTensor);

            VecStrain.times(dV);
            VecStress.times(dV);

            for ( j = 0; j < 6; j++ ) {
                index = Mask(j); //indexes in Mask from 1
                if ( index ) {
                    SumStrain(j) += VecStrain(index - 1);
                    SumStress(j) += VecStress(index - 1);
                }
            }

            //VecStrain.printYourself();
            //SumStrain.printYourself();
        }
    }

    //average
    SumStrain.times(1. / VolTot * scale);
    SumStress.times(1. / VolTot * scale);
    //SumStrain.printYourself();
    //SumStress.printYourself();
    answer.resize(6);
    answer = SumStress;
}

MaterialStatus *
MicroMaterial :: CreateStatus(GaussPoint *gp) const
{
    MicroMaterialStatus *status =
        new  MicroMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}