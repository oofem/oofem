/* $Header: /home/cvs/bp/oofem/sm/src/HOMExportModule.C,v 1.11.4.1 2004/04/05 15:19:47 bp Exp $ */
/*
 *
 *****    *****   ******  ******  ***   ***
 **   **  **   **  **      **      ** *** **
 **   **  **   **  ****    ****    **  *  **
 **   **  **   **  **      **      **     **
 **   **  **   **  **      **      **     **
 *****    *****   **      ******  **     **
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 2009 - 2009   Borek Patzak, Vit Smilauer
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

#include "homexportmodule.h"
#include "timestep.h"
#include "engngm.h"
#include "strreader.h"
#include "node.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "materialinterface.h"
#include "mathfem.h"
#include "oofem_limits.h"
#include "gausspnt.h"
#ifndef __MAKEDEPEND
 #include <vector>
#endif

namespace oofem {

//inherit LinearElasticMaterial for accessing stress/strain transformation functions
HOMExportModule :: HOMExportModule(EngngModel *e) : ExportModule(e), LinearElasticMaterial(0, NULL)
{
  this->matnum.resize(0);
}


HOMExportModule :: ~HOMExportModule()
{ }

IRResultType
HOMExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    IRResultType val;
    this->scale = 1.;
    ExportModule :: initializeFrom(ir);
    val = IR_GIVE_OPTIONAL_FIELD(ir, this->scale, IFT_HOMExportModule_scale, "scale"); // Macro
    if ( val == IRRT_NOTFOUND ) {
        this->scale = 1.;
    }

    val = IR_GIVE_OPTIONAL_FIELD(ir, this->matnum, IFT_HOMExportModule_matnum, "matnum"); // Macro
    return IRRT_OK;
}


void
HOMExportModule :: doOutput(TimeStep *tStep)
{
    int i, j, index, ielem;
    Element *elem;
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    Domain *d  = emodel->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    StructuralElement *structElem;
    //int nnodes = d->giveNumberOfDofManagers();
    double dV, VolTot = 0.;
    FloatArray VecStrain, VecStress, VecEigStrain, SumStrain(6), SumStress(6), SumEigStrain(6), tempFloatAr;
    FloatMatrix baseGCS;
    IntArray Mask;
    GaussPoint *gp;
    IntegrationRule *iRule;

    //stress and strain vectors are always in global c.s.
    SumStrain.zero(); //xx, yy, zz, yz, zx, xy
    SumStress.zero();
    SumEigStrain.zero();

    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        elem = d->giveElement(ielem);
        if( this->matnum.giveSize() == 0 || this->matnum.contains(elem->giveMaterial()->giveNumber()) ){
            iRule = elem->giveDefaultIntegrationRulePtr();
            for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
                gp  = iRule->getIntegrationPoint(i);
                structElem = ( StructuralElement * ) gp->giveElement();
                structElem->computeResultingIPEigenstrainAt(VecEigStrain, tStep, gp, VM_Incremental);
                //gp->giveCoordinate(1);
                dV  = elem->computeVolumeAround(gp);
                //OOFEM_LOG_INFO("Element %d GP %d Vol %f\n", elem->giveNumber(), gp->giveNumber(), dV);
                //fprintf(this->stream, "Element %d GP %d stress %f\n", elem->giveNumber(), gp->giveNumber(), 0.0);
                //((StructuralCrossSection*) gp->giveCrossSection())->giveFullCharacteristicVector(helpVec, gp, strainVector);
                elem->giveIPValue(VecStrain, gp, IST_StrainTensor, tStep);
                elem->giveIPValue(VecStress, gp, IST_StressTensor, tStep);
                elem->giveIntVarCompFullIndx(Mask, IST_StrainTensor);

                //truss element has strains and stresses in the first array so transform them to global coordinates
                if(elem->giveClassID() == Truss3dClass){
                  elem->giveLocalCoordinateSystem(baseGCS);
                  //this->giveStressVectorTranformationMtrx(transfMatrix,baseGCS,1);//from structuralMaterial
                  //baseGCS.printYourself();
                  VecStress.resize(6);
                  VecStrain.resize(6);
                  VecEigStrain.resize(6);
                  this->transformStressVectorTo( tempFloatAr, baseGCS, VecStress, 0);
                  VecStress = tempFloatAr;
                  this->transformStrainVectorTo( tempFloatAr, baseGCS, VecStrain, 0);
                  VecStrain = tempFloatAr;
                  tempFloatAr = VecEigStrain;
                  VecEigStrain = tempFloatAr;

                  for(j=1; j<=6; j++)
                        Mask.at(j)=j;
                  //VecStrain.printYourself();
                  //VecEigStrain.printYourself();
                  //VecStress.printYourself();
                }

                VolTot += dV;
                VecStrain.times(dV);
                VecStress.times(dV);
                VecEigStrain.times(dV);

                for ( j = 0; j < 6; j++ ) {
                    index = Mask(j); //indexes in Mask from 1
                    if ( index ) {
                        SumStrain(j) += VecStrain(index - 1);
                        if(VecEigStrain.giveSize())
                          SumEigStrain(j) += VecEigStrain(index - 1);
                        SumStress(j) += VecStress(index - 1);
                    }
                }

                //VecStrain.printYourself();
                //VecEigStrain.printYourself();
                //SumStrain.printYourself();
            }
        }
    }

    //average
    SumStrain.times(1. / VolTot * this->scale);
    SumEigStrain.times(1. / VolTot * this->scale);
    SumStress.times(1. / VolTot * this->scale);
    //SumStrain.printYourself();
    //SumEigStrain.printYourself();
    //SumStress.printYourself();
    fprintf( this->stream, "%f   ", tStep->giveTime() );
    for ( j = 0; j < 6; j++ ) { //strain
        fprintf( this->stream, "%06.5e ", SumStrain(j) );
    }

    fprintf(this->stream, "    ");
    for ( j = 0; j < 6; j++ ) { //stress
        fprintf( this->stream, "%06.5e ", SumStress(j) );
    }

    fprintf(this->stream, "    ");
    for ( j = 0; j < 6; j++ ) { //eigenstrain
        fprintf( this->stream, "%06.5e ", SumEigStrain(j) );
    }

    fprintf( this->stream, "Vol %06.5e ", VolTot );
    fprintf(this->stream, "\n");
    fflush(this->stream);
}

void
HOMExportModule :: initialize()
{
    char baseFileName [ MAX_FILENAME_LENGTH ];
    char fileName [ MAX_FILENAME_LENGTH ];

    emodel->giveOutputBaseFileName(baseFileName, MAX_FILENAME_LENGTH);
    sprintf(fileName, "%s.hom", baseFileName);
    if ( ( this->stream = fopen(fileName, "w") ) == NULL ) {
        OOFEM_ERROR2("HOMExportModule::giveOutputStream: failed to open file %s", fileName);
    }

    fprintf(this->stream, "#Time      Strain (xx, yy, zz, yz, zx, xy)                                              Stress (xx, yy, zz, yz, zx, xy)                                                Eigenstrain (xx, yy, zz, yz, zx, xy)\n");
}


void
HOMExportModule :: terminate()
{
    fclose(this->stream);
}

} // end namespace oofem
