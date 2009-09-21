/* $Header: /home/cvs/bp/oofem/oofemlib/src/micromaterial.C,v 1.9 2009/09/20 14:08:24 vs Exp $ */
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
#include "structuralmaterial.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "domain.h"
#include "dofmanager.h"
#include "usrdefsub.h"
#include "oofemtxtdatareader.h"
#include "util.h"

#ifndef __MAKEDEPEND
 #include <stdlib.h>
#endif


// constructor
//strainVector, tempStrainVector, stressVector, tempStressVector are defined on StructuralMaterialStatus
MicroMaterialStatus :: MicroMaterialStatus(int n, Domain *d, GaussPoint *gp) : StructuralMaterialStatus(n, d, gp) { }

MicroMaterialStatus :: ~MicroMaterialStatus() { }

void MicroMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void MicroMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
}

void MicroMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{ }

contextIOResultType MicroMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}

contextIOResultType MicroMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}

/// Constructor
MicroMaterial :: MicroMaterial(int n, Domain *d) : StructuralMaterial(n, d), UnknownNumberingScheme()
{
    isDefaultNumbering = true;
    //full stiffness matrix of microproblem
    SparseMtrxType sparseMtrxType = ( SparseMtrxType ) 0; //SMT_Skyline symmetric skyline
    stiffnessMatrixMicro = :: CreateUsrDefSparseMtrx(sparseMtrxType);
}

///destructor
MicroMaterial :: ~MicroMaterial() {
    if ( stiffnessMatrixMicro ) {
        delete stiffnessMatrixMicro;
    }
}

IRResultType MicroMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro


    IR_GIVE_FIELD2(ir, inputFileNameMicro, IFT_MicroMaterial_tmp, "file", MAX_FILENAME_LENGTH);

    OOFEM_LOG_INFO("** Instanciating microproblem with BC from file %s\n", inputFileNameMicro);
    OOFEMTXTDataReader drMicro(inputFileNameMicro);
    problemMicro = :: InstanciateProblem(& drMicro, _processor, 0); //0=contextFlag-store/resore
    drMicro.finish();
    OOFEM_LOG_INFO("** Microproblem at address %p instanciated\n", problemMicro);

    return IRRT_OK;
}




//original pure virtual function has to be declared here
//this function should not be used, internal forces are calculated based on reactions not stresses in GPs
void MicroMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime) {
    //perform average over microproblem
    int j, index, ielem;
    Element *elem;
    double dV, VolTot = 0.;
    double scale = 1.;
    FloatArray VecStrain, VecStress, SumStrain(6), SumStress(6);
    IntArray Mask;
    GaussPoint *gpL;
    IntegrationRule *iRule;
    Domain *microDomain = problemMicro->giveDomain(1); //from engngm.h
    EngngModel *microEngngModel = microDomain->giveEngngModel();
    StructuralMaterialStatus *status = ( StructuralMaterialStatus * ) this->giveStatus(gp);

    OOFEM_ERROR("\n MicroMaterial :: giveRealStressVector should not be called at all\n");

    this->initGpForNewStep(gp);
    int nelem = microDomain->giveNumberOfElements();
    //int nnodes = microDomain->giveNumberOfDofManagers();

    //need to update stress so change boundary conditions and recalculate
    OOFEM_LOG_INFO( "\n*** giveRealStress microproblem %p on macroElement %d GP %d, step %d, time %f\n", this, this->macroLSpaceElement->giveNumber(), gp->giveNumber(), microEngngModel->giveCurrentStep()->giveNumber(), microEngngModel->giveCurrentStep()->giveTime() );




    this->macroLSpaceElement->changeMicroBoundaryConditions(atTime);
    microEngngModel->solveYourselfAt( microEngngModel->giveCurrentStep() );
    microEngngModel->terminate( microEngngModel->giveCurrentStep() );
    OOFEM_LOG_INFO("\n*** giveRealStress microproblem %p done\n", this);

    for ( ielem = 1; ielem <= nelem; ielem++ ) { //return stress as average through all elements of the same MicroMaterial
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
    //answer.printYourself();

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}

MaterialStatus *
MicroMaterial :: CreateStatus(GaussPoint *gp) const
{
    MicroMaterialStatus *status =
        new  MicroMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}

//from UnknownNumberingScheme
void MicroMaterial :: init(void) {
    DofManager *DofMan;
    int i, j, counter = 0;
    IntArray ut, dofIDArry(3);
    for ( i = 1; i <= 3; i++ ) {
        dofIDArry.at(i) = i;
    }

    if ( problemMicro->giveNumberOfDomains() != 1 ) {
        OOFEM_ERROR("Number of domains on microproblem is greater than 1");
    }

    //the pointer to underlying problem is problemMicro
    for ( i = 1; i <= problemMicro->giveDomain(1)->giveNumberOfDofManagers(); i++ ) { //for each node
        DofMan = problemMicro->giveDomain(1)->giveDofManager(i);
        //printf("%d\n",DofMan->giveNumberOfPrimaryMasterDofs(dofIDArry));
        //DofMan->giveLocationArray(dofIDArry, ut, *this );
        DofMan->giveCompleteLocationArray(ut, * this);
        for ( j = 1; j <= DofMan->giveNumberOfDofs(); j++ ) {
            counter++;
        }
    }

    totalNumberOfDomainEquation = counter;
}

//answer must be of size 24x24 (linear brick 3*8=24)
void MicroMaterial :: giveCondensedStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep, CharType type, IntArray microNodes, IntArray microDOFs) {
    int i, j;
    Domain *microDomain = problemMicro->giveDomain(1);
    EngngModel *microEngngModel = microDomain->giveEngngModel();
    stiffnessMatrixMicro->zero();
    FloatMatrix stiffnessMatrixMicroFloat;

    answer.resize(24, 24);
    answer.zero();

    this->isDefaultNumbering = false; //total number of equations corresponds to all DOFs
    //stiffnessMatrixMicro->buildInternalStructure(microEngngModel, 1, EID_MomentumBalance, EModelDefaultEquationNumbering());

    stiffnessMatrixMicro->buildInternalStructure(microEngngModel, 1, EID_MomentumBalance, * this);
    stiffnessMatrixMicro->zero();
    //stiffnessMatrixMicro->printYourself();
    microEngngModel->assemble(stiffnessMatrixMicro, tStep, EID_MomentumBalance, type, * this, microDomain);
    this->isDefaultNumbering = true; //switch back
    stiffnessMatrixMicro->toFloatMatrix(stiffnessMatrixMicroFloat); //must be converted for condensation
    //now we have upper part of stiffness matrix. Need to symmetrize. It's ommited on macrolspace element;
    stiffnessMatrixMicroFloat.symmetrized();
    //stiffnessMatrixMicroFloat.printYourself();

    //perform condensation of the full stiffness matrix, MUST be FloatMatrix type
    //algorithm based on structuralelement.C, Rayleigh-Ritz method
    //zeroes on particular condensed places will appear in FloatMatrix
    IntArray what; //DOF array to be condensed out
    what.resize(0);
    for ( i = 1; i <= totalNumberOfDomainEquation; i++ ) {
        if ( !microDOFs.contains(i) ) {
            what.followedBy(i);
        }
    }

    int ii, k;
    int nkon = what.giveSize(); //how many DOF will be condensed out
    int size = stiffnessMatrixMicroFloat.giveNumberOfRows();
    int ndofs = totalNumberOfDomainEquation;
    long double coeff, dii;


    for ( i = 1; i <= nkon; i++ ) {
        ii  = what.at(i);
        if ( ( ii > ndofs ) || ( ii <= 0 ) ) {
            OOFEM_ERROR("condense: wrong DOF number");
        }

        dii = stiffnessMatrixMicroFloat.at(ii, ii);

        for ( j = 1; j <= size; j++ ) {
            coeff = -stiffnessMatrixMicroFloat.at(j, ii) / dii;
            if ( ii != j ) {
                for ( k = 1; k <= size; k++ ) {
                    stiffnessMatrixMicroFloat.at(j, k) += stiffnessMatrixMicroFloat.at(ii, k) * coeff;
                }
            }
        }

        for ( k = 1; k <= size; k++ ) {
            stiffnessMatrixMicroFloat.at(ii, k) = 0.;
            stiffnessMatrixMicroFloat.at(k, ii) = 0.;
        }
    }

    //    printf("\nstiffnessMatrixMicroFloat\n");
    //    for(i=1;i<=36;i++){
    //        printf("%d:%.3f  ", i,stiffnessMatrixMicroFloat.at(3,i));
    //    }

    //now stiffness matrix is condensed. Node 1 has 3 DOFs and corresponds to code numbers 1,2,3 etc. 0's are everywhere on condensed out rows and colums. It is necessary to leave them out and rearrange unknowns so they correspond to node numbering specified by keyword "micronodes" in macrolspace.C, 3x8=24 DOFs



    //   printf("\n");
    //microDOFs.printYourself();
    //   printf("\n");
    for ( i = 1; i <= microDOFs.giveSize(); i++ ) { //24 components
        for ( j = 1; j <= 24; j++ ) {
            answer.at(i, j) = stiffnessMatrixMicroFloat.at( microDOFs.at(i), microDOFs.at(j) ); //row copy
            answer.at(j, i) = stiffnessMatrixMicroFloat.at( microDOFs.at(j), microDOFs.at(i) ); //column copy
        }
    }

    //answer.printYourself();

    //stiffnessMatrixMicroFloat.printYourself();
}

void MicroMaterial :: setMacroProperties(Domain *macroDomain, MacroLSpace *macroLSpaceElement) {
    this->macroDomain = macroDomain;
    this->macroLSpaceElement = macroLSpaceElement;
}




//DOF number is from 1 to 3 for every node
int MicroMaterial :: giveDofEquationNumber(Dof *dof) const {
    int answer;

    if ( dof->giveDofManager()->giveNumberOfDofs() != 3 ) {
        OOFEM_ERROR("Needs three degrees of freedom at each node");
    }

    //number equation numbers consecutively as going through DoFs
    answer = 3 * ( dof->giveDofManNumber() - 1 ) + dof->giveNumber();
    return answer;
}


int MicroMaterial :: giveRequiredNumberOfDomainEquation() const {
    return totalNumberOfDomainEquation;
}
