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

#include "micromaterial.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "domain.h"
#include "dofmanager.h"
#include "classfactory.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "classfactory.h"
#include "node.h"
#include "engngm.h"

namespace oofem {
//valgrind --leak-check=full --show-reachable=no -v --log-file=valgr.txt ./oofem -f Macrolspace_1.in

//     FloatArray A;
//     FloatMatrix K, B;
//     stiffnessMatrix->toFloatMatrix(K);
//     K.symmetrized();
//     K.printYourself();
//     displacementVector.printYourself();
//     A.beProductOf(K, displacementVector);
//     A.printYourself();
//
//     B.beInverseOf(K);
//     A.beProductOf(B, displacementVector);
//     A.printYourself();

REGISTER_Material(MicroMaterial);

// constructor
//strainVector, tempStrainVector, stressVector, tempStressVector are defined on StructuralMaterialStatus
MicroMaterialStatus :: MicroMaterialStatus(GaussPoint *gp) : StructuralMaterialStatus(gp) { }


void MicroMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void MicroMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
}

void MicroMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{ }

void MicroMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
}

void MicroMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
}


MicroMaterial :: MicroMaterial(int n, Domain *d) : StructuralMaterial(n, d), UnknownNumberingScheme()
{}


void MicroMaterial :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->inputFileNameMicro, _IFT_MicroMaterial_fileName);

    OOFEM_LOG_INFO( "** Instanciating microproblem with BC from file %s\n", inputFileNameMicro.c_str() );
    OOFEMTXTDataReader drMicro( inputFileNameMicro.c_str() );
    this->problemMicro = InstanciateProblem(drMicro, _processor, 0); //0=contextFlag-store/resore
    drMicro.finish();
    OOFEM_LOG_INFO("Microproblem instanciated\n");
}


//original pure virtual function has to be declared here
//this function should not be used, internal forces are calculated based on reactions not stresses in GPs
FloatArrayF<6> MicroMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &totalStrain, GaussPoint *gp, TimeStep *tStep) const
{
    //perform average over microproblem
    //     int index;
    //     double dV, VolTot = 0.;
    //     double scale = 1.;
    //     FloatArray VecStrain, VecStress, SumStrain(6), SumStress(6);
    //     GaussPoint *gpL;
    //     Domain *microDomain = problemMicro->giveDomain(1); //from engngm.h
    //     EngngModel *microEngngModel = microDomain->giveEngngModel();
    //     StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    OOFEM_ERROR("Should not be called, use giveInternalForcesVector instead");

    //     int nelem = microDomain->giveNumberOfElements();
    //     //int nnodes = microDomain->giveNumberOfDofManagers();
    //
    //     //need to update stress so change boundary conditions and recalculate
    //     OOFEM_LOG_INFO( "\n*** giveRealStress microproblem %p on macroElement %d GP %d, step %d, time %f\n", this, this->macroLSpaceElement->giveNumber(), gp->giveNumber(), microEngngModel->giveCurrentStep()->giveNumber(), microEngngModel->giveCurrentStep()->giveTime() );
    //
    //     this->macroLSpaceElement->changeMicroBoundaryConditions(tStep);
    //     microEngngModel->solveYourselfAt( microEngngModel->giveCurrentStep() );
    //     microEngngModel->terminate( microEngngModel->giveCurrentStep() );
    //     OOFEM_LOG_INFO("\n*** giveRealStress microproblem %p done\n", this);
    //
    //     for ( int ielem = 1; ielem <= nelem; ielem++ ) { //return stress as average through all elements of the same MicroMaterial
    //         Element *elem = microDomain->giveElement(ielem);
    //         for ( GaussPoint *gpL: *elem->giveDefaultIntegrationRulePtr() ) {
    //             dV  = elem->computeVolumeAround(gpL);
    //             VolTot += dV;
    //             //OOFEM_LOG_INFO("Element %d GP %d Vol %f\n", elem->giveNumber(), gp->giveNumber(), dV);
    //             //fprintf(this->stream, "Element %d GP %d stress %f\n", elem->giveNumber(), gp->giveNumber(), 0.0);
    //             //((StructuralCrossSection*) gp->giveCrossSection())->giveFullCharacteristicVector(helpVec, gp, strainVector);
    //             elem->giveIPValue(VecStrain, gpL, IST_StrainTensor, tStep);
    //             elem->giveIPValue(VecStress, gpL, IST_StressTensor, tStep);
    //             elem->giveIntVarCompFullIndx(Mask, IST_StrainTensor);
    //
    //             VecStrain.times(dV);
    //             VecStress.times(dV);
    //
    //             for ( int j = 1; j <= 6; j++ ) {
    //                 SumStrain.at(j) += VecStrain.at(j);
    //                 SumStress.at(j) += VecStress.at(j);
    //             }
    //
    //             //VecStrain.printYourself();
    //             //SumStrain.printYourself();
    //         }
    //     }
    //
    //     //average
    //     SumStrain.times(1. / VolTot * scale);
    //     SumStress.times(1. / VolTot * scale);
    //     //SumStrain.printYourself();
    //     //SumStress.printYourself();
    //     answer.resize(6);
    //     answer = SumStress;
    //     //answer.printYourself();

    // update gp
    //     status->letTempStrainVectorBe(totalStrain);
    //     status->letTempStressVectorBe(answer);
}

MaterialStatus *
MicroMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MicroMaterialStatus(gp);
}


//answer must be of size 24x24 (linear brick 3*8=24)
void MicroMaterial :: giveMacroStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep, MatResponseMode rMode, const IntArray &microMasterNodes, const IntArray &microBoundaryNodes)
{
    int row, col;
    double tmpDouble;
    Domain *microDomain = this->problemMicro->giveDomain(1);
    EngngModel *microEngngModel = microDomain->giveEngngModel();
    DofManager *DofMan;
    Dof *dof;

    FloatMatrix Kbb; //contains reduced problem with boundary nodes and without interior nodes
    FloatMatrix Kbi; //can be zero size if no internal DOFs
    FloatMatrix Kii1KbiT; //can be zero size if no internal DOFs
    FloatMatrix Kbb_1; //help matrix

    FloatMatrix slaveMasterOnBoundary; //transformation matrix representing displacements on boundary tied to master nodes

    SparseMtrxType sparseMtrxType = ( SparseMtrxType ) 0; //SMT_Skyline symmetric skyline
    std::unique_ptr<SparseMtrx> stiffnessMatrixMicro; //full stiffness matrix without any constraint
    std::unique_ptr<SparseMtrx> Kii; //submatrix of internal DOFs

    Kbb.resize(totalBoundaryDofs, totalBoundaryDofs);
    Kbb.zero();

    this->isDefaultNumbering = false; //total number of equations can be now set arbitrarily
    this->reqNumberOfDomainEquation = this->maxNumberOfDomainEquation;

    //assemble sparse matrix K_ii of DOFS of internal nodes to be condensed out
    //K_ii is generally large, inversion of FloatMatrix consumes a lot of time and memory, sparse matrix is used
    this->reqNumberOfDomainEquation = totalInternalDofs;
    printf( "Internal DOFs %d\n", this->giveRequiredNumberOfDomainEquation() );
    this->DofEquationNumbering = InteriorNodes;
    if ( totalInternalDofs ) {
        Kbi.resize(totalBoundaryDofs, totalInternalDofs);
        Kbi.zero();
        Kii1KbiT.resize(totalInternalDofs, totalBoundaryDofs);
        Kii1KbiT.zero();
        Kii = classFactory.createSparseMtrx(sparseMtrxType);
        Kii->buildInternalStructure(microEngngModel, 1, * this);
        microEngngModel->assemble(*Kii, tStep, TangentAssembler(rMode), * this, microDomain);
    }



    //build a full stiffness matrix for the extraction of submatrices
    this->reqNumberOfDomainEquation = this->maxNumberOfDomainEquation;
    this->DofEquationNumbering = AllNodes;

    stiffnessMatrixMicro = classFactory.createSparseMtrx(sparseMtrxType);
    stiffnessMatrixMicro->buildInternalStructure(microEngngModel, 1, * this);
    microEngngModel->assemble(*stiffnessMatrixMicro, tStep, TangentAssembler(rMode), * this, microDomain);


    for ( int i = 1; i <= totalBoundaryDofs; i++ ) {
        for ( int j = 1; j <= totalBoundaryDofs; j++ ) { //Kbb
            row = microBoundaryDofsArr.at(i);
            col = microBoundaryDofsArr.at(j);
            if ( stiffnessMatrixMicro->isAllocatedAt(row, col) ) {
                //printf("%d %d   ", row, col);
                Kbb.at(i, j) = stiffnessMatrixMicro->at(row, col);
            }
        }

        for ( int j = 1; j <= totalInternalDofs; j++ ) { //Kbi
            row = microBoundaryDofsArr.at(i);
            col = microInternalDofsArr.at(j);
            if ( stiffnessMatrixMicro->isAllocatedAt(row, col) ) {
                Kbi.at(i, j) = stiffnessMatrixMicro->at(row, col);
            }
        }
    }

    //Kbi.printYourself();

    if ( totalInternalDofs ) {
        FloatArray xVector( Kii->giveNumberOfColumns() );
        //Kii->printYourself();
        Kii->factorized();

        for ( int i = 1; i <= totalInternalDofs; i++ ) {
            xVector.zero();
            xVector.at(i) = 1.;
            Kii->backSubstitutionWith(xVector);
            //xVector.printYourself();
            for ( int b = 1; b <= totalBoundaryDofs; b++ ) { //do not store Kii^-1, it is a dense matrix, compute multiplication directly
                tmpDouble = 0.;
                for ( int j = 1; j <= totalInternalDofs; j++ ) {
                    tmpDouble += xVector.at(j) * Kbi.at(b, j);
                }

                Kii1KbiT.at(i, b) = tmpDouble;
            }

            //OOFEM_LOG_INFO("%d ", i);
        }

        OOFEM_LOG_INFO("\n");

        Kbb_1.beProductOf( Kbi, Kii1KbiT );
        for ( int i = 1; i <= totalBoundaryDofs; i++ ) {
            for ( int j = 1; j <= totalBoundaryDofs; j++ ) {
                Kbb.at(i, j) -= Kbb_1.at(i, j);
            }
        }
    }

    //Kbb_1-printYourself();
    //OOFEM_ERROR("Stop");
    this->isDefaultNumbering = true; //switch back to default numbering
    this->DofEquationNumbering = AllNodes;

    //     IntArray interiorDofNode; //equation numbers to be condensed out, sorted
    //     IntArray boundaryDofNode;//equation numbers in rows (or columns) of stiffness matrix without interior nodes, sorted
    //     interiorDofNode.clear();
    //     boundaryDofNode.clear();
    //     for ( i = 1; i <= microDomain->giveNumberOfDofManagers(); i++){
    //       DofMan = microDomain->giveDofManager(i);
    //       for ( Dof *dof: *DofMan ){
    //         eqNumber = giveDofEquationNumber(dof);
    //         if(microBoundaryNodes.contains( DofMan->giveGlobalNumber())){
    //           boundaryDofNode.followedBy(eqNumber);
    //         }
    //         else {
    //           interiorDofNode.followedBy(eqNumber);
    //         }
    //       }
    //     }
    //     //Tmp.beInverseOf(stiffnessMatrixMicroFloat);
    //     /*Perform static condensation of internal nodes, leave microBoundaryNodes which also contains microMasterNodes
    //       algorithm based on structuralelement.C, Rayleigh-Ritz method
    //       zeroes on particular rows and columns will appear in the FloatMatrix
    //     */
    //     int ii, k;
    //     int size = stiffnessMatrixMicroFloat.giveNumberOfRows();
    //     int ndofs = maxNumberOfDomainEquation;
    //     long double coeff, dii;
    //
    //     for ( i = 1; i <= interiorDofNode.giveSize(); i++ ) {//how many DOFs will be condensed out
    //         ii  = interiorDofNode.at(i);
    //         if ( ( ii > ndofs ) || ( ii <= 0 ) ) {
    //             OOFEM_ERROR("wrong DOF number");
    //         }
    //
    //         dii = stiffnessMatrixMicroFloat.at(ii, ii);
    //
    //         for ( j = 1; j <= size; j++ ) {
    //             coeff = -stiffnessMatrixMicroFloat.at(j, ii) / dii;
    //             if ( ii != j ) {
    //                 for ( k = 1; k <= size; k++ ) {
    //                     stiffnessMatrixMicroFloat.at(j, k) += stiffnessMatrixMicroFloat.at(ii, k) * coeff;
    //                 }
    //             }
    //         }
    //
    //         for ( k = 1; k <= size; k++ ) {
    //             stiffnessMatrixMicroFloat.at(ii, k) = 0.;
    //             stiffnessMatrixMicroFloat.at(k, ii) = 0.;
    //         }
    //     }
    //     //stiffnessMatrixMicroFloat.printYourself();
    //
    //     //copy non-zeroed rows and columns to reduced stiffness matrix
    //     stiffnessMatrixMicroReducedFloat.resize(boundaryDofNode.giveSize(), boundaryDofNode.giveSize());
    //     stiffnessMatrixMicroReducedFloat.zero();
    //     for ( i = 1; i <= boundaryDofNode.giveSize(); i++ ) {
    //       for ( j = 1; j <= boundaryDofNode.giveSize(); j++ ) {
    //         stiffnessMatrixMicroReducedFloat.at(i,j) = stiffnessMatrixMicroFloat.at( boundaryDofNode.at(i), boundaryDofNode.at(j) );
    //       }
    //     }
    //stiffnessMatrixMicroReducedFloat.printYourself();

    //build the transformation matrix of size (x,24) relating slave nodes on the boundary to master nodes on the boundary
    slaveMasterOnBoundary.resize(Kbb.giveNumberOfColumns(), 24);
    slaveMasterOnBoundary.zero();

    FloatArray n(8);
    int node, nodePos;

    //master nodes
    //     for ( i = 1; i <= microMasterNodes.giveSize(); i++ ) {//8 nodes
    //       node = microMasterNodes.at(i);
    //       DofMan = microDomain->giveDofManager(node);
    //       dof = DofMan->giveDofWithID(D_u);
    //       j = boundaryDofNode.findFirstIndexOf(giveDofEquationNumber(dof));
    //       if(!j)
    //         OOFEM_ERROR("Not found equation number %d in reduced stiffness matrix of node %d\n", giveDofEquationNumber(dof), DofMan->giveGlobalNumber());
    //       microMasterNodesLoc.followedBy(j);
    //     }


    //boundary nodes
    for ( int i = 1; i <= microBoundaryNodes.giveSize(); i++ ) {
        node = microBoundaryNodes.at(i);
        DofMan = microDomain->giveDofManager(node);
        dof = DofMan->giveDofWithID(D_u);
        nodePos = microBoundaryDofsArr.findFirstIndexOf( giveDofEquationNumber(dof) ); //row(column) of reduced stiffness matrix
        if ( !nodePos ) {
            OOFEM_ERROR("Not found equation number %d in reduced stiffness matrix of node %d\n", giveDofEquationNumber(dof), DofMan->giveGlobalNumber() );
        }

        this->macroLSpaceElement->evalInterpolation( n, this->microMasterCoords, DofMan->giveCoordinates() );

        //construct transformation matrix relating displacement of slave boundary nodes to macroelement nodes
        //the answer is returned to macroelement so the columns correspond to x,y,z DOFs of each macroelement node

        for ( int j = 1; j <= this->macroLSpaceElement->giveNumberOfNodes(); j++ ) { //linhex - 8 nodes
            for ( int k = 0; k < 3; k++ ) { //the same interpolation for x,y,z
                row = nodePos + k;
                col = 3 * j + k - 2;
                if ( row > slaveMasterOnBoundary.giveNumberOfRows() ) {
                    OOFEM_ERROR("Row is outside the reduced stiffness matrix ");
                }

                if ( col > slaveMasterOnBoundary.giveNumberOfColumns() ) {
                    OOFEM_ERROR("Column is outside the reduced stiffness matrix ");
                }

                slaveMasterOnBoundary.at(row, col) = n.at(j);
            }
        }
    }

    //slaveMasterOnBoundary.printYourself();
#  if 0
    //check of the transformation matrix - the sum of each third column must be either zero or one
    double sum;
    for ( int i = 1; i <= slaveMasterOnBoundary.giveNumberOfRows(); i++ ) {
        sum = 0;
        for ( int j = 1; j <= slaveMasterOnBoundary.giveNumberOfColumns(); j++ ) {
            if ( j % 3 == 0 ) {
                sum += slaveMasterOnBoundary.at(i, j);
            }
        }

        OOFEM_LOG_INFO("Sum of %i row of transformation matrix row %f\n", i, sum);
    }

#  endif

    //slaveMasterOnBoundary.printYourself();

    //perform K(24,24) = T^T * K * T
    FloatMatrix A;
    A.beProductOf(Kbb,  slaveMasterOnBoundary);
    //A.printYourself();
    //slaveMasterOnBoundary.printYourself();
    //answer.resize(24, 24);
    //answer.zero();
    //OOFEM_ERROR("Stop");
    answer.beTProductOf(slaveMasterOnBoundary, A);
}

//should be called before the calculation of micromaterial
void MicroMaterial :: setMacroProperties(Domain *macroDomain, MacroLSpace *macroLSpaceElement, const IntArray &microMasterNodes, const IntArray &microBoundaryNodes)
{
    Domain *microDomain = this->problemMicro->giveDomain(1);
    EngngModel *microEngngModel = microDomain->giveEngngModel();
    MetaStep *mstep;
    DofManager *DofMan;

    int numDofs, numDofMan;
    int counterDefault = 1, counterBoundary = 1, counterInternal = 1;

    this->macroDomain = macroDomain;
    this->macroLSpaceElement = macroLSpaceElement;

    this->microMasterCoords.clear();
    this->microMasterCoords.reserve(microMasterNodes.giveSize());
    for ( int nodeNum: microMasterNodes ) { //8 nodes
        this->microMasterCoords.push_back( microDomain->giveNode( nodeNum )->giveCoordinates() );
    }

    microEngngModel->giveNextStep(); //set the first time step
    mstep = microEngngModel->giveCurrentMetaStep();
    mstep->setNumberOfSteps(this->macroDomain->giveEngngModel()->giveMetaStep(1)->giveNumberOfSteps() + 1);

    //separate DOFs into boundary and internal
    this->NumberOfDofManagers = microDomain->giveNumberOfDofManagers();
    microBoundaryDofs.resize( this->NumberOfDofManagers );
    microInternalDofs.resize( this->NumberOfDofManagers );
    microDefaultDofs.resize( this->NumberOfDofManagers );
    for ( int i = 0; i < this->NumberOfDofManagers; i++ ) {
        microBoundaryDofs [ i ].resize(3);
        microInternalDofs [ i ].resize(3);
        microDefaultDofs [ i ].resize(3);
        microBoundaryDofs [ i ].zero();
        microInternalDofs [ i ].zero();
        microDefaultDofs [ i ].zero();
    }

    for ( int i = 0; i < this->NumberOfDofManagers; i++ ) {
        DofMan = microDomain->giveDofManager(i + 1);
        numDofMan = DofMan->giveGlobalNumber();
        numDofs = DofMan->giveNumberOfDofs();
        if ( numDofs != 3 ) {
            OOFEM_ERROR("Node %d does not have three degrees of freedom", DofMan->giveGlobalNumber() );
        }

        for ( int j = 0; j < numDofs; j++ ) {
            microDefaultDofs [ i ] [ j ] = counterDefault; //equation number
            if ( microBoundaryNodes.contains(numDofMan) ) { //boundary (and master) nodes
                microBoundaryDofs [ i ] [ j ] = counterBoundary;
                microBoundaryDofsArr.followedBy(counterDefault);
                counterBoundary++;
            } else { //internal nodes to be condensed out
                microInternalDofs [ i ] [ j ] = counterInternal;
                microInternalDofsArr.followedBy(counterDefault);
                counterInternal++;
            }

            counterDefault++;
        }
    }

    this->totalBoundaryDofs = counterBoundary - 1;
    this->totalInternalDofs = counterInternal - 1;

    //     //the pointer to underlying problem is problemMicro
    //   DofManager *DofMan;
    //   int i, j, counter = 0;
    //   IntArray ut, dofIDArry = {D_u, D_v, D_w};1
    //
    //   for ( i = 1; i <= problemMicro->giveDomain(1)->giveNumberOfDofManagers(); i++ ) { //for each node
    //     DofMan = problemMicro->giveDomain(1)->giveDofManager(i);
    //         //printf("%d\n",DofMan->giveNumberOfPrimaryMasterDofs(dofIDArry));
    //         //DofMan->giveLocationArray(dofIDArry, *this );
    //     DofMan->giveCompleteLocationArray(* this);
    //     for ( j = 1; j <= DofMan->giveNumberOfDofs(); j++ ) {
    //       counter++;
    //     }
    //   }

    this->maxNumberOfDomainEquation = counterBoundary - 1 + counterInternal - 1;
}



//from class UnknownNumberingScheme
void MicroMaterial :: init(void)
{
    if ( this->problemMicro->giveNumberOfDomains() != 1 ) {
        OOFEM_ERROR("Number of domains on microproblem is greater than 1");
    }

    microMatIsUsed = true;
}

//from class UnknownNumberingScheme
//Each node has three dofs (x,y,z direction)
int MicroMaterial :: giveDofEquationNumber(Dof *dof) const
{
    int numDofMan = dof->giveDofManNumber();
    int numDof = dof->giveDofID(); // Note: Relieas on D_u, D_v, D_w being 1, 2, 3

    //depending on the assembly of submatrix, swith to equation numbering
    switch ( DofEquationNumbering ) {
    case AllNodes: //the default
        return  microDefaultDofs [ numDofMan - 1 ] [ numDof - 1 ];
        break;
    case BoundaryNodes:
        return  microBoundaryDofs [ numDofMan - 1 ] [ numDof - 1 ];
        break;
    case InteriorNodes:
        return  microInternalDofs [ numDofMan - 1 ] [ numDof - 1 ];
        break;
    default:
        OOFEM_ERROR("Node numbering undefined");
    }
}

//from class UnknownNumberingScheme
int MicroMaterial :: giveRequiredNumberOfDomainEquation() const
{
    return this->reqNumberOfDomainEquation;
}
} // end namespace oofem
