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

#include "bilinearczmaterialFagerstrom.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "shell7base.h"

namespace oofem {

    REGISTER_Material( BilinearCZMaterialFagerstrom );

    BilinearCZMaterialFagerstrom :: BilinearCZMaterialFagerstrom(int n, Domain *d) : StructuralMaterial(n, d)
        //
        // constructor
        //
    {
    }


    BilinearCZMaterialFagerstrom :: ~BilinearCZMaterialFagerstrom()
        //
        // destructor
        //
    { }

    int
        BilinearCZMaterialFagerstrom :: hasMaterialModeCapability(MaterialMode mode)
    {
        // returns whether receiver supports given mode
        if ( mode == _3dInterface ) {
            return 1;
        } else {
            return 0;
        }
    }




    void
        BilinearCZMaterialFagerstrom :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
        const FloatArray &jumpVector,
        TimeStep *atTime)
        //
        // returns real stress vector in 3d stress space of receiver according to
        // previous level of stress and current
        // strain increment, the only way, how to correctly update gp records
        //
    {
        BilinearCZMaterialFagerstromStatus *status = static_cast< BilinearCZMaterialFagerstromStatus * >( this->giveStatus(gp) );

        this->initGpForNewStep(gp);

        //@todoMartin Insert code for cohesive zone law
        // jumpVector = [dx dy dz]

        answer.resize(3);
        answer.zero();

        FloatMatrix Finv(3,3), F(3,3);
        FloatArray d(3), dJ(3);


        d.at(1) = jumpVector.at(1);
        d.at(2) = jumpVector.at(2);
        d.at(3) = jumpVector.at(3);

        F.at(1,1) = jumpVector.at(4);
        F.at(2,2) = jumpVector.at(5);
        F.at(3,3) = jumpVector.at(6);
        F.at(2,3) = jumpVector.at(7);
        F.at(1,3) = jumpVector.at(8);
        F.at(1,2) = jumpVector.at(9);
        F.at(2,1) = jumpVector.at(10);
        F.at(3,1) = jumpVector.at(11);
        F.at(3,2) = jumpVector.at(12);

        double xi = jumpVector.at(13);  //local xi coord ///@todo remove as input argument (should not be needed)

        Finv.beInverseOf(F);

        status->letTempInverseDefGradBe(Finv);

        dJ.beProductOf(Finv,d);
        status->letTempMaterialJumpBe(dJ);

        double oldDamage = status->giveDamage();	
        double dAlpha = 0.0;
        FloatArray Qold(3), Qtemp(3);
        Qtemp.zero();



        // 	SUBROUTINE stress_damage_XFEM_direct_Mandel(dJ,N,Qold,old_alpha,Fci,Q,Ea,new_alpha,adtim,dalpha_new,dalpha_old,diss,sig_f,fall);
        if (oldDamage < 0.99) {

            FloatArray N;
            FloatMatrix Gcov;
            //ShellInterface *shell = gp->giveElement()->giveInterface(IT_ShellInterface);
            Shell7Base *shell = dynamic_cast< Shell7Base* >(gp->giveElement());
            if ( !shell ) {
                OOFEM_ERROR("BilinearCZMaterialFagerstrom :: giveRealStressVector - oh no wrong element type");
            }
            FloatArray lCoords(3);
            lCoords.at(1) = gp->giveCoordinate(1);
            lCoords.at(2) = gp->giveCoordinate(2);
            lCoords.at(3) = xi;
            shell->evalInitialCovarBaseVectorsAt(lCoords, Gcov);

            FloatArray G1, G2;
            Gcov.copyColumn(G1,1);
            Gcov.copyColumn(G2,2);
            N.beVectorProductOf(G1,G2);
            N.normalize();

            FloatArray Qtrial = status->giveEffectiveMandelTraction(); 

            FloatMatrix Rot(3,3);
            G1.normalize();
            G2.beVectorProductOf(N, G1);
            Rot.setColumn(G1,1);
            Rot.setColumn(G2,2);
            Rot.setColumn(N,3);		

            Qtrial.rotatedWith(Rot,'n');

            FloatMatrix Kstiff(3,3);
            FloatArray help;

            Kstiff.zero();
            Kstiff.at(1,1) = this->ks0;
            Kstiff.at(2,2) = this->ks0;
            Kstiff.at(3,3) = this->kn0;

            dJ.subtract(status->giveOldMaterialJump());	 //@todo Martin: check with Mikael/Jim
            dJ.rotatedWith(Rot,'n');

            help.beProductOf(Kstiff,dJ);
            Qtrial.add(help); 

            double Qn = Qtrial.at(3);
            FloatArray QN(3);
            QN.zero();
            QN.at(3) = Qn;

            FloatArray QtrialShear;
                        
            QtrialShear = Qtrial;
            QtrialShear.subtract(QN);

            double Qt = QtrialShear.computeNorm();


            double S = this->GIc/this->sigf;
            double sigf = this->sigf;
            double gamma = this->gamma;
            double gammaGf = this->GIIc/this->GIc;

            double Qn_M = 0.5*(Qn + fabs(Qn));
            double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;	//@todo Martin: no use of parameter mu!!!

            Qold = status->giveEffectiveMandelTraction();  //@todo Martin: check with Mikael/Jim
            Qtemp = Qtrial;

            Qold.rotatedWith(Rot,'n');

            double alphaOld = status->giveDamage();

            if (alphaOld>0.1) {
                double bbb=1;       ///@todo Should this be used for anything Martin? /JB
            }

            if (loadFun/sigf < 0.0000001) {
                dAlpha = 0.0;											// new_alpha=old_alpha
            } else {
                                                                // dalpha = datr
                double Qt1,Qt2;
                Qt1 = Qtemp.at(1);
                Qt2 = Qtemp.at(2);

                FloatArray M(3);										// M = (/2*Q_t/(sig_f*gamma**2), 2*Q_n/sig_f, 0.d0/)
                M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);					// Qt = sqrt(Qt1^2 + Qt2^2)
                M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
                M.at(3) = 2*Qn_M/sigf;

                FloatArray dJElastic;
                dJElastic = dJ;											// dJtn_e(1:2) = dJtn_v(1:2) - S*dalpha*M(1:2)
                help = M;
                help.times(S*dAlpha);
                dJElastic.subtract(help);

                FloatMatrix Smat(4,4), Smati(4,4);
                FloatArray R(4), inc(4);

                const double errorTol = 0.0001;
                for( int iter = 1; fabs(loadFun)/sigf > errorTol; iter++) {
                    if (iter>40) {
                        //OOFEM_ERROR("BilinearCZMaterialFagerstrom :: giveRealStressVector - no convergence in constitutive driver");
                        double a=1;   ///@todo Should this be used for anything Martin? /JB
                    }
                    Smat.zero();										// S_mat=0.d0

                    R.at(1) = dJElastic.at(1) - (dJ.at(1) - gammaGf*S*dAlpha*M.at(1));		// R(1:2) = dJtn_e(1:2) - (dJtn_v(1:2) - S*dalpha*M(1:2))
                    R.at(2) = dJElastic.at(2) - (dJ.at(2) - gammaGf*S*dAlpha*M.at(2));
                    R.at(3) = dJElastic.at(3) - (dJ.at(3) - S*dAlpha*M.at(3));
                    R.at(4) = loadFun;	 // R(3) = F/sig_f

                    // dMdJtn_e(1,1:2) = (/2/(sig_f*gamma**2),0.d0/)
                    // IF (Q_nM>0) THEN
                    // dMdJtn_e(2,1:2) = (/0.d0,2/sig_f/)
                    // ELSE
                    // dMdJtn_e(2,1:2) = (/0.d0,0.d0/)
                    // END IF
                                
                    // dMdJtn_e = MATMUL(dMdJtn_e,Keye3(1:2,1:2))

                    Smat.at(1,1) = 1.0 + gammaGf*dAlpha*S*2*Kstiff.at(1,1)/(pow(gamma,2)*sigf);		// S_mat(1:2,1:2) = eye3(1:2,1:2)+ dalpha*S*dMdJtn_e(1:2,1:2)
                    Smat.at(2,2) = 1.0 + gammaGf*dAlpha*S*2*Kstiff.at(2,2)/(pow(gamma,2)*sigf);				
                    //dJElastic.printYourself();
                    if (Qn_M>0) {
                        Smat.at(3,3) = 1.0 + dAlpha*S*2*Kstiff.at(3,3)/(sigf);
                    } else {
                        Smat.at(3,3) = 1.0;
                    }

                    Smat.at(1,4) = gammaGf*S*M.at(1);		// S_mat(1:2,3) = S*M(1:2)
                    Smat.at(2,4) = gammaGf*S*M.at(2);
                    Smat.at(3,4) = S*M.at(3);
                    
                    Smat.at(4,1) = M.at(1)*Kstiff.at(1,1);	// S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
                    Smat.at(4,2) = M.at(2)*Kstiff.at(2,2);
                    Smat.at(4,3) = M.at(3)*Kstiff.at(3,3);

                    //FloatArray inc;
                    //bool transpose = false;
                    //Smat.SolveforRhs(R, inc, transpose);

                    Smati.zero();
                    Smati.beInverseOf(Smat);

                    inc.beProductOf(Smati,R);

                    dJElastic.at(1) = dJElastic.at(1) - inc.at(1);
                    dJElastic.at(2) = dJElastic.at(2) - inc.at(2);
                    dJElastic.at(3) = dJElastic.at(3) - inc.at(3);
                    dAlpha = dAlpha - inc.at(4);

                    Qtemp.at(1) = Qold.at(1) + Kstiff.at(1,1)*dJElastic.at(1);
                    Qtemp.at(2) = Qold.at(2) + Kstiff.at(2,2)*dJElastic.at(2);
                    Qtemp.at(3) = Qold.at(3) + Kstiff.at(3,3)*dJElastic.at(3);

                    Qt1 = Qtemp.at(1);
                    Qt2 = Qtemp.at(2);
                    Qt = sqrt(pow(Qt1,2) + pow(Qt2,2));
                    Qn_M = 0.5*(Qtemp.at(3)+fabs(Qtemp.at(3)));

                    M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);	// Qt = sqrt(Qt1^2 + Qt2^2)
                    M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
                    M.at(3) = 2*Qn_M/sigf;	

                    loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;
                }
                
                FloatMatrix Iep(3,3);
                Iep = Smati;
                /*
                Iep.at(1,1) = Smati.at(1,1);
                Iep.at(1,2) = Smati.at(1,2);
                Iep.at(1,3) = Smati.at(1,3);
                Iep.at(2,1) = Smati.at(2,1);
                Iep.at(2,2) = Smati.at(2,2);
                Iep.at(2,3) = Smati.at(2,3);
                Iep.at(3,1) = Smati.at(3,1);
                Iep.at(3,2) = Smati.at(3,2);
                Iep.at(3,3) = Smati.at(3,3);
                */
                status->letTempIepBe(Iep);
                                                                
                FloatArray alpha_v(3);
                alpha_v.at(1) = Smati.at(4,1);      // alpha_v(1:2) = S_mati(3,1:2)
                alpha_v.at(2) = Smati.at(4,2);
                alpha_v.at(3) = Smati.at(4,3);
                status->letTempAlphavBe(alpha_v);


            }

            Qtemp.rotatedWith(Rot,'t');				// Q=Qe
            status->letTempRotationMatrix(Rot);
        } else {
            dAlpha = 1.0-oldDamage;
        }

        answer.beTProductOf(Finv,Qtemp);			// t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
        answer.times(1-oldDamage-dAlpha);			// t1_s = (1-al)*t_1_hat

        
        status->letTempDamageBe(oldDamage + dAlpha);
        status->letTempEffectiveMandelTractionBe(Qtemp);		// NEW!
        //printf("damage %e \n", oldDamage + dAlpha );
    }

void
BilinearCZMaterialFagerstrom :: giveStiffnessMatrix(FloatMatrix &answer,
                                                    MatResponseMode rMode,
                                                    GaussPoint *gp, 
                                                    TimeStep *atTime)
        //
        // Returns characteristic material stiffness matrix of the receiver
        //
    {
        MaterialMode mMode = gp->giveMaterialMode();
        switch ( mMode ) {
        case _3dInterface:
        case _3dMat:
            give3dInterfaceMaterialStiffnessMatrix(answer, rMode, gp, atTime);
            break;
        default:
            //StructuralMaterial :: giveCharacteristicMatrix(answer, rMode, gp, atTime);
            StructuralMaterial :: give3dMaterialStiffnessMatrix(answer, rMode, gp, atTime);
        }
    }


    void
        BilinearCZMaterialFagerstrom :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
        GaussPoint *gp, TimeStep *atTime)
    {
        BilinearCZMaterialFagerstromStatus *status = static_cast< BilinearCZMaterialFagerstromStatus * >( this->giveStatus(gp) );


        double damage = status->giveTempDamage();
        FloatMatrix Finv = status->giveTempInverseDefGrad();
        FloatMatrix help;
        FloatMatrix Kstiff(3,3);
        FloatMatrix Rot = status->giveTempRotationMatrix();
        Kstiff.zero();
        Kstiff.at(1,1) = this->ks0;
        Kstiff.at(2,2) = this->ks0;
        Kstiff.at(3,3) = this->kn0;
        Kstiff.rotatedWith(Rot);

        if (damage >= 1.0) {
            answer.resize(3,3);
            answer.zero();
        } else {
        
            if (status->giveTempDamage()-status->giveDamage()==0.0) {
                //Kstiff.printYourself();
                //Rot.printYourself();
                //Kstiff.rotatedWith(Rot);
                //Kstiff.printYourself();
                help.beProductOf(Kstiff,Finv);
                answer.beTProductOf(Finv,help);
                answer.times(1-damage);						// Ea=(1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Keye3,Fci))
            } else {

                FloatMatrix Iep = status->giveTempIep();

                answer.beProductOf(Kstiff, Iep);
                answer.rotatedWith(Rot);						// Ea_h = MATMUL(TRANSPOSE(Rot),MATMUL(Keye3,Iep))
                                                            // Ea_h = MATMUL(Ea_h,Rot)

                FloatArray alpha_v = status->giveTempAlphav();
                alpha_v.rotatedWith(Rot, 't');					// alpha_v = MATMUL(TRANSPOSE(Rot),alpha_v)

                FloatMatrix t1hatFinvOpen;
                FloatArray temp1, temp2, Qtemp;
                Qtemp = status->giveTempEffectiveMandelTraction();

                temp1.beTProductOf(Finv,Qtemp);				// CALL gmopen33(MATMUL(TRANSPOSE(Fci),Q),MATMUL(alpha_v,Fci),t1halFci_o)
                temp2.beTProductOf(Finv,alpha_v);

                t1hatFinvOpen.beDyadicProductOf(temp1,temp2);

                help.beProductOf(answer,Finv);			// Ea = (1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Ea_h,Fci)) -&
                answer.beTProductOf(Finv,help);			//			t1halFci_o
                answer.times(1-damage);
                answer.subtract(t1hatFinvOpen);
            }
        }
                                                            
        //@todoMartin Insert code for full compressive stiffness!!!
        //Finv.printYourself();
        //Kstiff.printYourself();
        //answer.printYourself();

    }





int
BilinearCZMaterialFagerstrom :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_DamageScalar ) {
        BilinearCZMaterialFagerstromStatus *status = static_cast< BilinearCZMaterialFagerstromStatus * >( this->giveStatus(aGaussPoint) );
        answer.resize(1);
        answer.at(1) = status->giveTempDamage();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }

}


    InternalStateValueType
        BilinearCZMaterialFagerstrom :: giveIPValueType(InternalStateType type)
    {
        //@todoMartin Insert code here if necessaryfor returning type of internal state
        return StructuralMaterial :: giveIPValueType(type);
    }


    int
        BilinearCZMaterialFagerstrom :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
    {
        //@todoMartin check if this needs to be updated
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }


    int
        BilinearCZMaterialFagerstrom :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
    {
        //@todoMartin Insert code for returning size of internal state
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }



    const double tolerance = 1.0e-12; // small number
    IRResultType
        BilinearCZMaterialFagerstrom :: initializeFrom(InputRecord *ir)
    {
        const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
        IRResultType result;                    // Required by IR_GIVE_FIELD macro

        IR_GIVE_FIELD(ir, kn0, _IFT_BilinearCZMaterialFagerstrom_kn);
        this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
        IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_BilinearCZMaterialFagerstrom_knc);

        this->ks0 = kn0;                        // Defaults to kn0
        IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_BilinearCZMaterialFagerstrom_ks);

        IR_GIVE_FIELD(ir, GIc, _IFT_BilinearCZMaterialFagerstrom_g1c);

        this->GIIc = GIc;						//Defaults to GIc
        IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_BilinearCZMaterialFagerstrom_g2c);

        IR_GIVE_FIELD(ir, sigf, _IFT_BilinearCZMaterialFagerstrom_sigf);

        IR_GIVE_FIELD(ir, mu, _IFT_BilinearCZMaterialFagerstrom_mu);

        IR_GIVE_FIELD(ir, gamma, _IFT_BilinearCZMaterialFagerstrom_gamma);

        this->checkConsistency();                                // check validity of the material paramters
        this->printYourself();
        return IRRT_OK;
    }

    int
        BilinearCZMaterialFagerstrom :: checkConsistency()
    {
        if ( this->kn0 < 0.0 ) {
            OOFEM_ERROR2("BilinearCZMaterialFagerstrom :: initializeFrom - stiffness kn0 is negative (%.2e)", this->kn0);
        } else if ( this->ks0 < 0.0 ) {
            OOFEM_ERROR2("BilinearCZMaterialFagerstrom :: initializeFrom - stiffness ks0 is negative (%.2e)", this->ks0);
        } else if ( this->GIc < 0.0 ) {
            OOFEM_ERROR2("BilinearCZMaterialFagerstrom :: initializeFrom - GIc is negative (%.2e)", this->GIc);
        } else if ( this->GIIc < 0.0 ) {
            OOFEM_ERROR2("BilinearCZMaterialFagerstrom :: initializeFrom - GIIc is negative (%.2e)", this->GIIc);
        } else if ( this->gamma < 0.0  ) { 
            OOFEM_ERROR2("BilinearCZMaterialFagerstrom :: initializeFrom - gamma (%.2e) is below zero which is unphysical" ,
                this->gamma);
        }
        return 1;
    }

    void
        BilinearCZMaterialFagerstrom  :: printYourself()
    {
        printf("Paramters for BilinearCZMaterial: \n");

        printf("-Strength paramters \n");
        printf("  sigf = %e \n", this->sigf);
        printf("  GIc   = %e \n", this->GIc);
        printf("  GIIc   = %e \n", this->GIIc);
        printf("  gamma = sigfs/sigfn = %e \n", this->gamma);
        printf("  mu   = %e \n", this->GIc);

        //printf("\n");

        printf("-Stiffness parameters \n");
        printf("  kn0   = %e \n", this->kn0);
        printf("  ks0   = %e \n", this->ks0);
        printf("  knc   = %e \n", this->knc);
        //printf("\n");

    }

    BilinearCZMaterialFagerstromStatus :: BilinearCZMaterialFagerstromStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
    {
        oldMaterialJump.resize(3);
        oldMaterialJump.zero();
        tempMaterialJump = oldMaterialJump;

        damage = tempDamage = 0.0;

        QEffective = oldMaterialJump;
        tempQEffective = oldMaterialJump;
        
        tempFInv.resize(3,3);
        tempFInv.zero();
        tempFInv.at(1,1)=1.;
        tempFInv.at(2,2)=1.;
        tempFInv.at(3,3)=1.;


        //@todo Martin: Very bad implementation of intialisation of Rot
        //*************************************************************
        FloatMatrix Gcov;
        FloatArray N;

        Shell7Base *shell = dynamic_cast< Shell7Base* >(gp->giveElement());
        if ( !shell ) {
            OOFEM_ERROR("BilinearCZMaterialFagerstrom :: giveRealStressVector - oh no wrong element type");
        }
        FloatArray lCoords(3);
        lCoords.at(1) = gp->giveCoordinate(1);
        lCoords.at(2) = gp->giveCoordinate(2);
        shell->evalInitialCovarBaseVectorsAt(lCoords, Gcov);

        FloatArray G1, G2;
        Gcov.copyColumn(G1,1);
        Gcov.copyColumn(G2,2);
        N.beVectorProductOf(G1,G2);
        N.normalize();

        tempRot.resize(3,3);
        G1.normalize();
        G2.beVectorProductOf(N, G1);
        tempRot.setColumn(G1,1);
        tempRot.setColumn(G2,2);
        tempRot.setColumn(N,3);		


        Iep = tempFInv;
        alphav = oldMaterialJump;	
    }


    BilinearCZMaterialFagerstromStatus :: ~BilinearCZMaterialFagerstromStatus()
    { }


    void
        BilinearCZMaterialFagerstromStatus :: printOutputAt(FILE *file, TimeStep *tStep)
    {
        //@todo Martin: kolla behovet av denna
        StructuralMaterialStatus :: printOutputAt(file, tStep);
        /*
        fprintf(file, "status { ");
        if ( this->damage > 0.0 ) {
        fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
        }

        fprintf(file, "}\n");
        */
    }


    void
        BilinearCZMaterialFagerstromStatus :: initTempStatus()
    {
        //@todo Martin: Is this really necessary in this particular case??
        StructuralMaterialStatus :: initTempStatus();

        tempMaterialJump = oldMaterialJump;
        tempDamage = damage;
        tempQEffective = QEffective;

        tempFInv.resize(3,3);
        tempFInv.zero();

        Iep = tempFInv;
        alphav = oldMaterialJump;	

        tempFInv.at(1,1)=1.;
        tempFInv.at(2,2)=1.;
        tempFInv.at(3,3)=1.;



    }

    void
        BilinearCZMaterialFagerstromStatus :: updateYourself(TimeStep *atTime)
    {
        //@todo Martin: kolla behovet av denna
        StructuralMaterialStatus :: updateYourself(atTime);
        damage = tempDamage;
        oldMaterialJump = tempMaterialJump;
        QEffective = tempQEffective;
    }


} // end namespace oofem
