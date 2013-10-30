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


#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "shell7base.h"
#include "intmatbilinczfagerstrom.h"
//#include "vld.h"

namespace oofem {

    REGISTER_Material( IntMatBilinearCZFagerstrom );


IntMatBilinearCZFagerstrom :: IntMatBilinearCZFagerstrom(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{
    // constructor
}


IntMatBilinearCZFagerstrom :: ~IntMatBilinearCZFagerstrom()
{
    // destructor
}

int
IntMatBilinearCZFagerstrom :: hasMaterialModeCapability(MaterialMode mode)
{
    // returns whether receiver supports given mode
    if ( mode == _3dInterface ) {
        return 1;
    } else {
        return 0;
    }
}




void 
IntMatBilinearCZFagerstrom :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &d,
                                                     const FloatMatrix &F, TimeStep *tStep)
{
    // returns real stress vector in 3d stress space of receiver according to
    // previous level of stress and current
    // strain increment, the only way, how to correctly update gp records
    //
    
    IntMatBilinearCZFagerstromStatus *status = static_cast< IntMatBilinearCZFagerstromStatus * >( this->giveStatus(gp) );

    this->initGpForNewStep(gp);

    FloatMatrix Finv(3,3);

    Finv.beInverseOf(F);
    status->letTempInverseDefGradBe(Finv);

    FloatArray dJ(3);
    dJ.beProductOf(Finv,d);
    status->letTempMaterialJumpBe(dJ);

    double oldDamage = status->giveDamage();	
    double dAlpha = 0.0;
    FloatArray Qold(3), Qtemp(3);
    Qtemp.zero();



    // 	SUBROUTINE stress_damage_XFEM_direct_Mandel(dJ,N,Qold,old_alpha,Fci,Q,Ea,new_alpha,adtim,dalpha_new,dalpha_old,diss,sig_f,fall);
    if (oldDamage < 0.99) {
        
        #if 0
        FloatArray N;
        FloatMatrix Gcov;

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
        #endif

        FloatArray Qtrial = status->giveEffectiveMandelTraction(); 

        FloatMatrix Kstiff(3,3);
        FloatArray help;

        Kstiff.zero();
        Kstiff.at(1,1) = this->ks0;
        Kstiff.at(2,2) = this->ks0;
        if (dJ.at(3)>=0) {
			Kstiff.at(3,3) = this->kn0;
		} else {
			Kstiff.at(3,3) = this->kn0/(1-oldDamage);
		}

        dJ.subtract(status->giveOldMaterialJump());	 ///@todo Martin: check with Mikael/Jim
        //dJ.rotatedWith(Rot,'n');

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
        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;	///@todo Martin: no use of parameter mu!!!

        Qold = status->giveEffectiveMandelTraction();  ///@todo Martin: check with Mikael/Jim
        Qtemp = Qtrial;

        //Qold.rotatedWith(Rot,'n');

        //double alphaOld = status->giveDamage();

        //if (alphaOld>0.1) {
        //    double bbb=1; ///@todo Should this be used for anything Martin? /JB
        //}

        if (loadFun/sigf < 0.0000001) {
            dAlpha = 0.0;	// new_alpha=old_alpha
		    status->letTempEffectiveMandelTractionBe(Qtemp);		
			Qtemp.times(1-oldDamage);
		} else {
            // dalpha = datr
            double Qt1,Qt2;
            Qt1 = Qtemp.at(1);
            Qt2 = Qtemp.at(2);

            FloatArray M(3);						// M = (/2*Q_t/(sig_f*gamma**2), 2*Q_n/sig_f, 0.d0/)
            M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);	// Qt = sqrt(Qt1^2 + Qt2^2)
            M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
            M.at(3) = 2*Qn_M/sigf;

            FloatArray dJElastic;
            dJElastic = dJ;							// dJtn_e(1:2) = dJtn_v(1:2) - S*dalpha*M(1:2)
            help = M;
            help.times(S*dAlpha);
            dJElastic.subtract(help);

            FloatMatrix Smat(4,4), Smati(4,4);
            FloatArray R(4), inc(4);

            const double errorTol = 0.0001;
            for( int iter = 1; fabs(loadFun)/sigf > errorTol; iter++) {
                if (iter>40) {
                    OOFEM_ERROR("BilinearCZMaterialFagerstrom :: giveRealStressVector - no convergence in constitutive driver");
                    }
                Smat.zero();	// S_mat=0.d0

                R.at(1) = dJElastic.at(1) - (dJ.at(1) - gammaGf*S*dAlpha*M.at(1));		// R(1:2) = dJtn_e(1:2) - (dJtn_v(1:2) - S*dalpha*M(1:2))
                R.at(2) = dJElastic.at(2) - (dJ.at(2) - gammaGf*S*dAlpha*M.at(2));
                R.at(3) = dJElastic.at(3) - (dJ.at(3) - S*dAlpha*M.at(3));
                R.at(4) = loadFun;	// R(3) = F/sig_f

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

                Smat.at(1,4) = gammaGf*S*M.at(1);			// S_mat(1:2,3) = S*M(1:2)
                Smat.at(2,4) = gammaGf*S*M.at(2);
                Smat.at(3,4) = S*M.at(3);
                    
                Smat.at(4,1) = M.at(1)*Kstiff.at(1,1);      // S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
                Smat.at(4,2) = M.at(2)*Kstiff.at(2,2);
                Smat.at(4,3) = M.at(3)*Kstiff.at(3,3);

                //FloatArray inc;
                //bool transpose = false;
                //Smat.SolveforRhs(R, inc, transpose);

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

                M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);    // Qt = sqrt(Qt1^2 + Qt2^2)
                M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
                M.at(3) = 2*Qn_M/sigf;	

                loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;
            }
                
            FloatMatrix Iep(3,3);
            Iep = Smati;
            status->letTempIepBe(Iep);
                                                                
            FloatArray alpha_v(3);
            alpha_v.at(1) = Smati.at(4,1);			// alpha_v(1:2) = S_mati(3,1:2)
            alpha_v.at(2) = Smati.at(4,2);
            alpha_v.at(3) = Smati.at(4,3);
            status->letTempAlphavBe(alpha_v);
			
			dJ = status->giveTempJump();

			if (dJ.at(3)>=0) {
			    status->letTempEffectiveMandelTractionBe(Qtemp);		
				Qtemp.times(1-oldDamage-dAlpha);
			} else {
				if (oldDamage + dAlpha<1) {
					Qtemp.at(3) = (1-oldDamage)/(1-oldDamage + dAlpha)*Qtemp.at(3);
					status->letTempEffectiveMandelTractionBe(Qtemp);
					Qtemp.times(1-oldDamage-dAlpha);
				} else {
					status->letTempEffectiveMandelTractionBe(Qtemp);						// SHOULD NEVER BE USED
					Qtemp.times(1-oldDamage-dAlpha);
					Qtemp.at(3) = (1-oldDamage)*(Qold.at(3) + Kstiff.at(3,3)*dJElastic.at(3));
				}
			}
		}
			
        

    //Qtemp.rotatedWith(Rot,'t');							// Q=Qe
    //status->letTempRotationMatrix(Rot);
    } else {
        dAlpha = 1.0 - oldDamage;
		dJ = status->giveTempJump();
	    status->letTempEffectiveMandelTractionBe(Qtemp);		// SHOULD NEVER BE USED!!
		if (dJ.at(3)<0) {
			Qtemp.at(3) = kn0*dJ.at(3);
		}
    }

    answer.beTProductOf(Finv,Qtemp);					// t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
//    answer.times(1-oldDamage-dAlpha);					// t1_s = (1-al)*t_1_hat

        
    status->letTempDamageBe(oldDamage + dAlpha);
//    status->letTempEffectiveMandelTractionBe(Qtemp);		// NEW!

    status->letTempJumpBe(d);
    status->letTempFirstPKTractionBe(answer);
    status->letTempFBe(F);
}


void

IntMatBilinearCZFagerstrom :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{

    IntMatBilinearCZFagerstromStatus *status = static_cast< IntMatBilinearCZFagerstromStatus * >( this->giveStatus(gp) );

    double damage = status->giveTempDamage();
    FloatMatrix Finv = status->giveTempInverseDefGrad();
    FloatMatrix help;
    FloatMatrix Kstiff(3,3);
	FloatArray J = status->giveTempJump();

    //FloatMatrix Rot = status->giveTempRotationMatrix();
    Kstiff.zero();
    Kstiff.at(1,1) = this->ks0;
    Kstiff.at(2,2) = this->ks0;
    Kstiff.at(3,3) = this->kn0;
    //Kstiff.rotatedWith(Rot);

    if (damage >= 1.0) {
        answer.resize(3,3);
        answer.zero();
		if (J.at(3)<0) {
			Kstiff.at(1,1) = 0.0;
			Kstiff.at(2,2) = 0.0;
			help.beProductOf(Kstiff,Finv);
            answer.beTProductOf(Finv,help);
			//answer.printYourself();
		}
    } else {
        
        if ( status->giveTempDamage() - status->giveDamage()==0.0 ) {
            //Kstiff.printYourself();
            //Rot.printYourself();
            //Kstiff.rotatedWith(Rot);
            //Kstiff.printYourself();
			if (J.at(3)<0) {
				Kstiff.at(3,3) = (this->kn0)/(1-damage);
			}
			help.beProductOf(Kstiff,Finv);
            answer.beTProductOf(Finv,help);
            answer.times(1-damage);						// Ea=(1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Keye3,Fci))
        } else {
            FloatMatrix Iep = status->giveTempIep();
			
			//Iep.printYourself();
			
			if (J.at(3)<0) {
				Kstiff.at(1,1) = (1-damage)*Kstiff.at(1,1);
				Kstiff.at(2,2) = (1-damage)*Kstiff.at(2,2);
				Kstiff.at(3,3) = Kstiff.at(3,3);
			} else {
				Kstiff.times((1-damage));
			}

            answer.beProductOf(Kstiff, Iep);
            //answer.rotatedWith(Rot);						// Ea_h = MATMUL(TRANSPOSE(Rot),MATMUL(Keye3,Iep))
                                                        // Ea_h = MATMUL(Ea_h,Rot)

            FloatArray alpha_v = status->giveTempAlphav();
            //alpha_v.rotatedWith(Rot, 't');					// alpha_v = MATMUL(TRANSPOSE(Rot),alpha_v)

            FloatMatrix t1hatFinvOpen;
            FloatArray temp1, temp2, Qtemp;
            Qtemp = status->giveTempEffectiveMandelTraction();

            temp1.beTProductOf(Finv,Qtemp);				// CALL gmopen33(MATMUL(TRANSPOSE(Fci),Q),MATMUL(alpha_v,Fci),t1halFci_o)
            temp2.beTProductOf(Finv,alpha_v);

            t1hatFinvOpen.beDyadicProductOf(temp1,temp2);

            help.beProductOf(answer,Finv);			// Ea = (1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Ea_h,Fci)) -&
            answer.beTProductOf(Finv,help);			//			t1halFci_o
//            answer.times(1-damage);
            answer.subtract(t1hatFinvOpen);
        }
    }
                                                            
    //Finv.printYourself();
    //Kstiff.printYourself();
    //answer.printYourself();

}





int
IntMatBilinearCZFagerstrom :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    IntMatBilinearCZFagerstromStatus *status = static_cast< IntMatBilinearCZFagerstromStatus * >( this->giveStatus(aGaussPoint) );
    if ( type == IST_DamageScalar ) {     
        answer.resize(1);
        answer.at(1) = status->giveTempDamage();
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }

}


InternalStateValueType
IntMatBilinearCZFagerstrom :: giveIPValueType(InternalStateType type)
{
    //@todoMartin Insert code here if necessaryfor returning type of internal state
    return StructuralInterfaceMaterial :: giveIPValueType(type);
}






const double tolerance = 1.0e-12; // small number
IRResultType
IntMatBilinearCZFagerstrom :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn0, _IFT_IntMatBilinearCZFagerstrom_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_IntMatBilinearCZFagerstrom_knc);

    this->ks0 = kn0;                        // Defaults to kn0
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_IntMatBilinearCZFagerstrom_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatBilinearCZFagerstrom_g1c);

    this->GIIc = GIc;						//Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_IntMatBilinearCZFagerstrom_g2c);

    IR_GIVE_FIELD(ir, sigf, _IFT_IntMatBilinearCZFagerstrom_sigf);

    IR_GIVE_FIELD(ir, mu, _IFT_IntMatBilinearCZFagerstrom_mu);

    IR_GIVE_FIELD(ir, gamma, _IFT_IntMatBilinearCZFagerstrom_gamma);

    this->checkConsistency();                                // check validity of the material paramters
    this->printYourself();
    return IRRT_OK;
}

void IntMatBilinearCZFagerstrom :: giveInputRecord(DynamicInputRecord &input)
{
	StructuralInterfaceMaterial::giveInputRecord(input);

	input.setField(kn0, _IFT_IntMatBilinearCZFagerstrom_kn);
	input.setField(knc, _IFT_IntMatBilinearCZFagerstrom_knc);
	input.setField(ks0, _IFT_IntMatBilinearCZFagerstrom_ks);

	input.setField(GIc, _IFT_IntMatBilinearCZFagerstrom_g1c);
	input.setField(GIIc, _IFT_IntMatBilinearCZFagerstrom_g2c);

	input.setField(sigf, _IFT_IntMatBilinearCZFagerstrom_sigf);
	input.setField(mu, _IFT_IntMatBilinearCZFagerstrom_mu);
	input.setField(gamma, _IFT_IntMatBilinearCZFagerstrom_gamma);

}

int
IntMatBilinearCZFagerstrom :: checkConsistency()
{
    if ( this->kn0 < 0.0 ) {
        OOFEM_ERROR2("IntMatBilinearCZFagerstrom :: initializeFrom - stiffness kn0 is negative (%.2e)", this->kn0);
    } else if ( this->ks0 < 0.0 ) {
        OOFEM_ERROR2("IntMatBilinearCZFagerstrom :: initializeFrom - stiffness ks0 is negative (%.2e)", this->ks0);
    } else if ( this->GIc < 0.0 ) {
        OOFEM_ERROR2("IntMatBilinearCZFagerstrom :: initializeFrom - GIc is negative (%.2e)", this->GIc);
    } else if ( this->GIIc < 0.0 ) {
        OOFEM_ERROR2("IntMatBilinearCZFagerstrom :: initializeFrom - GIIc is negative (%.2e)", this->GIIc);
    } else if ( this->gamma < 0.0  ) { 
        OOFEM_ERROR2("IntMatBilinearCZFagerstrom :: initializeFrom - gamma (%.2e) is below zero which is unphysical" ,
            this->gamma);
    }
    return 1;
}

void
IntMatBilinearCZFagerstrom  :: printYourself()
{
    printf("Paramters for BilinearCZMaterial: \n");

    printf("-Strength paramters \n");
    printf("  sigf  = %e \n", this->sigf);
    printf("  GIc   = %e \n", this->GIc);
    printf("  GIIc  = %e \n", this->GIIc);
    printf("  gamma = sigfs/sigfn = %e \n", this->gamma);
    printf("  mu    = %e \n", this->GIc);

    printf("-Stiffness parameters \n");
    printf("  kn0  = %e \n", this->kn0);
    printf("  ks0  = %e \n", this->ks0);
    printf("  knc  = %e \n", this->knc);

}


//IntMatBilinearCZFagerstromStatus :: IntMatBilinearCZFagerstromStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
IntMatBilinearCZFagerstromStatus :: IntMatBilinearCZFagerstromStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
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


	#if 0
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
#endif

    Iep = tempFInv;
    alphav = oldMaterialJump;	
}


IntMatBilinearCZFagerstromStatus :: ~IntMatBilinearCZFagerstromStatus()
{ }


void
    IntMatBilinearCZFagerstromStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    ///@todo Martin: check need of this
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    /*
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
    fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
    */
}


void
    IntMatBilinearCZFagerstromStatus :: initTempStatus()
{
    ///@todo Martin: Is this really necessary in this particular case??
    StructuralInterfaceMaterialStatus :: initTempStatus();

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
    IntMatBilinearCZFagerstromStatus :: updateYourself(TimeStep *atTime)
{
    ///@todo Martin: kolla behovet av denna
    StructuralInterfaceMaterialStatus :: updateYourself(atTime);
    damage = tempDamage;
    oldMaterialJump = tempMaterialJump;
    QEffective = tempQEffective;
}


} // end namespace oofem
