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
#include "intmatbilinczjansson.h"

namespace oofem {

REGISTER_Material( IntMatBilinearCZJansson );


IntMatBilinearCZJansson :: IntMatBilinearCZJansson(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }


IntMatBilinearCZJansson :: ~IntMatBilinearCZJansson() { }


void 
IntMatBilinearCZJansson :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &d,
                                                     const FloatMatrix &F, TimeStep *tStep)
{
    // returns real stress vector in 3d stress space of receiver according to
    // previous level of stress and current
    // strain increment, the only way, how to correctly update gp records
    
    IntMatBilinearCZJanssonStatus *status = static_cast< IntMatBilinearCZJanssonStatus * >( this->giveStatus(gp) );

    FloatMatrix Finv;

    Finv.beInverseOf(F);
    status->letTempInverseDefGradBe(Finv);

    FloatArray dJ;
    dJ.beProductOf(Finv,d);


    status->letTempMaterialJumpBe(dJ);

    double oldDamage = status->giveDamage();
    double dAlpha = 0.0;
    FloatArray Qold(3), Qtemp(3), Qtemp_comp(3);
    Qtemp.zero();
    Qtemp_comp.zero();

    if ( dJ.at(3) < 0 ) {
        Qtemp_comp.at(3) = this->knc*dJ.at(3);
    }


    if (oldDamage < 0.99) {
        
        FloatArray Qtrial = status->giveEffectiveMandelTraction(); 

        FloatMatrix Kstiff(3,3);
        FloatArray help;

        Kstiff.zero();
        Kstiff.at(1,1) = this->ks0;
        Kstiff.at(2,2) = this->ks0;
        Kstiff.at(3,3) = this->kn0;
        //} else {
        //    Kstiff.at(3,3) = this->kn0/(1-oldDamage);
        //}

        dJ.subtract(status->giveOldMaterialJump());

        help.beProductOf(Kstiff, dJ);
        Qtrial.add(help); 

        double Qn = Qtrial.at(3);
        FloatArray QN(3);
        QN.zero();
        QN.at(3) = Qn;

        FloatArray QtrialShear;
                        
        QtrialShear = Qtrial;
        QtrialShear.subtract(QN);

        double Qt = QtrialShear.computeNorm();


//        double S = this->GIc/this->sigf;
        double sigf = this->sigf;
        double gamma = this->gamma;
//        double gammaGf = this->GIIc/this->GIc;
        double mu = this->mu;
        double c = mu/gamma; 

        double Qn_M = 0.5*(Qn + fabs(Qn));
//        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;    ///@todo Martin: no use of parameter mu!!!

        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*(1-c)*pow((Qn_M/sigf),2) + sigf*c*(Qn/sigf) - sigf;

        Qold = status->giveEffectiveMandelTraction();  
        Qtemp = Qtrial;

        //Qold.rotatedWith(Rot,'n');

    
        if (loadFun/sigf < 0.0000001) {
            dAlpha = 0.0;   // new_alpha=old_alpha
            status->letTempEffectiveMandelTractionBe(Qtemp);
            Qtemp.times(1-oldDamage);
        } else {
            status->letTempDamageDevBe(true);

             // dalpha = datr
            double C1,C2;
            
            C1 = (pow((Qt/gamma),2)+(1-c)*pow(Qn_M,2))/(pow(sigf,2));
            C2 = c*Qn_M/sigf;

            //double xi = (-C2 + sqrt(pow(C2,2)+(1-c*Qn_M/sigf)*4*C1))/(2*C1);
            double xi = 0.0;

            if (Qn >=0) {
                xi = (-C2 + sqrt(pow(C2,2)+4*C1))/(2*C1);
            } else {
                if (1-c*Qn/sigf>0) {
                    xi = (sigf*gamma/Qt)*sqrt(1-c*Qn/sigf);
                } else {
                    OOFEM_ERROR("Inconsistent cohesive model specification, 1-c*Qn/sigf =  %e", 1 - c*Qn / sigf);
                }
            }

            Qt = xi*Qt;
            Qn = 0.5*((1+xi)-(1-xi)*sgn(Qn))*Qn;
            Qn_M = 0.5*(Qn + fabs(Qn));
            
            double beta = pow(Qt,2)*Kstiff.at(3,3)/(pow(Qn_M,2)*Kstiff.at(1,1) + pow(Qt,2)*Kstiff.at(3,3));

            double G_beta = beta*(this->GIIc - pow(gamma*sigf,2)/this->ks0) + (1-beta)*(this->GIc - pow(sigf,2)/this->kn0); // assuming linear interpolation between mode I and II
            
            
            double eta = (pow(Qn_M,2) + pow(Qt,2)*Kstiff.at(3,3)/Kstiff.at(1,1))/(G_beta*sigf);
            
            dAlpha = (1/xi-1)*sigf/(2*Kstiff.at(3,3))*eta;

            if ( oldDamage + dAlpha > 1 ) {
                dAlpha = 1-oldDamage;
            }
            
            double Qt_trial = QtrialShear.computeNorm();

            double Qt1,Qt2;

            if ( Qt_trial > 0 ) {
                Qt1 = Qt*QtrialShear.at(1)/Qt_trial;
                Qt2 = Qt*QtrialShear.at(2)/Qt_trial;
            } else {
                Qt1 = 0.0;
                Qt2 = 0.0;
            }


            FloatArray Mstar(3), M(3);

            Mstar.at(1) = 2*Qt1*this->kn0/Kstiff.at(1,1)/(sigf);    // Qt = sqrt(Qt1^2 + Qt2^2)
            Mstar.at(2) = 2*Qt2*this->kn0/Kstiff.at(2,2)/(sigf);
            Mstar.at(3) = 2*Qn_M/sigf;

            M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);
            M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
            M.at(3) = 2*(1-c)/(sigf)*Qn_M + c;

            FloatMatrix Smat(4,4), Smati(4,4);
            Smat.zero();    // S_mat=0.d0

            Smat.at(1,1) = 1.0 + dAlpha*(1/eta)*2*this->kn0/(sigf);     // S_mat(1:2,1:2) = eye3(1:2,1:2)+ dalpha*S*dMdJtn_e(1:2,1:2)
            Smat.at(2,2) = 1.0 + dAlpha*(1/eta)*2*this->kn0/(sigf);
            
            if ( Qn_M > 0 ) {
                Smat.at(3,3) = 1.0 + dAlpha*(1/eta)*2*Kstiff.at(3,3)/(sigf);
            } else {
                Smat.at(3,3) = 1.0;
            }

            Smat.at(1,4) = (1/eta)*Mstar.at(1);         // S_mat(1:2,3) = S*M(1:2)
            Smat.at(2,4) = (1/eta)*Mstar.at(2);
            Smat.at(3,4) = (1/eta)*Mstar.at(3);
                    
            Smat.at(4,1) = M.at(1)*Kstiff.at(1,1);      // S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
            Smat.at(4,2) = M.at(2)*Kstiff.at(2,2);
            Smat.at(4,3) = M.at(3)*Kstiff.at(3,3);
            //Smat.printYourself();
            Smati.beInverseOf(Smat);

            Qtemp.at(1) = Qt1;
            Qtemp.at(2) = Qt2;
            Qtemp.at(3) = Qn;

            IntArray Indx(3);
            Indx.at(1) = 1;
            Indx.at(2) = 2;
            Indx.at(3) = 3;

            FloatMatrix Iep(3,3);
            Iep.beSubMatrixOf(Smati,Indx,Indx);

            status->letTempIepBe(Iep);

            FloatArray alpha_v(3);
            alpha_v.at(1) = Smati.at(4,1);              // alpha_v(1:2) = S_mati(3,1:2)
            alpha_v.at(2) = Smati.at(4,2);
            alpha_v.at(3) = Smati.at(4,3);
            status->letTempAlphavBe(alpha_v);

            dJ = status->giveTempJump();
            status->letTempEffectiveMandelTractionBe(Qtemp);

            Qtemp.times(1-oldDamage-dAlpha);



#if 0
            if (dJ.at(3)>=0) {
                status->letTempEffectiveMandelTractionBe(Qtemp);
                Qtemp.times(1-oldDamage-dAlpha);
            } else {
                if (oldDamage + dAlpha<1) {
                    Qtemp.at(3) = (1-oldDamage)/(1-oldDamage + dAlpha)*Qtemp.at(3);
                    status->letTempEffectiveMandelTractionBe(Qtemp);
                    Qtemp.times(1-oldDamage-dAlpha);
                } else {
                    status->letTempEffectiveMandelTractionBe(Qtemp);        // SHOULD NEVER BE USED
                    Qtemp.times(1-oldDamage-dAlpha);
                    Qtemp.at(3) = (1-oldDamage)*(Qold.at(3) + Kstiff.at(3,3)*dJElastic.at(3)); // dJElastic.at(3) is always equal to dJ.at(3) if dJ.at(3)<0
                }
            }
#endif
        }
    } else {
        dAlpha = 1.0 - oldDamage;
        dJ = status->giveTempJump();
        status->letTempEffectiveMandelTractionBe(Qtemp);    // SHOULD NEVER BE USED!!
        //if (dJ.at(3)<0) {
        //    Qtemp.at(3) = kn0*dJ.at(3);
        //}
    }

    Qtemp.add(Qtemp_comp);

    answer.beTProductOf(Finv,Qtemp);                    // t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
//    answer.times(1-oldDamage-dAlpha);                 // t1_s = (1-al)*t_1_hat

        
    status->letTempDamageBe(oldDamage + dAlpha);
//    status->letTempEffectiveMandelTractionBe(Qtemp);  // NEW!

    status->letTempJumpBe(d);
    status->letTempFirstPKTractionBe(answer);
    status->letTempFBe(F);
}


void

IntMatBilinearCZJansson :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{

    answer.resize(3,3);
    answer.zero();
    //this->give3dStiffnessMatrix_dTdj_num(answer, rMode, gp, tStep);
    //OOFEM_WARNING("numerical tangent");
    //answer.printYourself();

    IntMatBilinearCZJanssonStatus *status = static_cast< IntMatBilinearCZJanssonStatus * >( this->giveStatus(gp) );

    if (status->giveOldDamageDev()) {
        answer = status->giveOlddTdJ();
        //answer.printYourself();
        status->letOldDamageDevBe(false);
    } else {
    
        double damage = status->giveTempDamage();
        const FloatMatrix &Finv = status->giveTempInverseDefGrad();
        FloatMatrix help;
        FloatMatrix Kstiff(3,3);
        const FloatArray &J = status->giveTempJump();
        answer.zero();

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
                Kstiff.at(3,3) = this->knc;

                help.beProductOf(Kstiff, Finv);
                answer.beTProductOf(Finv, help);
                //answer.printYourself();
            }
        } else {
            if ( status->giveTempDamage() - status->giveDamage()==0.0 ) {
            
                if ( J.at(3) < 0 ) {
                    Kstiff.at(3,3) = Kstiff.at(3,3) + (this->knc)/(1-damage);
                }

                help.beProductOf(Kstiff, Finv);
                answer.beTProductOf(Finv, help);
                answer.times(1-damage);  // Ea=(1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Keye3,Fci))
            } else {
                FloatMatrix Iep = status->giveTempIep();

                //Iep.printYourself();

                //if (J.at(3)<0) {
                //    Kstiff.at(1,1) = (1-damage)*Kstiff.at(1,1);
                //    Kstiff.at(2,2) = (1-damage)*Kstiff.at(2,2);
                //    Kstiff.at(3,3) = Kstiff.at(3,3);
                //} else {
                //    Kstiff.times((1-damage));
                //}
                Kstiff.times(1-damage);

                answer.beProductOf(Kstiff, Iep);
                //answer.rotatedWith(Rot);                  // Ea_h = MATMUL(TRANSPOSE(Rot),MATMUL(Keye3,Iep))
                                                            // Ea_h = MATMUL(Ea_h,Rot)

                FloatArray alpha_v = status->giveTempAlphav();
                //alpha_v.rotatedWith(Rot, 't');            // alpha_v = MATMUL(TRANSPOSE(Rot),alpha_v)

                FloatMatrix t1hatFinvOpen;
                FloatArray temp1, temp2, Qtemp;
                Qtemp = status->giveTempEffectiveMandelTraction();

                temp1.beTProductOf(Finv,Qtemp);             // CALL gmopen33(MATMUL(TRANSPOSE(Fci),Q),MATMUL(alpha_v,Fci),t1halFci_o)
                temp2.beTProductOf(Finv,alpha_v);

                t1hatFinvOpen.beDyadicProductOf(temp1,temp2);

                if ( J.at(3) < 0 ) {
                    answer.at(3,3) += this->knc;
                }

                help.beProductOf(answer,Finv);              // Ea = (1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Ea_h,Fci)) -&
                answer.beTProductOf(Finv,help);             // t1halFci_o
    //            answer.times(1-damage);
                answer.subtract(t1hatFinvOpen);
            }
        }
    }
    status->letTempdTdJBe(answer);

    //Finv.printYourself();
    //Kstiff.printYourself();
    //printf("analytical tangent \n");
    //answer.printYourself();

}


int
IntMatBilinearCZJansson :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    IntMatBilinearCZJanssonStatus *status = static_cast< IntMatBilinearCZJanssonStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {     
        answer.resize(1);
        answer.at(1) = status->giveTempDamage();
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, atTime);
    }

}






const double tolerance = 1.0e-12; // small number
IRResultType
IntMatBilinearCZJansson :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn0, _IFT_IntMatBilinearCZJansson_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_IntMatBilinearCZJansson_knc);

    this->ks0 = kn0;                        // Defaults to kn0
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_IntMatBilinearCZJansson_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatBilinearCZJansson_g1c);

    this->GIIc = GIc;                       //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_IntMatBilinearCZJansson_g2c);

    IR_GIVE_FIELD(ir, sigf, _IFT_IntMatBilinearCZJansson_sigf);

    IR_GIVE_FIELD(ir, mu, _IFT_IntMatBilinearCZJansson_mu);

    IR_GIVE_FIELD(ir, gamma, _IFT_IntMatBilinearCZJansson_gamma);

    // check validity of the material paramters
    if ( this->kn0 < 0.0 ) {
        OOFEM_WARNING("stiffness kn0 is negative (%.2e)", this->kn0);
        return IRRT_BAD_FORMAT;
    } else if ( this->ks0 < 0.0 ) {
        OOFEM_WARNING("stiffness ks0 is negative (%.2e)", this->ks0);
        return IRRT_BAD_FORMAT;
    } else if ( this->GIc < 0.0 ) {
        OOFEM_WARNING("GIc is negative (%.2e)", this->GIc);
        return IRRT_BAD_FORMAT;
    } else if ( this->GIIc < 0.0 ) {
        OOFEM_WARNING("GIIc is negative (%.2e)", this->GIIc);
        return IRRT_BAD_FORMAT;
    } else if ( this->gamma < 0.0  ) { 
        OOFEM_WARNING("gamma (%.2e) is below zero which is unphysical",  this->gamma);
        return IRRT_BAD_FORMAT;
    }
    return IRRT_OK;
}

int
IntMatBilinearCZJansson :: checkConsistency()
{
    return 1;
}

void
IntMatBilinearCZJansson  :: printYourself()
{
    printf("Paramters for BilinearCZMaterial: \n");

    printf("-Strength paramters \n");
    printf("  sigf  = %e \n", this->sigf);
    printf("  GIc   = %e \n", this->GIc);
    printf("  GIIc  = %e \n", this->GIIc);
    printf("  gamma = sigfs/sigfn = %e \n", this->gamma);
    printf("  mu    = %e \n", this->mu);

    printf("-Stiffness parameters \n");
    printf("  kn0  = %e \n", this->kn0);
    printf("  ks0  = %e \n", this->ks0);
    printf("  knc  = %e \n", this->knc);

}


//IntMatBilinearCZJanssonStatus :: IntMatBilinearCZJanssonStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
IntMatBilinearCZJanssonStatus :: IntMatBilinearCZJanssonStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
{
    oldMaterialJump.resize(3);
    oldMaterialJump.zero();
    tempMaterialJump = oldMaterialJump;

    damage = tempDamage = 0.0;

    QEffective = oldMaterialJump;
    tempQEffective = oldMaterialJump;

        
    tempFInv.resize(3,3);
    tempFInv.beUnitMatrix();

    old_dTdJ.resize(3,3);
    old_dTdJ.zero();

    oldDamageDev = false;


#if 0
    ///@todo Martin: Very bad implementation of intialisation of Rot
    //*************************************************************
    FloatMatrix Gcov;
    FloatArray N;

    Shell7Base *shell = dynamic_cast< Shell7Base* >(gp->giveElement());
    if ( !shell ) {
        OOFEM_ERROR("BilinearCZMaterialJansson :: giveRealStressVector - oh no wrong element type");
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


IntMatBilinearCZJanssonStatus :: ~IntMatBilinearCZJanssonStatus()
{ }


void
    IntMatBilinearCZJanssonStatus :: printOutputAt(FILE *file, TimeStep *tStep)
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
IntMatBilinearCZJanssonStatus :: initTempStatus()
{
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

    tempDamageDev = false;
}

void
IntMatBilinearCZJanssonStatus :: updateYourself(TimeStep *atTime)
{
    StructuralInterfaceMaterialStatus :: updateYourself(atTime);
    damage = tempDamage;
    oldMaterialJump = tempMaterialJump;
    QEffective = tempQEffective;

    old_dTdJ = temp_dTdJ;
    oldDamageDev = tempDamageDev;
}


} // end namespace oofem
