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
#include "intmatbilinczfagerstrom.h"
#include "dynamicinputrecord.h"


namespace oofem {

REGISTER_Material( IntMatBilinearCZFagerstrom );


IntMatBilinearCZFagerstrom :: IntMatBilinearCZFagerstrom(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }


FloatArrayF<3>
IntMatBilinearCZFagerstrom :: giveFirstPKTraction_3d(const FloatArrayF<3> &d, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const
{
    ///@todo - need to rearrange traction and stiffness matrix so the first component is normal
    // returns vector in 3d stress space of receiver according to
    // previous level of stress and current
    // strain increment, the only way, how to correctly update gp records

    IntMatBilinearCZFagerstromStatus *status = static_cast< IntMatBilinearCZFagerstromStatus * >( this->giveStatus(gp) );

    auto Finv = inv(F);
    status->letTempInverseDefGradBe(Finv);

    auto dJ = dot(Finv, d);
    status->letTempMaterialJumpBe(dJ);

    double oldDamage = status->giveDamage();
    double dAlpha = 0.0;
    FloatArrayF<3> Qold, Qtemp, Qtemp_comp;

    if ( dJ.at(3) < 0 ) {
        Qtemp_comp.at(3) = this->knc*dJ.at(3);
    } 


    // SUBROUTINE stress_damage_XFEM_direct_Mandel(dJ,N,Qold,old_alpha,Fci,Q,Ea,new_alpha,adtim,dalpha_new,dalpha_old,diss,sig_f,fall);
    if ( oldDamage < 0.99 ) {

        auto Kstiff = diag<3>({this->ks0, this->ks0, this->kn0});

        dJ -= status->giveOldMaterialJump(); ///@todo Martin: check with Mikael/Jim
        auto Qtrial = status->giveEffectiveMandelTraction() - dot(Kstiff, dJ); 

        double Qn = Qtrial.at(3);
        FloatArrayF<2> QtrialShear = {Qtrial.at(1), Qtrial.at(2)};

        double Qt = norm(QtrialShear);
        double sigf = this->sigf;
        double gamma = this->gamma;
        double mu = this->mu;

        double S = (this->GIc - pow(sigf,2)/(2*this->kn0))/sigf;
        double gammaGf = gamma*(this->GIIc - pow(gamma*sigf,2)/(2*this->ks0))/(this->GIc - pow(sigf,2)/(2*this->kn0));

        double Qn_M = 0.5*(Qn + fabs(Qn));
        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + (sigf/gamma)*(gamma-mu)*pow((Qn_M/sigf),2) - 1/gamma*(gamma*sigf-mu*Qn); // Added support for parameter mu

        Qold = status->giveEffectiveMandelTraction();  ///@todo Martin: check with Mikael/Jim
        Qtemp = Qtrial;

        if ( loadFun/sigf < 1e-7 ) {
            dAlpha = 0.0;  // new_alpha=old_alpha
            status->letTempEffectiveMandelTractionBe(Qtemp);
            Qtemp *= 1 - oldDamage;
        } else {
            // dalpha = datr
            double Qt1 = Qtemp.at(1);
            double Qt2 = Qtemp.at(2);

            FloatArrayF<3> M = {                        // M = (/2*Q_t/(sig_f*gamma**2), 2*Q_n/sig_f, 0.d0/)
                2*Qt1/(pow(gamma,2)*sigf),    // Qt = sqrt(Qt1^2 + Qt2^2)
                2*Qt2/(pow(gamma,2)*sigf),
                2*Qn_M/sigf,
            };

            auto dJElastic = dJ - S * dAlpha * M;                         // dJtn_e(1:2) = dJtn_v(1:2) - S*dalpha*M(1:2)

            FloatMatrixF<4,4> Smati;

            const double errorTol = 0.0001;
            for( int iter = 1; fabs(loadFun)/sigf > errorTol; iter++) {
                //printf("loadfun = %e \n",loadFun);
                if (iter>40) {
                    OOFEM_ERROR("BilinearCZMaterialFagerstrom :: giveRealStressVector - no convergence in constitutive driver");
                }
                status->letTempDamageDevBe(true);

                FloatArrayF<4> R;
                R.at(1) = dJElastic.at(1) - (dJ.at(1) - gammaGf*S*dAlpha*M.at(1));              // R(1:2) = dJtn_e(1:2) - (dJtn_v(1:2) - S*dalpha*M(1:2))
                R.at(2) = dJElastic.at(2) - (dJ.at(2) - gammaGf*S*dAlpha*M.at(2));
                R.at(3) = dJElastic.at(3) - (dJ.at(3) - S*dAlpha*M.at(3));
                R.at(4) = loadFun;  // R(3) = F/sig_f

                FloatMatrixF<4,4> Smat;    // S_mat=0.d0
                Smat.at(1,1) = 1.0 + gammaGf*dAlpha*S*2*Kstiff.at(1,1)/(pow(gamma,2)*sigf);     // S_mat(1:2,1:2) = eye3(1:2,1:2)+ dalpha*S*dMdJtn_e(1:2,1:2)
                Smat.at(2,2) = 1.0 + gammaGf*dAlpha*S*2*Kstiff.at(2,2)/(pow(gamma,2)*sigf);
                //dJElastic.printYourself();
                if (Qn_M>0) {
                    Smat.at(3,3) = 1.0 + dAlpha*S*2*Kstiff.at(3,3)/(sigf);
                } else {
                    Smat.at(3,3) = 1.0;
                }

                Smat.at(1,4) = gammaGf*S*M.at(1);       // S_mat(1:2,3) = S*M(1:2)
                Smat.at(2,4) = gammaGf*S*M.at(2);
                Smat.at(3,4) = S*M.at(3);

                Smat.at(4,1) = M.at(1)*Kstiff.at(1,1);      // S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
                Smat.at(4,2) = M.at(2)*Kstiff.at(2,2);
                //Smat.at(4,3) = M.at(3)*Kstiff.at(3,3);
                Smat.at(4,3) = (2*(gamma-mu)/(gamma*sigf)*Qn_M + mu/gamma)*Kstiff.at(3,3);      //Martin: included parameter mu

                //FloatArray inc;
                //bool transpose = false;
                //Smat.SolveforRhs(R, inc, transpose);

                Smati = inv(Smat);

                auto inc = dot(Smati, R);

                dJElastic.at(1) -= inc.at(1);
                dJElastic.at(2) -= inc.at(2);
                dJElastic.at(3) -= inc.at(3);
                dAlpha = dAlpha - inc.at(4);

                Qtemp = Qold + dot(Kstiff, dJElastic);

                double Qt1 = Qtemp.at(1);
                double Qt2 = Qtemp.at(2);
                Qt = sqrt(Qt1*Qt1 + Qt2*Qt2);
                Qn = Qtemp.at(3);                       // Martin: included parameter mu
                Qn_M = 0.5*(Qtemp.at(3)+fabs(Qtemp.at(3)));

                M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);    // Qt = sqrt(Qt1^2 + Qt2^2)
                M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
                M.at(3) = 2*Qn_M/sigf;

                //loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;
                loadFun = sigf*pow(Qt/(gamma*sigf),2) + (sigf/gamma)*(gamma-mu)*pow((Qn_M/sigf),2) - 1/gamma*(gamma*sigf-mu*Qn); //Martin: included parameter mu
            }
            //printf("converged, loadfun = %e, oldDamage = %e \n",loadFun, oldDamage);

            if (oldDamage + dAlpha > 1) {
              dAlpha = 1. - oldDamage;
            }

            auto Iep = Smati({0,1,2},{0,1,2});
            status->letTempIepBe(Iep);

            FloatArrayF<3> alpha_v = {
                Smati.at(4,1),          // alpha_v(1:2) = S_mati(3,1:2)
                Smati.at(4,2),
                Smati.at(4,3),
            };
            status->letTempAlphavBe(alpha_v);
            status->letTempEffectiveMandelTractionBe(Qtemp);
            Qtemp *= 1 - oldDamage - dAlpha;

        }

        //Qtemp.rotatedWith(Rot,'t');               // Q=Qe
        //status->letTempRotationMatrix(Rot);
    } else {
        dAlpha = 1.0 - oldDamage;
        dJ = status->giveTempJump();
        status->letTempEffectiveMandelTractionBe(Qtemp);    // SHOULD NEVER BE USED!!
        //if (dJ.at(3)<0) {
        //  Qtemp.at(3) = kn0*dJ.at(3);
        //}
    }

    Qtemp += Qtemp_comp;

    auto answer = Tdot(Finv, Qtemp);                // t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
    //* (1-oldDamage-dAlpha);             // t1_s = (1-al)*t_1_hat

    status->letTempDamageBe(oldDamage + dAlpha);
    //status->letTempEffectiveMandelTractionBe(Qtemp);  // NEW!

    status->letTempJumpBe(d);
    status->letTempFirstPKTractionBe(answer);
    status->letTempFBe(F);

    return answer;
}


FloatMatrixF<3,3>
IntMatBilinearCZFagerstrom :: give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    IntMatBilinearCZFagerstromStatus *status = static_cast< IntMatBilinearCZFagerstromStatus * >( this->giveStatus(gp) );

    FloatMatrixF<3,3> answer;
    if ( status->giveOldDamageDev() ) {
        answer = status->giveOlddTdJ();
        //answer.printYourself();
        status->letOldDamageDevBe(false);
    } else {

        double damage = status->giveTempDamage();
        const auto &Finv = status->giveTempInverseDefGrad();
        const auto &J = status->giveTempJump();

        FloatMatrixF<3,3> Kstiff = diag<3>({this->ks0, this->ks0, this->kn0});
        //FloatMatrix Rot = status->giveTempRotationMatrix();
        //Kstiff = rotate(Kstiff, Rot);

        if ( damage >= 1.0 ) {
            if ( J.at(3) < 0 ) {
                Kstiff = diag<3>({0.0, 0.0, this->knc});
                answer = rotate(Kstiff, Finv);
                //printf("fully damaged");
                //answer.printYourself();
            }
        } else {
            if ( status->giveTempDamage() - status->giveDamage()==0.0 ) {

                if ( J.at(3) < 0 ) {
                    Kstiff.at(3,3) = Kstiff.at(3,3) + (this->knc)/(1-damage);
                }

                // Ea=(1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Keye3,Fci))
                answer = (1-damage) * rotate(Kstiff, Finv);
                //printf("elastic step");
            } else {

                const auto &Iep = status->giveTempIep();
                //Iep.ourselfourself();

                //if (J.at(3)<0) {
                //Kstiff.at(1,1) = (1-damage)*Kstiff.at(1,1);
                //Kstiff.at(2,2) = (1-damage)*Kstiff.at(2,2);
                //Kstiff.at(3,3) = Kstiff.at(3,3);
                //} else {
                //Kstiff.times((1-damage));
                //}

                answer = (1 - damage) * dot(Kstiff, Iep);
                //answer.rotatedWith(Rot);
                // Ea_h = MATMUL(TRANSPOSE(Rot),MATMUL(Keye3,Iep))
                // Ea_h = MATMUL(Ea_h,Rot)

                const auto &alpha_v = status->giveTempAlphav();
                //alpha_v = Tdot(Rot, alpha_v);    // alpha_v = MATMUL(TRANSPOSE(Rot),alpha_v)

                const auto &Qtemp = status->giveTempEffectiveMandelTraction();

                auto temp1 = Tdot(Finv, Qtemp);     // CALL gmopen33(MATMUL(TRANSPOSE(Fci),Q),MATMUL(alpha_v,Fci),t1halFci_o)
                auto temp2 = Tdot(Finv, alpha_v);

                auto t1hatFinvOpen = dyad(temp1, temp2);

                if ( J.at(3) < 0 ) {
                    answer.at(3,3) = answer.at(3,3) + this->knc;
                }

                //printf("plastic step");
                // Ea = (1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Ea_h,Fci)) - t1halFci_o
                answer = rotate(answer, Finv) - t1hatFinvOpen; 
            }
        }
    }
    status->letTempdTdJBe(answer);

    return answer;
}


//const double tolerance = 1.0e-12; // small number
void
IntMatBilinearCZFagerstrom :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, kn0, _IFT_IntMatBilinearCZFagerstrom_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_IntMatBilinearCZFagerstrom_knc);

    this->ks0 = kn0;                        // Defaults to kn0
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_IntMatBilinearCZFagerstrom_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatBilinearCZFagerstrom_g1c);

    this->GIIc = GIc;                       //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_IntMatBilinearCZFagerstrom_g2c);

    IR_GIVE_FIELD(ir, sigf, _IFT_IntMatBilinearCZFagerstrom_sigf);

    IR_GIVE_FIELD(ir, mu, _IFT_IntMatBilinearCZFagerstrom_mu);

    IR_GIVE_FIELD(ir, gamma, _IFT_IntMatBilinearCZFagerstrom_gamma);

    StructuralInterfaceMaterial ::initializeFrom(ir);

    // check validity of the material paramters
    if ( this->kn0 < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZFagerstrom_kn, "must be positive");
    } else if ( this->ks0 < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZFagerstrom_ks, "must be positive");
    } else if ( this->GIc < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZFagerstrom_g2c, "must be positive");
    } else if ( this->GIIc < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZFagerstrom_g2c, "must be positive");
    } else if ( this->gamma < 0.0  ) { 
        throw ValueInputException(ir, _IFT_IntMatBilinearCZFagerstrom_gamma, "must be positive");
    }
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
    return 1;
}

void
IntMatBilinearCZFagerstrom  :: printYourself()
{
    printf("Parameters for BilinearCZMaterial: \n");

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


IntMatBilinearCZFagerstromStatus :: IntMatBilinearCZFagerstromStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{
    tempFInv = eye<3>();
    Iep = tempFInv;
}


void
IntMatBilinearCZFagerstromStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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
    StructuralInterfaceMaterialStatus :: initTempStatus();

    tempMaterialJump = oldMaterialJump;
    tempDamage = damage;
    tempQEffective = QEffective;

    tempFInv = eye<3>();
    Iep = tempFInv;
    alphav = oldMaterialJump;

    tempDamageDev = false;
}

void
IntMatBilinearCZFagerstromStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    damage = tempDamage;
    oldMaterialJump = tempMaterialJump;
    QEffective = tempQEffective;
    old_dTdJ = temp_dTdJ;
    oldDamageDev = tempDamageDev;
}


} // end namespace oofem
