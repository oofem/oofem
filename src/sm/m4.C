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

#include "m4.h"
#include "microplane.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"

namespace oofem {
M4Material :: M4Material(int n, Domain *d) :
    MicroplaneMaterial_Bazant(n, d)
{ }


inline double
M4Material  :: macbra(double x) /* Macauley bracket = positive part of x */
{
    return ( max(x, 0.0) );
}

inline double
M4Material :: FVplus(double ev, double k1, double c13, double c14, double c15, double Ev)
/*positive volumetric boundary */
{
    return ( Ev * k1 * c13 / ( 1 + ( c14 / k1 ) * macbra(ev - c13 * c15 * k1) ) );
}

inline double
M4Material :: FVminus(double ev, double k1, double k3, double k4, double E)
/*negative volumetric boundary */
{
    return ( -E * k1 * k3 * exp( -ev / ( k1 * k4 ) ) );
}


inline double
M4Material :: FDminus(double ed, double k1, double c7, double c8, double c9,
                      double E)
/*negative deviatoric boundary */
{
    double a;
    a = macbra(-ed - c8 * c9 * k1) / ( k1 * c7 );
    return ( -E * k1 * c8 / ( 1 + a * a ) );
}

inline double
M4Material  :: FDplus(double ed, double k1, double c5, double c6, double c7, double c20,
                      double E)
/*positive deviatoric bondary */
{
    double a;
    a = ( macbra(ed - c5 * c6 * k1) / ( k1 * c20 * c7 ) );
    return ( E * k1 * c5 / ( 1 + a * a ) );
}

inline double
M4Material :: FN(double en, double sv, double k1, double c1, double c2, double c3, double c4,
                 double E, double Ev)
/*normal boundary */
{
    return ( E * k1 * c1 * exp( -macbra(en - c1 * c2 * k1) / ( k1 * c3 + macbra( -c4 * ( sv / Ev ) ) ) ) );
}

inline double
M4Material ::  FT(double sn, double ev, double k1, double k2, double c10,
                  double c11, double c12, double Et)
/*shear boundary */
{
    double a, sn0;
    //oprava sn0 fabs<=>macbra
    sn0 = Et * k1 * c11 / ( 1 + c12 * macbra(ev) );
    a = macbra(-sn + sn0);
    return ( Et * k1 * k2 * c10 * a / ( Et * k1 * k2 + c10 * a ) );
}


void
M4Material :: giveCharacteristicMatrix(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    // elastic stiffness matrix
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = E / ( 2. + 2. * nu );
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = E * ( 1. - nu ) / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
    answer.at(1, 2) = answer.at(2, 1) = answer.at(1, 3) = answer.at(3, 1) =
                                                              answer.at(2, 3) = answer.at(3, 2) = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
}


void
M4Material :: giveRealMicroplaneStressVector(FloatArray &answer,
                                             Microplane *mplane,
                                             const FloatArray &strain,
                                             TimeStep *tStep)
{
    M4MaterialStatus *status =  ( M4MaterialStatus * ) this->giveMicroplaneStatus(mplane);
    FloatArray previousStress, previousStrain;
    FloatArray stressIncrement, strainIncrement;
    double EpsN, DEpsN, DEpsL, DEpsM;
    double EpsV, DEpsV;
    double SEV, SVdash, EpsD, DEpsD, SED, SEM, SEL, SD;
    double SNdash, F;
    double CV, CD;


    // size answer
    answer.resize(4);
    answer.zero();

    // ask status for tempVH parameter
    previousStress = status->giveStressVector();
    previousStrain = status->giveStrainVector();
    strainIncrement.beDifferenceOf(strain, status->giveStrainVector() );
    if ( !previousStress.isNotEmpty() ) {
        previousStress.resize(4);
        previousStress.zero();
    }

    EpsV =  strain.at(1);
    EpsN =  strain.at(2);
    DEpsV = strainIncrement.at(1);
    DEpsN = strainIncrement.at(2);
    DEpsL = strainIncrement.at(3);
    DEpsM = strainIncrement.at(4);
    EpsD = EpsN - EpsV;
    DEpsD = DEpsN - DEpsV;

    /* previousStress.at(i)...Vector of history parameters
     * (1)...previous volumetric stress
     * (2)..previous normal stress for  microplane
     * (3)..previous l-shear stress for  microplane
     * (4)..previous m-shear stress for  microplane
     *
     */


    // SVsum=0.0;
    // novy koncept s odtizenim
    CV = EV;
    CD = ED;
    /*  EpsV0 = EpsV - DEpsV;
     * EpsD0 = EpsD - DEpsD;
     * SigV0 = previousStress.at(1);
     * SigD0 = previousStress.at(2)-SigV0;
     *
     * if ((SigV0<0.0)&&(EpsV0<0.0)&&(DEpsV>0.0))
     *  CV = EV*(c11/(c11-EpsV0)+SigV0*EpsV0/(c11*c12*EV));
     * if ((SigD0>0.0)&&(EpsD0>0.0)&&(DEpsD<0.0)) CD = (SigD0/EpsD0);
     * if ((SigD0<0.0)&&(EpsD0<0.0)&&(DEpsD>0.0)) CD = (SigD0/EpsD0);
     * // if (SigD0*DEpsD<0.0) CD = fabs(SigD0/EpsD0);
     */

    SEV = previousStress.at(1) + CV * DEpsV;

    SVdash = min( max( SEV, this->FVminus(EpsV, k1, k3, k4, E) ), this->FVplus(EpsV, k1, c13, c14, c15, EV) );

    SED = previousStress.at(2) - previousStress.at(1) + CD * DEpsD;
    SD = min( max( SED, this->FDminus(EpsD, k1, c7, c8, c9, E) ),
             this->FDplus(EpsD, k1, c5, c6, c7, c20, E) );

    SNdash = SVdash + SD;
    answer.at(2) = min( SNdash, this->FN(EpsN, previousStress.at(1), k1, c1, c2, c3, c4, E, EV) );

    SEM = previousStress.at(4) + ET * DEpsM;
    SEL = previousStress.at(3) + ET * DEpsL;

    F = this->FT(answer.at(2), EpsV, k1, k2, c10, c11, c12, ET);

    if ( SEL > F ) {
        answer.at(3) = F;
    } else if ( SEL < -F ) {
        answer.at(3) = -F;
    } else {
        answer.at(3) = SEL;
    }

    if ( SEM > F ) {
        answer.at(4) = F;
    } else if ( SEM < -F ) {
        answer.at(4) = -F;
    } else {
        answer.at(4) = SEM;
    }

    answer.at(1) = SVdash;

    // uncomment this
    //stressIncrement = answer;
    //stressIncrement.subtract (previousStress);
    //status->letStressIncrementVectorBe (stressIncrement);
    //status->letStrainIncrementVectorBe (strainIncrement);

    // update gp
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
}


IRResultType
M4Material :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // je-li funkce prepsana, je treba zavolat predchazejici
    MicroplaneMaterial_Bazant :: initializeFrom(ir);

    c1 = 6.20e-1;
    c2 = 2.76;
    // c3 = 4.0;
    IR_GIVE_FIELD(ir, c3, IFT_M4Material_c3, "c3"); // Macro
    c4 = 70.;
    IR_GIVE_OPTIONAL_FIELD(ir, c4, IFT_M4Material_c4, "c4"); // Macro
    c5 = 2.50;
    c6 = 1.30;
    c7 = 50.;
    c8 = 8.0;
    c9 = 1.30;
    c10 = 0.73;
    c11 = 0.2;
    c12 = 7000.;
    c13 = 0.23;
    c14 = 0.8;
    c15 = 1.0;
    c16 = 0.02;
    c17 = 0.01;
    c18 = 1.0;
    c19 = 0.4;
    //c20 = 14.0e-2;
    IR_GIVE_FIELD(ir, c20, IFT_M4Material_c20, "c20"); // Macro

    IR_GIVE_FIELD(ir, k1, IFT_M4Material_k1, "k1"); // Macro
    IR_GIVE_FIELD(ir, k2, IFT_M4Material_k2, "k2"); // Macro
    IR_GIVE_FIELD(ir, k3, IFT_M4Material_k3, "k3"); // Macro
    IR_GIVE_FIELD(ir, k4, IFT_M4Material_k4, "k4"); // Macro
    //k5   = this->readDouble (initString,"k5");
    // k6   = this->readDouble (initString,"k6");
    IR_GIVE_FIELD(ir, E, IFT_M4Material_e, "e"); // Macro
    IR_GIVE_FIELD(ir, nu, IFT_M4Material_n, "n"); // Macro
    IR_GIVE_FIELD(ir, talpha, IFT_M4Material_talpha, "talpha"); // Macro
    //mu   = this->readDouble (initString,"mu");
    mu = 1.0;
    EV = E / ( 1 - 2 * nu );
    ED = 5 * E / ( 2 + 3 * mu ) / ( 1 + nu );
    ET = mu * ED;

    return IRRT_OK;
}

int
M4Material :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode ==  _3dMat ) {
        return 1;
    }

    return 0;
}


void
M4Material :: updateVolumetricStressTo(Microplane *mPlane, double sigv)
{
    //FloatArray stressIncrement;
    FloatArray stress;
    M4MaterialStatus *status =  ( M4MaterialStatus * ) this->giveStatus(mPlane);
    //stressIncrement = status->giveStressIncrementVector();
    //stressIncrement.at(1) = sigv - status->giveStressVector().at(1);
    //status->letStressIncrementVectorBe (stressIncrement);
    stress = status->giveTempStressVector();
    stress.at(1) = sigv;
    status->letTempStressVectorBe(stress);
}

int
M4Material :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
{
    if ( mode == _3dMicroplane ) {
        return 4;
    }

    return MicroplaneMaterial_Bazant :: giveSizeOfReducedStressStrainVector(mode);
}

void
M4Material :: giveThermalDilatationVector(FloatArray &answer,
                                          GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(6);
    answer.zero();
    answer.at(1) = talpha;
    answer.at(2) = talpha;
    answer.at(3) = talpha;
}


////////////////////////////////////////////////////////////////////////////

M4MaterialStatus :: M4MaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g)
{ }


M4MaterialStatus :: ~M4MaterialStatus()
{ }

void
M4MaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void
M4MaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
}

contextIOResultType
M4MaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    return StructuralMaterialStatus :: saveContext(stream, mode, obj);
}

contextIOResultType
M4MaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    return StructuralMaterialStatus :: restoreContext(stream, mode, obj);
}
} // end namespace oofem
