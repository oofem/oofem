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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "m4.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(M4Material);

M4Material :: M4Material(int n, Domain *d) :
    MicroplaneMaterial_Bazant(n, d)
{ }

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
M4Material :: FDplus(double ed, double k1, double c5, double c6, double c7, double c20,
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
M4Material :: FT(double sn, double ev, double k1, double k2, double c10,
                 double c11, double c12, double Et)
/*shear boundary */
{
    double a, sn0;
    //oprava sn0 fabs<=>macbra
    sn0 = Et * k1 * c11 / ( 1 + c12 * macbra(ev) );
    a = macbra(-sn + sn0);
    return ( Et * k1 * k2 * c10 * a / ( Et * k1 * k2 + c10 * a ) );
}


MicroplaneState
M4Material :: giveRealMicroplaneStressVector(GaussPoint *gp, int mnumber,
                                             const MicroplaneState &strain,
                                             TimeStep *tStep) const
{
    M4MaterialStatus *status = static_cast< M4MaterialStatus * >( this->giveStatus(gp) );
    double SEV, SVdash, SED, SEM, SEL, SD;
    double SNdash, F;
    double CV, CD;

    // ask status for tempVH parameter
    auto &prevStrain = status->giveMicroplaneStrain(mnumber);
    auto &previousStress = status->giveMicroplaneStress(mnumber);

    double EpsV =  strain.v;
    double EpsN =  strain.n;
    double DEpsV = strain.v - prevStrain.v;
    double DEpsN = strain.n - prevStrain.n;
    double DEpsL = strain.l - prevStrain.l;
    double DEpsM = strain.m - prevStrain.m;
    double EpsD = EpsN - EpsV;
    double DEpsD = DEpsN - DEpsV;

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

    MicroplaneState answer;

    SEV = previousStress.v + CV * DEpsV;

    SVdash = min( max( SEV, this->FVminus(EpsV, k1, k3, k4, E) ), this->FVplus(EpsV, k1, c13, c14, c15, EV) );

    SED = previousStress.n - previousStress.v + CD * DEpsD;
    SD = min( max( SED, this->FDminus(EpsD, k1, c7, c8, c9, E) ),
              this->FDplus(EpsD, k1, c5, c6, c7, c20, E) );

    SNdash = SVdash + SD;
    answer.n = min( SNdash, this->FN(EpsN, previousStress.v, k1, c1, c2, c3, c4, E, EV) );

    SEM = previousStress.m + ET * DEpsM;
    SEL = previousStress.l + ET * DEpsL;

    F = this->FT(answer.n, EpsV, k1, k2, c10, c11, c12, ET);

    if ( SEL > F ) {
        answer.l = F;
    } else if ( SEL < -F ) {
        answer.l = -F;
    } else {
        answer.l = SEL;
    }

    if ( SEM > F ) {
        answer.m = F;
    } else if ( SEM < -F ) {
        answer.m = -F;
    } else {
        answer.m = SEM;
    }

    answer.v = SVdash;

    // uncomment this
    //stressIncrement = answer;
    //stressIncrement.subtract (previousStress);

    // update gp
    status->letTempMicroplaneStrainBe(mnumber, strain);
    status->letTempMicroplaneStressBe(mnumber, answer);
    return answer;
}


void
M4Material :: initializeFrom(InputRecord &ir)
{
    MicroplaneMaterial_Bazant :: initializeFrom(ir);

    c1 = 6.20e-1;
    c2 = 2.76;
    // c3 = 4.0;
    IR_GIVE_FIELD(ir, c3, _IFT_M4Material_c3);
    c4 = 70.;
    IR_GIVE_OPTIONAL_FIELD(ir, c4, _IFT_M4Material_c4);
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
    IR_GIVE_FIELD(ir, c20, _IFT_M4Material_c20);

    IR_GIVE_FIELD(ir, k1, _IFT_M4Material_k1);
    IR_GIVE_FIELD(ir, k2, _IFT_M4Material_k2);
    IR_GIVE_FIELD(ir, k3, _IFT_M4Material_k3);
    IR_GIVE_FIELD(ir, k4, _IFT_M4Material_k4);
    IR_GIVE_FIELD(ir, talpha, _IFT_M4Material_talpha);
    mu = 1.0;
    EV = E / ( 1 - 2 * nu );
    ED = 5 * E / ( 2 + 3 * mu ) / ( 1 + nu );
    ET = mu * ED;
}


void
M4Material :: updateVolumetricStressTo(GaussPoint *gp, int mnumber, double sigv) const
{
    //FloatArray stressIncrement;
    M4MaterialStatus *status = static_cast< M4MaterialStatus * >( this->giveStatus(gp) );
    //stressIncrement = status->giveStressIncrementVector();
    //stressIncrement.at(1) = sigv - status->giveStressVector().at(1);
    //status->letStressIncrementVectorBe (stressIncrement);
    auto stress = status->giveTempMicroplaneStress(mnumber);
    stress.v = sigv;
    status->letTempMicroplaneStressBe(mnumber, stress);
}


FloatArrayF<6>
M4Material :: giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    return {
        talpha,
        talpha,
        talpha,
        0.,
        0.,
        0.,
    };
}


////////////////////////////////////////////////////////////////////////////

M4MaterialStatus :: M4MaterialStatus(GaussPoint *g, int nplanes) :
    StructuralMaterialStatus(g),
    microplaneStrain(nplanes), tempMicroplaneStrain(nplanes),
    microplaneStress(nplanes), tempMicroplaneStress(nplanes)
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
    microplaneStrain = tempMicroplaneStrain;
    microplaneStress = tempMicroplaneStress;
}

void
M4MaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: saveContext(stream, mode);
    /// @todo, save microplanestress etc.
}

void
M4MaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);
    /// @todo, save microplanestress etc.
}
} // end namespace oofem
