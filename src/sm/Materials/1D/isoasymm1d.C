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

#include "isoasymm1d.h"
#include "sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_Material(IsotropicAsymmetric1DMaterial);

IsotropicAsymmetric1DMaterial :: IsotropicAsymmetric1DMaterial(int n, Domain *d) :
    StructuralMaterial(n, d)
{ }

IsotropicAsymmetric1DMaterial :: IsotropicAsymmetric1DMaterial(int n, Domain *d,
                                                                 double _Ec, double _Et) :
    StructuralMaterial(n, d),
    Et(_Et),
    Ec(_Ec),
    a(0.),
    m(15.)
{}


void
IsotropicAsymmetric1DMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, Ec, _IFT_IsotropicAsymmetric1DMaterial_ec);
    IR_GIVE_FIELD(ir, Et, _IFT_IsotropicAsymmetric1DMaterial_et);
    IR_GIVE_FIELD(ir, a, _IFT_IsotropicAsymmetric1DMaterial_talpha);
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_IsotropicAsymmetric1DMaterial_m);

}


void
IsotropicAsymmetric1DMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->Ec, _IFT_IsotropicAsymmetric1DMaterial_ec);
    input.setField(this->Et, _IFT_IsotropicAsymmetric1DMaterial_et);
    input.setField(this->a, _IFT_IsotropicAsymmetric1DMaterial_talpha);
    input.setField(this->m, _IFT_IsotropicAsymmetric1DMaterial_m);

}


void IsotropicAsymmetric1DMaterial :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterial :: saveContext(stream, mode);

    if ( ( mode & CM_Definition ) ) {
        if ( !stream.write(Ec) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(Et) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(a) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(m) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void IsotropicAsymmetric1DMaterial :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterial :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.read(Ec) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(Et) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(a) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(m) ) {
            THROW_CIOERR(CIO_IOERR);
        }
 
    }
}



double
IsotropicAsymmetric1DMaterial :: give(int aProperty, GaussPoint *gp) const
{
    return this->StructuralMaterial :: give(aProperty, gp);
}


bool
IsotropicAsymmetric1DMaterial:: hasMaterialModeCapability(MaterialMode mode) const
{
    if (mode == _1dMat) 
        return true;
    else 
        return false;
}


FloatMatrixF<1,1>
IsotropicAsymmetric1DMaterial :: give1dStressStiffMtrx(MatResponseMode mode,
                                                        GaussPoint *gp,
                                                        TimeStep *tStep) const
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double eps = status->giveTempStrainVector().at(1);
    double E = this->Ec+0.5*(1+tanh(this->m*eps))*(this->Et-this->Ec);

    return {E};
}

FloatArrayF< 1 >
IsotropicAsymmetric1DMaterial::giveRealStressVector_1d(const FloatArrayF< 1 > &reducedE, GaussPoint *gp, TimeStep *tStep) const
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double eps = reducedE.at(1);
    double s=(0.5*this->Et-0.5*this->Ec)*((log(cosh(this->m*eps)))/(this->m))+eps*(0.5*this->Ec+0.5*this->Et);

    status->letTempStressVectorBe({s});
    status->letTempStrainVectorBe(reducedE);

    return {s};
}

} // end namespace oofem
