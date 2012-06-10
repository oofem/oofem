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

#include "stressstrainbasevector.h"
#include "intarray.h"
#include "flotmtrx.h"
#include "error.h"
#include "datastream.h"
#include "materialmode.h"
#include "matresponseform.h"

namespace oofem {
StressStrainBaseVector :: StressStrainBaseVector(MaterialMode m) : FloatArray()
{
    this->resize( this->giveReducedSize(m) );
    this->zero();
    this->mode = m;
}

StressStrainBaseVector :: StressStrainBaseVector(const FloatArray &src, MaterialMode m) : FloatArray(src)
{
    if ( this->giveReducedSize(m) != src.giveSize() ) {
        OOFEM_ERROR("StressStrainBaseVector::StressStrainBaseVector: size mismatch");
    }

    this->mode = m;
}

StressStrainBaseVector &
StressStrainBaseVector :: operator = ( const StressStrainBaseVector & src )
{
    // assignment: cleanup and copy
    double *srcVal;

    if ( this != & src ) { // beware of s=s;
        this->resize(src.size);

        srcVal = src.givePointer();
        for ( int i = 0; i < size; i++ ) {
            this->values [ i ] = srcVal [ i ];
        }
    }

    this->mode = src.mode;
    return * this;
}

void
StressStrainBaseVector :: convertToFullForm(FloatArray &answer) const
{
    IntArray indx;
    int i, j, answerSize = 6;

    if ( mode == _3dMat ) {
        answer = * this;
        return;
    }

    answer.resize(answerSize);
    answer.zero();
    if ( mode == _3dRotContinuum ) {
        OOFEM_ERROR("StressStrainBaseVector::convertToFullForm: Fullform not available for 3dRotContinuum");
    }

    this->giveStressStrainMask(indx, ReducedForm, ( MaterialMode ) mode);
    for ( i = 1; i <= indx.giveSize(); i++ ) {
        if ( ( j = indx.at(i) ) ) {
            answer.at(j) = this->at(i);
        }
    }

    return;
}

void
StressStrainBaseVector :: convertFromFullForm(const FloatArray &vector, MaterialMode mode)
{
    IntArray indx;
    int i, j;

    if ( mode == _3dMat ) {
        if ( size != 6 ) {
            OOFEM_ERROR("convertFromFullForm - full vector size mismatch");
        }

        this->resize(6);
        for ( i = 1; i <= 6; i++ ) {
            this->at(i) = vector.at(i);
        }

        return;
    } else {
        this->giveStressStrainMask(indx, ReducedForm, mode);
        this->resize( giveReducedSize(mode) );
        this->zero();

        for ( i = 1; i <= indx.giveSize(); i++ ) {
            if ( ( j = indx.at(i) ) ) {
                this->at(i) = vector.at(j);
            }
        }

        return;
    }
}


contextIOResultType
StressStrainBaseVector :: storeYourself(DataStream *stream, ContextMode mode)
{
    contextIOResultType iores;
    if ( ( iores = FloatArray :: storeYourself(stream, mode) ) != CIO_OK ) {
        return CIO_OK;
    }

    // write material mode
    if ( !stream->write(& mode, 1) ) {
        return CIO_IOERR;
    }

    return CIO_OK;
}

contextIOResultType
StressStrainBaseVector :: restoreYourself(DataStream *stream, ContextMode mode)
{
    contextIOResultType iores;
    if ( ( iores = FloatArray :: restoreYourself(stream, mode) ) != CIO_OK ) {
        return iores;
    }

    // read material mode
    if ( !stream->read(& mode, 1) ) {
        return CIO_IOERR;
    }

    return CIO_OK;
}


int
StressStrainBaseVector :: giveReducedSize(MaterialMode mode)
//
// returns the size of reduced stress-strain vector
// acording to mode given by gp.
//
{
    switch ( mode ) {
    case _3dMat:
        return 6;

    case _PlaneStress:
        return 3;

    case _PlaneStrain:
        return 4;

    case _1dMat:
        return 1;

    case _3dRotContinuum:
        return 4;

    case _Unknown:
        return 0;

    case _3dMatGrad:
        return 7;

    case _PlaneStressGrad:
        return 4;

    case _PlaneStrainGrad:
        return 5;

    case _1dMatGrad:
        return 2;

    default:
        OOFEM_ERROR2( "StressStrainBaseVector::giveReducedSize : unknown mode (%s)", __MaterialModeToString(mode) );
    }

    return 0;
}


void
StressStrainBaseVector :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                               MaterialMode mmode) const
//
// this function returns mask of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain vector in full or
// reduced StressStrainVector
// acording to stressStrain mode of given gp.
//
//
// mask has size of reduced or full StressStrain Vector and  i-th component
// is index to full or reduced StressStrainVector where corresponding
// stressStrain resides.
//
// Reduced form is sub-vector (of stress or strain components),
// where components corresponding to imposed zero stress (plane stress,...)
// are not included. On the other hand, if zero strain component is imposed
// (Plane strain, ..) this condition must be taken into account in geometrical
// relations, and corresponding component is included in reduced vector.
//
{
    int i;

    if ( form == ReducedForm ) {
        switch ( mmode ) {
        case _3dMat:
            answer.resize(6);
            for ( i = 1; i <= 6; i++ ) {
                answer.at(i) = i;
            }

            break;
        case _PlaneStress:
            answer.resize(3);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 6;
            break;
        case _PlaneStrain:
            answer.resize(4);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(4) = 6;
            break;
        case _1dMat:
            answer.resize(1);
            answer.at(1) = 1;
            break;
        case _3dRotContinuum:
            OOFEM_ERROR("giveStressMask :  No Mask for 3dRotContinuum");
            break;
        default:
            OOFEM_ERROR2( "giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
        }
    } else if ( form == FullForm ) {
        switch ( mmode ) {
        case _3dMat:
            answer.resize(6);
            answer.zero();
            for ( i = 1; i <= 6; i++ ) {
                answer.at(i) = i;
            }

            break;
        case _PlaneStress:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
            break;
        case _PlaneStrain:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(6) = 4;
            break;
        case _1dMat:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            break;
        case _3dRotContinuum:
            OOFEM_ERROR("giveStressMask :  No Mask for 3dRotContinuum");
            break;
        default:
            OOFEM_ERROR2( "giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
        }
    } else {
        OOFEM_ERROR("giveStressStrainMask : unknown form mode");
    }

    return;
}

void
StressStrainBaseVector :: letStressStrainModeBe(const MaterialMode newMode)
{
    this->mode = ( StressStrainMatMode ) newMode;
    this->resize( this->giveReducedSize(newMode) );
    this->zero();
}

void
StressStrainBaseVector :: transformTo(StressStrainBaseVector &answer, const FloatMatrix &base, int transpose) const
//
//
// performs transformation of given vector to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation. These vectors must
// be expressed in the same coordinate system as strainVector
// If transpose == 1 we transpose base matrix before transforming
//

{
    FloatMatrix tt;
    FloatArray fullReceiver(6), fullAnswer(6);

    this->giveTranformationMtrx(tt, base, transpose);
    // convert receiver to full mode
    this->convertToFullForm(fullReceiver);
    fullAnswer.beProductOf(tt, fullReceiver);
    // convert back to reduced form
    answer.convertFromFullForm( fullAnswer, this->giveStressStrainMode() );
}
} // end namespace oofem
