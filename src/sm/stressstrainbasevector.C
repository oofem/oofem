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

#include "stressstrainbasevector.h"
#include "intarray.h"
#include "floatmatrix.h"
#include "error.h"
#include "datastream.h"
#include "materialmode.h"
#include "../sm/Materials/structuralmaterial.h"
#include "datastream.h"

namespace oofem {
StressStrainBaseVector :: StressStrainBaseVector(MaterialMode m) : FloatArray()
{
    this->resize( StructuralMaterial :: giveSizeOfVoigtSymVector(m) );
    this->zero();
    this->mode = m;
}

StressStrainBaseVector :: StressStrainBaseVector(const FloatArray &src, MaterialMode m) : FloatArray(src)
{
    if ( StructuralMaterial :: giveSizeOfVoigtSymVector(m) != src.giveSize() ) {
        OOFEM_ERROR( "The source has size %d and a new MaterialMode %s has reduced size %d", src.giveSize(), __MaterialModeToString(m), StructuralMaterial :: giveSizeOfVoigtSymVector(m) );
    }

    this->mode = m;
}

StressStrainBaseVector &
StressStrainBaseVector :: operator = ( const StressStrainBaseVector & src )
{
    if ( this != & src ) { // beware of s=s;
        this->values = src.values;
    }

    this->mode = src.mode;
    return * this;
}

void
StressStrainBaseVector :: convertToFullForm(FloatArray &answer) const
{
    IntArray indx;
    int answerSize = 6;

    if ( mode == _3dMat ) {
        answer = * this;
        return;
    }

    answer.resize(answerSize);
    answer.zero();

    StructuralMaterial :: giveVoigtSymVectorMask(indx, ( MaterialMode ) mode);
    answer.assemble(* this, indx);
}

void
StressStrainBaseVector :: convertFromFullForm(const FloatArray &vector, MaterialMode mode)
{
    IntArray indx;

    if ( mode == _3dMat ) {
        if ( this->giveSize() != 6 ) {
            OOFEM_ERROR("full vector size mismatch");
        }

        this->resize(6);
        for ( int i = 1; i <= 6; i++ ) {
            this->at(i) = vector.at(i);
        }
    } else {
        StructuralMaterial :: giveVoigtSymVectorMask(indx, mode);
        this->resize( StructuralMaterial :: giveSizeOfVoigtSymVector(mode) );
        this->zero();

        for ( int i = 1; i <= indx.giveSize(); i++ ) {
            int j;
            if ( ( j = indx.at(i) ) ) {
                this->at(i) = vector.at(j);
            }
        }
    }
}


contextIOResultType
StressStrainBaseVector :: storeYourself(DataStream &stream)
{
    contextIOResultType iores;
    if ( ( iores = FloatArray :: storeYourself(stream) ) != CIO_OK ) {
        return CIO_OK;
    }

    // write material mode
    if ( !stream.write((int)mode) ) {
        return CIO_IOERR;
    }

    return CIO_OK;
}

contextIOResultType
StressStrainBaseVector :: restoreYourself(DataStream &stream)
{
    contextIOResultType iores;
    if ( ( iores = FloatArray :: restoreYourself(stream) ) != CIO_OK ) {
        return iores;
    }

    // read material mode
    int val;
    if ( !stream.read(val) ) {
        return CIO_IOERR;
    }
    mode = (StressStrainMatMode)val;

    return CIO_OK;
}


void
StressStrainBaseVector :: letStressStrainModeBe(const MaterialMode newMode)
{
    this->mode = ( StressStrainMatMode ) newMode;
    this->resize( StructuralMaterial :: giveSizeOfVoigtSymVector(newMode) );
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
    FloatArray fullReceiver, fullAnswer;

    this->giveTranformationMtrx(tt, base, transpose);
    // convert receiver to full mode
    this->convertToFullForm(fullReceiver);
    fullAnswer.beProductOf(tt, fullReceiver);
    // convert back to reduced form
    answer.convertFromFullForm( fullAnswer, this->giveStressStrainMode() );
}

double
StressStrainBaseVector :: computeVolumetricPart() const
{
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        // 1D model
        OOFEM_ERROR("No Split for 1D!");
        return 0.0;
    } else if ( myMode == _PlaneStress ) {
        // plane stress problem
        OOFEM_ERROR("No Split for plane stress!");
        return 0.0;
    } else {
        // 3d, plane strain or axisymmetric problem
        return ( this->at(1) + this->at(2) + this->at(3) ) / 3.0;
    }
}
} // end namespace oofem
