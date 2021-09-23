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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "prescribeddispsliphomogenization.h"
#include "domain.h"
#include "dynamicinputrecord.h"
#include "set.h"
#include "feinterpol.h"
#include "element.h"
#include "mathfem.h"

namespace oofem {

void PrescribedDispSlipHomogenization::initializeFrom(InputRecord &ir)
{
    mCenterCoord.resize( dispGradient.giveNumberOfColumns() );
    mCenterCoord.zero();
    dispField.resize( dispGradient.giveNumberOfColumns() );
    dispField.zero(); // prescribed displacement field results only in rigid body motion, no reason to require this input

    IR_GIVE_OPTIONAL_FIELD(ir, dispField, _IFT_PrescribedDispSlipHomogenization_dispField);
    IR_GIVE_OPTIONAL_FIELD(ir, dispGradient, _IFT_PrescribedDispSlipHomogenization_dispGrad);
    IR_GIVE_OPTIONAL_FIELD(ir, slipField, _IFT_PrescribedDispSlipHomogenization_slipField);
    IR_GIVE_OPTIONAL_FIELD(ir, slipGradient, _IFT_PrescribedDispSlipHomogenization_slipGrad);
    IR_GIVE_OPTIONAL_FIELD(ir, mCenterCoord, _IFT_PrescribedDispSlipHomogenization_centercoords)
}

void PrescribedDispSlipHomogenization::giveInputRecord(DynamicInputRecord &input)
{
    input.setField(dispField, _IFT_PrescribedDispSlipHomogenization_dispField);
    input.setField(slipField, _IFT_PrescribedDispSlipHomogenization_slipField);
    input.setField(dispGradient, _IFT_PrescribedDispSlipHomogenization_dispGrad);
    input.setField(slipGradient, _IFT_PrescribedDispSlipHomogenization_slipGrad);
    input.setField(mCenterCoord, _IFT_PrescribedDispSlipHomogenization_centercoords);
}


void PrescribedDispSlipHomogenization::setDispField( const FloatArray &t )
{
    int n = t.giveSize();
    if ( n == 2 ) { // Then 2D
        this->dispField.resize(2);
        this->dispField.at(1) = t.at(1);
        this->dispField.at(2) = t.at(2);
    } else {
        OOFEM_ERROR("Field is in strange format. Should be 2 or 3.");
    }
}


void PrescribedDispSlipHomogenization::setSlipField( const FloatArray &t )
{
    int n = t.giveSize();
    if ( n == 2 ) { // Then 2D
        this->slipField.resize(2);
        this->slipField.at(1) = t.at(1);
        this->slipField.at(2) = t.at(2);
    } else {
        OOFEM_ERROR("Field is in strange format. Should be 2 or 3.");
    }
}


void PrescribedDispSlipHomogenization::setDispGradient(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 3 ) { // Then 2D
        this->dispGradient.resize(2, 2);
        this->dispGradient.at(1, 1) = t.at(1);
        this->dispGradient.at(2, 2) = t.at(2);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->dispGradient.at(1, 2) = this->dispGradient.at(2, 1) = t.at(3) * 0.5;
    } else {
        OOFEM_ERROR("Tensor is in strange voigt format. Should be 3.");
    }
}


void PrescribedDispSlipHomogenization::setSlipGradient( const FloatArray &t )
{
    int n = t.giveSize();
    if ( n == 4 ) { // Then 2D
        this->slipGradient.resize(2, 2);
        this->slipGradient.at(1, 1) = t.at(1);
        this->slipGradient.at(2, 2) = t.at(2);
        this->slipGradient.at(1, 2) = t.at(3);
        this->slipGradient.at(2,1) = t.at(4);
    } else {
        OOFEM_ERROR("Tensor is in strange voigt format. Should be 4.");
    }

}


void PrescribedDispSlipHomogenization::giveDispField( FloatArray &oField ) const
{
    int size = dispField.giveSize();
    if ( size == 2 ) {
        oField = { dispField.at(1), dispField.at(2) };
    } else {
        OOFEM_ERROR("PrescribedDispSlipHomogenization :: giveDispField not implemented for 3D.\n");
    }
}


void PrescribedDispSlipHomogenization::giveSlipField( FloatArray &oField ) const
{
    int size = slipField.giveSize();
    if ( size == 2 ) {
        oField = { slipField.at(1), slipField.at(2) };
    } else {
        OOFEM_ERROR("PrescribedDispSlipHomogenization :: giveSlipField not implemented for 3D.\n");
    }
}


void PrescribedDispSlipHomogenization::giveDispGradient( FloatArray &oGradient ) const
{
    int numRows = dispGradient.giveNumberOfRows();
    if ( numRows == 2 ) {
        oGradient = { dispGradient.at(1, 1), dispGradient.at(2, 2), dispGradient.at(1, 2), dispGradient.at(2, 1) };
    } else {
        OOFEM_ERROR("PrescribedDispSlipHomogenization :: giveDispGradient not implemented for 3D.\n");
    }
}


void PrescribedDispSlipHomogenization::giveSlipGradient( FloatArray &oGradient ) const
{
    int numRows = slipGradient.giveNumberOfRows();
    if ( numRows == 2 ) {
        oGradient = { slipGradient.at(1, 1), slipGradient.at(2, 2), slipGradient.at(1, 2), slipGradient.at(2, 1) };
    } else {
        OOFEM_ERROR("PrescribedDispSlipHomogenization :: giveSlipGradient not implemented for 3D.\n");
    }
}


double PrescribedDispSlipHomogenization::domainSize(Domain *d, int setNum)
{
    int nsd = d->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = d->giveSet(setNum);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = d->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
}
} /* namespace oofem */
