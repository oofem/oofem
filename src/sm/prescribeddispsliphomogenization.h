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

#ifndef prescribeddispsliphomogenization_h
#define prescribeddispsliphomogenization_h

#include "inputrecord.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include "error.h"


///@name Input fields for PrescribedDispSlipHomogenization
//@{
#define _IFT_PrescribedDispSlipHomogenization_dispField "disp"
#define _IFT_PrescribedDispSlipHomogenization_dispGrad "dispgrad"
#define _IFT_PrescribedDispSlipHomogenization_slipField "slip"
#define _IFT_PrescribedDispSlipHomogenization_slipGrad "slipgrad"
#define _IFT_PrescribedDispSlipHomogenization_centercoords "ccoord"
//@}

namespace oofem {
class TimeStep;
class DynamicInputRecord;
class Domain;

/**
 * Class for homogenization of multiple applied gradients and fields.
 * Specifically, applied to boundary conditions in multiscale analysis of reinfoced concrete structure
 * with displacement and reinforcement slip fields at the macroscale.
 * Currently implemented only for 2D problems.
 * 
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedDispSlipHomogenization
{
protected:
    /// Prescribed gradients
    FloatMatrix dispGradient;
    FloatMatrix slipGradient;

    /// Prescribed fields
    FloatArray dispField;
    FloatArray slipField;

    /// Center coordinates
    FloatArray mCenterCoord;

public:
    PrescribedDispSlipHomogenization() { }
    virtual ~PrescribedDispSlipHomogenization() { }

    virtual void initializeFrom(InputRecord &ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    // Methods for field homogenization implemented by respective BCs
    virtual void computeStress(FloatArray &stress, TimeStep *tStep) = 0;
    virtual void computeTransferStress(FloatArray &bStress, TimeStep *tStep) = 0;
    virtual void computeReinfStress(FloatArray &rStress, TimeStep *tStep) = 0;

    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep) = 0;

    virtual void setDispField(const FloatArray &t);
    virtual void setSlipField(const FloatArray &t);
    virtual void setDispGradient(const FloatArray &t);
    virtual void setSlipGradient(const FloatArray &t);

    void giveDispField(FloatArray &oField) const;
    void giveSlipField(FloatArray &oField) const;
    void giveDispGradient(FloatArray &oGradient) const;
    void giveSlipGradient(FloatArray &oGradient) const;

    virtual void setCenterCoordinate(FloatArray &x) { mCenterCoord = x; }
    FloatArray &giveCenterCoordinate() { return mCenterCoord; }

    virtual double domainSize(Domain *d, int set);
};
} // end namespace oofem

#endif // prescribeddispsliphomogenization_h
