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


#ifndef PRESCRIBEDGRADIENTBC_H_
#define PRESCRIBEDGRADIENTBC_H_

#include "activebc.h"
#include "valuemodetype.h"
#include "floatmatrix.h"
#include "floatarray.h"

#define _IFT_PrescribedGradientBC_centercoords "ccoord"
#define _IFT_PrescribedGradientBC_gradient "gradient"

namespace oofem {
/**
 * Base class for boundary conditions that prescribe a displacement gradient
 * (in a strong or weak sense) on the boundary of a RVE.
 *
 * Under development.
 *
 * @author Erik Svenning
 * @date Mar 5, 2014
 */

class OOFEM_EXPORT PrescribedGradientBC : public ActiveBoundaryCondition
{
protected:
    /// Prescribed gradient @f$ d_{ij} @f$
    FloatMatrix mGradient;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray mCenterCoord;

    double domainSize();

public:
    PrescribedGradientBC(int n, Domain *d);
    virtual ~PrescribedGradientBC();

    virtual bcType giveType() const { return UnknownBT; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void scale(double s) { mGradient.times(s); }

    void giveGradientVoigt(FloatArray &oGradient) const;
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBC_H_ */
