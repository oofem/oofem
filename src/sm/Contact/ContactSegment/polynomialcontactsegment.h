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



#pragma once
#include "classfactory.h"
#include "nrfunctioncontactsegment2d.h"

#define _IFT_PolynomialContactSegment_Name "polynomialcontactsegment"
#define _IFT_PolynomialContactSegment_coeffs "coeffs"

namespace oofem {
    class PolynomialContactSegment : public NRFunctionContactSegment2D
    {
    public:
        PolynomialContactSegment(int n, Domain *aDomain) : NRFunctionContactSegment2D(n, aDomain) { ; }
        ~PolynomialContactSegment() {};

        IRResultType initializeFrom(InputRecord * ir) override;

        const char *giveClassName() const override { return "Polynomialcontactsegment"; }
        const char *giveInputRecordName() const override { return _IFT_PolynomialContactSegment_Name; }

    private:
        int order;
        FloatArray coeffs;

    protected:

        double functionValue(const double x) const override;
        double derivativeValue(const double x) const override;
        double doubleDerivativeValue(const double x) const override;
    };

}