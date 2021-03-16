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
#include "functioncontactsegment.h"
#include "floatarray.h"

#define NRFunctionContact_Maxiter 40
#define NRFunctionContact_Tolerance 1.e-6

namespace oofem {
    //virtual class implementing Newton-Rhapson iteration for analytical function contact segments
    //children need only implement the function, derivative and double derivative functions
    class NRFunctionContactSegment2D : public FunctionContactSegment
    {
    public:
        NRFunctionContactSegment2D(int n, Domain *aDomain) : FunctionContactSegment(n, aDomain) { ; }
        ~NRFunctionContactSegment2D() {};

    protected:

        virtual void computeContactPoint(FloatArray& answer, FloatArray& normal, const FloatArray& nodeCoords) override;

        virtual double functionValue(const double x) const = 0;
        virtual double derivativeValue(const double x) const = 0;
        virtual double doubleDerivativeValue(const double x) const = 0;
    };

}