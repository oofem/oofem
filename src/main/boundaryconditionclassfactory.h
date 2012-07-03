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

#include "generalbc.h"
#include "deadwght.h"
#include "nodload.h"
#include "linearedgeload.h"
#include "constantedgeload.h"
#include "constantsurfaceload.h"
#include "pointload.h"
#include "prescribedgradient.h"
#include "mixedgradientpressurebc.h"
#include "surfacetensionbc.h"

#ifdef __SM_MODULE
 #include "structtemperatureload.h"
 #include "structeigenstrainload.h"
 #include "usrdeftempfield.h"
 #include "tf1.h"
 #include "rotatingboundary.h"
#endif //__SM_MODULE

#ifdef __FM_MODULE
 #include "tractionpressurebc.h"
#endif //__FM_MODULE

REGISTER_CLASS(BoundaryCondition, "boundarycondition", BoundaryConditionClass)
REGISTER_CLASS(DeadWeight, "deadweight", DeadWeightClass)
REGISTER_CLASS(NodalLoad, "nodalload", NodalLoadClass)

REGISTER_CLASS(PrescribedGradient, "prescribedgradient", PrescribedGradientClass)
REGISTER_CLASS(MixedGradientPressureBC, "mixedgradientpressure", MixedGradientClass)
REGISTER_CLASS(LinearEdgeLoad, "linearedgeload", LinearEdgeLoadClass)
REGISTER_CLASS(ConstantEdgeLoad, "constantedgeload", ConstantEdgeLoadClass)
REGISTER_CLASS(ConstantSurfaceLoad, "constantsurfaceload", ConstantSurfaceLoadClass)
REGISTER_CLASS(PointLoad, "pointload", PointLoadClass)
REGISTER_CLASS(SurfaceTensionBoundaryCondition, "surfacetension", SurfaceTensionBoundaryConditionClass)

#ifdef __SM_MODULE
REGISTER_CLASS(StructuralTemperatureLoad, "structtemperatureload", StructuralTemperatureLoadClass)
REGISTER_CLASS(StructuralEigenstrainLoad, "structeigenstrainload", StructuralEigenstrainLoadClass)
REGISTER_CLASS(UserDefinedTemperatureField, "usrdeftempfield", UserDefinedTemperatureFieldClass)
REGISTER_CLASS(RotatingBoundary, "rotatingboundary", RotatingBoundaryClass)
REGISTER_CLASS(TF1, "tf1", TF1Class)
#endif //__SM_MODULE
#ifdef __FM_MODULE
REGISTER_CLASS(TractionPressureBC, "prescribedtractionpressurebc", TractionPressureBCClass)
#endif //__FM_MODULE

