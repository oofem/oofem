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

/*
 * matforceevaluator.h
 *
 *  Created on: Nov 12, 2014
 *      Author: svennine
 */

#ifndef MATFORCEEVALUATOR_H_
#define MATFORCEEVALUATOR_H_

namespace oofem {

class TipInfo;
class Domain;
class FloatArray;
class TimeStep;

/**
 * Evaluates material forces.
 *
 * Under development. Currently, only elastic material and traction free cracks are considered.
 *
 * @author Erik Svenning
 * @date Nov 12, 2014
 */
class MaterialForceEvaluator {
public:
    MaterialForceEvaluator();
    virtual ~MaterialForceEvaluator();

    void computeMaterialForce(FloatArray &oMatForce, Domain &iDomain, const TipInfo &iTipInfo, TimeStep *tStep, const double &iRadius);

    double computeWeightFunctionInPoint(const FloatArray &iCoord, const FloatArray &iTipCoord, const double &iRadius) const;
};

} /* namespace oofem */

#endif /* MATFORCEEVALUATOR_H_ */
