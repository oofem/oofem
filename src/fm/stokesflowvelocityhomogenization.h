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


#ifndef stokesflowvelocityhomogenization_h
#define stokesflowvelocityhomogenization_h

#include "stokesflow.h"
#include "engngm.h"
#include "sparselinsystemnm.h"
#include "classfactory.h"

#define _IFT_StokesFlowVelocityHomogenization_Name "stokesflowvelocityhomogenization"

namespace oofem
{
/**
 * Class for using the stokes flow class as an rve/constitutive model.
 *
 * @author Carl Sandström
 *
 */
class StokesFlowVelocityHomogenization : public StokesFlow
{
public:
    StokesFlowVelocityHomogenization(int i, EngngModel * _master = NULL);
    virtual ~StokesFlowVelocityHomogenization();

    /** Compute area of domain (includes holes)*/
    double giveAreaOfRVE();

    virtual const char *giveClassName() const { return "StokesFlowVelocityHomogenization"; }
    virtual const char *giveInputRecordName() const { return _IFT_StokesFlowVelocityHomogenization_Name; }


    void computeTangent(FloatMatrix &answer, TimeStep *tStep);
    /** Computes the mean velocity and pressure gradient */
    void computeSeepage(FloatArray &v, TimeStep *tStep);
    void applyPressureGradient(const FloatArray &grad);

private:
    void integrateNMatrix(FloatMatrix &N, Element &elem, TimeStep *tStep);
};
}
#endif // stokesflowvelocityhomogenization_h
