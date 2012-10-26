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


#ifndef stokesflowvelocityhomogenization_h
#define stokesflowvelocityhomogenization_h

#include "stokesflow.h"
#include "rveengngmodel.h"
#include "engngm.h"
#include "sparselinsystemnm.h"
#include "usrdefsub.h"

namespace oofem
{
/**
 * Class for using the stokes flow class as an rve/constitutive model.
 *
 * @author Carl Sandström
 *
 */
class StokesFlowVelocityHomogenization : public StokesFlow, public rveEngngModel
{
protected:
    double areaOfDomain;
    double areaOfRVE;

public:
    StokesFlowVelocityHomogenization(int i, EngngModel *_master = NULL);
    virtual ~StokesFlowVelocityHomogenization();

    virtual void solveYourselfAt(TimeStep *tStep);

    void handlePrescribedValues();

    /** Compute area of domain (excludes holes)*/
    double giveAreaOfDomain();

    /** Compute area of domain (includes holes)*/
    double giveAreaOfRVE();

    virtual const char *giveClassName() const { return "StokesFlowVelocityHomogenization"; }
    virtual classType giveClassID() const { return StokesFlowVelocityHomogenizationClass; }

    void updateC();

    void computeTangent(FloatMatrix &answer, TimeStep *atTime);

    virtual void rveSetBoundaryConditions(int type, FloatArray value);
    virtual void rveGiveCharacteristicData(int type, void *value, void *answer, TimeStep *atTime);

private:
    /** Computes the mean velocity and pressure gradient */
    void getMeans(FloatArray &gradP, FloatArray &v, TimeStep *atTime);
};
}
#endif // stokesflowvelocityhomogenization_h
