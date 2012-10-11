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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#ifndef reinforcement_h
#define reinforcement_h

#include "bodyload.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

namespace oofem {
/**
 * This class implements an influence of reinforcement into flow problems, especially concrete (binhamfluid).
 * It is modeled as a special body load acting on elements in area, where reinfocement is placed. The inherited attribute 'componentArray' contains the components of an
 * loading prescribed per unit volume.
 *
 * Its task is to return the body force @f$ \rho a @f$
 */
class Reinforcement : public BodyLoad
{
 protected:
  double porosity;
  double shapefactor;
  FloatArray permeability;


public:
    /// Constructor
    Reinforcement(int i, Domain *d) : BodyLoad(i, d) { }
    /**
     * Computes components values of deadweight field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value while respecting load response mode.
     * @param answer Component values at given point and time.
     * @param atTime Time step.
     * @param coords Global coordinates, which are used to evaluate components values.
     * @param mode Determines response mode-
     */
    void computeValueAt(FloatArray &answer, TimeStep *atTime, FloatArray &coords, ValueModeType mode)
    { computeComponentArrayAt(answer, atTime, mode); }

    bcValType giveBCValType() const { return ReinforceBVT; }
    bcGeomType giveBCGeoType() const { return BodyLoadBGT; }
    
    IRResultType initializeFrom(InputRecord *ir);
    classType giveClassID() const { return ReinforcementClass; }
    const char *giveClassName() const { return "Reinforcement"; }
    const char *giveInputRecordName() const { return "Reinforcement"; }

    ///Accessor
    double givePorosity(){ return porosity; }
    double giveshapefactor() { return shapefactor; }
    FloatArray* givePermeability() { return &permeability; }
};
} // end namespace oofem
#endif // deadwght_h






