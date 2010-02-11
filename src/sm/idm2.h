/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.h,v 1.9.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *           Copyright (C) 2010 Christian Hoover, Vit Smilauer
 *
 *
 *
 *   Czech Technical University, Faculty of Civil Engineering,
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

//   ******************************************************************************
//   *** CLASS NONLOCAL ISOTROPIC DAMAGE MODEL 2 FOR CONCRETE IN TENSION ************
//   ******************************************************************************

#ifndef idm2_h
#define idm2_h

#include "idm1.h"

namespace oofem {
/**
 * This class implements a simple local isotropic damage model for concrete in tension.
 * The material model operates on an equivalent strain(s) and fracture energy
 * Linear, bi-linear and exponential softening curves are implemented
 */

class IsotropicDamageMaterial2 : public IsotropicDamageMaterial1
{
protected:

    /// This is a switch that will determine which softenign law to use
    int SofteningType;
    /// Max effective strain at peak
    double eps0;
    /// Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area under the stress/strain curve.  For an exponential law, it is the area bounded by the elastic range and a tangent to the softening part of the curve at the peak stress. For a bilinear law, Gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress
    double gf;
    /// Determines the softening for the bilinear law -> corresponds to the strain at the knee point
    double epsk;
    /// Determines the softening for the bilinear law -> corresponds to the total energy
    double gft;

public:

    /// Constructor
    IsotropicDamageMaterial2(int n, Domain * d);
    /// Destructor
    ~IsotropicDamageMaterial2();

    // identification and auxiliary functions
    const char *giveClassName() const { return "IsotropicDamageMaterial2"; }
    classType giveClassID()         const { return IsotropicDamageMaterial2Class; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "idm2"; }

    /**
     * Initializes receiver acording to object description stored in input record.
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);
    virtual void initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp);
};
} // end namespace oofem
#endif // idm2_h
