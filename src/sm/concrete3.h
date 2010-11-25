/* $Header: /home/cvs/bp/oofem/sm/src/concrete3.h,v 1.5 2003/04/06 14:08:30 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *********************************************************************************
//   *** CLASS CONCRETE1
//   ********************************************************************************

#ifndef concrete3_h
#define concrete3_h

#include "rcm2.h"
#include "isolinearelasticmaterial.h"

namespace oofem {
class Concrete3 : public RCM2Material
{
    /*
     * This class implements a Concrete3 material in a finite element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     * DESCRIPTION
     * Concrete3 is NonLinear elasto-plastic material model of concrete without hardening
     * in compression.
     * Cracking is described using Rotating Crack Model based on fracture energy
     * criterion. Softening is linear or exponencial, with user defined unloading (parameter beta).
     *
     * TASK
     *
     */
    enum Concrete3_softeningMode { linearSoftening, exponentialSoftening };
private:

    //double shearRetFactor; // shearRetentionFactor
    Concrete3_softeningMode softeningMode;

public:

    Concrete3(int n, Domain *d);
    ~Concrete3() { delete linearElasticMaterial; }
    // identification and auxiliary functions

    IRResultType initializeFrom(InputRecord *ir);
    int       hasNonLinearBehaviour()     { return 1; }
    const char *giveClassName() const { return "Concrete3"; }
    classType giveClassID()          const { return Concrete3Class; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:

    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                       double crackStrain, int i);
    //virtual     double giveShearRetentionFactor(GaussPoint* gp, double eps_cr, int i);
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i);
    virtual double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i);
    //virtual     void   updateStatusForNewCrack( GaussPoint*, int, double);
    virtual double computeStrength(GaussPoint *, double);
    virtual int    checkSizeLimit(GaussPoint *gp, double);
};
} // end namespace oofem
#endif // concrete3_h
