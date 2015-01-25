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

#ifndef concrete3_h
#define concrete3_h

#include "Materials/rcm2.h"
#include "Materials/isolinearelasticmaterial.h"

///@name Input fields for Concrete3
//@{
#define _IFT_Concrete3_Name "concrete3"
#define _IFT_Concrete3_exp_soft "exp_soft"
//@}

namespace oofem {
/**
 * This class implements a Concrete3 material in a finite element problem.
 *
 * Concrete3 is NonLinear elasto-plastic material model of concrete without hardening
 * in compression.
 * Cracking is described using Rotating Crack Model based on fracture energy
 * criterion. Softening is linear or exponential, with user defined unloading (parameter beta).
 */
class Concrete3 : public RCM2Material
{
    enum Concrete3_softeningMode { linearSoftening, exponentialSoftening };

private:
    //double shearRetFactor; // shearRetentionFactor
    Concrete3_softeningMode softeningMode;

public:
    Concrete3(int n, Domain * d);
    virtual ~Concrete3() {
        delete linearElasticMaterial;
    }

    // identification and auxiliary functions
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveClassName() const { return "Concrete3"; }
    virtual const char *giveInputRecordName() const { return _IFT_Concrete3_Name; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                       double crackStrain, int i);
    //virtual double giveShearRetentionFactor(GaussPoint* gp, double eps_cr, int i);
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i);
    virtual double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i);
    //virtual void updateStatusForNewCrack( GaussPoint*, int, double);
    virtual double computeStrength(GaussPoint *, double);
    virtual int checkSizeLimit(GaussPoint *gp, double);
};
} // end namespace oofem
#endif // concrete3_h
