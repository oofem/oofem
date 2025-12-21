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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "sm/Materials/rcm2.h"
#include "sm/Materials/isolinearelasticmaterial.h"

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
    enum Concrete3_softeningMode { linearSoftening, exponentialSoftening, hordijkSoftening };

private:
    //double shearRetFactor; // shearRetentionFactor
    Concrete3_softeningMode softeningMode = linearSoftening;

public:
    Concrete3(int n, Domain * d);
    virtual ~Concrete3() {
        delete linearElasticMaterial;
    }

     void initializeFrom(InputRecord &ir) override;
     const char *giveClassName() const override { return "Concrete3"; }
     const char *giveInputRecordName() const override { return _IFT_Concrete3_Name; }

     std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

protected:
    double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                               double crackStrain, int i) const override;
    //double giveShearRetentionFactor(GaussPoint* gp, double eps_cr, int i) override;
    double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i) const override;
    double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i) const override;
    //void updateStatusForNewCrack( GaussPoint*, int, double) override;
    double computeStrength(GaussPoint *, double) const override;
    int checkSizeLimit(GaussPoint *gp, double) const override;
};
} // end namespace oofem
#endif // concrete3_h
