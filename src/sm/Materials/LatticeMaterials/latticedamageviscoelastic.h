
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

#ifndef latticedamageviscoelastic_h
#define latticedamageviscoelastic_h

#include "latticedamage.h"
#include "../RheoChainMaterials/rheoChM.h"

///@name Input fields for LatticeDamage
//@{
#define _IFT_LatticeDamageViscoelastic_Name "latticedamageviscoelastic"
#define _IFT_LatticeDamageViscoelastic_viscoMat "viscomat"
#define _IFT_LatticeDamageViscoelastic_timeFactor "timefactor"

//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeDamageViscoelastic.
 * @author: Petr Havlasek
 */
class LatticeDamageViscoelasticStatus : public LatticeDamageStatus

{
protected:
    std :: unique_ptr< GaussPoint >slaveGpVisco;

public:

    /// Constructor
    LatticeDamageViscoelasticStatus(GaussPoint *g);

    /// Prints the receiver state to given stream
    void printOutputAt(FILE *file, TimeStep *tStep) const override;


    const char *giveClassName() const override { return "LatticeDamageViscoelasticStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override; // update after new equilibrium state reached


    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

    GaussPoint *giveSlaveGaussPointVisco() const { return this->slaveGpVisco.get(); }
};


/**
 * This class implements a local viscoelastic model for concrete in tension for 3D lattice elements.
 * @author: Petr Havlasek
 */
class LatticeDamageViscoelastic : public LatticeDamage

{
protected:
    /// 'slave' (= viscoelastic) material model number.
    int viscoMat = 0;


public:

    /// Constructor
    LatticeDamageViscoelastic(int n, Domain *d);

    const char *giveInputRecordName() const override { return _IFT_LatticeDamageViscoelastic_Name; }
    const char *giveClassName() const override { return "LatticeDamageViscoelastic"; }

    void initializeFrom(InputRecord &ir) override;


    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rmode,
                                                     GaussPoint *gp,
                                                     TimeStep *atTime) const override;

    FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                        GaussPoint *gp,
                                        TimeStep *tStep) override;


    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    RheoChainMaterial *giveViscoelasticMaterial();

protected:

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;

    int checkConsistency(void) override;
};
} // end namespace oofem


#endif
