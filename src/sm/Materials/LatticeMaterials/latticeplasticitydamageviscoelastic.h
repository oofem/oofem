
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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#ifndef latticeplasticitydamageviscoelastic_h
#define latticeplasticitydamageviscoelastic_h

#include "latticeplasticitydamage.h"
#include "../RheoChainMaterials/rheoChM.h"

///@name Input fields for LatticePlasticityDamageViscoelastic
//@{
#define _IFT_LatticePlasticityDamageViscoelastic_Name "latticeplasticitydamageviscoelastic"
#define _IFT_LatticePlasticityDamageViscoelastic_viscoMat "viscomat"
#define _IFT_LatticePlasticityDamageViscoelastic_timeFactor "timefactor"

//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticePlasticityDamageViscoelastic.
 * @author: Petr Havalasek
 */
class LatticePlasticityDamageViscoelasticStatus : public LatticePlasticityDamageStatus

{
protected:
    std :: unique_ptr< GaussPoint >slaveGpVisco;

public:

    /// Constructor
    LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *g);

    /// Prints the receiver state to given stream
    void printOutputAt(FILE *file, TimeStep *tStep) const override;


    const char *giveClassName() const override { return "LatticePlasticityDamageViscoelasticStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;


    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

    GaussPoint *giveSlaveGaussPointVisco() const { return this->slaveGpVisco.get(); }
};



/**
 * This class implements a local viscoelastic model for concrete in tension for 3D lattice elements.
 * @author: Petr Havlasek
 */
class LatticePlasticityDamageViscoelastic : public LatticePlasticityDamage

{
protected:
    /// 'slave' (= viscoelastic) material model number.
    int viscoMat = 0;

public:

    /// Constructor
    LatticePlasticityDamageViscoelastic(int n, Domain *d);

    const char *giveInputRecordName() const override { return _IFT_LatticePlasticityDamageViscoelastic_Name; }
    const char *giveClassName() const override { return "LatticePlasticityDamageViscoelastic"; }

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                        GaussPoint *gp,
                                        TimeStep *tStep) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;

    int checkConsistency(void) override;
};
} // end namespace oofem


#endif
