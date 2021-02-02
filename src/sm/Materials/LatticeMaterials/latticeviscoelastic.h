
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

#ifndef latticeviscoelastic_h
#define latticeviscoelastic_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"
#include "../RheoChainMaterials/rheoChM.h"

///@name Input fields for LatticeDamage
//@{
#define _IFT_LatticeViscoelastic_Name "latticeviscoelastic"
#define _IFT_LatticeViscoelastic_viscoMat "viscomat"

//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeViscoelastic.
 */
class LatticeViscoelasticStatus : public LatticeMaterialStatus

{
protected:
    std :: unique_ptr< GaussPoint >slaveGpVisco;

public:

    /// Constructor
    LatticeViscoelasticStatus(GaussPoint *g);

    /// Prints the receiver state to given stream
    void printOutputAt(FILE *file, TimeStep *tStep) const override;


    const char *giveClassName() const override { return "LatticeViscoelasticStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;


    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

    GaussPoint *giveSlaveGaussPointVisco() const { return this->slaveGpVisco.get(); }
};





/**
 * This class implements a local viscoelastic model for concrete in tension for 3D lattice elements.
 */
class LatticeViscoelastic : public LatticeLinearElastic

{
protected:
    /// 'slave' material model number.
    int viscoMat = 0;


public:

    /// Constructor
    LatticeViscoelastic(int n, Domain *d);


    const char *giveInputRecordName() const override { return _IFT_LatticeViscoelastic_Name; }
    const char *giveClassName() const override { return "LatticeViscoelastic"; }

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rmode,
                                                     GaussPoint *gp,
                                                     TimeStep *atTime) const override;

    FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rmode,
                                                     GaussPoint *gp,
                                                     TimeStep *atTime) const override;


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
