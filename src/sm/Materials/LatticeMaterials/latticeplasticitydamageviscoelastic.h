
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
#define _IFT_LatticePlasticityDamageViscoelastic_slaveMat "slavemat"

//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticePlasticityDamageViscoelastic.
 * @author: Petr Havalasek
 */
class LatticePlasticityDamageViscoelasticStatus : public LatticePlasticityDamageStatus

{
protected:
    GaussPoint *viscoelasticGP;
    /// 'slave' material model number.
    int slaveMat;

public:

    /// Constructor
    LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *g, int s);
    /// Destructor
    ~LatticePlasticityDamageViscoelasticStatus() {}

    /// Prints the receiver state to given stream
    void   printOutputAt(FILE *file, TimeStep *tStep);


    const char *giveClassName() const { return "LatticePlasticityDamageViscoelasticStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached


    virtual void saveContext(DataStream &stream, ContextMode mode);

    virtual void restoreContext(DataStream &stream, ContextMode mode);

    GaussPoint *giveViscoelasticGaussPoint() { return viscoelasticGP; }

    MaterialStatus *giveViscoelasticMatStatus();
};





/**
 * This class implements a local viscoelastic model for concrete in tension for 3D lattice elements.
 * @author: Petr Havlasek
 */
class LatticePlasticityDamageViscoelastic : public LatticePlasticityDamage

{
protected:
    /// 'slave' (= viscoelastic) material model number.
    int slaveMat;


public:

    /// Constructor
    LatticePlasticityDamageViscoelastic(int n, Domain *d);

    /// Destructor
    virtual ~LatticePlasticityDamageViscoelastic();

    virtual const char *giveInputRecordName() const { return _IFT_LatticePlasticityDamageViscoelastic_Name; }
    virtual const char *giveClassName() const { return "LatticePlasticityDamageViscoelastic"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }


    virtual void give2dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode rmode,
                                        GaussPoint *gp,
                                        TimeStep *atTime);

    virtual void give3dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode rmode,
                                        GaussPoint *gp,
                                        TimeStep *atTime);


    virtual int hasMaterialModeCapability(MaterialMode mode);


    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);


    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    RheoChainMaterial *giveViscoelasticMaterial();

    virtual void giveReducedStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);
};
} // end namespace oofem


#endif
