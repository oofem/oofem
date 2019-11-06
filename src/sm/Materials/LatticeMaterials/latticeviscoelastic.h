
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
#define _IFT_LatticeViscoelastic_slaveMat "slavemat"

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
    /// Destructor
    ~LatticeViscoelasticStatus() {}

    /// Prints the receiver state to given stream
    void   printOutputAt(FILE *file, TimeStep *tStep);


    const char *giveClassName() const { return "LatticeViscoelasticStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached


    virtual void saveContext(DataStream &stream, ContextMode mode);

    virtual void restoreContext(DataStream &stream, ContextMode mode);

    GaussPoint *giveSlaveGaussPointVisco() { return this->slaveGpVisco.get(); }

    MaterialStatus *giveViscoelasticMatStatus();
};





/**
 * This class implements a local viscoelastic model for concrete in tension for 3D lattice elements.
 */
class LatticeViscoelastic : public LatticeLinearElastic

{
protected:
    /// 'slave' material model number.
    int slaveMat;


public:

    /// Constructor
    LatticeViscoelastic(int n, Domain *d);

    /// Destructor
    virtual ~LatticeViscoelastic();

    virtual const char *giveInputRecordName() const { return _IFT_LatticeViscoelastic_Name; }
    virtual const char *giveClassName() const { return "LatticeViscoelastic"; }

    virtual void initializeFrom(InputRecord &ir);

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

protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);
};
} // end namespace oofem


#endif
