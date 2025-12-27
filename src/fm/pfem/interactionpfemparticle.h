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

#ifndef interactionpfemparticle_h
#define interactionpfemparticle_h

#include "pfemparticle.h"

///@name Input fields for Pfemparticle
//@{
#define _IFT_InteractionPFEMParticle_Name "interactionpfemparticle"
//@}

namespace oofem {
class FloatArray;
class IntArray;
class StructuralEngngModel;
class FluidStructureProblem;
class ParamKey;

/**
 * This class represents a fluid particle attached to a node on the structural part
 * of the interface. The linked structural node is attribute of this class.
 */
class OOFEM_EXPORT InteractionPFEMParticle : public PFEMParticle
{
protected:
    int coupledNode;

    static ParamKey IPK_InteractionPFEMParticle_CoupledNode;
public:
    /**
     * Constructor. Creates a particle  with number n, belonging to aDomain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    InteractionPFEMParticle(int n, Domain *aDomain);
    /**
     * Destructor.
     */
    virtual ~InteractionPFEMParticle(void) { }

    void initializeFrom(InputRecord &ir, int priority) override;
    int checkConsistency() override;

    void updateYourself(TimeStep *tStep) override;

    void givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                     ValueModeType mode, TimeStep *stepN) override;

    void giveCoupledVelocities(FloatArray &answer, TimeStep *stepN);

    const char *giveClassName() const override { return "InteractionPFEMParticle"; }
    const char *giveInputRecordName() const override { return _IFT_InteractionPFEMParticle_Name; }

#ifdef __OOFEG
    void drawScalar(oofegGraphicContext &gc) override;
#endif

private:
    StructuralEngngModel *giveStructuralProblem();
    FluidStructureProblem *giveFluidStructureMasterProblem();
};
} // end namespace oofem
#endif // interactionpfemparticle_h
