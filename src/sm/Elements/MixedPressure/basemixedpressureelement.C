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


#include "../sm/Elements/MixedPressure/basemixedpressureelement.h"
#include "../sm/Materials/MixedPressure/mixedpressurematerialextensioninterface.h"


#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Materials/structuralms.h"

#include "material.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "domain.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "unknownnumberingscheme.h"


#include <cstdio>

namespace oofem {
BaseMixedPressureElement :: BaseMixedPressureElement()
{}



void
BaseMixedPressureElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_p, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u, const IntArray &dofIdArray_p)
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_p.clear();
    NLStructuralElement *el = this->giveElement();
    int k = 0;
    IntArray nodalArray;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = el->giveDofManager(i);
        int itt = 1;
        for ( int j = 1; j <= dofIdArray_u.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_u.at(j) ) ) {
                //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
                locationArray_u.followedBy(k + itt);
            }
            itt++;
        }
        for ( int j = 1; j <= dofIdArray_p.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_p.at(j) ) ) {
                //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
                locationArray_p.followedBy(k + itt);
            }
            itt++;
        }
        k += dMan->giveNumberOfDofs();
    }
}


void
BaseMixedPressureElement :: computeStressVector(FloatArray &stress, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    MixedPressureMaterialExtensionInterface *mixedPressureMat = static_cast< MixedPressureMaterialExtensionInterface * >( cs->giveMaterialInterface(MixedPressureMaterialExtensionInterfaceType, gp) );
    if ( !mixedPressureMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }

    double pressure;

    if ( elem->giveGeometryMode() == 0 ) {
        FloatArray strain;
        this->computeStrainVector(strain, gp, tStep);
        this->computePressure(pressure, gp, tStep);
        mixedPressureMat->giveRealStressVector(stress, gp, strain, pressure, tStep);
    } else {
        /*      FloatArray devF;
         * elem->computeDeviatoricDeformationGradientVector(devF, gp, tStep);
         * this->giveDofManDofIDMask_p( IdMask_p );
         * this->computePressure(pressure,IdMask_p, gp, tStep);
         * mixedPressureMat->giveFiniteStrainStressVectors(stress, gp, devF, pressure, tStep);
         */
        OOFEM_ERROR("Large-strain formulaton is not available yet")
    }
}


void
BaseMixedPressureElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray d_u;
    FloatMatrix B;
    NLStructuralElement *elem = this->giveElement();
    IntArray IdMask_u;
    this->giveDofManDofIDMask_u(IdMask_u);
    this->giveElement()->computeVectorOf(IdMask_u, VM_Total, tStep, d_u);
    elem->computeBmatrixAt(gp, B);
    answer.beProductOf(B, d_u);
}


void
BaseMixedPressureElement :: computePressure(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_p;
    FloatArray d_p, N_p;

    this->giveDofManDofIDMask_p(IdMask_p);
    this->giveElement()->computeVectorOf(IdMask_p, VM_Total, tStep, d_p);
    this->computePressureNMatrixAt(gp, N_p);
    answer = N_p.dotProduct(d_p);
}



void
BaseMixedPressureElement :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatArray BS, vStress, s, M;
    FloatMatrix B;
    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        if ( useUpdatedGpRecord == 1 ) {
            if ( this->giveElement()->giveGeometryMode() == 0 ) {
                vStress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();
            } else {
                vStress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveTempPVector();
            }
        } else {
            this->computeStressVector(vStress, gp, tStep);
        }

        MixedPressureMaterialExtensionInterface *mixedPressureMat = dynamic_cast< MixedPressureMaterialExtensionInterface * >( cs->giveMaterialInterface(MixedPressureMaterialExtensionInterfaceType, gp) );
        if ( !mixedPressureMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }
        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = elem->computeVolumeAround(gp);
        elem->computeBmatrixAt(gp, B);
        BS.beTProductOf(B, vStress);
        answer.add(dV, BS);
    }
}


void
BaseMixedPressureElement :: giveInternalForcesVector_p(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    double pressure, kappa, factor;
    IntArray IdMask_u;
    FloatArray N_p, Bvol, d_u;
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();


    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        MixedPressureMaterialExtensionInterface *mixedPressureMat = dynamic_cast< MixedPressureMaterialExtensionInterface * >( cs->giveMaterialInterface(MixedPressureMaterialExtensionInterfaceType, gp) );
        if ( !mixedPressureMat ) {
            OOFEM_ERROR("Material doesn't implement the required Mixed pressure interface!");
        }
        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = elem->computeVolumeAround(gp);
        this->computePressure(pressure, gp, tStep);
        //      this->computePressureNMatrixAt(gp, N_p);
        this->computeVolumetricBmatrixAt(gp, Bvol, elem);
        this->giveDofManDofIDMask_u(IdMask_u);
        elem->computeVectorOf(IdMask_u, VM_Total, tStep, d_u);

        double eps_V = Bvol.dotProduct(d_u);
        mixedPressureMat->giveInverseOfBulkModulus(kappa, TangentStiffness, gp, tStep);
        factor = eps_V + pressure * kappa;
        answer.times(dV * factor);
    }
}

void
BaseMixedPressureElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.resize( this->giveNumberOfDofs() );
    answer.zero();

    FloatArray answer_u( this->giveNumberOfDisplacementDofs() );
    answer_u.zero();
    FloatArray answer_p( this->giveNumberOfPressureDofs() );
    answer_p.zero();


    this->giveInternalForcesVector_u(answer_u, tStep, 0);
    this->giveInternalForcesVector_p(answer_p, tStep, 0);
    answer.assemble(answer_u, locationArray_u);
    answer.assemble(answer_p, locationArray_p);
}

void
BaseMixedPressureElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    FloatArray localForces( this->giveNumberOfDisplacementDofs() );
    answer.resize( this->giveNumberOfDofs() );
    this->computeLocForceLoadVector(localForces, tStep, mode);
    answer.assemble(localForces, locationArray_u);
}


/************************************************************************/
void
BaseMixedPressureElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveElement();
    //@todo check this
    //    elem->computeLocalForceLoadVector(answer, tStep, mode);

    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    } else {
        answer.resize( this->giveNumberOfDisplacementDofs() );
        answer.zero();
    }
}


void
BaseMixedPressureElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    answer.resize( this->giveNumberOfDofs(), this->giveNumberOfDofs() );
    answer.zero();

    FloatMatrix Kuu, Kup, Kpu, Kpp;
    this->computeStiffnessMatrix_uu(Kuu, rMode, tStep);
    this->computeStiffnessMatrix_up(Kup, rMode, tStep);
    Kpu.beTranspositionOf(Kup);
    this->computeStiffnessMatrix_pp(Kpp, rMode, tStep);

    answer.assemble(Kuu, locationArray_u);
    answer.assemble(Kup, locationArray_u, locationArray_p);
    answer.assemble(Kpu, locationArray_p, locationArray_u);
    answer.assemble(Kpp, locationArray_p);
}


void
BaseMixedPressureElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, D, DB;
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();
    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        MixedPressureMaterialExtensionInterface *mixedPressureMat = dynamic_cast< MixedPressureMaterialExtensionInterface * >(
            cs->giveMaterialInterface(MixedPressureMaterialExtensionInterfaceType, gp) );
        if ( !mixedPressureMat ) {
            OOFEM_ERROR("Material doesn't implement the required Mixed Pressure Material interface!");
        }


        mixedPressureMat->giveDeviatoricConstitutiveMatrix(D, rMode, gp, tStep);
        elem->computeBmatrixAt(gp, B);
        double dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);

        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}



void
BaseMixedPressureElement :: computeStiffnessMatrix_up(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double dV;
    FloatArray N_p, Bvol;
    NLStructuralElement *elem = this->giveElement();

    answer.clear();

    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        this->computeVolumetricBmatrixAt(gp, Bvol, elem);
        this->computePressureNMatrixAt(gp, N_p);
        FloatMatrix mNp=FloatMatrix::fromArray(N_p, true);
        FloatMatrix mBvol=FloatMatrix::fromArray(Bvol, true);

        dV = elem->computeVolumeAround(gp);
        answer.plusProductUnsym(mBvol, mNp, -dV);
    }
}

void
BaseMixedPressureElement :: computeStiffnessMatrix_pp(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double dV, kappa;
    FloatArray N_p;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    answer.clear();

    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
        MixedPressureMaterialExtensionInterface *mixedPressureMat = dynamic_cast< MixedPressureMaterialExtensionInterface * >( cs->giveMaterialInterface(MixedPressureMaterialExtensionInterfaceType, gp) );
        if ( !mixedPressureMat ) {
            OOFEM_ERROR("Material doesn't implement the required Mixed Pressure Material interface!");
        }

        mixedPressureMat->giveInverseOfBulkModulus(kappa, rMode, gp, tStep);
        this->computePressureNMatrixAt(gp, N_p);
        FloatMatrix mN_p=FloatMatrix::fromArray(N_p, true);
        dV = elem->computeVolumeAround(gp);
        answer.plusProductUnsym(mN_p, mN_p, -dV * kappa);
    }
}



void
BaseMixedPressureElement :: initializeFrom(InputRecord &ir)
{
    // @todo Is this function necessary???
}

void
BaseMixedPressureElement :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray stress, strain;
    /*
     * // force updating strains & stresses
     * for ( auto &iRule: integrationRulesArray ) {
     *  for ( GaussPoint *gp: *iRule ) {
     *      this->computeStrainVector(strain, gp, tStep);
     *      this->computeStressVector(stress, strain, gp, tStep);
     *  }
     * }
     */
}


void
BaseMixedPressureElement :: postInitialize()
{
    IntArray IdMask_u, IdMask_p;
    this->giveDofManDofIDMask_u(IdMask_u);
    this->giveDofManDofIDMask_p(IdMask_p);
    this->giveLocationArrayOfDofIDs(locationArray_u, locationArray_p, EModelDefaultEquationNumbering(), IdMask_u, IdMask_p);
}
} // end namespace oofem
