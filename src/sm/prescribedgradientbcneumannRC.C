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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "prescribedgradientbcneumannRC.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"
#include "fei2dquadlin.h"


namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCNeumannRC);

PrescribedGradientBCNeumannRC :: PrescribedGradientBCNeumannRC(int n, Domain *d) :
    PrescribedGradientBCNeumann(n, d)
{
}

PrescribedGradientBCNeumannRC :: ~PrescribedGradientBCNeumannRC()
{
}


double PrescribedGradientBCNeumannRC::domainSize( Domain *d, int set )
{
    double omegaBox = PrescribedGradientHomogenization::domainSize(d, this->giveSetNumber());

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant in 2D
        Element *e = this->giveDomain()->giveElement( this->giveDomain()->giveSet( this->giveSetNumber() )->giveBoundaryList().at(1) );
        std::unique_ptr<IntegrationRule> ir = e->giveInterpolation()->giveIntegrationRule( e->giveInterpolation()->giveInterpolationOrder() );
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(0);
        double thickness = cs->give(CS_Thickness, gp);
        return omegaBox * thickness;
    } else {
        return omegaBox;
    }
}


void PrescribedGradientBCNeumannRC::computeField( FloatArray &sigma, TimeStep *tStep )
{
    double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
    double areaRVE = PrescribedGradientHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());
    double thick = volRVE/areaRVE;
    mpSigmaHom->giveUnknownVector(sigma, mSigmaIds, VM_Total, tStep);
    sigma.times(1/thick);
}


void PrescribedGradientBCNeumannRC::assembleVector( FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorm, void *lock )
{
    IntArray sigma_loc;  // For the displacements and stress respectively
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for sigma, the load is simply the prescribed gradient.
        double rve_size = PrescribedGradientHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());
        FloatArray stressLoad;
        FloatArray gradVoigt;
        giveGradientVoigt(gradVoigt);

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        stressLoad.beScaled(-rve_size*loadLevel, gradVoigt);

        answer.assemble(stressLoad, sigma_loc);

    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ke;
        FloatArray fe_v, fe_s;
        FloatArray sigmaHom, e_u;
        IntArray loc, masterDofIDs, sigmaMasterDofIDs;

        // Fetch the current values of the stress;
        mpSigmaHom->giveUnknownVector(sigmaHom, mSigmaIds, mode, tStep);
        // and the master dof ids for sigmadev used for the internal norms
        mpSigmaHom->giveMasterDofIDArray(mSigmaIds, sigmaMasterDofIDs);

        // Assemble
        Set *setPointer = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = setPointer->giveBoundaryList();
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            e->giveLocationArray(loc, s, & masterDofIDs);
            e->computeVectorOf(mode, tStep, e_u);
            this->integrateTangent(Ke, e, boundary);

            // We just use the tangent, less duplicated code (the addition of sigmaDev is linear).
            fe_v.beProductOf(Ke, e_u);
            fe_s.beTProductOf(Ke, sigmaHom);

            // temporary hack for inconsistent normal vectors
            auto interp = dynamic_cast< FEI2dQuadLin * >( e->giveInterpolation());
            if ( interp != nullptr ) {
            } else {
                fe_v.negated();
                fe_s.negated();
            }


            answer.assemble(fe_s, loc); // Contributions to delta_v equations
            answer.assemble(fe_v, sigma_loc); // Contribution to delta_s_i equations
            if ( eNorm != NULL ) {
                eNorm->assembleSquared(fe_s, masterDofIDs);
                eNorm->assembleSquared(fe_v, sigmaMasterDofIDs);
            }
        }
    }
}


void PrescribedGradientBCNeumannRC::assemble( SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale, void *lock )
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;
        Set *set = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = set->giveBoundaryList();

        // Fetch the columns/rows for the stress contributions;
        mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_r, r_s);
        mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_c, c_s);

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            e->giveLocationArray(loc_r, r_s);
            e->giveLocationArray(loc_c, c_s);

            this->integrateTangent(Ke, e, boundary);

            //temporary hack for inconsistent normal vectors
            auto interp = dynamic_cast< FEI2dQuadLin * >( e->giveInterpolation() );
            if ( interp != nullptr ) {
            } else {
                Ke.negated();
            }

            Ke.times(scale);
            KeT.beTranspositionOf(Ke);

            answer.assemble(sigma_loc_r, loc_c, Ke); // Contribution to delta_s_i equations
            answer.assemble(loc_r, sigma_loc_c, KeT); // Contributions to delta_v equations

        }
    } else   {
        OOFEM_LOG_DEBUG("Skipping assembly in PrescribedGradientBCNeumann::assemble().");
    }

}


} /* namespace oofem */
