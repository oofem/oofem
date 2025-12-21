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

#include "prescribeddispslipbcdirichletrc.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "function.h"
#include "engngm.h"
#include "set.h"
#include "node.h"
#include "element.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "assemblercallback.h"
#include "mathfem.h"
#include "floatarrayf.h"
#include "crosssection.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedDispSlipBCDirichletRC);

double PrescribedDispSlipBCDirichletRC :: give(Dof *dof, ValueModeType mode, double time)
{
    DofIDItem id = dof->giveDofID();
    const auto &coords = dof->giveDofManager()->giveCoordinates();

    double factor = 0;
    if ( mode == VM_Total ) {
        factor = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        factor = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        factor = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }

    FloatArray dx;
    dx.beDifferenceOf(coords, this->mCenterCoord);

    dispGradient.resizeWithData(coords.giveSize(), coords.giveSize());
    dispField.resizeWithValues(coords.giveSize());
    slipField.resizeWithValues(coords.giveSize());
    slipGradient.resizeWithData(coords.giveSize(), coords.giveSize());

    FloatArray u;
    u.beProductOf(dispGradient, dx);
    u.at(1) += dispField.at(1);
    u.at(2) += dispField.at(2);
    u.times( factor );


    FloatArray us;
    us.beProductOf(slipGradient, dx);
    us.at(1) += slipField.at(1);
    us.at(2) += slipField.at(2);
    us.times( factor );

    int pos = this->dofs.findFirstIndexOf(id);

    bool onConcrete = true;
    bool onSteel = false;

    if ( reinfXBound && reinfYBound ) {
        Set* setx = domain->giveSet(reinfXBound);
        Set* sety = domain->giveSet(reinfYBound);

        if ( setx->giveNodeList().contains( dof->giveDofManGlobalNumber() ) || sety->giveNodeList().contains( dof->giveDofManGlobalNumber() ) ) {
            onConcrete = false;
            onSteel = true;
        }
    }

    if ( onConcrete ) { //node on concrete
        return this->giveOnConcrete(dof, pos, u);
    } else if ( onSteel ) { //node on steel
        return this->giveOnSteel(dof, pos, u, us);
    } else {
        return 0.0;
    }

}


double PrescribedDispSlipBCDirichletRC::giveOnConcrete( Dof *dof, int pos, const FloatArray u )
{
    if ( pos > 0 && pos <= u.giveSize() ) {
        return u.at(pos);
    } else {
        return 0.0;
    }

}


double PrescribedDispSlipBCDirichletRC::giveOnSteel( Dof *dof, int pos, const FloatArray u, const FloatArray us )
{
    bool isXReinf = false;
    bool isYReinf = false;

    domain->giveSet(reinfXBound)->giveNodeList().contains( dof->giveDofManager()->giveGlobalNumber() ) ? isXReinf = true : isXReinf = false;
    domain->giveSet(reinfYBound)->giveNodeList().contains( dof->giveDofManager()->giveGlobalNumber() ) ? isYReinf = true : isYReinf = false;

    if ( isXReinf ) {
        if ( pos == 1 ) {
            return u.at(pos) + us.at(pos);
        } else if ( pos == 2 ) {
            return u.at(pos);
        } else if ( pos == 3 ) {
            return -dispGradient.at(2,1);
        } else {
            return 0.0;
        }
    } else if ( isYReinf ) {
        if ( pos == 1 ) {
            return u.at(pos);
        } else if ( pos == 2 ) {
            return u.at(pos) + us.at(pos);
        } else if ( pos == 3 ) {
            return dispGradient.at(2,1);
        } else {
            return 0.0;
        }
    } else {
        return 0.0;
    }
}


void PrescribedDispSlipBCDirichletRC::updateCoefficientMatrix( FloatMatrix &C )
{
// This is written in a very general way, supporting both fm and sm problems.
// v_prescribed = C.d = (x-xbar).d;
// Modified by ES.
// C = [x 0 0 y]
//     [0 y x 0]
//     [ ... ] in 2D, voigt form [d_11, d_22, d_12 d_21]
//Modified by AS:
//Include end moments from the reinforcement.
//\sum (R_L e_l + R_perp e_{\perp} ) \outerp (x-\bar{x}) already included in C^T.R_c (in computeField)
//Added term \sum R_M e_{\perp} \outerp e_l
// C = [x                0              0               y  ]
//     [0                y              y               0  ]
//     [ePerp1*eL1   ePerp2*eL2    ePerp1*eL2    ePerp2*eL1]
//  for DofManagers with rotational degrees of freedom.
    Domain *domain = this->giveDomain();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    C.resize(npeq, nsd * nsd);
    C.zero();

    FloatArray &cCoords = this->giveCenterCoordinate();
    double xbar = cCoords.at(1), ybar = cCoords.at(2);

    for ( auto &n : domain->giveDofManagers() ) {
        const auto &coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID( this->dofs[0] );
        Dof *d2 = n->giveDofWithID( this->dofs[1] );
        int k1 = d1->__givePrescribedEquationNumber();
        int k2 = d2->__givePrescribedEquationNumber();
        int k3 = 0;
        FloatArrayF<2> eL = { 0., 0. };
        FloatArrayF<2> ePerp = { 0., 0. };
        if ( (reinfXBound && reinfYBound) && n->giveNumberOfDofs() == 3 ) {
            Dof *d3 = n->giveDofWithID( this->dofs[2] );
            k3 = d3->__givePrescribedEquationNumber();
            Set* setX = domain->giveSet(reinfXBound);
            Set* setY = domain->giveSet(reinfYBound);
            if ( setX->giveNodeList().contains( n->giveGlobalNumber() ) ) {
                eL = { 1., 0.};
                ePerp = { 0., 1.};
            } else if ( setY->giveNodeList().contains( n->giveGlobalNumber() ) ) {
                eL = { 0., 1.};
                ePerp = { -1., 0.};
            }
        }
        if ( nsd == 2 ) {
            if ( k1 ) {
                C.at(k1, 1) = coords.at(1) - xbar;
                C.at(k1, 4) = coords.at(2) - ybar;
            }

            if ( k2 ) {
                C.at(k2, 2) = coords.at(2) - ybar;
                C.at(k2, 3) = coords.at(1) - xbar;
            }

            if ( k3 ) {
                C.at(k3, 1) = ePerp.at(1) * eL.at(1);
                C.at(k3, 2) = ePerp.at(2) * eL.at(2);
                C.at(k3, 3) = ePerp.at(1) * eL.at(2);
                C.at(k3, 4) = ePerp.at(2) * eL.at(1);
            }
        } else {
            OOFEM_ERROR("Only 2D mode supported.\n");
        }
    }

}



void PrescribedDispSlipBCDirichletRC :: computeStress(FloatArray &sigma, TimeStep *tStep)
{
    if ( dispGradON ) {
        EngngModel *emodel = this->domain->giveEngngModel();
        int npeq           = emodel->giveNumberOfDomainEquations( this->giveDomain()->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
        FloatArray R_c( npeq ), R_ext( npeq );

        R_c.zero();
        R_ext.zero();
        emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
        emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
        R_c.subtract( R_ext );

        // Condense it;
        FloatMatrix C;
        this->updateCoefficientMatrix( C );
        sigma.beTProductOf( C, R_c );
        double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
        sigma.times( 1. / volRVE );
    }
}

void PrescribedDispSlipBCDirichletRC::computeTransferStress( FloatArray &bStress, TimeStep *tStep )
{
    //According to 1/(\Omega_\Box) * \sum R_L \be{l}
    // [ tau_bxx, tau_byy ]
    // only compute when applied to reinforcement

    if ( slipON ) {
        EngngModel *emodel = this->domain->giveEngngModel();
        Domain *domain = this->giveDomain();

        int nsd = domain->giveNumberOfSpatialDimensions();
        int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );

        FloatArray R_c(npeq), R_ext(npeq);

        R_c.zero();
        R_ext.zero();
        emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
                                EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
        //R_c contains reactions at the boundary - coming from concrete tractions, beam end forces (normal and shear) and end moments

        emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
                                EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );

        R_c.subtract(R_ext);

        //taking contribution only from end nodes of the reinforcement
        IntArray nodesX = domain->giveSet(reinfXBound)->giveNodeList();
        IntArray nodesY = domain->giveSet(reinfYBound)->giveNodeList();
        FloatArray eL;
        double R_L;
        eL.resize(nsd); eL.zero();
        bStress.resize(nsd); bStress.zero();

        for (int i=1; i<=nodesX.giveSize(); i++) {
            DofManager *d1 = domain->giveDofManager( nodesX.at(i) );
            eL.at(1) = 1;
            eL.at(2) = 0;
            R_L = R_c.at( d1->giveDofWithID(1)->__givePrescribedEquationNumber() );
            bStress.at(1) += R_L*eL.at(1);
            bStress.at(2) += R_L*eL.at(2);
        }

        for (int i=1; i<=nodesY.giveSize(); i++) {
            DofManager *d2 = domain->giveDofManager( nodesY.at(i) );
            eL.at(1) = 0;
            eL.at(2) = 1;
            R_L = R_c.at( d2->giveDofWithID(2)->__givePrescribedEquationNumber() );
            bStress.at(1) += R_L*eL.at(1);
            bStress.at(2) += R_L*eL.at(2);
        }

        double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
        bStress.times( 1. / volRVE );
    }
}


void PrescribedDispSlipBCDirichletRC::computeReinfStress( FloatArray &rStress, TimeStep *tStep )
{
    //According to 1/(\Omega_\Box) * \sum R_L \be{l} \outerp [x - \bar{x} ]
    //order: sxx, syy, sxy, syx
    // only compute when applied to reinforcement

    if ( slipGradON ) {
        EngngModel *emodel = this->domain->giveEngngModel();
        Domain *domain = this->giveDomain();

        int nsd = domain->giveNumberOfSpatialDimensions();
        int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );

        FloatArray R_c(npeq), R_ext(npeq);

        R_c.zero();
        R_ext.zero();
        emodel->assembleVector( R_c, tStep, InternalForceAssembler(), VM_Total,
                                EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
        //R_c contains reactions at the boundary - coming from concrete tractions, beam end forces (normal and shear) and end moments

        emodel->assembleVector( R_ext, tStep, ExternalForceAssembler(), VM_Total,
                                EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );

        R_c.subtract(R_ext);

        //taking contribution only from end nodes of the reinforcement
        IntArray nodesX = domain->giveSet(reinfXBound)->giveNodeList();
        IntArray nodesY = domain->giveSet(reinfYBound)->giveNodeList();
        FloatArray eL, cCoords = this->giveCenterCoordinate();
        double xbar=cCoords.at(1), ybar=cCoords.at(2), R_L;
        eL.resize(nsd); eL.zero();
        rStress.resize(nsd*2); rStress.zero();

        for (int i=1; i<=nodesX.giveSize(); i++) {
            DofManager *d1 = domain->giveDofManager( nodesX.at(i) );
            const auto &coords = d1->giveCoordinates();
            eL.at(1) = 1;
            eL.at(2) = 0;
            R_L = R_c.at( d1->giveDofWithID(1)->__givePrescribedEquationNumber() );
            rStress.at(1) += R_L*eL.at(1) * (coords.at(1) - xbar);
            rStress.at(2) += R_L*eL.at(2) * (coords.at(2) - ybar);
            rStress.at(3) += R_L*eL.at(1) * (coords.at(2) - ybar);
            rStress.at(4) += R_L*eL.at(2) * (coords.at(1) - xbar);
        }

        for (int i=1; i<=nodesY.giveSize(); i++) {
            DofManager *d2 = domain->giveDofManager( nodesY.at(i) );
            const auto &coords = d2->giveCoordinates();
            eL.at(1) = 0;
            eL.at(2) = 1;
            R_L = R_c.at( d2->giveDofWithID(2)->__givePrescribedEquationNumber() );
            rStress.at(1) += R_L*eL.at(1) * (coords.at(1) - xbar);
            rStress.at(2) += R_L*eL.at(2) * (coords.at(2) - ybar);
            rStress.at(3) += R_L*eL.at(1) * (coords.at(2) - ybar);
            rStress.at(4) += R_L*eL.at(2) * (coords.at(1) - xbar);
        }

        double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
        rStress.times( 1. / volRVE );
    }
}


void PrescribedDispSlipBCDirichletRC :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
// a = [a_c; a_f];
// K.a = [R_c,0];
// [K_cc, K_cf; K_fc, K_ff].[a_c;a_f] = [R_c; 0];
// a_c = d.[x-x_b] = [x-x_b].d = C.d
// E = C'.(K_cc - K_cf.K_ff^(-1).K_fc).C
//   = C'.(K_cc.C - K_cf.(K_ff^(-1).(K_fc.C)))
//   = C'.(K_cc.C - K_cf.a)
//   = C'.X
{
    if ( dispGradON && !slipON && !slipGradON ) {
        // Fetch some information from the engineering model
        EngngModel *rve = this->giveDomain()->giveEngngModel();
        std ::unique_ptr<SparseLinearSystemNM> solver( classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
        SparseMtrxType stype = solver->giveRecommendedMatrix( true );
        EModelDefaultEquationNumbering fnum;
        EModelDefaultPrescribedEquationNumbering pnum;

        // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
        std ::unique_ptr<SparseMtrx> Kff( classFactory.createSparseMtrx( stype ) );
        std ::unique_ptr<SparseMtrx> Kfp( classFactory.createSparseMtrx( stype ) );
        std ::unique_ptr<SparseMtrx> Kpf( classFactory.createSparseMtrx( stype ) );
        std ::unique_ptr<SparseMtrx> Kpp( classFactory.createSparseMtrx( stype ) );
        if ( !Kff ) {
            OOFEM_ERROR( "Couldn't create sparse matrix of type %d\n", stype );
        }
        Kff->buildInternalStructure( rve, 1, fnum );
        Kfp->buildInternalStructure( rve, 1, fnum, pnum );
        Kpf->buildInternalStructure( rve, 1, pnum, fnum );
        Kpp->buildInternalStructure( rve, 1, pnum );
        rve->assemble( *Kff, tStep, TangentAssembler( TangentStiffness ), fnum, this->domain );
        rve->assemble( *Kfp, tStep, TangentAssembler( TangentStiffness ), fnum, pnum, this->domain );
        rve->assemble( *Kpf, tStep, TangentAssembler( TangentStiffness ), pnum, fnum, this->domain );
        rve->assemble( *Kpp, tStep, TangentAssembler( TangentStiffness ), pnum, this->domain );

        FloatMatrix C, X, Kpfa, KfpC, a;

        this->updateCoefficientMatrix( C );
        Kpf->timesT( C, KfpC );
        solver->solve( *Kff, KfpC, a );
        Kpp->times( C, X );
        Kpf->times( a, Kpfa );
        X.subtract( Kpfa );
        tangent.beTProductOf( C, X );
        tangent.times( 1. / this->domainSize( this->giveDomain(), this->giveSetNumber() ) );
    }
}

void PrescribedDispSlipBCDirichletRC :: initializeFrom(InputRecord &ir)
{
    GeneralBoundaryCondition :: initializeFrom(ir);
    PrescribedDispSlipHomogenization::initializeFrom(ir);
    IR_GIVE_FIELD(ir, conBoundSet, _IFT_PrescribedDispSlipBCDirichletRC_ConcreteBoundary);

    if ( this->dispGradient.isNotEmpty() ) {
        this->dispGradON = true;
    }

    if ( this->slipField.isNotEmpty() ) {
        this->slipON = true;
    }

    if ( this->slipGradient.isNotEmpty() ) {
        this->slipGradON = true;
    }

    if ( slipON || slipGradON) {
        IR_GIVE_OPTIONAL_FIELD( ir, reinfXBound, _IFT_PrescribedDispSlipBCDirichletRC_ReinfXBound );
        IR_GIVE_OPTIONAL_FIELD( ir, reinfYBound, _IFT_PrescribedDispSlipBCDirichletRC_ReinfYBound );
    }
}


void PrescribedDispSlipBCDirichletRC :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);
    PrescribedDispSlipHomogenization :: giveInputRecord(input);
}


void PrescribedDispSlipBCDirichletRC::scale( double s )
{
    dispField.times(s);
    slipField.times(s);
    dispGradient.times(s);
    slipGradient.times(s);
}

double PrescribedDispSlipBCDirichletRC::domainSize( Domain *d, int set )
{
    double omegaBox = PrescribedDispSlipHomogenization::domainSize(d, conBoundSet);

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant in 2D
        Element *e = this->giveDomain()->giveElement( this->giveDomain()->giveSet(conBoundSet)->giveBoundaryList().at(1) );
        std::unique_ptr<IntegrationRule> ir = e->giveInterpolation()->giveIntegrationRule( e->giveInterpolation()->giveInterpolationOrder(), e->giveGeometryType() );
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(0);
        double thickness = cs->give(CS_Thickness, gp);
        return omegaBox * thickness;
    } else {
        return omegaBox;
    }

}


} // end namespace oofem
