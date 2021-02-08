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

#include "prescribedgradientbcdirichletRC.h"
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
#include "crosssection.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCDirichletRC);

double PrescribedGradientBCDirichletRC :: give(Dof *dof, ValueModeType mode, double time)
{
    DofIDItem id = dof->giveDofID();
    int pos = this->dofs.findFirstIndexOf(id);

    if  ( reinfYBound && reinfYBound && ( pos == 3 ) ) { //reinforcement ends - rotation
        //Prescribing rotations at beams' ends. Defining sets and checking whether the node belongs to one of them
        //Works only in 2d, i.e. for a reinforcement with 3 degrees of freedom
        bool isXReinf = false;
        bool isYReinf = false;

        domain->giveSet(reinfXBound)->giveNodeList().contains( dof->giveDofManager()->giveGlobalNumber() ) ? isXReinf = true : isXReinf = false;
        domain->giveSet(reinfYBound)->giveNodeList().contains( dof->giveDofManager()->giveGlobalNumber() ) ? isYReinf = true : isYReinf = false;

        int pos = this->dofs.findFirstIndexOf(dof->giveDofID());
        if (pos == 3 && isXReinf) {
            return -mGradient.at(2,1);
        } else if (pos ==3 && isYReinf) {
            return mGradient.at(2,1);
        } else {
            return 0.;
        }
    } else {
        return PrescribedGradient::give(dof, mode, time);
    }
}


void PrescribedGradientBCDirichletRC :: updateCoefficientMatrix(FloatMatrix &C)
//Modified by AS:
//Include end moments from the reinforcement.
//\sum (R_L e_l + R_perp e_{\perp} ) \outerp (x-\bar{x}) already included in C^T.R_c (in computeField)
//Added term \sum R_M e_{\perp} \outerp e_l
// C = [x                0              0               y  ]
//     [0                y              y               0  ]
//     [ePerp1*eL1   ePerp2*eL2    ePerp1*eL2    ePerp2*eL1]
//  for DofManagers with rotational degrees of freedom.
{
    PrescribedGradient::updateCoefficientMatrix(C);

    Domain *domain = this->giveDomain();
    int nsd = domain->giveNumberOfSpatialDimensions();

    if ( ( reinfXBound && reinfYBound && ( nsd == 2 ) ) ) {
        Set* setX = domain->giveSet(reinfXBound);
        Set* setY = domain->giveSet(reinfYBound);

        FloatArrayF<2> eLX, ePerpX, eLY, ePerpY; //tangential and normal basis vectors, hardcoded for horizontal and vertical reinforcement
        eLX = {1., 0.};
        ePerpX = {0., 1.};
        eLY = {0., 1.};
        ePerpY = {-1., 0.};

        for (auto &n : domain->giveDofManagers() ) {
            if ( n->giveNumberOfDofs() == 3 ) { //beam element node in 2d
                Dof *drot = n->giveDofWithID( this->dofs[2] );
                int krot = drot->__givePrescribedEquationNumber();

                if ( setX->giveNodeList().contains( n->giveGlobalNumber() ) ) {
                    C.at( krot, 1 ) = ePerpX.at( 1 ) * eLX.at( 1 );
                    C.at( krot, 2 ) = ePerpX.at( 2 ) * eLX.at( 2 );
                    C.at( krot, 3 ) = ePerpX.at( 1 ) * eLX.at( 2 );
                    C.at( krot, 4 ) = ePerpX.at( 2 ) * eLX.at( 1 );
                } else if ( setY->giveNodeList().contains( n->giveGlobalNumber() ) ) {
                    C.at( krot, 1 ) = ePerpY.at( 1 ) * eLY.at( 1 );
                    C.at( krot, 2 ) = ePerpY.at( 2 ) * eLY.at( 2 );
                    C.at( krot, 3 ) = ePerpY.at( 1 ) * eLY.at( 2 );
                    C.at( krot, 4 ) = ePerpY.at( 2 ) * eLY.at( 1 );
                }
            } else {
                continue;
            }
        }
    }
}


double PrescribedGradientBCDirichletRC::domainSize( Domain *d, int set )
{
    double omegaBox = PrescribedGradientHomogenization::domainSize(d, conBoundSet);

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant in 2D
        Element *e = this->giveDomain()->giveElement( this->giveDomain()->giveSet(conBoundSet)->giveBoundaryList().at(1) );
        std::unique_ptr<IntegrationRule> ir = e->giveInterpolation()->giveIntegrationRule( e->giveInterpolation()->giveInterpolationOrder() );
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(1);
        double thickness = cs->give(CS_Thickness, gp);
        return omegaBox * thickness;
    } else {
        return omegaBox;
    }
}


void PrescribedGradientBCDirichletRC :: initializeFrom(InputRecord &ir)
{
    PrescribedGradient :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, conBoundSet, _IFT_PrescribedGradientBCDirichletRC_ConcreteBoundary);
    IR_GIVE_OPTIONAL_FIELD(ir, reinfXBound , _IFT_PrescribedGradientBCDirichletRC_ReinfXBound);
    IR_GIVE_OPTIONAL_FIELD(ir, reinfYBound , _IFT_PrescribedGradientBCDirichletRC_ReinfYBound);
}


void PrescribedGradientBCDirichletRC :: giveInputRecord(DynamicInputRecord &input)
{
    PrescribedGradient :: giveInputRecord(input);
    input.setField(conBoundSet, _IFT_PrescribedGradientBCDirichletRC_ConcreteBoundary);
    input.setField(reinfXBound, _IFT_PrescribedGradientBCDirichletRC_ReinfXBound);
    input.setField(reinfYBound, _IFT_PrescribedGradientBCDirichletRC_ReinfYBound);
}
} // end namespace oofem
