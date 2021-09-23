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

#include "latticeframeelastic.h"
#include "latticematstatus.h"
#include "latticestructuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeFrameElastic);

bool
LatticeFrameElastic::hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeFrameElastic::initializeFrom(InputRecord &ir)
{
    LatticeStructuralMaterial::initializeFrom(ir);

    //Young's modulus of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->e, _IFT_LatticeFrameElastic_e); // Macro

    //Poisson's ratio of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->nu, _IFT_LatticeFrameElastic_n); // Macro
}

MaterialStatus *
LatticeFrameElastic::CreateStatus(GaussPoint *gp) const
{
    return new LatticeMaterialStatus(gp);
}

MaterialStatus *
LatticeFrameElastic::giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( !status ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status ) {
            gp->setMaterialStatus(status);
        }
    }

    return status;
}


FloatArrayF< 6 >
LatticeFrameElastic::giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);


    return {
        alpha, 0., 0., 0., 0., 0.
    };
}


FloatArrayF< 6 >
LatticeFrameElastic::giveFrameForces3d(const FloatArrayF< 6 > &strain,
                                       GaussPoint *gp,
                                       TimeStep *tStep)
{
    auto status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    auto stiffnessMatrix = LatticeFrameElastic::give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep);
    auto stress = dot(stiffnessMatrix, strain);
    status->letTempLatticeStrainBe(strain);
    status->letTempLatticeStressBe(stress);

    return stress;
}


Interface *
LatticeFrameElastic::giveInterface(InterfaceType type)
{
    return nullptr;
}


FloatMatrixF< 6, 6 >
LatticeFrameElastic::give3dFrameStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
    static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    double g = this->e / ( 2. * ( 1. + this->nu ) );

    const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
    const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
    const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();
    const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();
    const double shearareay = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaY();
    const double shearareaz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaZ();

    //Peter: You need to put here the correct values. Please check this.
    FloatArrayF< 6 >d = {
        this->e * area,
        g *shearareay,
        g *shearareaz,
        g *ik,
        this->e * iy,
        this->e * iz
    };

    return diag(d);
}
}
