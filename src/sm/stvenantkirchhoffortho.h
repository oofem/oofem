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

/*
 * stvenantkirchhoffortho.h
 *
 *		Description: Material model for transversely isotropic St.Venant-Kirchhoff elasticity
 *
 *		Based on:
 *		J. Bonet and A.J. Burton: A simple orthotropic, transversely isotropic
 *		hyperelastic constitutive equation for large strain computations.
 *		Comput. Methods Appl. Mech. Engrg. 162 (1998) pp 151-164
 *
 *		Note the following two errors in that paper:
 *		Eqn 55g reads n = EA/E. Should (in my opinion) read n = E/EA.
 *		Eqn 52a reads ... -alpha*A. Should (in my opinion) read -2*alpha*A.
 *
 *  	Created on: May 28, 2013
 *      Author: Erik Svenning
 *
 */

#ifndef STVENANTKIRCHHOFFORTHO_H_
#define STVENANTKIRCHHOFFORTHO_H_

#include "structuralmaterial.h"
#include "structuralms.h"

///@name Input fields for StVenantKirchhoffOrtho
//@{
#define _IFT_StVenantKirchhoffOrtho_Name "stvenantkirchhoffortho"
#define _IFT_StVenantKirchhoffOrtho_ex "ex"
#define _IFT_StVenantKirchhoffOrtho_ey "ey"
#define _IFT_StVenantKirchhoffOrtho_ez "ez"
#define _IFT_StVenantKirchhoffOrtho_E_i "ei"
#define _IFT_StVenantKirchhoffOrtho_nu_i "nui"
#define _IFT_StVenantKirchhoffOrtho_EA "ea"
#define _IFT_StVenantKirchhoffOrtho_GA "ga"
//@}

namespace oofem {

/**
 * This class implements associated MaterialStatus for StVenantKirchhoffOrtho.
 */
class StVenantKirchhoffOrthoMaterialStatus : public StructuralMaterialStatus
{
public:
    /// Constructor
    StVenantKirchhoffOrthoMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~StVenantKirchhoffOrthoMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "StVenantKirchhoffOrthoMaterialStatus"; }
    virtual classType giveClassID() const { return StVenantKirchhoffOrthoMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
};

class StVenantKirchhoffOrtho : public StructuralMaterial
{
protected:
    FloatArray VecA; // Fiber direction
    double E_i, nu_i, EA, GA;

public:
    StVenantKirchhoffOrtho(int n, Domain *d);
    virtual ~StVenantKirchhoffOrtho();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int hasMaterialModeCapability(MaterialMode);
    virtual const char *giveInputRecordName() const { return _IFT_StVenantKirchhoffOrtho_Name; }
    virtual const char *giveClassName() const { return "StVenantKirchhoffOrtho"; }
    virtual classType giveClassID() const { return StVenantKirchhoffOrthoMaterialClass; }

};

} // end namespace oofem

#endif /* STVENANTKIRCHHOFFORTHO_H_ */
