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



#pragma once
#include "Contact/fecontactsurface.h"
#include "Contact/ContactElement/structuralcontactelement.h"
#include "set.h"


#define _IFT_StructuralFEContactSurface_Name "structuralfecontactsurface"
#define _IFT_StructuralFEContactSurface_contactElementSetNumber "sce_set"


namespace oofem {
class StructuralFEContactSurface : public FEContactSurface

/**
 * Structural Contact Surface made of Finite Elements
 *
 * Tasks:
 * keep the set of contact elements
 */

{
protected:
int sce_set;

public:
     StructuralFEContactSurface(int n, Domain *aDomain) : FEContactSurface(n, aDomain) {; }
    ~StructuralFEContactSurface() {};
    const char *giveClassName() const override { return "StructuralFEContactSurface"; }
    const char *giveInputRecordName() const override {return _IFT_StructuralFEContactSurface_Name;}

  //comment
protected:
  
  
};
} //end namespace oofem
