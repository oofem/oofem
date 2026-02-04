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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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
#include "Contact/ContactElements/thermalcontactelement.h"
#include "set.h"

#define _IFT_ThermalFEContactSurface_Name "thermalfecontactsurface"
#define _IFT_ThermalFEContactSurface_contactElementSetNumber "tce_set"

namespace oofem {
class ThermalFEContactSurface : public FEContactSurface
{
protected:
  int tce_set;

  
public:
    ThermalFEContactSurface(int n, Domain *aDomain) : FEContactSurface(n, aDomain) {; }
    ~ThermalFEContactSurface() {};

  const char *giveClassName() const override { return "ThermalFEContactSurface"; }
  const char *giveInputRecordName() const override {return _IFT_ThermalFEContactSurface_Name;}
protected:
  std::unique_ptr<ThermalContactElement> createContactElement(int iElem, Domain *d, IntArray nodes, FEInterpolation &interpolation, int boundaryNumber);

  
};
} //end namespace oofem
