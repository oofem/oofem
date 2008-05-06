/* $Header: /home/cvs/bp/oofem/oofemlib/src/structtemperatureload.h,v 1.10 2003/04/06 14:08:26 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   *****************************************
//   *** CLASS STRUCTURAL TEMPERATURE LOAD ***
//   *****************************************


#ifndef structtemperatureload_h
#define structtemperatureload_h


#include "load.h"


class Element;
class TimeStep;


class StructuralTemperatureLoad : public Load
{
    /*
     * This class implements temperature (constant) load over the element
     * component array contains one or two numbers;
     * componentArray->at(1) contains increment of temperature in mid-surface
     * componentArray->at(2) contains increment of gradient of temperature
     *                     over the thickness of element (optional)
     *
     *
     *
     */
public:
    StructuralTemperatureLoad(int i, Domain *d) : Load(i, d) { }    // constructor

    /**
     * Computes components values of temperature field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value respecting load response mode.
     * @param answer component values at given point and time
     * @param stepN time step representing time
     * @param coords gp global coordinates, which are used to evaluate components values.
     * @param mode determines response mode.
     */
    virtual void         computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode)
    { this->computeComponentArrayAt(answer, tStep, mode); }

    classType    giveClassID() const { return StructuralTemperatureLoadClass; }
    bcValType    giveBCValType() const { return TemperatureBVT; }
    bcGeomType   giveBCGeoType() const { return BodyLoadBGT; }
    const char *giveClassName() const { return "StructuralTemperatureLoad"; }
};

#endif // structtemperatureload_h








