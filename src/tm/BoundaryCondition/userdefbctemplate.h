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

#include "activebc.h"
#include "floatarray.h"

#include <memory>

#define _IFT_UserDefinedBC_Name "usrdefbc" // change to preferred keyword
#define _IFT_UserDefinedBC_nlocations "nlocations"
#define _IFT_UserDefinedBC_location "location"

namespace oofem {
class Node;
class Element;

/**
 * Template illustrating implementation of user defined boundary condition. The boundary condition can contibute to RHS as well as LHS terms.
 * 
 * @note Experimental code
 * 
 * @author Borek Patzak
 */
class OOFEM_EXPORT UserDefinedBC : public ActiveBoundaryCondition //, public Homogenization
{
public:
    UserDefinedBC(int n, Domain *d);

    // initialize receiver from input record
    void initializeFrom(InputRecord &ir) override;
    // returns input record of receiver
    void giveInputRecord(DynamicInputRecord &input) override;
    void postInitialize() override;

    bcType giveType() const override { return UnknownBT; }

    // Assembles B.C. contributions to specified vector (typically rhs term)
    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr, void*lock=nullptr) override;

    // Assembles B.C. contributions to specified matrix (typically lhs contribution)
    void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s,
                  const UnknownNumberingScheme &c_s, double scale=1.0, void*lock=nullptr) override;
    // Gives a list of location arrays that will be assembled.
    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                            const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

    const char *giveClassName() const override { return "UserDefinedBC"; }
    const char *giveInputRecordName() const override { return _IFT_UserDefinedBC_Name; }

protected:
    int nlocations; // number of sampling locations
    std::vector <FloatArray> locations;
};
} /* namespace oofem */
