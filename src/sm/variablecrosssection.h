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

#ifndef variablecrosssection_h
#define variablecrosssection_h

#include "simplecrosssection.h"
#include "structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "parser.h"
#include <string>

///@name Input fields for SimpleCrossSection
//@{
#define _IFT_VariableCrossSection_Name "variablecs"
//@}

namespace oofem {

/**
 * Class implementing cross section model in finite element problem.
 * A cross section properties are not a constant values, but they can be set as 
 * expressions as a fuction of the spatial position in terms of x,y,z variables.
 * example "0.1+x*0.2+y*0.0", where (x,y,z) are global (or optionally local element) coordinates.
 * This cross section model is a integral model, so it does not perform any integration over cross-section volume.
 *
 * Note: varible cross section model uses the same keywords as simple cs model. As it is derived from simple cross section class,
 * the initializeFrom method does not calls the parent method, because the same keywords are used to read 
 * variables of different type.
 */
class OOFEM_EXPORT VariableCrossSection : public SimpleCrossSection
{
 protected:
  /// Expression parser
  Parser exprParser;
  /// Expression for cross section thickness
  std::string thicknessExpr;
  /// Expression for cross section width
  std::string widthExpr;
  /// Expression for cross section area
  std::string areaExpr;
  /// Expression for cross section inertia moment $I_y$
  std::string iyExpr;
  /// Expression for cross section inertia moment $I_z$
  std::string izExpr;
  /// Expression for cross section torsion moment $I_x$
  std::string ixExpr;
  /// Expression for cross section beam shear area $A_y$
  std::string shearAreayExpr;
  /// Expression for cross section beam shear area $A_z$
  std::string shearAreazExpr;
  /// Expression for cross section beam drilling stiffness
  std::string drillingStiffnessExpr;

  /// if set to true, all expressions are in element local cs, otherwise are expressed in global cs
  bool localFormulationFlag;
  
 public:
    /** 
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
     VariableCrossSection(int n, Domain *d) : SimpleCrossSection(n, d) {localFormulationFlag = false;}

    /**
     * Initializes receiver acording to object description stored in input record.
     * Calls CrossSection initializeFrom service and reads the values of
     * - 'thick' thickness
     * - 'width' width
     * - 'area' area
     * - 'iy' Moment of inertia around y
     * - 'iz' Moment of inertia around z
     * - 'ik' Torsion moment around x
     * - 'beamshearcoeff' Beam shear coefficient
     * @param ir Record to read off.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "VariableCrossSection"; }
    virtual const char *giveInputRecordName() const { return _IFT_VariableCrossSection_Name; }

    virtual double give(CrossSectionProperty a, GaussPoint *gp);
    virtual double give(CrossSectionProperty a, const FloatArray* coords, Element* elem, bool local);

 protected:
    std::string giveExpression (CrossSectionProperty aProperty);

};
} // end namespace oofem
#endif // variablecrosssection_h
