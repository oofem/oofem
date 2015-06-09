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

#ifndef NEUMANNMOMENTLOAD_H
#define NEUMANNMOMENTLOAD_H

#include "boundaryload.h"
#include "element.h"

#define _IFT_NeumannMomentLoad_Name "momentload"
#define _IFT_NeumannMomentLoad_Gradient "gradient"
#define _IFT_NeumannMomentLoad_Constant "constant"
#define _IFT_NeumannMomentLoad_CenterSet "cset"

/**
 * This class contains a Neumann type boundary condition given as
 * @f[ t=p+g\cdot[x-\bar{x}]\otimes n @f]
 * where @f$ p @f$ is a prescribed constant (eg pressure), @f$ g @f$ is the gradient (pressure gradient), @f$ x @f$ is the coordinate, @f$ \bar{x} @f$ is the
 * centre of the structure and $n$ is the outward pointing normal.
 */
namespace oofem {
class OOFEM_EXPORT NeumannMomentLoad : public BoundaryLoad
{
private:
    /// Center of structure
    FloatArray xbar;
    /// Set containing elements used to calculate xbar
    int cset;
    /// Array containing elements elements in set cset
    IntArray celements;
    /// Gradient
    FloatArray g;
    /// Constant
    double p;

    /// Compute centre of mass for set cset
    void computeXbar();

    /// Compute normal at center of element
    void computeNormal(FloatArray &answer, Element *e, int side);
public:
    NeumannMomentLoad(int i, Domain * d) : BoundaryLoad(i, d) { }

    // Overloaded methods:
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);
    virtual int giveApproxOrder() { return 0; }

    FormulationType giveFormulationType() { return FT_Global; }

    virtual void computeValueAtBoundary(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode, Element *e, int boundary);
    /**
     * Sets a new load vector.
     * @param newValue New load.
     */
    void updateLoad(const FloatArray &newValue) { componentArray = newValue; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual bcGeomType giveBCGeoType() const { return SurfaceLoadBGT; }

    virtual const char *giveClassName() const { return "NeumannMomentLoad"; }
    virtual const char *giveInputRecordName() const { return _IFT_NeumannMomentLoad_Name; }

private:
    virtual void computeNArray(FloatArray &answer, const FloatArray &coords) const { answer.clear(); }
};
} // end namespace oofem


#endif // NEUMANNMOMENTLOAD_H
