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

#ifndef libeam3dboundaryvoigt_h
#define libeam3dboundaryvoigt_h

#include "sm/Elements/Beams/libeam3dboundary.h"

///@name Input fields for LIBeam3dBoundaryVoigt
//@{
#define _IFT_LIBeam3dBoundaryVoigt_Name "libeam3dboundaryvoigt"
//@}

namespace oofem {
/**
 * This class implements a boundary version of the 3-dimensional mindlin theory Linear Isoparametric
 * beam element, with reduced integration. Useful for prescribing periodicity in multiscale analyses.
 * MACROSCOPIC INPUT: STRAIN TENSOR (3D, 6 COMPONENTS, VOIGT NOTATION: Exx Eyy Ezz Gyz Gxz Gxy)
 *
 * @author: Adam Sciegaj
 */
class LIBeam3dBoundaryVoigt : public LIBeam3dBoundary
{
public:
    LIBeam3dBoundaryVoigt(int n, Domain *d);
    virtual ~LIBeam3dBoundaryVoigt() { }

    void initializeFrom(InputRecord &ir) override;

    int computeNumberOfDofs() override { return 18; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LIBeam3dBoundaryVoigt_Name; }
    const char *giveClassName() const override { return "LIBeam3dBoundaryVoigt"; }

protected:
    void computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // libeam3dboundaryvoigt_h
