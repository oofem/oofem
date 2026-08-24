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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#ifndef lattice3dboundary_mt_h
#define lattice3dboundary_mt_h

#include "lattice3d_mt.h"
#define _MYPI 3.14159265359

///@name Input fields for lattice3dboundary_mt
//@{
#define _IFT_Lattice3dboundary_mt_Name "latticemt3dboundary"
#define _IFT_Lattice3dboundary_mt_location "location"
//@}

namespace oofem {
/**
 * This class implements a 3-dimensional lattice mass transport element
 */

class Lattice3dboundary_mt : public Lattice3d_mt
{
protected:
    IntArray location;
    int numberOfLocations;

public:

    // constructor
    Lattice3dboundary_mt(int, Domain *, ElementMode em = HeatTransferEM);
    ~Lattice3dboundary_mt();                       // destructor

    void computeBCSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx);

    void computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint);

    double computeVolumeAround(GaussPoint *);

    void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rmode, TimeStep *tStep);

    void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep);

    ///a function that gives the location array of the boundary element
    void giveLocationArray(IntArray &locArray) { locArray = location; }

    void giveVTKCoordinates(int nodeNumber, FloatArray &coords);

    virtual void computeInternalForcesVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);

    virtual void computeHomogenisedInternalForcesVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, FloatArray &unknowns);

    virtual void computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_Lattice3dboundary_mt_Name; }
    const char *giveClassName() const { return "Lattice3dboundary_mt"; }

    virtual void computeGeometryProperties();

    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

    virtual int            computeNumberOfDofs() { return 5; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;

#ifdef __OOFEG

    virtual void  drawRawGeometry(oofegGraphicContext &, TimeStep *tStep);

#endif

protected:

    void giveSwitches(IntArray &answer, int location);
};
} // end namespace oofem
#endif
