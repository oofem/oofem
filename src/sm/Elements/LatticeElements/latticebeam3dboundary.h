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

#ifndef latticebeam3dboundary_h
#define latticebeam3dboundary_h

#include "latticebeam3d.h"

///@name Input fields for Lattice2d
//@{
#define _IFT_LatticeBeam3dBoundary_Name "latticebeam3dboundary"
#define _IFT_LatticeBeam3dBoundary_location "location"
//@}

namespace oofem {
class ParamKey;
/**
 * This class implements a 3-dimensional lattice element for the boundaries of a periodic cell.
 */

class LatticeBeam3dBoundary : public LatticeBeam3d
{
protected:
    IntArray location;
    static ParamKey IPK_LatticeBeam3dBoundary_location;
public:
    LatticeBeam3dBoundary(int n, Domain *);
    virtual ~LatticeBeam3dBoundary();


    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 18; }

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;

    double computeVolumeAround(GaussPoint *) override;

    void giveDofManDofIDMask(int inode, IntArray &) const override;

    virtual void  giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0) override;
    virtual void computeGeometryProperties() override;

    virtual const char *giveInputRecordName() const override { return _IFT_LatticeBeam3dBoundary_Name; }
    virtual const char *giveClassName() const override { return "LatticeBeam3dBoundary"; }
    void initializeFrom(InputRecord &ir, int priority) override;
    void postInitialize() override;
    
    virtual Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;

    void giveVTKCoordinates(int nodeNumber, FloatArray &coords);

#ifdef __OOFEG
    virtual void drawYourself(oofegGraphicContext &context, TimeStep *tStep) override;
    virtual void drawRawGeometry(oofegGraphicContext &, TimeStep *tStep) override;
    virtual void drawDeformedGeometry(oofegGraphicContext &, TimeStep *tStep, UnknownType) override;
#endif

protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override {};
    bool computeGtoLRotationMatrix(FloatMatrix &) override;
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    void giveSwitches(IntArray &answer, int location);
};
} // end namespace oofem
#endif
