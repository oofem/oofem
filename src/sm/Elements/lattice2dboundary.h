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


#ifndef lattice2dboundary_h
#define lattice2dboundary_h

#include "../sm/Elements/lattice2d.h"

///@name Input fields for Lattice2d
//@{
#define _IFT_Lattice2dBoundary_Name "latticeboundary2d"
#define _IFT_Lattice2dBoundary_location "location"
#define TOL 1.e-8
//@}

namespace oofem {
class Lattice2dBoundary : public Lattice2d
{
    /*
     *   This class implements a 2-dimensional lattice element
     */
protected:
    int location;

public:
    Lattice2dBoundary(int, Domain *);                     // constructor
    ~Lattice2dBoundary();                                 // destructor

    virtual void giveInternalForcesVector(FloatArray &answer,
                                           TimeStep *, int useUpdatedGpRecord = 0) override;

    int giveLocation() { return this->location; }

    virtual int            computeNumberOfDofs() override { return 9; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const override;
    virtual double        computeVolumeAround(GaussPoint *) override;

    //
    // definition & identification
    //
    virtual const char *giveInputRecordName() const override { return _IFT_Lattice2dBoundary_Name; }
    virtual const char *giveClassName() const override { return "Lattice2dBoundary"; }
    virtual IRResultType initializeFrom(InputRecord *ir) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    
#ifdef __OOFEG
    void          drawYourself(oofegGraphicContext &context,TimeStep *tStep);
    void          drawRawGeometry(oofegGraphicContext &,TimeStep *tStep);
    void          drawRawCrossSections(oofegGraphicContext &,TimeStep *tStep);
    void drawDeformedGeometry(oofegGraphicContext &,  TimeStep *tStep, UnknownType);
    void drawSpecial(oofegGraphicContext &gc,TimeStep *tStep);
    virtual void giveCrossSectionCoordinates(FloatArray &coords);
#endif


protected:
    void          computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS) override;
    virtual bool           computeGtoLRotationMatrix(FloatMatrix &) override;
    
    void  computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN) override;

    virtual double        giveLength() override;
    double        givePitch();
    void giveSwitches(FloatArray &answer);

};
} // end namespace oofem

#endif
