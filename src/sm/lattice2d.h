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

#ifndef lattice2d_h
#define lattice2d_h

#include "latticestructuralelement.h"

///@name Input fields for Lattice2d
//@{
#define _IFT_Lattice2d_Name "lattice2d"
#define _IFT_Lattice2d_thick "thick"
#define _IFT_Lattice2d_width "width"
#define _IFT_Lattice2d_gpcoords "gpcoords"
#define _IFT_Lattice2d_couplingflag "couplingflag"
#define _IFT_Lattice2d_couplingnumber "couplingnumber"
//@}

namespace oofem {
/**
 * This class implements a 2-dimensional lattice element
 */

class Lattice2d : public LatticeStructuralElement
{
protected:
    double kappa, pitch, length;

    double width, thickness;
    FloatArray gpCoords;
    int couplingFlag, couplingNumber;

public:
    Lattice2d(int n, Domain * d);
    virtual ~Lattice2d();

    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);

    /**
     * This function is different from the standard computeGlobalCorrdinates
     * function as it returns the global coordinates of the gausspoint
     * independent to the value of the lcoords.
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    /**
     * This function returns the length of the element
     * independent of the FloatArray.
     */
    virtual double giveCharacteristicLength(const FloatArray &)
    { return this->giveLength(); }

    virtual double giveLength();

    virtual double giveNormalStress();

    virtual double giveArea() { return this->width * this->thickness; }

    virtual int computeNumberOfDofs() { return 6; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int giveCrackFlag();

    virtual double giveCrackWidth();
    virtual double giveDissipation();
    virtual double giveDeltaDissipation();

    virtual int giveCouplingFlag() { return couplingFlag; }

    virtual int giveCouplingNumber() { return couplingNumber; }
    //
    // definition & identification
    //
    virtual const char *giveInputRecordName() const { return _IFT_Lattice2d_Name; }
    virtual const char *giveClassName() const { return "Lattice2d"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context);
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawSpecial(oofegGraphicContext &gc);
    void drawRawCrossSections(oofegGraphicContext &gc);
    void giveCrossSectionCoordinates(FloatArray &coords);
    void drawDeformedCrossSections(oofegGraphicContext &, UnknownType);
#endif

protected:

    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual int giveNumberOfCrossSectionNodes() { return 2; }
    double givePitch();
    virtual void computeGaussPoints();
    virtual integrationDomain giveIntegrationDomain() const { return _Line; }
    virtual void  giveGpCoordinates(FloatArray &coords);
};
} // end namespace oofem
#endif
