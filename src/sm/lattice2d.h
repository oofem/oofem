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

#ifndef lattice2d_h
#define lattice2d_h

#include "latticestructuralelement.h"

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

public:
    Lattice2d(int n, Domain *d);
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
    virtual double giveCharacteristicLenght(GaussPoint *, const FloatArray &)
    { return this->giveLength(); }

    virtual double giveLength();

    virtual double giveArea() { return this->width * this->thickness; }

    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int giveCrackFlag();

    virtual double giveCrackWidth();
    virtual double giveDissipation();
    virtual double giveDeltaDissipation();
    //
    // definition & identification
    //
    virtual const char *giveClassName() const { return "Lattice2d"; }
    virtual classType giveClassID() const { return Lattice2dClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context);
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawSpecial(oofegGraphicContext &gc);
#endif

protected:

    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual int giveNumberOfCrossSectionNodes() { return 2; }
    double givePitch();
    virtual void computeGaussPoints();
    virtual integrationDomain giveIntegrationDomain() { return _Line; }
};
} // end namespace oofem
#endif
