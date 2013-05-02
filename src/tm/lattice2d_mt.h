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

#ifndef lattice2d_mt_h
#define lattice2d_mt_h

#include "latticetransportelement.h"
#include "spatiallocalizer.h"

///@name Input fields for Lattice2d_mt
//@{
#define _IFT_Lattice2DMT_dim "dim"
#define _IFT_Lattice2DMT_thick "thick"
#define _IFT_Lattice2DMT_width "width"
#define _IFT_Lattice2DMT_gpcoords "gpcoords"
#define _IFT_Lattice2DMT_crackwidth "crackwidth"
#define _IFT_Lattice2DMT_couplingflag "couplingflag"
#define _IFT_Lattice2DMT_couplingnumber "couplingnumber"
//@}


namespace oofem {
/**
 * This class implements a 2-dimensional lattice mass transport element
 */

class Lattice2d_mt : public LatticeTransportElement
{
protected:
    int numberOfGaussPoints;
    double area;
    double length;


    double dimension, width, thickness;
    FloatArray gpCoords;

    double crackWidth;

    int couplingFlag, couplingNumber;

public:

    // constructor
    Lattice2d_mt(int, Domain *, ElementMode em = HeatTransferEM);
    ~Lattice2d_mt();                       // destructor

    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void          computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    double                computeVolumeAround(GaussPoint *);

    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);

    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep);

    const char *giveClassName() const { return "Lattice2d_mtElement"; }
    classType                giveClassID() const { return Lattice2d_mtClass; }

    Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

    virtual int            computeNumberOfDofs(EquationID ut) { return 2; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    IRResultType           initializeFrom(InputRecord *ir);

    void updateInternalState(TimeStep *stepN);

#ifdef __OOFEG
    //
    // Graphics output
    //

    void drawYourself(oofegGraphicContext &gc);

    void  drawRawGeometry(oofegGraphicContext &);

    void drawRawCrossSections(oofegGraphicContext &gc);

    void  giveCrossSectionCoordinates(FloatArray &coords);

#endif

protected:
    void                  computeGaussPoints();

    virtual void  computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *);
    virtual void  computeNmatrixAt(FloatMatrix &n, FloatArray *);

    virtual double giveCrackFactor();

    virtual double givePressure();

    virtual double giveMass();

    virtual int giveCouplingFlag() { return couplingFlag; }

    virtual int giveCouplingNumber() { return couplingNumber; }

    /* computes the submatrix of interpolation matrix cooresponding to single unknown.
     */
    virtual void  computeNSubMatrixAt(FloatMatrix &n, FloatArray *);

    virtual double giveLength();

    virtual double giveArea() { return width * thickness; }

    virtual double giveWidth() { return width; }

    virtual double giveCrackWidth() { return this->crackWidth; }

    virtual void  giveGpCoordinates(FloatArray &coords);

    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) { return 0; }

    int giveApproxOrder(int unknownIndx) { return 1; }
};
} // end namespace oofem
#endif
