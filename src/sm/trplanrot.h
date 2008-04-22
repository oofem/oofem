/* $Header: /home/cvs/bp/oofem/sm/src/trplanrot.h,v 1.4 2003/04/06 14:08:32 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *********************************************************
//   *** CLASS PLANE STRAIN WITH INDEPENDENT ROTATION FIELD***
//   *********************************************************
// 5.5.1995
//
#ifndef trplanrot_h
#define trplanrot_h


#include "trplanstrss.h"

class TrPlaneStrRot : public TrPlaneStress2d
{
    /*
     * This class implements an triangular three-node  plane-
     * stress elasticity finite element with independent rotation field.
     * Each node has 3 degrees of freedom.
     *
     * DESCRIPTION :
     *
     * TASKS :
     *
     * - calculating its B,D,N matrices and dV.
     */

protected:

    int numberOfRotGaussPoints;

public:

    TrPlaneStrRot(int, Domain *);            // constructor
    ~TrPlaneStrRot();                        // destructor


    FloatArray *GivePitch();
    FloatArray *GiveDerivativeUX(GaussPoint *);
    FloatArray *GiveDerivativeVX(GaussPoint *);
    FloatArray *GiveDerivativeUY(GaussPoint *);
    FloatArray *GiveDerivativeVY(GaussPoint *);
    void           computeStrainVector(FloatArray &answer, GaussPoint *, TimeStep *);

    virtual int            computeNumberOfDofs(EquationID ut) { return 9; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    // characteristic length in gp (for some material models)
    double        giveCharacteristicLenght(GaussPoint *, const FloatArray &) { return 0.; }

    virtual int testElementExtension(ElementExtension ext) { return 0; }
    //int    hasEdgeLoadSupport () {return 0;}


    // definition & identification

    MaterialMode          giveMaterialMode() { return _Unknown; }
    const char *giveClassName() const { return "TrPlaneStrRot"; }
    classType      giveClassID() const { return TrPlaneStrRotClass; }
    IRResultType initializeFrom(InputRecord *ir);

protected:
    void           computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    void           computeNmatrixAt(GaussPoint *, FloatMatrix &);
    void           computeGaussPoints();
    integrationDomain  giveIntegrationDomain() { return _Triangle; }
};

#endif //  trplanrot_h
