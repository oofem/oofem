/* $Header: /home/cvs/bp/oofem/sm/src/truss1d.h,v 1.8 2003/04/06 14:08:32 bp Exp $ */
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

//   *********************
//   *** CLASS TRUSS2D ***
//   *********************


#ifndef qtruss1dgrad_h
#define qtruss1dgrad_h

#include "fei1dlin.h"
#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "graddpelement.h"
#include "qtruss1d.h"

namespace oofem {
  class QTruss1dGrad : public QTruss1d, public GradDpElement
{
    /*
     * This class implements a two-node truss bar element for two-dimensional
     * analysis.
     * DESCRIPTION :
     * A truss bar element is characterized by its 'length' and its 'pitch'. The
     * pitch is the angle in radians between the X-axis anf the axis of the
     * element (oriented node1 to node2).
     * The 'rotationMatrix' R is such that u{loc}=R*u{glob}.
     * Note: element is formulated in global c.s.
     * TASKS :
     * - calculating its Gauss points ;
     * - calculating its B,D,N matrices and dV ;
     * - expressing M,K,f,etc, in global axes. Methods like 'computeStiffness-
     *   Matrix' of class Element are here overloaded in order to account for
     *   rotational effects.
     */

protected:
    double length;
    static FEI1dLin interpolation;
    int numberOfGaussPoints;
    // FloatMatrix*  rotationMatrix ;

public:
    QTruss1dGrad(int, Domain *);                         // constructor
    ~QTruss1dGrad()   { }                                // destructor

  
 
    const char *giveClassName() const { return "QTruss1dGrad"; }
    classType            giveClassID() const { return QTruss1dGradClass; }
    //IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }
    integrationDomain  giveIntegrationDomain() { return _Line; }
    MaterialMode          giveMaterialMode()  { return _1dMatGrad; }
    IRResultType initializeFrom(InputRecord *ir);
   
     int getNprimNodes(){return 3;}
     int getNprimVars(){return 1;}
     int getNsecNodes(){return 2;}
     int getNsecVars(){return 1;}
protected:
    virtual void          computeBkappaMatrixAt(GaussPoint *, FloatMatrix &);
    virtual void          computeNkappaMatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    virtual void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void          computeGaussPoints();
    void          giveDofManDofIDMask(int inode, EquationID,IntArray &) const;
    //double        giveLength();
    //double        givePitch();
    StructuralElement* giveStructuralElement(){return this;}
    //    int           giveApproxOrder() { return 1; }
};
} // end namespace oofem
#endif // truss1d_h
