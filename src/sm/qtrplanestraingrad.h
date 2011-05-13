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


#ifndef qtrplanestraingrad_h
#define qtrplanestraingrad_h
#include "structuralelement.h"
#include "gaussintegrationrule.h"
#include "graddpelement.h"
#include "qtrplanestrain.h"

namespace oofem {
  class QTrPlaneStrainGrad : public QTrPlaneStrain,public GradDpElement
{
 
 protected:
 
  int numberOfGaussPoints;
  // FloatMatrix*  rotationMatrix ;

public:
    QTrPlaneStrainGrad(int, Domain *);                         // constructor
    ~QTrPlaneStrainGrad()   { }                                // destructor

  
 
    const char *giveClassName() const { return "QTrPlaneStrainGrad"; }
    classType            giveClassID() const { return QTrPlaneStrainGradClass; }
    //IRResultType initializeFrom(InputRecord *ir);
    IRResultType initializeFrom(InputRecord *ir);
  
protected:
    virtual void          computeBkappaMatrixAt(GaussPoint *, FloatMatrix &);
    virtual void          computeNkappaMatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep){GradDpElement ::computeStiffnessMatrix(answer, rMode,tStep);}
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) {GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode){GradDpElement :: computeForceLoadVector(answer, stepN,mode);}
    virtual void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode){GradDpElement :: computeNonForceLoadVector(answer,stepN,mode);}
    virtual int  computeNumberOfDofs(EquationID ut) { return 15; }	
    void          computeGaussPoints();
    void          giveDofManDofIDMask(int inode, EquationID,IntArray &) const;
    StructuralElement* giveStructuralElement(){return this;}

};
} // end namespace oofem
#endif // qtrplanestraingrad_h
