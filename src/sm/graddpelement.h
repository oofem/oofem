/* $Header: /home/cvs/bp/oofem/sm/src/QSpace.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   ************************************
//   *** Abstract class for gradient formulation of coupled damage-plasticity model(GradDp). Yield function is formulated in the effective stress space and damage is driven by the nonlocal(over-nonlocal) cumulated plastic strain. The new nonlocal degrees of freedom (with the meaning of the nonlocal cumulated plastic strain) are introduced with lower order of approximation functions than the displacement field to avoid spurious stress oscillations ***
//   ************************************

#ifndef graddpelement_h
#define graddpelement_h

#include "structuralelement.h"



namespace oofem {

class GradDpElement 
{
 protected:
  int numberOfGaussPoints;
  int nPrimNodes, nPrimVars,nSecNodes,nSecVars;
  IntArray locU,locK;
  int totalSize,nlSize,locSize;
  

 public:
  GradDpElement () ;     // constructor
  ~GradDpElement ()  {}             // destructor

  //IRResultType initializeFrom (InputRecord* ir);
  //virtual void giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const;
  //
  // definition & identification
  //
  const char* giveClassName () const { return "GradDpElement"; }
  classType   giveClassID   () const { return GradDpElementClass; }
  /***********************Predelat************************************/
  virtual int getNprimNodes(){return 0;}
  virtual int getNprimVars(){return 0;}
  virtual int getNsecNodes(){return 0;}
  virtual int getNsecVars(){return 0;}

  /************************************************************/
  //  virtual int  computeNumberOfDofs (EquationID ut){return (nPrimNodes*nPrimVars+nSecNodes*nSecVars);};
  //virtual int    checkConsistency();

 protected:
  virtual StructuralElement* giveStructuralElement() = 0;
  virtual void computeNkappaMatrixAt(GaussPoint*,FloatMatrix&) = 0;
  virtual void computeBkappaMatrixAt(GaussPoint*,FloatMatrix&) = 0;
  //void initialize();
  void setDisplacementLocationArray(IntArray& answer,int nPrimNodes,int nPrimVars,int nSecNodes,int nSecVars);
  void setNonlocalLocationArray(IntArray& answer,int nPrimNodes,int nPrimVars,int nSecNodes,int nSecVars);

  void computeStiffnessMatrix(FloatMatrix&, MatResponseMode, TimeStep*);
  void computeStiffnessMatrix_uu(FloatMatrix&,MatResponseMode,TimeStep*);
  void computeStiffnessMatrix_uk(FloatMatrix&,MatResponseMode,TimeStep*);
  void computeStiffnessMatrix_kk(FloatMatrix&,MatResponseMode,TimeStep*);
  void computeStiffnessMatrix_ku(FloatMatrix&,MatResponseMode,TimeStep*);

  void computeDisplacementDegreesOfFreedom(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
  void computeNonlocalDegreesOfFreedom(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);

  void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
  void computeLocalStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
  void computeNonlocalCumPlasticStrain(double &answer, GaussPoint *gp, TimeStep *stepN);
  
  void giveInternalForcesVector(FloatArray &answer,TimeStep *tStep, int useUpdatedGpRecord);
  void giveLocalInternalForcesVector(FloatArray &answer,TimeStep *tStep, int useUpdatedGpRecord);
  void giveNonlocalInternalForcesVector(FloatArray &answer,TimeStep *tStep, int useUpdatedGpRecord);
  
  void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
  void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
  void computeLocForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
  void computeLocNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
  void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
///////////////////////////////////////////////////////////////////////////////
};
 
} // end namespace oofem

#endif
