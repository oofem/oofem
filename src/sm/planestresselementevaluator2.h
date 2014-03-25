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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef planestresselementevaluator2_h
#define planestresselementevaluator2_h

#include "structuralelementevaluator2.h"

namespace oofem {


class ThreeDimensionalStructuralElementEvaluator : public StructuralElementEvaluator2
{
public:
  ThreeDimensionalStructuralElementEvaluator(int);
  ThreeDimensionalStructuralElementEvaluator();
private:
  IntArray *locationArray;
protected:

   /**
     * Computes derivative of the interpolation matrix for element unknowns.
     * @param gp Integration point for which answer is assembled.
     * @param answer Matrix containing derivatives of the interpolation matrix evaluated at gp.
     */
  void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,int lowerIndx = 1, int upperIndx = ALL_STRAINS);

  virtual double computeVolumeAround(GaussPoint *gp,ElementGeometry* elem);
  virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
  {
	  answer.setValues(3, D_u, D_v,D_w);
  }

  virtual int giveNumberOfDofMans();
  

}; // end of ThreeDimensionalStructuralElementEvaluator definition


class AxisymmetricStructuralElementEvaluator : public StructuralElementEvaluator2
{
public:
  AxisymmetricStructuralElementEvaluator(int);
  AxisymmetricStructuralElementEvaluator();
private:
  IntArray *locationArray;
protected:

   /**
     * Computes derivative of the interpolation matrix for element unknowns.
     * @param gp Integration point for which answer is assembled.
     * @param answer Matrix containing derivatives of the interpolation matrix evaluated at gp.
     */
  void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,int lowerIndx = 1, int upperIndx = ALL_STRAINS);

  virtual double computeVolumeAround(GaussPoint *gp,ElementGeometry* elem);
  virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
  {
    answer.setValues(2, D_u, D_v);    
  }
  
  virtual int giveNumberOfDofMans();

 

}; // end of ThreeDimensionalStructuralElementEvaluator definition

/**
 * General purpose Plane stress structural element evaluator.
 */
class PlaneStressStructuralElementEvaluator2 : public StructuralElementEvaluator2
{
public:
  PlaneStressStructuralElementEvaluator2(int);
  PlaneStressStructuralElementEvaluator2();
private:
  IntArray *locationArray;
protected:

   /**
     * Computes derivative of the interpolation matrix for element unknowns.
     * @param gp Integration point for which answer is assembled.
     * @param answer Matrix containing derivatives of the interpolation matrix evaluated at gp.
     */
  void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,int lowerIndx = 1, int upperIndx = ALL_STRAINS);

  virtual double computeVolumeAround(GaussPoint *gp,ElementGeometry* elem);
  virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
  {
    answer.setValues(2, D_u, D_v);    

    }
  virtual int giveNumberOfDofMans();
 
 
}; // end of PlaneStressStructuralElementEvaluator definition



class PlaneStrainStructuralElementEvaluator : public StructuralElementEvaluator2
{
public:
  PlaneStrainStructuralElementEvaluator(int);
  PlaneStrainStructuralElementEvaluator();
private:
  IntArray *locationArray;
protected:

   /**
     * Computes derivative of the interpolation matrix for element unknowns.
     * @param gp Integration point for which answer is assembled.
     * @param answer Matrix containing derivatives of the interpolation matrix evaluated at gp.
     */
  void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,int lowerIndx = 1, int upperIndx = ALL_STRAINS);

  virtual double computeVolumeAround(GaussPoint *gp,ElementGeometry* elem);
  virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
  {
	  answer.setValues(2, D_u, D_v);    
  }
  virtual int giveNumberOfDofMans();
}; // end of PlaneStrainStructuralElementEvaluator definition


class OneDimensionalStructuralElementEvaluator : public StructuralElementEvaluator2
{
public:
  OneDimensionalStructuralElementEvaluator(int);
  OneDimensionalStructuralElementEvaluator();
private:
  IntArray *locationArray;
protected:

   /**
     * Computes derivative of the interpolation matrix for element unknowns.
     * @param gp Integration point for which answer is assembled.
     * @param answer Matrix containing derivatives of the interpolation matrix evaluated at gp.
     */
  void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,int lowerIndx = 1, int upperIndx = ALL_STRAINS);

  virtual double computeVolumeAround(GaussPoint *gp,ElementGeometry* elem);
  virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
  {
     answer.setValues(1, D_u);    
  }
 

}; // end of OneDimensionalStructuralElementEvaluator definition










} // end namespace oofem
#endif //planestresselementevaluator2_h
