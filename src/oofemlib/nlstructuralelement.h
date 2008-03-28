/* $Header: /home/cvs/bp/oofem/oofemlib/src/nlstructuralelement.h,v 1.13 2003/04/06 14:08:25 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/


//   *******************************************
//   *** CLASS NON LINEAR STRUCTURAL ELEMENT *** 
//   *******************************************
 

#ifndef nlstructuralelement_h
#define nlstructuralelement_h

#include "structuralelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"

class TimeStep ; class Node ; class Material ; class GaussPoint ;
class FloatMatrix ; class FloatArray ; class IntArray ;


/**
 Abstract base class for "structural" finite elements with geometrical nonlinearities.
 The base \ref StructuralElement supports nonlinearities at constitutive level.
 The new implementation of relevant services is extended to incorporate general geometric
 nonlinerarity. The services for computing nonlinear parts of geometrical equations for
 particular strain components are provided. The updated services for stiffness, 
 strain and internal forces computations are formulated.

 The activation of non-linear effects can be generally controlled on element level
 using \ref nlGeometry, allowing elements to be used also in linear computations.
*/
class NLStructuralElement : public StructuralElement
{
/*
 IMPORTANT :

   methods ComputeStressVector() and ComputeStrainVector() as they are implemented here
   assume Total Lagrange approach.


 TASKS :
   -defining itself :
   -calculating its contribution to the problem :
     .calculating its mass matrix M, its stiffness matrix K, its load vector
      f, its location array ;
     .calculating its contribution to the LHS and RHS of the linear system,
      using Static,Newmark,etc, formula. These contributions are usually
      combinations of M,K,f.
   -performing end-of-step operations :
     .calculating the strains and stresses at its Gauss points ;
     .printing its output in the data file and updating itself ;
*/
protected:
  /// flag indicating if geometrical nonlineraities apply.
    int nlGeometry;
public:
   /**
    Constructor. Creates element with given number, belonging to given domain.
    @param n element number
    @param d domain to which new material will belong
    */
   NLStructuralElement (int,Domain*) ;              // constructors
     /// Destructor
     ~NLStructuralElement () {}                       // destructor


      // characteristic  matrix

      /**
    Computes the stiffness matrix of receiver.
    The response is evaluated using \f$\int (B_1+B_2(r))^TD(B_1+B_2(r)) dv\f$, where 
    \f$B_2\f$ is nonlinear contribution evaluated using computeNLBMatrixAt service for each strain component 
    (\f$B_2(i) = \Delta r^T A(i)\f$)
    Necessary transformations and reduced integration are taken into account.
    @param answer computed stiffness matrix.
    @param rMode response mode
    @param tStep time step
   */
      virtual void  computeStiffnessMatrix (FloatMatrix& answer, 
                      MatResponseMode rMode, TimeStep* tStep) ;

    // stress equivalent vector (vector of internal forces) - for nonLinear Analysis.
    /**
   Evaluates nodal representation of real internal forces.
   The response is evaluated using $F = \int (B+B_2)^T\sigma dV$ formula, where
   B is linear strain-displacement contribution, $B_2$ is nonlinear contribution evaluated using
   computeNLBMatrixAt service for each strain component ($B_2(i) = \Delta r^T A(i)$).  
   Necessary transformations are taken into account.
   @param answer equivalent nodal forces vector
   @param tStep time step
   @param useUpdatedGpRecord if equal to zero, the stresses in integration points are computed (slow but safe), else if
   */
    virtual void giveInternalForcesVector (FloatArray& answer, 
                      TimeStep*, int useUpdatedGpRecord = 0) ;

    /**
   Compute strain vector of receiver evaluated at given integration point at time
    step stepN from element displacement vector. 
    Total (green) strains are computed using following scheme
   \f$\varepsilon_G = Br+B_2(r)\f$, where \f$B_2\f$ is obtained using 
   computeNLBMatrixAt service for each strain component (\f$B_2(i) = \Delta r^T A(i)\f$)
   and \f$B\f$ is usual linear strain-displacement matrix.
   Necessary transformation are taken into account.
   @param answer element strain vector
   @param gp integration point
   @param stepN time step
   */
   void   computeStrainVector (FloatArray& answer, GaussPoint*,TimeStep*) ; 
    
  // data management
     /**
    Initializes receiver acording to object description stored in input record.
    The nlgeo flag is read and StructuralElement instanciateFromString service called.
    @param corresponding data record of receiver
   */
     IRResultType initializeFrom (InputRecord* ir);

      // time step termination
//    virtual void          printOutputAt (TimeStep*) ;
//    void                  updateYourself (TimeStep*) ;

      // definition
   /// Returns "NLStructuralElement" - class name of the receiver.
   const char* giveClassName () const { return "NLStructuralElement" ;}
    /// Returns NLStructuralElementClass - class name of the receiver.
  classType                giveClassID () const
   { return NLStructuralElementClass; }
/*      int                   giveNumber ()
          { return FEMComponent::giveNumber() ;} */

protected:
  // void   computeStressVector (FloatArray& answer, GaussPoint*,TimeStep*) ;   // real one
  /**
  Computes the nonlinear part of strain-displacement (geometrical) equation 
  related to i-th component of strain vector.
  @param answer returned nonlinear strain vector component contribution
  @param gp integration point
  @param i determines the component of strain vector for which contribution is assembled
  @see computeStrainVector service
  */
  virtual void computeNLBMatrixAt (FloatMatrix& answer, GaussPoint*, int ) {
  answer.resize (0,0); return ;}

} ;


#endif // nlstructuralelement_h








