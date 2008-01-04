/* $Header: /home/cvs/bp/oofem/oofemlib/src/ppdeelement.h,v 1.12 2003/04/06 14:08:25 bp Exp $ */
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

//   ************************************************
//   *** CLASS PARABOLIC PARTIAL EQUATION ELEMENT ***
//   ************************************************
 

#ifndef ppdeelement_h


#include "element.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"

class TimeStep ; class Node ; class Material ; class GaussPoint ;
class FloatMatrix ; class FloatArray ; class IntArray ;


class PPdeElement : public Element
{
/*
   This abstract class is the most important class of the program. It is the
   superclass of all classes implementing structural finite elements 
   (bar, shell, etc).
   An element is an attribute of a domain.
 DESCRIPTION :
   The basic data of an element are the numbers of its 'numberOfNodes' nodes,
   stored in 'nodeArray', of its 'material', of its body loads (eg, the dead
   weight) stored in 'loadArray'. These data are obtained from the domain.
   The element contains logical reference to Volume object (
   which stores material,geometrical characteristics, gaussPoints, state of 
   material and so on)
   The calculated data of an element are its 'massMatrix', its 'stiffnessMa-
   trix', its 'locationArray'. Since the load vector is recalculated at every
   time step, it is not given the status of attribute.
 TASKS :
   -defining itself :
     .typing itself (methods 'typed' and 'ofType'). When the domain creates
      an element, it actually creates a temporary instance of Element, then
      asks this element to transform itself into an element of the right type
      (PlaneStrain, Truss2D, etc) ;
     .obtaining its basic data : the element reads in the data file the num-
      ber of these objects, then obtains the data from the domain (methods
      'giveNode', 'giveMaterial',etc) ;
   -calculating its contribution to the problem :
     .calculating its conductivity matrix K, its capacity matrix C, its load vector
      f, its location array ;
     .calculating its contribution to the LHS and RHS of the linear system,
      using Static,Newmark,etc, formula. These contributions are usually
      combinations of K,C,f.
   -performing end-of-step operations :
     .calculating the strains and stresses at its Gauss points ;
     .printing its output in the data file and updating itself ;
*/
protected:

  FloatMatrix*  rotationMatrix ;
  int           rotationMatrixDefined ;
  CharType      unknownType ;
  
public:
  PPdeElement (int,Domain*) ;              // constructors
  ~PPdeElement () ;                        // destructor

  // characteristic  matrix
  void  giveCharacteristicMatrix (FloatMatrix& answer, CharType, TimeStep *) ;
  void  giveCharacteristicVector (FloatArray& answer, CharType, ValueModeType, TimeStep *) ;
  
  
  virtual void          computeCapacityMatrix (FloatMatrix& answer, TimeStep*) ;
  virtual void          computeConductivityMatrix (FloatMatrix& answer, 
                                                   MatResponseMode rMode, TimeStep* tStep) ;
  virtual void          computeConstitutiveMatrixAt (FloatMatrix& answer, 
                                                     MatResponseMode rMode,GaussPoint*,
                                                     TimeStep* tStep);
  
  // load vector
  void                  computeRhsVectorAt (FloatArray& answer, TimeStep*, ValueModeType mode) ;
  
  virtual int  giveLocalCoordinateSystem(FloatMatrix& answer) {answer.beEmptyMtrx();return 0;}
  // returns a unit vectors of local coordinate system at element
  // stored rowwise (mainly used by some materials with ortho and anisotrophy)
  
  // data management
  
  // time step termination
  void                  printOutputAt (FILE *, TimeStep*) ;
  void                  updateYourself (TimeStep*) ;
  
  virtual int    checkConsistency ();
  
  // definition
  //       Element*              typed () ;
  //       Element*              ofType (char*) ; 
  const char* giveClassName () const { return "PPdeElement" ;}
  classType                giveClassID () const
    { return PPdeElementClass; }
  /*      int                   giveNumber ()
          { return FEMComponent::giveNumber() ;} */
  
#ifdef __OOFEG
  //
  // Graphics output 
  //
  void          drawYourself (oofegGraphicContext&);
  virtual void  drawRawGeometry (oofegGraphicContext&) {}
  virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
  /*
     returns internal state variable (like stress,strain) at node of element,
     the way how is obtained is eleemnt dependent.
     returns zero if element is unable to respont to DrawMode.
     */
#endif
  
protected:
  virtual void  computeBmatrixAt (GaussPoint*, FloatMatrix& answer)  =0;
  virtual void  computeNmatrixAt (GaussPoint*, FloatMatrix& n)  =0;
  /*
     updates rotation matrix r(l)=T r(g) between  local and global coordinate sysytem
     taking into account also possible local - coordinate system in some elements
     nodes
     @return nonzero if transformation necessary.
     */
  int          updateRotationMatrix () ;      // 
  // give Transformation matrix from global coord. sysyt. to element-local c.s
  virtual int  computeGtoLRotationMatrix (FloatMatrix& answer) {answer.beEmptyMtrx(); return 0;}
  // give Transformation matrix from global coord. syst. to local coordinate system in nodes. 
  virtual int  computeGNDofRotationMatrix (FloatMatrix& answer, DofManTrasfType mode) ;
  
  void                  computeBcRhsVectorAt (FloatArray& answer, TimeStep*, ValueModeType mode) ;
  virtual void          computeGeneratorRhsVectorAt (FloatArray& answer, Load*,TimeStep*, ValueModeType mode) ;
  virtual void            computeEdgeRhsVectorAt (FloatArray& answer, Load*,TimeStep*, ValueModeType mode) ;
  virtual void            computeSurfaceRhsVectorAt (FloatArray& answer, Load*, TimeStep*, ValueModeType mode) ;
  virtual void            computeBcLhsDueConvection (FloatMatrix& answer, TimeStep*);
     
} ;

#define ppdeelement_h
#endif








