/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

//   **************************************************************
//   *** CLASS GENERAL ELEMENT CLASS FOR FLUID DYNAMIC PROBLEMS ***
//   **************************************************************
 
#ifndef fluiddynamicelement_h


#include "element.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "primaryfield.h"

class TimeStep ; class Node ; class Material ; class GaussPoint ;
class FloatMatrix ; class FloatArray ; class IntArray ;

/**
   This abstract class represent a general base element class for 
   fluid dynamic problems.
*/
class FluidDynamicElement : public Element
{
public:
  enum ElementMode {IncompressibleEM, CompressibleEM};
protected:
  ElementMode emode;

public:
  // constructor
  FluidDynamicElement (int,Domain*, ElementMode em = IncompressibleEM) ;
  ~FluidDynamicElement () ;                        // destructor
  
  // characteristic  matrix
  void  giveCharacteristicMatrix (FloatMatrix& answer, CharType, TimeStep *) ;
  void  giveCharacteristicVector (FloatArray& answer, CharType, ValueModeType, TimeStep *) ;
 
  /** Computes the mass matrix M of the receiver */
  virtual void          computeMassMatrix (FloatMatrix& answer, TimeStep*) = 0;
  /** Computes the diffusion matrix K of the receiver */
  virtual void          computeDiffusionMatrix (FloatMatrix& answer, TimeStep* tStep) = 0;
  /** Computes the advection matrix N(u) of the receiver */
  virtual void          computeAdvectionMatrix (FloatMatrix& answer, TimeStep* tStep) = 0;
  /** Computes the coupling u-p matrix C of the receiver */
  virtual void          computeCouplingUPMatrix (FloatMatrix& answer,TimeStep* tStep) = 0;
  /** Computes the RHS contribution to balance equation(s) due to boundary conditions */
  virtual void          computeBCVectorAt (FloatArray& answer, TimeStep*, ValueModeType mode) = 0;
  
  // time step termination
  /**
     Updates element state corresponding to newly reached solution.
     It computes stress vector in each element integration point (to ensure that data in integration point's
     statuses are valid). 
     @param tStep finished time step
  */
  //void                  updateInternalState (TimeStep*) ;
  void                  printOutputAt (FILE *, TimeStep*) ;
  virtual int           checkConsistency ();
  
  // definition
  const char* giveClassName () const { return "FluidDynamicElement" ;}
  classType                giveClassID () const { return FluidDynamicElementClass; }
  
  virtual void           giveElementDofIDMask  (IntArray& answer) const = 0;
 
 /**
  Returns the mask of reduced indexes of Internal Variable component .
  @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
  @param type determines the internal variable requested (physical meaning)
  @returns nonzero if ok or error is generated for unknown mat mode.
  */
 virtual int giveIntVarCompFullIndx (IntArray& answer, InternalStateType type);

#ifdef __OOFEG
int giveInternalStateAtNode (FloatArray& answer, InternalStateType type, InternalStateMode mode, 
                             int node, TimeStep* atTime);
 //
 // Graphics output 
 //
 //void          drawYourself (oofegGraphicContext&);
 //virtual void  drawRawGeometry (oofegGraphicContext&) {}
 //virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif
 
protected:
} ;

#define fluiddynamicelement_h
#endif







