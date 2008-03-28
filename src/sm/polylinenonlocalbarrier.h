/* $Header: /home/cvs/bp/oofem/sm/src/polylinenonlocalbarrier.h,v 1.2.4.1 2004/04/05 15:19:47 bp Exp $ */
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


//
// class PolylineNonlocalBarrier
//
#ifndef polylinenonlocalbarrier_h
#define polylinenonlocalbarrier_h

#include "nonlocalbarrier.h"
#include "domain.h"
#include "cltypes.h"
#include "flotarry.h"
/**
 Implemantation of polyline nonlocal barier.
 It is a composite one-dimensional cell consisting of one or more connected lines. 
 The polyline is defined by an ordered list of n+1 vertices (nodes), 
 where n is the number of lines in the polyline. 
 Each pair of points (i,i+1) defines a line.

 The purpose of this class is to 
 model barrier for nonlocal averaging process (visibility criterion).
 Usually, the given remote integration point influences to the source point
 nonlocal average if the averaging function at source point and evaluated for
 remote point has nonzero value. The barrier allows to exclude additional points,
 which may be close enough, but due to several reasons there is no influence
 between these points (for example, they can be  separated by a notch).

 @see NonlocalBarrier class.
 */
class PolylineNonlocalBarrier : public NonlocalBarrier
{
 
protected:
 /// local coordinate indexes
 int localXCoordIndx, localYCoordIndx;
 /// list of polyline vertices
 IntArray vertexNodes;
 
public:
 /**
  Constructor. Creates an element with number n belonging to domain aDomain.
  @param n Element's number
  @param aDomain Pointer to the domain to which element belongs.
  */
 PolylineNonlocalBarrier (int n,Domain* aDomain) ;    // constructors
 /// Virtual destructor.
  virtual ~PolylineNonlocalBarrier () ;                // destructor

 /**
  Returns true if the barrier is activated 
  by interaction of two given points. In this case the nonlocal influence
  is not considered. Otherwise returns false.
  @param c1 coordinates of first point
  @param c2 coordinates of second point
  @return true if barrier is activated, false otherwise
  */
 virtual bool isActivated (const FloatArray& c1, const FloatArray& c2);
 /**
    Abstract method modifying the integration weight between master (c1) and source (c2) point.
    @param c1 coordinates of master point
    @param c2 coordinates of source point
    @param weight original integration weight; on output modified weight
    @param shieldFlag (output param; set to true if shielding is activated)
  */
 virtual void applyConstraint (const FloatArray& c1, const FloatArray& c2, double &weight, 
                               bool &shieldFlag, NonlocalMaterialExtensionInterface* nei);
 /** Initializes receiver acording to object description stored in input record.
  This function is called immediately after creating object using
  constructor. Input record can be imagined as data record in component database
  belonging to receiver. Receiver may use value-name extracting functions 
  to extract particular field from record. 
  @see readInteger, readDouble and similar functions */
 virtual IRResultType initializeFrom (InputRecord* ir);

 
 /// Returns class name of the receiver.
 const char* giveClassName () const { return "PolylineNonlocalBarrier" ;}
 /** Returns classType id of receiver.
  @see FEMComponent::giveClassID 
  */
 classType                giveClassID () const { return PolylineNonlocalBarrierClass; }

 

};

#endif // polylinenonlocalbarrier_h

