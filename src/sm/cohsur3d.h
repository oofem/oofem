/*
class CohesiveSurface3d added by Milan Jirasek on 1 Feb 2010

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       



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

//   *******************************
//   *** CLASS CohesiveSurface3d ***
//   *******************************

#ifndef cohsur3d_h
#define cohsur3d_h

#include "structuralelement.h"

namespace oofem {

class CohesiveSurface3d : public  StructuralElement
{
/*
   This class implements a cohesive surface element used by the 
   cohesive particle model.
*/
protected:
  double area, length;
  FloatArray center; // coordinates of the center of the cohesive surface
  FloatMatrix lcs; // matrix defining the local coordinate system
  
  // shift constants of periodic particles (near boundary of periodic cell)
  int kx, ky, kz;
  double kxa, kyb, kzc;

public:
  CohesiveSurface3d (int,Domain*); // constructor
  ~CohesiveSurface3d () {};        // destructor

  void          computeBmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer, int li, int ui);
  double        computeVolumeAround (GaussPoint*) ;
  virtual int   computeNumberOfDofs (EquationID ut) {return 6*giveNumberOfNodes();}
  virtual void  giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double        giveLength ();
  virtual void  computeNmatrixAt(GaussPoint*, FloatMatrix&){};
  virtual int   computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords);  
  // definition & identification

  const char* giveClassName () const { return "CohesiveSurface3d" ;}
  classType            giveClassID () const { return CohesiveSurface3dClass; }

  // input and output 

  IRResultType initializeFrom (InputRecord* ir);
  void  printOutputAt (FILE* , TimeStep* );

  
#ifdef __OOFEG
  void drawRawGeometry (oofegGraphicContext&);
  void drawDeformedGeometry (oofegGraphicContext&, UnknownType);
  void drawScalar(oofegGraphicContext& context);
#endif
  

protected:
  void computeGaussPoints ();
  void evaluateCenter ();
  void evaluateLocalCoordinateSystem ();
public:
  integrationDomain  giveIntegrationDomain () {return _Line;}
  MaterialMode          giveMaterialMode()  {return _3dInterface;}
  
} ;
} // namespace oofem
#endif
