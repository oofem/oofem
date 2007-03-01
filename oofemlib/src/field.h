/* $Header: /home/cvs/bp/oofem/oofemlib/src/field.h,v 1.1.4.1 2004/04/05 15:19:43 bp Exp $ */
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

#ifndef field_h

#include "domain.h"


enum FieldType
{
 FT_Unknown,
 FT_Temperature,
 FT_HumidityConcentration,
 FT_Velocity,
};


enum FieldBaseID
{
  FBID_FluxField,
  FBID_PressureField,
  FBID_VelocityField,
  FBID_VelocityPressureField,
};

/**
 Abstract class representing field. Field represent the spatial distribution of certain variable.
 Field is able to evaluate its value at any point of interest. The field is usually associated to 
 specific domain. 
*/
class Field 
{
protected:
 //Domain* domain;
 FieldBaseID type;
 

public:
 /** Constructor. Creates a field of given type associated to given domain.
  */
 Field (FieldBaseID b) {type = b;}
 virtual ~Field() {}
 /** Evaluates the field at given point
  @param coords coordinates of the point of interest
  @return error code (0-ok, 1-point not found in domain)
  */
 virtual int evaluateAt (FloatArray& answer, FloatArray& coords, IntArray& dofId, 
                         ValueModeType mode, TimeStep* atTime) = 0;
 /// Returns the type of receiver
 FieldBaseID giveType () {return type;}


 /** Stores receiver state to output stream. 
  Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
  @param stream output stream 
  @return contextIOResultType
  @exception throws an ContextIOERR exception if error encountered
  */
 virtual contextIOResultType    saveContext (FILE* stream) = 0;
 /** Restores the receiver state previously written in stream.
  Reads the FEMComponent class-id in order to allow test consistency.
  @see saveContext member function.
  @return contextIOResultType
  @exception throws an ContextIOERR exception if error encountered
  */
 virtual contextIOResultType    restoreContext(FILE* stream) = 0;
 

 /** Returns class name of the receiver.
  */
 virtual const char*  giveClassName () const {return "Field" ;}
  /**@name error and warning reporting methods
    These methods will print error (or warning) message using oofem default loggers.
    Do not use these methods directly, to avoid specify file and line parameters.
    More preferably, use these methods via corresponding OOFEM_CLASS_ERROR and OOFEM_CLASS_WARNING macros,
    that will include file and line parameters automatically.
    
    Uses variable number of arguments, so a format string followed by optional argumens is expected 
    (according to printf conventions).
    @param file  source file name, where error encountered (where error* function called)
    @param line  source file line number, where error encountered
  */
  //@{
  /// prints error message and exits
  void error (char* file, int line, char *format, ...) const ;
  /// prints warning message
  void warning (char* file, int line, char *format, ...) const ;
 //@}
};

#define field_h
#endif
