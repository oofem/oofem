/* $Header: /home/cvs/bp/oofem/oofemlib/src/fieldmanager.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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

#ifndef fieldmanager_h
#define fieldmanager_h

#include "field.h"
#ifndef __MAKEDEPEND
#include <map>
#endif

/**
*/
class FieldManager 
{
protected:

 /** Field container. Stores only pointers to objects (not object themselves) 
  to avoid copiing elements and to preserve the use of polymorphic types.
  This implies that the objects in container must be explicitly deleted.
  */
 std::map<FieldType, Field*> externalFields;


public:

 FieldManager () {}
 /**
  Registers the given field (the receiver is not assumed to own given field).
  The field is registered under given key. Using this key, it can be later accessed.
  */
 void registerField (Field* eField, FieldType key);
 /** Returns true if field is registered under key */
 bool isFieldRegistered (FieldType key);
 /**
  Returns the previously registered field under given key; NULL otherwise
  */
 Field* giveField (FieldType key) ;
 /**
  Unregisters (deletes) the field registered under given key.
  */
 void unregisterField (FieldType key);

};

#endif // fieldmanager_h
