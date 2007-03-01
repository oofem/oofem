/* $Header: /home/cvs/bp/oofem/oofemlib/src/inputrecord.h,v 1.4 2003/05/19 13:03:57 bp Exp $ */
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

#ifndef oofemtxtinputrecord_h

#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif
#include "inputrecord.h"
#include "intarray.h"
#include "flotarry.h"
#include "dictionr.h"
#include "dynalist.h"
#include "range.h"

#define IR_MAX_ERROR_LENGTH 100

/** Define this directive to use tokenizer to parse records.
  This also enables to perform additional check for input records, since
  unread fields can be detected
*/
#define IR_USE_TOKENIZER

#ifdef IR_USE_TOKENIZER 
#include "tokenizer.h"
#endif




/**
 Class representing the Input Record for OOFEM txt input file format. The input record is represented as string consisting of several fields.
*/
class OOFEMTXTInputRecord : public InputRecord
{
 protected:
#ifdef IR_USE_TOKENIZER 
  Tokenizer tokenizer;
  bool readFlag[MAX_TOKENS];
#endif
  // record representation
  char record[OOFEM_MAX_LINE_LENGTH];

public:
  /// Constructor. Creates an empty input record.
  OOFEMTXTInputRecord ();
  /// Constructor. Creates the input record corresponding to given string
  OOFEMTXTInputRecord (char* source) ;
  /// Copy constructor
  OOFEMTXTInputRecord (const OOFEMTXTInputRecord&); 
  /// Destructor
  ~OOFEMTXTInputRecord () {}
  /// Assingnment operator
  OOFEMTXTInputRecord&  operator=  (const OOFEMTXTInputRecord&);

  /** Creates a newly allocated copy of the receiver */
  virtual InputRecord* GiveCopy () {return new OOFEMTXTInputRecord (*this);}

 public:
 /// Sets the record string
 void setRecordString (char*);
 /// Returns record string
 char* giveRecordAsString () {return this->record;}

 /** terminates the current record session and if flag is true warnin is printed for unscanned tokens */
 void finish (bool wrn = true);



 public:
 /**@name Compulsory field extraction methods
  Reads the field value identified by keyword
  @param answer contains result
  @param idString field keyword
  @return IRResultType 
  */
  //@{
 /// Reads the record id field  (type of record) and its corresponding number
 virtual IRResultType giveRecordKeywordField(char* answer, int& value, int maxchar) ;
 /// Reads the record id field  (type of record) 
 virtual IRResultType giveRecordKeywordField(char* answer, int maxchar);
 /// Reads the integer field value
 virtual IRResultType giveField (int& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the double field value
 virtual IRResultType giveField (double& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the char* field value
 virtual IRResultType giveField (char* answer, int maxchar, const InputFieldType fieldI, const char* idString);
 /// Reads the FloatArray field value
 virtual IRResultType giveField (FloatArray& answer, const InputFieldType fieldI, const char* idString);
 /// Reads the IntArray field value
 virtual IRResultType giveField (IntArray& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the Dictionary field value
 virtual IRResultType giveField (Dictionary& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the dynaList<Range> field value
 virtual IRResultType giveField (dynaList<Range> &answer, const InputFieldType fieldID, const char* idString);
 //@}

 /**@name Optional field extraction methods
  Reads the field value identified by keyword
  @param answer contains result
  @param idString field keyword
  @return IRResultType 
  */
 //@{
 /// Reads the integer field value
 virtual IRResultType giveOptionalField (int& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the double field value
 virtual IRResultType giveOptionalField (double& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the char* field value
 virtual IRResultType giveOptionalField (char* answer, int maxchar, const InputFieldType fieldID, const char* idString);
 /// Reads the FloatArray field value
 virtual IRResultType giveOptionalField (FloatArray& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the IntArray field value
 virtual IRResultType giveOptionalField (IntArray& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the Dictionary field value
 virtual IRResultType giveOptionalField (Dictionary& answer, const InputFieldType fieldID, const char* idString);
 /// Reads the dynaList<Range> field value
 virtual IRResultType giveOptionalField (dynaList<Range> &answer, const InputFieldType fieldID, const char* idString);
 //@}

 /// Returns true if record contains field identified by idString keyword
 virtual bool         hasField  (const InputFieldType fieldID, const char* idString);

protected:

#ifdef IR_USE_TOKENIZER 
  int giveKeywordIndx (const char* kwd);
  int scanInteger (const char* source, int& value);
  int scanDouble  (const char* source, double& value);
  void setReadFlag (int itok) {readFlag[itok-1]=true;}
 
#endif

#ifndef IR_USE_TOKENIZER 

  char* getPosAfter (char*, const char *) ;
  char* scanInteger (char* source, int* value);
  char* scanDouble  (char* source, double* value);
  char* skipNextWord (char* src) ;

  char*  readSimpleString (char* source, char* simpleString, int maxchar, char** remain); 
  char*  readKeyAndVal    (char* source, char* key, int* val, int maxchar, char** remain);
  char*  readKeyAndVal    (char* source, char* key, double* val, int maxchar, char** remain);
#endif
 /**
  Reads single range record from input record represented by *helpSource  string.
  @param helpSource pointer to current string possition, on return helpSource points
  to next charcter after reading range record.
  @param li starting range index
  @param hi end range index
  @return on success nonzero valur returned
  */
 int    readRange (const char** helpSource, int& li, int& hi);

};

#define oofemtxtinputrecord_h
#endif
