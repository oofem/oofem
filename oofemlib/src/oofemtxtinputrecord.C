/* $Header: /home/cvs/bp/oofem/oofemlib/src/inputrecord.C,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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
// file inputrecord.cc
//

#include "oofemtxtinputrecord.h"
#ifndef __MAKEDEPEND
#include <ctype.h>
#endif

#ifdef IR_USE_TOKENIZER 

OOFEMTXTInputRecord::OOFEMTXTInputRecord () : InputRecord(), tokenizer() {
  record[0]='\0';
}

OOFEMTXTInputRecord::OOFEMTXTInputRecord (const OOFEMTXTInputRecord &src) : InputRecord (src), tokenizer() {
  strncpy (this->record, src.record, OOFEM_MAX_LINE_LENGTH);
  this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
  tokenizer.tokenizeLine(this->record);
  int ntok = tokenizer.giveNumberOfTokens();
  for (int i=0; i<ntok; i++) readFlag[i]=src.readFlag[i];
}

OOFEMTXTInputRecord::OOFEMTXTInputRecord (char* source) : InputRecord(), tokenizer()
{
 strncpy (this->record, source, OOFEM_MAX_LINE_LENGTH);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
 tokenizer.tokenizeLine(this->record);
 int ntok = tokenizer.giveNumberOfTokens();
 for (int i=0; i<ntok; i++) readFlag[i]=false;
}

OOFEMTXTInputRecord&
OOFEMTXTInputRecord::operator= (const OOFEMTXTInputRecord& src)
{
 strncpy (this->record, src.record, OOFEM_MAX_LINE_LENGTH);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
 tokenizer.tokenizeLine(this->record);
 int ntok = tokenizer.giveNumberOfTokens();
 for (int i=0; i<ntok; i++) readFlag[i]=src.readFlag[i];
 return *this;
}

void
OOFEMTXTInputRecord::setRecordString (char* newRec) 
{
 strncpy (this->record, newRec, OOFEM_MAX_LINE_LENGTH-1);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
 tokenizer.tokenizeLine(this->record);
 int ntok = tokenizer.giveNumberOfTokens();
 for (int i=0; i<ntok; i++) readFlag[i]=false;
}

IRResultType 
OOFEMTXTInputRecord::giveRecordKeywordField (char* answer, int& value, int maxchar)
{
  if (tokenizer.giveNumberOfTokens()) {
    strncpy (answer, tokenizer.giveToken (1), maxchar);
    answer[maxchar-1] = '\0'; setReadFlag(1);
    if(scanInteger(tokenizer.giveToken (2), value) == 0)  return IRRT_BAD_FORMAT;
    setReadFlag(2);
    
    return IRRT_OK;
  } else {
    return IRRT_NOTFOUND;
  }
}

IRResultType 
OOFEMTXTInputRecord::giveRecordKeywordField (char* answer, int maxchar)
{
  if (tokenizer.giveNumberOfTokens()) {
    strncpy (answer, tokenizer.giveToken (1), maxchar);
    answer[maxchar-1] = '\0'; setReadFlag(1);

    return IRRT_OK;
  } else {
    return IRRT_NOTFOUND;
  }
}

IRResultType 
OOFEMTXTInputRecord::giveField (int& answer, const InputFieldType fieldID, const char* idString)
{
  int indx = this->giveKeywordIndx (idString);
  if (indx) {
    if(scanInteger(tokenizer.giveToken (indx+1), answer) == 0)  return IRRT_BAD_FORMAT;

    setReadFlag(indx);
    setReadFlag(indx+1);
    return IRRT_OK;
  } else return IRRT_NOTFOUND;
}

IRResultType
OOFEMTXTInputRecord::giveField(double& answer, const InputFieldType fieldID, const char* idString)
{
  int indx = this->giveKeywordIndx (idString);
  if (indx) {
    if(scanDouble(tokenizer.giveToken (indx+1), answer) == 0)  return IRRT_BAD_FORMAT;
    setReadFlag(indx);
    setReadFlag(indx+1);
    return IRRT_OK;
  } else return IRRT_NOTFOUND;
 
}

IRResultType 
OOFEMTXTInputRecord::giveField(char* answer, int maxchar, const InputFieldType fieldID, const char* idString)
{
  int indx = 0;
  if (idString) {
    if ((indx = this->giveKeywordIndx (idString)) == 0) return IRRT_NOTFOUND;
    setReadFlag(indx);
    indx++;
  } else indx = 1;
  
  strncpy (answer, tokenizer.giveToken (indx), maxchar);
  answer[maxchar-1] = '\0'; setReadFlag(indx);
  return IRRT_OK;
}

IRResultType 
OOFEMTXTInputRecord :: giveField (IntArray& answer, const InputFieldType fieldID, const char* idString)
{
  int value,size ;
  int indx = this->giveKeywordIndx (idString);
  if (indx) {
    setReadFlag(indx);
    if(scanInteger(tokenizer.giveToken (++indx), size) == 0)  return IRRT_BAD_FORMAT;
    answer.resize(size); setReadFlag(indx);
    
    for (int i = 1 ; i<=size; i++){
      if(scanInteger(tokenizer.giveToken (indx+i), value) == 0)  return IRRT_BAD_FORMAT;
      answer.at(i) = value; setReadFlag(indx+i);
    }
    return IRRT_OK;
  } else return IRRT_NOTFOUND;
}


IRResultType 
OOFEMTXTInputRecord::giveField(FloatArray& answer, const InputFieldType fieldID, const char* idString)
{
  double value; int size ;
  int indx = this->giveKeywordIndx (idString);
  if (indx) {
    setReadFlag(indx);
    if(scanInteger(tokenizer.giveToken (++indx), size) == 0)  return IRRT_BAD_FORMAT;
    answer.resize(size); setReadFlag(indx);
    
    for (int i = 1 ; i<=size; i++){
      if(scanDouble(tokenizer.giveToken (indx+i), value) == 0)  return IRRT_BAD_FORMAT;
      answer.at(i) = value; setReadFlag(indx+i);
    }
    return IRRT_OK;
  } else return IRRT_NOTFOUND;
}


IRResultType 
OOFEMTXTInputRecord::giveField (Dictionary& answer, const InputFieldType fieldID, const char* idString)
{
  double value; int size ;
  char key;
  int indx = this->giveKeywordIndx (idString);
  if (indx) {
    setReadFlag(indx);
    if(scanInteger(tokenizer.giveToken (++indx), size) == 0)  return IRRT_BAD_FORMAT;
    setReadFlag(indx);

    answer.clear();
    for (int i = 1 ; i<=size; i++){
      key = (tokenizer.giveToken(++indx))[0]; setReadFlag(indx);
      if(scanDouble(tokenizer.giveToken (++indx), value) == 0)  return IRRT_BAD_FORMAT;
      setReadFlag(indx);
      answer.add(key, value);
    }
    return IRRT_OK;
  } else return IRRT_NOTFOUND;
}

IRResultType  
OOFEMTXTInputRecord:: giveField (dynaList<Range> &list, const InputFieldType fieldID, const char* idString)
{
  int li, hi;
  const char *rec;
  int indx = this->giveKeywordIndx (idString);
  if (indx) {
    setReadFlag(indx);
    rec = tokenizer.giveToken (++indx);
    if (*rec != '{') {
      OOFEM_WARNING ("OOFEMTXTInputRecord::readRangeList: parse error - missing left '{'");
      list.clear();
      return IRRT_BAD_FORMAT;
    }
    setReadFlag(indx);
    rec ++;
    // read ranges
    while (readRange (&rec, li, hi)) {
      Range range (li, hi);
      list.pushBack (range);
    }
    // skip whitespaces after last range
    while (isspace(*rec)) rec ++;
    // test for enclosing bracket
    if (*rec != '}') {
      OOFEM_WARNING ("OOFEMTXTInputRecord::readRangeList: parse error - missing end '}'");
      list.clear();
      return IRRT_BAD_FORMAT;
    }
    return IRRT_OK;
  } else return IRRT_NOTFOUND;
}

bool
OOFEMTXTInputRecord :: hasField (const InputFieldType fieldID, const char* idString)
{
  //returns nonzero if idString is present in source
  int indx = this->giveKeywordIndx (idString);
  if (indx) setReadFlag(indx);
  return (indx>0)?true:false;
}

int
OOFEMTXTInputRecord :: scanInteger (const char* source, int& value)
{
// 
// reads integer value from source, returns nonzero if o.k.
//
  int i;
  if (source == NULL) {value = 0; return 0;}

  i = sscanf (source,"%d", &value);
  if ((i == EOF)||(i == 0)){ value =0 ; return 0 ;}
  return 1;
}
  

int
OOFEMTXTInputRecord :: scanDouble (const char* source, double& value)
{
// 
// reads integer value from source, returns pointer to char after this number
//
  int i;

  if (source == NULL) {value = 0.0; return 0;}

  i = sscanf (source,"%lf", &value);
  if ((i == EOF)||(i == 0)){ value = 0.0 ; return 0;}
  return 1;
}

int
OOFEMTXTInputRecord :: giveKeywordIndx (const char* kwd)
{
  int ntokens = tokenizer.giveNumberOfTokens();
  for (int i=1; i<=ntokens; i++)
    if (strcmp (kwd, tokenizer.giveToken(i)) == 0) return i;
  return 0;
}

void
OOFEMTXTInputRecord :: finish (bool wrn)
{
  if (!wrn) return;

  char buff[MAX_ERROR_MSG_LENGTH];
  int pos=0;
  int j, pf = 1, wf = 0, ntokens = tokenizer.giveNumberOfTokens();
  for (int i=0; i<ntokens; i++) {
    //fprintf (stderr, "[%s] ", tokenizer.giveToken(i+1));
    if (!readFlag[i]) {
      if (pf) {
        pos=sprintf (buff, "Unread token(s) detected in the following record\n\"");
        for (j=0; j<40; j++) {
          if (record[j]=='\n' || record[j]=='\0') break; else buff[pos++]=record[j];
        }
        if (j == 40) pos+=sprintf (buff+pos, " ...\":\n"); else pos+=sprintf (buff+pos, "\":\n");
        pf = false; wf = 1;
      }
      pos+=sprintf (buff+pos, "[%s] ", tokenizer.giveToken(i+1));
    }
  }
  if (wf) {
    if (pos >=MAX_ERROR_MSG_LENGTH) OOFEM_ERROR3 ("OOFEMTXTInputRecord::finish : print buffer too small (%d,%d)",
                                                  pos, MAX_ERROR_MSG_LENGTH);

    OOFEM_WARNING (buff);
  }
}

#endif // #ifdef IR_USE_TOKENIZER 








#ifndef IR_USE_TOKENIZER 

OOFEMTXTInputRecord::OOFEMTXTInputRecord () : InputRecord() {
  record[0]='\0';
}

OOFEMTXTInputRecord::OOFEMTXTInputRecord (char* source)
{
 strncpy (this->record, source, OOFEM_MAX_LINE_LENGTH);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
}

OOFEMTXTInputRecord::OOFEMTXTInputRecord (const OOFEMTXTInputRecord& src) : InputRecord (src)
{
 strncpy (this->record, src.record, OOFEM_MAX_LINE_LENGTH);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
}



OOFEMTXTInputRecord&
OOFEMTXTInputRecord::operator= (const OOFEMTXTInputRecord& src)
{
 strncpy (this->record, src.record, OOFEM_MAX_LINE_LENGTH);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
 return *this;
}


void
OOFEMTXTInputRecord::setRecordString (char* newRec) 
{
 strncpy (this->record, newRec, OOFEM_MAX_LINE_LENGTH-1);
 this->record[OOFEM_MAX_LINE_LENGTH-1]='\0';
}

IRResultType 
OOFEMTXTInputRecord::giveRecordKeywordField (char* answer, int& value, int maxchar)
{
 char* curr = this->record; // read record keyword;
 int len = 0;

 if (!curr) { return IRRT_NOTFOUND;}
 // skip whitespaces
 while (isspace(*curr) || !*curr) curr++;
 if (!curr) return IRRT_BAD_FORMAT; //{fprintf (stderr,"End-of-line encountered\a\n"); exit(1);} 
 if (*curr == '"') { // read quoted string
   curr++;
   while (!(*curr == '"')) {
     if ((*curr == '\n') || !*curr) return IRRT_BAD_FORMAT; //{ fprintf (stderr,"readQuotedString: final \" expected\a\n"); exit(1);}
     if (++len == (maxchar)) { *answer = 0; break;}
     *answer++ = *curr++;
   }
   *answer = 0;
 } else { // read simple string
   char *s;
   readSimpleString (curr, answer, maxchar, &s);
 } 
 curr = scanInteger(curr, &value);
 if (curr == NULL) return IRRT_BAD_FORMAT;
 return IRRT_OK;
}


IRResultType 
OOFEMTXTInputRecord::giveRecordKeywordField (char* answer, int maxchar)
{
 char* curr = this->record; // read record keyword;
 int len = 0;

 if (!curr) { return IRRT_NOTFOUND;}
 // skip whitespaces
 while (isspace(*curr) || !*curr) curr++;
 if (!curr) return IRRT_BAD_FORMAT; //{fprintf (stderr,"End-of-line encountered\a\n"); exit(1);} 
 if (*curr == '"') { // read quoted string
   curr++;
   while (!(*curr == '"')) {
     if ((*curr == '\n') || !*curr) return IRRT_BAD_FORMAT; //{ fprintf (stderr,"readQuotedString: final \" expected\a\n"); exit(1);}
     if (++len == (maxchar)) { *answer = 0; break;}
     *answer++ = *curr++;
   }
   *answer = 0;
 } else { // read simple string
   char *s;
   readSimpleString (curr, answer, maxchar, &s);
 } 

 return IRRT_OK;
}


IRResultType 
OOFEMTXTInputRecord::giveField (int& answer, const InputFieldType fieldID, const char* idString)
{
  char* str = getPosAfter(this->record,idString);
  if (str == NULL) return IRRT_NOTFOUND;
 answer = atoi(str);
 return IRRT_OK;
}

IRResultType
OOFEMTXTInputRecord::giveField(double& answer, const InputFieldType fieldID, const char* idString)
{
  char* str = getPosAfter(this->record,idString);
  
  if (str == NULL) return IRRT_NOTFOUND;
  answer = strtod (str, NULL);
 return IRRT_OK;
 
}

IRResultType 
OOFEMTXTInputRecord::giveField(char* answer, int maxchar, const InputFieldType fieldID, const char* idString)
{
 char* curr;
 int len = 0;

 if (idString)  curr = getPosAfter(this->record,idString) ;
 else curr = this->record; // read record keyword

 if (!curr) { return IRRT_NOTFOUND;}
 // skip whitespaces
 while (isspace(*curr) || !*curr) curr++;
  if (!curr) return IRRT_BAD_FORMAT; //{fprintf (stderr,"End-of-line encountered\a\n"); exit(1);} 
 if (*curr == '"') { // read quoted string
  curr++;
  while (!(*curr == '"')) {
   if ((*curr == '\n') || !*curr) return IRRT_BAD_FORMAT; //{ fprintf (stderr,"readQuotedString: final \" expected\a\n"); exit(1);}
   if (++len == (maxchar)) { *answer = 0; break;}
   *answer++ = *curr++;
  }
  *answer = 0;
 } else { // read simple string
  char *s;
  readSimpleString (curr, answer, maxchar, &s);
 } 
 return IRRT_OK;
}


IRResultType 
OOFEMTXTInputRecord :: giveField (IntArray& answer, const InputFieldType fieldID, const char* idString)
{
  char *str1; int value,size ;

  str1 = getPosAfter(this->record,idString) ;
 if (str1 == NULL) return IRRT_NOTFOUND;
  str1 = scanInteger(str1,&size) ;
 if(str1 == NULL)  return IRRT_BAD_FORMAT;
  answer.resize(size);
 
  for (int i = 1 ; i<=size; i++){
    str1 = scanInteger (str1,&value);
  if (str1)
   answer.at(i) = value;
  else return IRRT_BAD_FORMAT;
  }
  return IRRT_OK;
}


IRResultType 
OOFEMTXTInputRecord::giveField(FloatArray& answer, const InputFieldType fieldID, const char* idString)
{
  char *str1; double  value; int size ;

  str1 = getPosAfter(this->record,idString) ;
 if (str1 == NULL) return IRRT_NOTFOUND;
  str1 = scanInteger(str1,&size) ;
  if(str1 == NULL)  return IRRT_BAD_FORMAT;
  answer.resize(size);
  for (int i = 1 ; i<=size; i++){
    str1 = scanDouble (str1,&value);
  if (str1)
   answer.at(i) = value;
  else return IRRT_BAD_FORMAT;
  }
  return IRRT_OK;
}


IRResultType 
OOFEMTXTInputRecord::giveField (Dictionary& answer, const InputFieldType fieldID, const char* idString)
{
  char *str1; double  value; int size ;
  char key [MAX_NAME_LENGTH+1]; // 'key' is eventually of size 1, but some words that are
                                // read in the meantime in the data file can be larger !

  str1 = getPosAfter(this->record,idString) ;  // move after idString
 if (str1 == NULL) return IRRT_NOTFOUND;
  str1 = scanInteger(str1,&size) ;       // read number of conditions
 if(str1 == NULL)  return IRRT_BAD_FORMAT;
 answer.clear();

  for (int i = 1 ; i<=size; i++){
    readSimpleString (str1,key,MAX_NAME_LENGTH+1,&str1);
    str1 = scanDouble(str1,&value);
  if (str1 == NULL) return IRRT_BAD_FORMAT;
    answer.add(key[0], value);
  }
  return IRRT_OK;
}

IRResultType  
OOFEMTXTInputRecord:: giveField (dynaList<Range> &list, const InputFieldType fieldID, const char* idString)
{
 int li, hi;
 char* str1;
 const char *helpSource = this->record ;
 // Range* range;

 // find first valid occurence of idString
 int len = strlen(idString);
 do {
  if ((str1 =strstr(helpSource,idString))==NULL) return IRRT_NOTFOUND;
  helpSource = str1+1;
 } while (! ((*(helpSource+len-1) == ' ') || (*(helpSource+len-1) == '\t') || (*(helpSource+len-1) == '\n')) );
 

 helpSource = str1 + len;
 // find first nonwhitespace character
 // skip whitespaces
 while (isspace(*helpSource)) helpSource ++;
 // test if list left bracked found
 if (*helpSource != '{') {
   OOFEM_WARNING ("OOFEMTXTInputRecord::readRangeList: parse error - missing left '{'");
   list.clear();
   return IRRT_BAD_FORMAT;
 }
 helpSource ++;
 // read ranges
 while (readRange (&helpSource, li, hi)) {
  Range range (li, hi);
  list.pushBack (range);
 }
 // skip whitespaces after last range
 while (isspace(*helpSource)) helpSource ++;
 // test for enclosing bracket
 if (*helpSource != '}') {
   OOFEM_WARNING ("OOFEMTXTInputRecord::readRangeList: parse error - missing end '}'");
   list.clear();
   return IRRT_BAD_FORMAT;
 }
 return IRRT_OK;
}

bool
OOFEMTXTInputRecord :: hasField (const InputFieldType fieldID, const char* idString)
{
  //returns nonzero if idString is present in source
  char* str = strstr(this->record,idString);
  if (str == NULL) return false;
  return true;
}

char* OOFEMTXTInputRecord :: readSimpleString (char* source, char* simpleString, int maxchar, char** remain)
// reads Simple string from source according to following rules:
// at begining skips whitespace (blank, tab)
// read string terminated by whitespace or end-of-line
// remain is unread remain of source string.
// maximum of maxchar (including terminating '\0') is copyied into simpleString.
{
  char *curr = source;
  char *ss = simpleString ;
 int count = 0;
  
  if (source == NULL) {*remain = NULL; return NULL;}

  while (isspace(*curr) || !*curr) curr++;
  if (!curr) { OOFEM_ERROR ("OOFEMTXTInputRecord::readSimpleString : unexpected end-of-line"); }
  while ((!(isspace(*curr) || !*curr)) && (++count < maxchar)) 
    *ss++ = *curr++;

  *ss = '\0' ;
  *remain = curr;
  return simpleString;
}
 

char*  OOFEMTXTInputRecord :: readKeyAndVal (char* source, char* key, int* val, int maxchar, char** remain)
//
// 
//
{

  key = readSimpleString (source,key,maxchar,remain);
  *remain = scanInteger(*remain,val);
  return *remain;
}


char*  
OOFEMTXTInputRecord :: readKeyAndVal (char* source, char* key, double* val, int maxchar, char** remain)
//
// 
//
{
  key = readSimpleString (source,key,maxchar,remain);
  *remain = scanDouble(*remain,val);
  return *remain;
}

char* OOFEMTXTInputRecord :: getPosAfter (char* source, const char* idString)
// 
// returns possition of substring idString in source
// return value pointer at the end of occurence idString in source
// (idString must be separated from rest by blank or by tabulator
// if string not found, returns NULL
//
{
 char* str1, *helpSource = source;
 int len = strlen(idString);
 int whitespaceBefore, whitespaceAfter;

 do {
  if ((str1 =strstr(helpSource,idString))==NULL) return NULL;
  helpSource = str1+1;
  whitespaceAfter = isspace(*(helpSource+len-1));
  if (str1 == source) whitespaceBefore = 1;
  else whitespaceBefore = isspace(*(str1-1));
 } while (! (whitespaceBefore&&whitespaceAfter));
 
 return str1+len;
}

char* OOFEMTXTInputRecord :: skipNextWord (char*src)
//
// skips next word in src ; returns pointer after it
//
{
  
    while (isspace(*src) || !*src) src++;     
    // skips whitespaces if any
    while (!(isspace(*src) || !*src)) src ++;
    // skips one word
    return src ;
  }
    

char* OOFEMTXTInputRecord :: scanInteger (char* source, int* value)
{
// 
// reads integer value from source, returns pointer to char after this number
//
  int i;

  if (source == NULL) {*value = 0; return NULL;}

  i = sscanf (source,"%d",value);
  if (i == EOF ){ *value =0 ; return NULL ;}
  return skipNextWord(source);
}
  

char* OOFEMTXTInputRecord :: scanDouble (char* source, double* value)
{
// 
// reads integer value from source, returns pointer to char after this number
//
  int i;

  if (source == NULL) {*value = 0; return NULL;}

  i = sscanf (source,"%lf",value);
  if (i == EOF ){ *value =0 ; return NULL ;}
  return skipNextWord(source);
}


void
OOFEMTXTInputRecord :: finish (bool wrn)
{}

#endif // #ifndef IR_USE_TOKENIZER 


IRResultType 
OOFEMTXTInputRecord::giveOptionalField (int& answer, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}

IRResultType 
OOFEMTXTInputRecord::giveOptionalField (double& answer, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}

IRResultType
OOFEMTXTInputRecord::giveOptionalField (char* answer, int maxchar, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, maxchar, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}

IRResultType
OOFEMTXTInputRecord::giveOptionalField (FloatArray& answer, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}

IRResultType
OOFEMTXTInputRecord::giveOptionalField (IntArray& answer, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}

IRResultType
OOFEMTXTInputRecord::giveOptionalField (Dictionary& answer, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}

IRResultType
OOFEMTXTInputRecord::giveOptionalField (dynaList<Range> &answer, const InputFieldType fieldID, const char* idString)
{
 IRResultType r = this->giveField (answer, fieldID, idString);
 if (r == IRRT_NOTFOUND) return IRRT_OK;
 else return r;
}


int
OOFEMTXTInputRecord :: readRange (const char** helpSource, int& li, int& hi)
{
 char* endptr;
 // skip whitespaces
 while (isspace(**helpSource)) (*helpSource) ++;
 // test if character is digit 
 if (isdigit (**helpSource)) {
  // digit character - read one value range
  li = hi = strtol (*helpSource, &endptr, 10);
  *helpSource = endptr;
  return 1;
 } else if (**helpSource == '(') {
  // range left parenthesis found
  (*helpSource)++;
  // read lower index
  li = strtol (*helpSource, &endptr, 10);
  *helpSource = endptr;
  // test whitespaces
  if (**helpSource != ' ' && **helpSource!='\t') {
   OOFEM_WARNING ("OOFEMTXTInputRecord::readRange: unexpected token while reading range value");
   return 0;
  }
  // read end index
  hi = strtol (*helpSource, &endptr, 10);
  *helpSource = endptr;
  // skip whitespaces
  while (isspace(**helpSource)) (*helpSource) ++;
  // test for enclosing bracket
  if (**helpSource == ')') {(*helpSource)++; return 1;}
  else {
   OOFEM_WARNING ("OOFEMTXTInputRecord::readRange: end ')' missing while parsing range value");
   return 0;
  }
 } 
 return 0;
}

