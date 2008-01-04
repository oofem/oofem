/* $Header: /home/cvs/bp/oofem/oofemlib/src/util.C,v 1.13.4.1 2004/04/05 15:19:44 bp Exp $ */
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

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#include <ctype.h> 
#endif
#include "cltypes.h"
#include "engngm.h"
#include "usrdefsub.h"
#include "util.h"
#include "inputrecord.h"


char* giveLineFromInput (FILE* inputStream, char* line, int len)
//
// reads one line from inputStream - for private use only.
//
{
 char *ptr;

 giveRawLineFromInput (inputStream, line, len);
 // convert line to lowercase
 for (ptr=line; (*ptr = tolower (*ptr)); ptr++);
  return line;
}

char* giveRawLineFromInput (FILE* inputStream, char* line, int len)
//
// reads one line from inputStream - for private use only.
//
{
  char* _res;
  do {
    _res = fgets(line,len,inputStream);
    if (_res == NULL) {
      OOFEM_ERROR ("giveRawLineFromInput : End of file encountered");
    }
  } while (*line == '#');   // skip comments
  return line;
}


char*  giveInputDataFileName (char* dataInputFileName, int maxlen)

{
 int len;
 // Returns the name of the file containing the data of the problem.
 // char s[MAX_FILENAME_LENGTH] ;
 
 printf ("please enter the name of the input data file : \n") ;
 //gets (s) ;
 //strcpy (dataInputFileName,s) ;
 char* _res = fgets(dataInputFileName, maxlen, stdin);
 if (_res == NULL) OOFEM_ERROR ("giveInputDataFileName: reading error or EOF encountered");
  // test if last read character is newline
 if ((len = strlen(dataInputFileName))) 
  if (dataInputFileName[len-1] == '\n')    // if yes remove it
   dataInputFileName[len-1] = '\0';
 
 return dataInputFileName ;
}


EngngModel* InstanciateProblem (DataReader* dr, problemMode mode, int contextFlag, EngngModel* _master)
{
  const char *__keyword, *__proc = "InstanciateProblem"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                                   // Required by IR_GIVE_FIELD macro
  EngngModel* problem;
  char desc [OOFEM_MAX_LINE_LENGTH+1], problemName[MAX_NAME_LENGTH];
  char dataOutputFileName[MAX_FILENAME_LENGTH];
  
  InputRecord* ir = dr->giveInputRecord (DataReader::IR_outFileRec, 1);
  __keyword = NULL; result = ir->giveField(dataOutputFileName, MAX_FILENAME_LENGTH, IFT_EngngModel_outfile, __keyword);
  if (result != IRRT_OK) IR_IOERR ("", __proc, IFT_EngngModel_outfile, "Output file record", ir, result);
  ir->finish();
  
  ir = dr->giveInputRecord (DataReader::IR_jobRec, 1); 
  __keyword = NULL; result = ir->giveField(desc, OOFEM_MAX_LINE_LENGTH, IFT_EngngModel_probdescription, __keyword);
  
  /* here we need copy of input record. The pointer returned by dr->giveInputRecord can (and will) 
     be updated as reading e-model components (nodes, etc). But we need this record being available
     through the whole e-model instanciation 
  */
  InputRecord* emodelir = dr->giveInputRecord (DataReader::IR_emodelRec, 1) -> GiveCopy();
  result = emodelir->giveRecordKeywordField(problemName, MAX_NAME_LENGTH);
  //result = IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);
  if (result != IRRT_OK) IR_IOERR ("", __proc, IFT_EngngModel_probname, "", emodelir, result);
  
  problem = CreateUsrDefEngngModelOfType (problemName,1,_master) ;
  problem->setProblemMode (mode);
  
  if (contextFlag) problem-> setContextOutputMode(ALWAYS);
  problem -> instanciateYourself (dr, emodelir, dataOutputFileName, desc) ;
  //emodelir.finish();
  delete (emodelir);
  
  return problem;
}
