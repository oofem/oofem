#include "stdio.h"
#include "string.h"

#define MAX_LINE_LENGTH 1024
#define MAX_TOKEN_LENGTH 80
#define MAX_TOKENS 20

typedef char tokenType[MAX_TOKENS][MAX_TOKEN_LENGTH] ;

char* giveLineFromInput (char* line) ;
int tokenizeLine (char separator, char* line, tokenType* tokens);


class Tokenizer {
private:
 tokenType tokens;
 FILE* inputStream;
 char separator;

 int  currentTokens, isEOFflag;
 char currentLine[MAX_LINE_LENGTH];
 int tokenizeLine ();

public:
 Tokenizer (FILE* inFile, char separator = 0);
 char* giveLineFromInput ();
 int   giveNumberOfTokens ();
 int   isEOF() {return isEOFflag;} 
 char* giveToken (int);
 char* giveCurrentLine () {return currentLine;}
 // return token number if token present; zero otherwise
 int   hasToken (const char* str);
 
};


