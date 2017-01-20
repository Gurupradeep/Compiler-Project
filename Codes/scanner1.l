%{
	#include <stdlib.h>
	int linecount = 1;
%}

VARNAME [_a-zA-Z][a-zA-Z0-9]* 
FUNCTION [_a-zA-Z][a-zA-Z0-9]*\([ ]*int[ ]+{VARNAME}[ ]*,[ ]*int[ ]+{VARNAME}[ ]*\) 
DIGIT [0-9]
KEYWORD int|main|return|for|do" "while|while|if|else
PREPROCESSOR #[a-zA-Z]+
HEADER <[a-zA-Z]+\.h>
STRING \".+\"
SINGLE_LINE_COMMENT \/\/.*
START_MULTI_LINE \/\*
%s startComment

%%

{START_MULTI_LINE} 	{BEGIN startComment; printf("START COMMENT at line %d\n",linecount);}
<startComment>\n 	{linecount++;}
<startComment>\*\/	{BEGIN 0; printf("END COMMENT at line %d\n",linecount);}
<startComment>[^(\*\/|\n)]*	{;}
{FUNCTION} printf("%s FUNCTION\n",yytext);
{STRING} printf("%s : STRING\n",yytext);
{PREPROCESSOR} printf("%s : PREPROCESSOR DIRECTIVE at line %d\n",yytext,linecount);
{SINGLE_LINE_COMMENT} printf("%s : COMMENTS at line %d\n",yytext,linecount);
"=" printf("%s : ASSIGNMENT at line %d\n ",yytext,linecount);
"+" printf("%s : ADD at line %d\n",yytext,linecount);
"-" printf("%s : SUBTRACT at line %d\n",yytext,linecount);
"*" printf("%s : MULTIPLY at line %d\n",yytext,linecount);
"/" printf("%s : DIVIDE at line %d\n",yytext,linecount);
"%" printf("%s : MODULO at line %d\n",yytext,linecount);
"++" printf("%s : INCREMENT at line %d\n",yytext,linecount);
"--" printf("%s : DECREMENT at line %d\n",yytext,linecount);
"+=" printf("%s : ADD_ASSIGN at line %d\n",yytext,linecount);
"-=" printf("%s : SUB_ASSIGN at line %d\n",yytext,linecount);
"*=" printf("%s : MULT_ASSIGN at line %d\n",yytext,linecount);
"/=" printf("%s : DIV_ASSIGN at line %d\n",yytext,linecount);
"%=" printf("%s : MOD_ASSIGN at line %d\n",yytext,linecount);
"==" printf("%s : EQUAL at line %d\n",yytext,linecount);
">=" printf("%s : GREAT_EQ at line %d\n",yytext,linecount);
"<=" printf("%s : LESS_EQ at line %d\n",yytext,linecount);
">" printf("%s : GREAT at line %d\n",yytext,linecount);
"<" printf("%s : LESS at line %d\n",yytext,linecount);
"?" printf("%s : CONDITIONAL at line %d\n",yytext,linecount);
"&" printf("%s : REFERENCE/BIT_AND at line %d\n",yytext,linecount);
"|" printf("%s : BIT_OR at line %d\n",yytext,linecount);
"||" printf("%s : LOGIC_OR at line %d\n",yytext,linecount);
"&&" printf("%s : LOGIC_AND at line %d\n",yytext,linecount);
"!" printf("%s : NEGATE at line %d\n",yytext,linecount);
":" printf("%s : COLON at line %d\n",yytext,linecount);
";" printf("%s : SEMICOLON at line %d\n",yytext,linecount);
"," printf("%s : COMMA at line %d\n",yytext,linecount);
"{" printf("%s : OPEN_BRACE at line %d\n",yytext,linecount);
"}" printf("%s : CLOSE_BRACE at line %d\n",yytext,linecount);
"(" printf("%s : OPEN_PAREN at line %d\n",yytext,linecount);
")" printf("%s : CLOSE_PAREN at line %d\n",yytext,linecount);
"do" printf("%s : DO at line %d\n",yytext,linecount);
"else" printf("%s : ELSE at line %d\n",yytext,linecount);
"for" printf("%s : FOR at line %d\n",yytext,linecount);
"if" printf("%s : IF at line %d\n",yytext,linecount);
"int" printf("%s : INTEGER at line %d\n",yytext,linecount);
"main" printf("%s : MAIN at line %d\n",yytext,linecount);
"printf" printf("%s : print at line %d\n", yytext, linecount);
"return" printf("%s : RETURN at line %d\n",yytext,linecount);
"scanf" printf("%s : SCAN at line %d\n",yytext, linecount);
"while" printf("%s : WHILE at line %d\n",yytext,linecount);
{HEADER} printf("%s : HEADER at line %d\n",yytext,linecount);
{VARNAME} printf("%s : VARIABLE NAME at line %d\n",yytext,linecount);
{DIGIT}+ printf("%s : NUMBER at line %d\n",yytext,linecount);
\n linecount++;
[ \t]+ ;
. printf("Unexpected character at line %d\n",linecount);

%%

int main(int argc, char *argv[])
{
	FILE *fd;
	if (argc == 2 && (fd = fopen(argv[1], "r")))
        	yyin = fd;
	else if(fd == NULL)
	{
		printf("File not found\n");
		exit(1);
	}
	yylex();
	printf("\nNumber of lines :  %d\n",linecount); 
	return 0;
}
	
