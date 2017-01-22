%{
	#include <stdlib.h>
	int linecount = 1;
	FILE* LEXOUT;
%}

DIGIT [0-9]
FUNCTION [_a-zA-Z][a-zA-Z0-9]*\([ ]*int[ ]+{VARNAME}[ ]*,[ ]*int[ ]+{VARNAME}[ ]*\) 
INVALID [^(' '|\t)]*[^[_a-zA-Z0-9=+?\-*/%!:;'"><\|]]+[^(' '|\t)]* 
HEADER <[a-zA-Z]+\.h>
PREPROCESSOR #[a-zA-Z]+
SINGLE_LINE_COMMENT \/\/.*
START_MULTI_LINE \/\*
STRING \".+\"
VARNAME [_a-zA-Z][_a-zA-Z0-9]* 

%s startComment

%%

{SINGLE_LINE_COMMENT} printf("%s : COMMENTS at line %d\n",yytext,linecount);
{START_MULTI_LINE} 	{BEGIN startComment; printf("START COMMENT at line %d\n",linecount);}
<startComment>\n 	{linecount++;}
<startComment>\*\/	{BEGIN 0; printf("END COMMENT at line %d\n",linecount);}
<startComment>[^(*/|\n)]*	;

{FUNCTION}	{
				printf("%s FUNCTION\n",yytext); 
				int i = 0; char func[100]; 
			 	while(yytext[i] != '(')
			 	{
			 		func[i] = yytext[i];
			 		i++;
			 	}
			 	func[i] = '\0';
			 	fprintf(LEXOUT,"@%s\n",func);
			}

{PREPROCESSOR} printf("%s : PREPROCESSOR DIRECTIVE at line %d\n",yytext,linecount);
{STRING} printf("%s : STRING\n",yytext);

"="	{
		printf("%s : ASSIGNMENT at line %d\n ",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"+"	{
		printf("%s : ADD at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"-"	{
		printf("%s : SUBTRACT at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"*" {
		printf("%s : MULTIPLY at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"/" {
		printf("%s : DIVIDE at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"%" {
		printf("%s : MODULO at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"++"	{
			printf("%s : INCREMENT at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"--"	{
			printf("%s : DECREMENT at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}


"+=" 	{
			printf("%s : ADD_ASSIGN at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"-="	{
			printf("%s : SUB_ASSIGN at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"*=" 	{
			printf("%s : MULT_ASSIGN at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"/="	{
			printf("%s : DIV_ASSIGN at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"%="	{
			printf("%s : MOD_ASSIGN at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"=="	{
			printf("%s : EQUAL at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

">="	{
			printf("%s : GREAT_EQ at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"<="	{
			printf("%s : LESS_EQ at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}
">"	{
		printf("%s : GREAT at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"<"	{
		printf("%s : LESS at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"?"	{
		printf("%s : CONDITIONAL at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"&"	{
		printf("%s : REFERENCE/BIT_AND at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"|"	{
		printf("%s : BIT_OR at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"||"	{
			printf("%s : LOGIC_OR at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"&&"	{
			printf("%s : LOGIC_AND at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"!"	{
		printf("%s : NEGATE at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

":"	{
		printf("%s : COLON at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

";"	printf("%s : SEMICOLON at line %d\n",yytext,linecount);
","	printf("%s : COMMA at line %d\n",yytext,linecount);

"{"	{
		printf("%s : OPEN_BRACE at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"}"	{
		printf("%s : CLOSE_BRACE at line %d\n",yytext,linecount);
		fprintf(LEXOUT,"%s\n",yytext);
	}

"("	printf("%s : OPEN_PAREN at line %d\n",yytext,linecount);
")" printf("%s : CLOSE_PAREN at line %d\n",yytext,linecount);
"do" printf("%s : DO at line %d\n",yytext,linecount);
"else" printf("%s : ELSE at line %d\n",yytext,linecount);
"for" printf("%s : FOR at line %d\n",yytext,linecount);
"if" printf("%s : IF at line %d\n",yytext,linecount);

"int"	{
			printf("%s : INT at line %d\n",yytext,linecount);
			fprintf(LEXOUT,"%s\n",yytext);
		}

"int main"	{
				printf("%s : INT at line %d\n","int",linecount);
				printf("%s : MAIN at line %d\n","main",linecount);
				fprintf(LEXOUT,"#%s\n","main");
			}

"printf" printf("%s : PRINT at line %d\n",yytext,linecount);
"return" printf("%s : RETURN at line %d\n",yytext,linecount);
"scanf"	printf("%s : SCAN at line %d\n",yytext,linecount);
"while" printf("%s : WHILE at line %d\n",yytext,linecount);

{HEADER} printf("%s : HEADER at line %d\n",yytext,linecount);

{VARNAME}	{
				printf("%s : VARIABLE NAME at line %d\n",yytext,linecount);
				fprintf(LEXOUT,"%s\n",yytext);
			}

{DIGIT}+	{
				printf("%s : NUMBER at line %d\n",yytext,linecount);
				fprintf(LEXOUT,"$%s\n",yytext);
			}

\n linecount++;
[ \t]+ ;
{INVALID}	printf("Unexpected character/s at line %d\n",linecount);

%%

int main(int argc, char *argv[])
{
	FILE* SRC;
	if (argc == 2 && (SRC = fopen(argv[1],"r")))
        	yyin = SRC;
	else if(SRC == NULL)
	{
		printf("File not found\n");
		exit(1);
	}
	if(!(LEXOUT = fopen("tokens.txt","w")))
	{
		printf("File not found\n");
		exit(2);
	}
	yylex();
	printf("\nNumber of lines :  %d\n",linecount); 
	return 0;
}
	
