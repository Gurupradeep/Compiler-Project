#include <bits/stdc++.h>
using namespace std;

string token,func = "NULL"; 
stack < map<string,string> > s;

// NOTE : include "{ }" for every new block for it to be recognized as a separate block
// Need to fix scope of for loop iterators and functions' formal parameters

// Get token type 
int getType()
{
	if(token == "int")
		return 1;
	if(token == "+" || token == "-" || token == "*" || token == "/" || token == "%" || token == "&" || token == "|" || token == "||" || token == "&&" || token == ">" || token == "<" || token == ">=" || token == "<=" || token == "!" || token == "==" || token == "+=" || token == "-=" || token == "*=" || token == "/=" || token == "%=" || token == "=" || token == "?" || token == ":" || token == "++" || token == "--")
		return 2;
	if(token[0] == '$')
		return 3;
	if(token[0] == '@')
		return 4;
	if(token[0] == '#')
		return 5;
	return 6;
}

// Insert token
int insertToken(int TOKENID, map <string,string> &table)
{
	// Search for redeclaration
	/*if(table.find(token)!=table.end())
		return 1;*/

	// Insert if there are no conflicts
	switch(TOKENID)
	{
		case 1 : table[token] = "KEYWORD"; break; 
		case 2 : table[token] = "OPERATOR"; break;
		case 3 : table[token.substr(1,token.length())] = "CONSTANT"; break;
		case 4 : table[token.substr(1,token.length())] = "FUNCTION"; break;
		case 5 : table[token.substr(1,token.length())] = "MAIN_FUNCTION"; break;
		default : table[token] = "INTEGER_VARIABLE"; 
	}
	// table[token]=type;

	// Update the symbol table for the current block by poping the existing top and pushing the modified table
	s.pop();
	s.push(table);
	return 0;
}

// To be defined later
int lookupTable();

// Display the contents of the block
void dispBlock(map <string,string> table)
{
	cout << "\nNEW BLOCK" << endl;
	for(map<string,string>::iterator it = table.begin(); it != table.end(); it++)
		cout << "< " << it->first << " , " << it->second << " >" << endl;
}

// Load block into new symbol table 
int loadBlock(ifstream &FP)
{
	map <string,string> table;

	// Empty table
	s.push(table);
	
	// Hack for inserting function names 
	if(func != "NULL")
	{
		token = func; func = "NULL";
		int TOKENID = getType();
		int error = insertToken(TOKENID,table);
		if(error)
			return error;
	}
	while(true)
	{
		FP >> token;
		if(FP.eof())
			break;
		// EOF bit is set only after EOF character is read from the file. If kept in while loop, last character repeats
		// Link : http://stackoverflow.com/questions/21647/reading-from-text-file-until-eof-repeats-last-line
		
		// Nested block
		if(token == "{")
		{
			int error = loadBlock(FP);
			if(error)
				return error;
		}

		// Block goes out of scope
		else if(token == "}")
		{
			dispBlock(s.top());
			s.pop();
			break; 
		}

		// Insert into table
		else 
		{
			int TOKENID = getType();
			int error = insertToken(TOKENID,table);
			if(error)
				return error;
		}
	}
	return 0;
}

// Global scope
int loadSymbolTable(ifstream &FP)
{
	map <string,string> table;
	
	// Empty table
	s.push(table);

	while(true)
	{
		FP >> token;
		// EOF bit is set only after EOF character is read from the file. If kept in while loop, last character repeats
		// Link : http://stackoverflow.com/questions/21647/reading-from-text-file-until-eof-repeats-last-line
		if(FP.eof())
			break;

		// Hack for functions
		if(token[0] == '@' || token[0] == '#')
		{
			func = token;
			FP >> token;
		}

		// New block 
		if(token == "{")
		{
			int error = loadBlock(FP);
			if(error)
				return error;
		}

		// Insert into table
		else 
		{
			int TOKENID = getType();
			int error = insertToken(TOKENID,table);
			if(error)
				return error;
		}
	}

	// Destroy table
	dispBlock(s.top());
	s.pop();
	if(s.empty())
		cout << "\nStack Empty!\n";
	
	return 0;
}

int main()
{
	ifstream FP;
	FP.open("tokens.txt");
	printf("\nSYMBOL TABLE\n------------\n");
	int error = loadSymbolTable(FP);
	switch(error)
	{
		case 1 : cout << "Redeclaration of" << token << endl; exit(error); break;  
		case 2 : cout << "Undeclared variable" << token << endl; exit(error); break;  
	}
	FP.close();
	cout << endl;
	return 0;
}