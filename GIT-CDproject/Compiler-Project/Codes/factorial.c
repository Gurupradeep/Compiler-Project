#include <stdio.h>

int main()
{
	int n;
	long long fact = 1;
	scanf("%d",&n);
	//Single line comment
	for(int i = 1; i < n + 1; i++)
	{
		int x;
		fact *= -i;
	}
	/*M
	ul
	ti
	line
	comment*/
	printf("Factorial : %lld",fact);
	return 0;
}