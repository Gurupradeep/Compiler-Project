#include <stdio.h>

int func(int a,int b)
{
	return 0;
}

int main()
{
	// Single Line Comment
	int var;
	for(var = 10; var <= 100; var ++)
	{
		int temp = 0;
		while(temp < var)
		{
			temp ++ ;
		}

		/* Multi-line
		Comment */

		do
		{
			temp += 10 * var;
			if(temp > 100)
				printf("%d", temp);
			else
				printf("%d", var); 
		} while(temp < 10000);
	}
	return 0;
}