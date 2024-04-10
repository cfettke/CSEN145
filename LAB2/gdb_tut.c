#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv){
	float sum = 0;
	int i = 0;
	for(i = 0; i < argc; i++){
		sum += atof(argv[i]);
		//printf("%f", (float)(*argv[i]));
	}
	printf("Sum: %f\n", sum);
	return 0;
}
