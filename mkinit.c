// Anthony DiGirolamo
// Random Initial Condition Generator

#define INITIAL "initial_condition"
#define LIVE '1'
#define DEAD '0'

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[] ) {
	FILE *file;
	int r;	// number of rows
	int c;	// number of columns
	int row = 0, col = 0;	// array iterators

	if (argc != 4) {
		fputs("usage: [rows] [columns] [file name]\n", stderr);
		exit(EXIT_FAILURE);
	}
	r = (int) strtol(argv[1], NULL, 10);
	c = (int) strtol(argv[2], NULL, 10);

	file = fopen(argv[3], "w");
	if (file == NULL) {
		fputs(INITIAL, stderr);
		fputs(" could not be opened.\n", stderr);
		exit(EXIT_FAILURE);		
	}
	
	srand(time(0));
	
	// Make a text image
	for (row=0; row<r; row++) {
		for (col=0; col<c; col++) {
			if (rand()%2 == 0)
				fputc(DEAD, file);
			else
				fputc(LIVE, file);
		}
		//fputc('\n', file);
	}
	
	fclose(file);		
	
	return EXIT_SUCCESS;
}
