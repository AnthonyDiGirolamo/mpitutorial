// Game Of Life Training Example
// Anthony DiGirolamo - anthony.d@asu.edu

/* This file is part of an MPI training course.
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2009 Anthony DiGirolamo
 */

#define INITIAL "initial_condition"
#define CHECKPOINT "checkpoint"
#define FINAL_RESULTS "final_results"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

#define CHECKMALLOC(var) if((var) == NULL) { printf("ERROR: malloc()\n"); abort(); }

void printhelp() {
	printf("Usage: life -i [ITERATIONS] -r [ROWS] -c [COLS]\n");
	printf("  -h, --help        Display this message\n");
	printf("Required: \n");
	printf("  -r  --rows        [integer]    Number of rows in the grid\n");
	printf("  -c  --columns     [integer]    Number of columns in the grid\n");
	printf("  -i  --iterations  [integer]    Number of iterations to perform\n");
	printf("Optional: \n");
	printf("  -C  --checkpoint  [integer]    Checkpoint after a number of iterations\n");
	printf("  -R  --resume      [file name]  Resume from a given checkpoint\n");
	printf("  -I  --initial     [file name]  Use a specific initial condition file\n");
	printf("  -O  --output      [file name]  Final output file name\n");
	printf("Default file names:\n");
	printf("  Initial condition:  %s\n", INITIAL);
	printf("  Checkpoints:        %s-[ITERATION]\n", CHECKPOINT);
	printf("  Final results:      %s\n", FINAL_RESULTS);
}


int exists(const char *filename) {
	return !access(filename, F_OK);
}


int liveOrDie(int cell, int neighbors) {
	cell -= 48; neighbors -= 384; // convert values from chars to ints, (48 = 0, 49 = 1)
	int newCell = 0;

	if (cell == 0) {
		if (neighbors == 3)
			newCell = 1;
	}
	else {
		if (neighbors == 2 || neighbors == 3)
			newCell = 1;
		else if (neighbors < 2 || neighbors > 3)
			newCell = 0;
	}

	return newCell+48; // back into a character
}


void gameOfLife(char **array1, char **array2, int rows, int cols) {
	int r, c;
	for (r = 1; r < rows+1; r++) {
		for (c = 1; c < cols+1; c++) {
			array2[r][c] = liveOrDie( array1[r][c],
					array1[r-1][c-1] + array1[r-1][c] + array1[r-1][c+1] +
					array1[r][c-1]   +                  array1[r][c+1]   +
					array1[r+1][c+1] + array1[r+1][c] + array1[r+1][c+1] );
		}
	}
}


int main(int argc, char *argv[] ) {

	// Initialize MPI
	int rank, processors;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &processors);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char *iterations_string = NULL;
	char *checkpoints_string = NULL;

	int making_checkpoints = 0; // are we making checkpoints?
	int checkpoint_resume = 0;	// resume from last checkpoint?

	char *rows_string = NULL;
	char *columns_string = NULL;
	char *initial_condition_name = NULL;
	char *checkpoint_name = NULL;
	char *final_output_name = NULL;

	// Parse the Command Line
	int opt = 0;
	while (opt != -1) {
		static struct option long_options[] = {
			{"help",		no_argument, 0, 'h'},
			{"rows",		required_argument, 0, 'r'},
			{"columns",		required_argument, 0, 'c'},
			{"iterations",	required_argument, 0, 'i'},
			{"checkpoint",	required_argument, 0, 'C'},
			{"resume",		required_argument, 0, 'R'},
			{"initial",		required_argument, 0, 'I'},
			{"output",		required_argument, 0, 'O'},
			{0, 0, 0, 0}
		};
		int option_index = 0;
		opt = getopt_long (argc, argv, "hr:c:i:C:R:I:O:", long_options, &option_index);
		if (opt == -1) { break; }

		switch (opt) {
		case 'h':
			if (rank==0) printhelp();
			MPI_Finalize();
			exit(EXIT_FAILURE);
			break;
		case 'i':
			iterations_string = optarg;
			break;
		case 'r':
			rows_string = optarg;
			break;
		case 'c':
			columns_string = optarg;
			break;
		case 'C':
			checkpoints_string = optarg;
			making_checkpoints = 1;
			break;
		case 'R':
			checkpoint_name = optarg;
			checkpoint_resume = 1;
			break;
		case 'I':
			initial_condition_name = optarg;
			break;
		case 'O':
			final_output_name = optarg;
			break;
		default:
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}

	if (iterations_string == NULL || rows_string == NULL || columns_string == NULL) {
		if (rank == 0)
			printf("Iterations, Rows, and Columns must be specified on the command line.\nRun life -h for more information.");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	if (initial_condition_name) {
		if (!exists(initial_condition_name)) {
			if (rank == 0)
				printf("The file \"%s\" cannot be found.\nRun life -h for more information.", initial_condition_name);
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}
	else {
		if (!exists(INITIAL)) {
			if (rank == 0)
				printf("The file \"%s\" cannot be found.\nRun life -h for more information.", INITIAL);
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}

	if (checkpoint_resume && checkpoint_name) {
		if (!exists(checkpoint_name)) {
			if (rank == 0)
				printf("The file \"%s\" cannot be found.\nRun life -h for more information.", checkpoint_name);
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}

	int iterations;				// number of iterations
	int checkpoint_after;		// create a checkpoint after this many iterations
	int checkpoint_counter;
	char *checkpoint_fullname;
	checkpoint_fullname = (char*) malloc(128 * sizeof(char));
	CHECKMALLOC(checkpoint_fullname);
	int i = 0;
	//int c, row=0, col=0;

	int grid_size[2] 		= {0,0}; // total size of the distributed array
	int proc_size[2] 		= {0,0}; // size of the communicator grid
	int local_size[2] 		= {0,0}; // size of the local chunk, no ghost rows
	int remainder_size[2] 	= {0,0}; // if the proc array doesn't divide the grid_size evently
	int coords[2] 			= {0,0}; // this ranks location in the proc grid
	int start_indices[2] 	= {0,0}; // this ranks chunk start indicies in the grid_size
	int periods[2] 			= {1,1}; // periodic proc grid, donut world!
	int mem_size[2]			= {0,0}; // local_size + ghost rows

	// This proc's local arrays. grid1 and grid2 for "double buffering"
	char **grid1;			// subscriptable array handle
	char *grid1_pointer;	// head of the array, for MPI functions
	char **grid2;
	char *grid2_pointer;
	char **temp1, *temp2;

	MPI_Status status;
	MPI_Datatype filetype, memtype;
	MPI_File file_handle;
	int file_open_error;

  	// Determine iterations
	iterations = (int) strtol(iterations_string, NULL, 10);
	if(rank==0) printf("SETUP: %d iterations\n", iterations);

	if (making_checkpoints) {
		checkpoint_after = (int) strtol(checkpoints_string, NULL, 10);
		if(rank==0) printf("SETUP: Checkpoint every %d iterations\n", checkpoint_after);
	}

	// Get the total grid size
	grid_size[0] = (int) strtol(rows_string, NULL, 10);
	grid_size[1] = (int) strtol(columns_string, NULL, 10);
	if (rank==0) printf("SETUP: grid_size: %d, %d\n", grid_size[0], grid_size[1]);

	// Create the communicator size
	MPI_Dims_create(processors, 2, proc_size);
	if (rank==0) printf("SETUP: proc_size: %d, %d\n", proc_size[0], proc_size[1]);

	// Determine the local chunk size
	local_size[0] = grid_size[0] / proc_size[0];
	local_size[1] = grid_size[1] / proc_size[1];
	if (rank==0) printf("SETUP: local_size: %d, %d\n", local_size[0], local_size[1]);

	// Check to see that grid_size divides evenly, otherwise exit
	remainder_size[0] = grid_size[0] % proc_size[0];
	remainder_size[1] = grid_size[1] % proc_size[1];
	if (rank==0) printf("SETUP: remainder_size: %d, %d\n", remainder_size[0], remainder_size[1]);
	if (remainder_size[0] != 0 || remainder_size[1] != 0) {
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	// Setup the communicator
	MPI_Comm comm;
	MPI_Cart_create(MPI_COMM_WORLD, 2, proc_size, periods, 0, &comm);
	MPI_Comm_rank(comm, &rank);
	MPI_Cart_coords(comm, rank, 2, coords);

	// Determine each procs starting index withing the global grid
	start_indices[0] = coords[0] * local_size[0];
	start_indices[1] = coords[1] * local_size[1];

	printf("RANK%d - Proc Location: (%d, %d) - Start Index: (%d, %d)\n",
		rank, coords[0], coords[1], start_indices[0], start_indices[1]);
	fflush(stdout);

	// Create an MPI Datatype for IO
	MPI_Type_create_subarray(2, grid_size, local_size, start_indices, MPI_ORDER_C, MPI_CHAR, &filetype);
	MPI_Type_commit(&filetype);

	// Find each rank's neighbor in the communicator
	int row_minus=0, row_plus=0, col_minus=0, col_plus=0, neighbors[4] = {0,0,0,0};

	MPI_Cart_shift(comm, 0, 1, &row_minus, &row_plus); // shift row +1
	MPI_Cart_shift(comm, 1, 1, &col_minus, &col_plus); // shift col +1

	enum neighbor { north, west, east, south };
	neighbors[north] = row_plus;
	neighbors[west] = col_plus;
	neighbors[east] = col_minus;
	neighbors[south] = row_minus;

	// Now let's set up an array with ghost rows
	mem_size[0] = local_size[0] + 2;
	mem_size[1] = local_size[1] + 2;
	start_indices[0] = start_indices[1] = 1;

	MPI_Type_create_subarray(2, mem_size, local_size, start_indices, MPI_ORDER_C, MPI_CHAR, &memtype);
	MPI_Type_commit(&memtype);

	// Open the checkpoint
	if (checkpoint_resume) {
		if (checkpoint_name) {
			file_open_error = MPI_File_open(MPI_COMM_WORLD, checkpoint_name,
				MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file_handle);
		}
	}
	else { // Open the init file
		if (initial_condition_name) {
			file_open_error = MPI_File_open(MPI_COMM_WORLD, initial_condition_name,
				MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file_handle);
		}
		else {
			file_open_error = MPI_File_open(MPI_COMM_WORLD, INITIAL,
				MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file_handle);
		}
	}

	if (file_open_error != MPI_SUCCESS) {
		printf("MPI_File_open: The initial condition could not be opened.");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	MPI_File_set_view(file_handle,0, MPI_CHAR, filetype, "native", MPI_INFO_NULL);

	// Allocate the local arrays
	grid2_pointer = (char *)  malloc(mem_size[0] * mem_size[1] * sizeof(char));
	CHECKMALLOC(grid2_pointer);
	grid2         = (char **) malloc(mem_size[0] * sizeof(char*));
	CHECKMALLOC(grid2);
	for(i = 0; i < mem_size[0]; i++)
		grid2[i] = &grid2_pointer[i * mem_size[1]];

	grid1_pointer = (char *)  malloc(mem_size[0] * mem_size[1] * sizeof(char));
	CHECKMALLOC(grid1_pointer);
	grid1         = (char **) malloc(mem_size[0] * sizeof(char*));
	CHECKMALLOC(grid1);
	for(i = 0; i < mem_size[0]; i++)
		grid1[i] = &grid1_pointer[i * mem_size[1]];

	// Populate the local arrays

	// for(int row=0; row<mem_size[0]; row++) {
	// 	for(int col=0; col<mem_size[1]; col++) {
	// 		if (rand()%2 == 0)
	// 			c = '0';
	// 		else
	// 			c = '1';
	// 		grid1[row][col] = grid2[row][col] = c;
	// 	}
	// }

	MPI_File_read_all(file_handle, grid1_pointer, 1, memtype, &status);
	MPI_File_close(&file_handle);

	// Construct a column data type
	MPI_Datatype col_type;
	MPI_Type_vector(local_size[0], 1, mem_size[1], MPI_CHAR, &col_type);
	MPI_Type_commit(&col_type);

	MPI_Barrier(comm); // everyone done?

	// start the iteration loop

	srand((unsigned int)time(NULL));

	if (making_checkpoints)
		checkpoint_counter = checkpoint_after;

	for (i = 0; i < iterations; i++) {

		// send/recv top row
		MPI_Sendrecv(&grid1[1][1], local_size[1], MPI_CHAR, neighbors[north], 0,
			&grid1[local_size[0]+1][1], local_size[1], MPI_CHAR, neighbors[south], 0, comm, &status);

		// send/recv bottom row
		MPI_Sendrecv(&grid1[local_size[0]][1], local_size[1], MPI_CHAR, neighbors[south], 0,
			&grid1[0][1], local_size[1], MPI_CHAR, neighbors[north], 0, comm, &status);

		// send/recv left column
		MPI_Sendrecv(&grid1[1][1], 1, col_type, neighbors[west], 0,
			&grid1[1][local_size[1]+1], 1, col_type, neighbors[east], 0, comm, &status);

		// send/recv right column
		MPI_Sendrecv(&grid1[1][local_size[1]], 1, col_type, neighbors[east], 0,
			&grid1[1][0], 1, col_type, neighbors[west], 0, comm, &status);

		// run the game of life

		gameOfLife(grid1, grid2, local_size[0], local_size[1]);

		// swap the arrays
		temp1 = grid1;
		grid1 = grid2;
		grid2 = temp1;

		temp2 = grid1_pointer;
		grid1_pointer = grid2_pointer;
		grid2_pointer = temp2;

		if (making_checkpoints) {
			checkpoint_counter--;
			// check to see if this iteration needs a checkpoint
			if (checkpoint_counter == 0) {
				checkpoint_counter = checkpoint_after;

				sprintf(checkpoint_fullname, "%s-%d", CHECKPOINT, i);
				MPI_File_open(MPI_COMM_WORLD, checkpoint_fullname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handle);

				MPI_File_set_view(file_handle, 0, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
				MPI_File_write_all(file_handle, grid1_pointer, 1, memtype, &status);
				MPI_File_close(&file_handle);

				if (rank == 0)
					printf("Iteration %d - Checkpoint made: %s\n", i, checkpoint_fullname);
			}
		} // end making_checkpoints?
	} // end iteration loop

	// All done! write out final result
	if (final_output_name)
		MPI_File_open(MPI_COMM_WORLD, final_output_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handle);
	else
		MPI_File_open(MPI_COMM_WORLD, FINAL_RESULTS, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_handle);

	MPI_File_set_view(file_handle, 0, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
	MPI_File_write_all(file_handle, grid1_pointer, 1, memtype, &status);
	MPI_File_close(&file_handle);

	if (rank == 0)
		printf("Final Results made: Iteration %d\n", i-1);

	free(grid1);
	free(grid1_pointer);
	free(grid2);
	free(grid2_pointer);

	MPI_Finalize();
	return EXIT_SUCCESS;
}
