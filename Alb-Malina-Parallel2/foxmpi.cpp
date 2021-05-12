
//Libraries included
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <ctime>
#include<math.h>
#include"mpi.h"

//for iostream
using namespace std;

//structure for grid
typedef struct {
    int proc_count;
    int dim;
    int row;
    int col;
    int rank;
    //mpi grid,row and column variables
    MPI_Comm grid_comm;
    MPI_Comm row_comm;
    MPI_Comm col_comm;
} GridStructure;

//grid function 
void GridSetup(GridStructure *grid) {
    //required arrays for grid setup
	int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    int free_coords[2];
    int world_rank;

  //mpi size and rank for each processes
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->proc_count));
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   //grid and its dimension
    grid->dim = (int) sqrt((double) grid->proc_count);
    dimensions[0] = dimensions[1] = grid->dim;

    wrap_around[0] = 0;
    wrap_around[1] = 1;
    
    //mpi cart creating grid and its dimensions
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &(grid->grid_comm));
    MPI_Comm_rank(grid->grid_comm, &(grid->rank));
    MPI_Cart_coords(grid->grid_comm, grid->rank, 2, coordinates);
    grid->row = coordinates[0];
    grid->col = coordinates[1];

    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(grid->grid_comm, free_coords, &(grid->row_comm));

    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid->grid_comm, free_coords, &(grid->col_comm));
}

//Function for array multiplication
void MultiplyLocal(int **a, int **b, int **c, int size) {
    int temp = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            temp = 0;
            for (int k = 0; k < size; k++)
                temp += (a[i][k] * b[k][j]);
            c[i][j] += temp;
        }
}
//Function for unpack matrix and assign buffer value to array
void UnpackMatrix(int *buff, int **a, int size) {
    int k = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            a[i][j] = buff[k];
            k++;
        }
}

//Function for pack matrix and assign value to buffer
void PackMatrix(int *buff, int **a, int size) {
    int k = 0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            buff[k] = a[i][j];
            k++;
        }
}
//Function for printing the matrix
void PrintMatrix(int **matrix, int size) {
    for (int i = 0; i < size; i++) {
        std::cout << "|";
        for (int j = 0; j < size; j++) {
            int el = matrix[i][j];
            if (el < 10)
                std::cout << " ";
            std::cout << el;
            std::cout << "|";
        }
        std::cout << std::endl;
    }
}
//This function will print the pack matrix
void PrintPackedMatrix(int *matrix, int size) {
    for (int i = 0; i < size; i++) {
        std::cout << "|";
        for (int j = 0; j < size; j++) {
            int el = matrix[i * size + j];
            if (el < 10)
                std::cout << " ";
            std::cout << el;
            std::cout << "|";
        }
        std::cout << std::endl;
    }
}

//The function for generating matrix from taking elements from files
void GenerateMatrices(int **first_matrix,int **second_matrix,int MATRIX_SIZE) {
    
    //This is a object for accessing files 
	ofstream fout;
	
	//Generating File1 and input random numbers in it
	fout.open("File1.txt");
	
	//generating random numbers and putting in file1
	for(int i = 0; i < MATRIX_SIZE*MATRIX_SIZE; i++){
		int n = rand()%50+1;
		fout<<n<<endl;
	}
	//this function will close the file1
	fout.close();
	
	//Generating file2 and input random numbers in it
	fout.open("File2.txt");
    //generating random numbers and putting the numbers in array
	for(int i = 0; i < MATRIX_SIZE*MATRIX_SIZE; i++){
		int n = rand()%100+1;
		fout<<n<<endl;
	}
	//this function will close the file2
	fout.close();

    //object for opening the file
	ifstream f;

    //Opening a file1 named File1.txt to read its content
    f.open("File1.txt", ios::in);
	
	//these variables helps in getting numbers from files
	char ch;
	int n,i,j=0; 
    //Taking numbers from file1 and Initilizing first matrix
    for(i=0;i<MATRIX_SIZE;i++){
	
        for(j=0;j<MATRIX_SIZE;j++)
        {
        	f>>n;
        	first_matrix[i][j] = n;
		}
    }
	//This function will close the file1
	f.close();
    
    //Taking numbers from file2 and Initilizing second matrix 
    f.open("File2.txt", ios::in);
    n=0;
	j=0;
	
    for(i=0;i<MATRIX_SIZE;i++){
	
        for(j=0;j<MATRIX_SIZE;j++)
        {
        	f>>n;
        	second_matrix[i][j] = n;
		}
    }
    //This function will close the file2
    f.close();

}

//function for multily elements of array in fox algorithm
void FoxMultiply(int n, GridStructure *grid, int **a, int **b, int **c) {
    int **temp_a, *buff, stage, root, submat_dim, src, dst;
    MPI_Status status;

    submat_dim = n / grid->dim;

    temp_a = new int*[submat_dim];
    for(int i = 0; i < submat_dim; ++i)
        temp_a[i] = new int[submat_dim];
    for (int i = 0; i < submat_dim; i++)
        for (int j = 0; j < submat_dim; j++)
            temp_a[i][j] = 0;

    buff = new int[submat_dim*submat_dim];
    for (int i = 0; i < submat_dim * submat_dim; i++)
        buff[i] = 0;

    src = (grid->row + 1) % grid->dim;
    dst = (grid->row + grid->dim - 1) % grid->dim;
    //this loop will perform some operation in pack matrix in mpi bcast
    for (stage = 0; stage < grid->dim; stage++) {
        root = (grid->row + stage) % grid->dim;
        if (root == grid->col) {
            PackMatrix(buff, a, submat_dim);
            MPI_Bcast(buff, submat_dim * submat_dim, MPI_INT, root, grid->row_comm);
            UnpackMatrix(buff, a, submat_dim);
            MultiplyLocal(a, b, c, submat_dim);
        } else {
            PackMatrix(buff, temp_a, submat_dim);
            MPI_Bcast(buff, submat_dim * submat_dim, MPI_INT, root, grid->row_comm);
            UnpackMatrix(buff, temp_a, submat_dim);
            MultiplyLocal(temp_a, b, c, submat_dim);
        }
        //sendig buffer to processes
        PackMatrix(buff, b, submat_dim);
        MPI_Sendrecv_replace(buff, submat_dim * submat_dim, MPI_INT, dst, 0, src, 0, grid->col_comm, &status);
        UnpackMatrix(buff, b, submat_dim);
    }
}


//this function is for testing purpose
int** TestSingleThread(int **first_matrix,int **second_matrix,int MATRIX_SIZE)
{
    //initlizating some variables   
    int **test_result, **test_a, **test_b;
    test_a = new int*[MATRIX_SIZE];
    test_b = new int*[MATRIX_SIZE];
    test_result = new int*[MATRIX_SIZE];
    //loop for testing purpose
    for(int i = 0; i < MATRIX_SIZE; ++i) {
        test_a[i] = new int[MATRIX_SIZE];
        test_b[i] = new int[MATRIX_SIZE];
        test_result[i] = new int[MATRIX_SIZE];
    }
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++) {
            test_a[i][j] = first_matrix[i][j];
            test_b[i][j] = second_matrix[i][j];
            test_result[i][j] = 0;
        }
        //multiply local rest result
    MultiplyLocal(test_a, test_b, test_result, MATRIX_SIZE);

    std::cout << "\nLocal result:" << std::endl;
    //This function calls for printing the matrix
	PrintMatrix(test_result, MATRIX_SIZE);
    return test_result;
}


//this is main function that will start and all initlization of all variables and arrays and mpi will take place here
int main(int argc, char *argv[]) {
    
    //this will take any number from user and assign it to matrix_size
    int N1;
	cout<<"\nEnter Size of Matrix :";
	cin>>N1;
	int MATRIX_SIZE = N1;
	
	//Two-D matrix initlization
	int **first_matrix = new int*[MATRIX_SIZE];
	int **second_matrix = new int*[MATRIX_SIZE];
	for(int i = 0; i < MATRIX_SIZE; ++i)
    {
    	//it will assign memory to two D matrix
        first_matrix[i] = new int[MATRIX_SIZE];
        second_matrix[i] = new int[MATRIX_SIZE];
    
    }
    
    //these are some local arrays that will help at different process to multiply matrix
	int block_size;
    int **local_a, **local_b, **local_c;
    //this is for time initlization
	clock_t start,finish;
    //mpi initlization
    MPI_Init(&argc, &argv);
   //grid initlization
    GridStructure grid;
    GridSetup(&grid);
    //call function for generating two matrix
    GenerateMatrices(first_matrix,second_matrix,N1);

    block_size = MATRIX_SIZE / grid.dim;
    int base_row = grid.row * block_size;
    int base_col = grid.col * block_size;

    //some local arrays help in individual
    local_a = new int*[MATRIX_SIZE];
    local_b = new int*[MATRIX_SIZE];
    local_c = new int*[MATRIX_SIZE];
    for(int i = 0; i < MATRIX_SIZE; ++i)
    {
    	//assigning memory to local arrays
        local_a[i] = new int[MATRIX_SIZE];
        local_b[i] = new int[MATRIX_SIZE];
        local_c[i] = new int[MATRIX_SIZE];
    }
    //assigning elements to local arrays from two 2-D matrix
    for (int i = base_row; i < base_row + block_size; i++)
        for (int j = base_col; j < base_col + block_size; j++) {
            local_a[i - (base_row)][j - (base_col)] = first_matrix[i][j];
            local_b[i - (base_row)][j - (base_col)] = second_matrix[i][j];
            local_c[i - (base_row)][j - (base_col)] = 0;
        }

    if (grid.rank == 0)
    {
        std::cout<<"Ready..." <<std::endl;
    }

    //Barrier for synchronization of all processes
    MPI_Barrier(grid.grid_comm);
    if (grid.rank == 0)
    {
        start = clock();
    }
    //call function for multiply
    FoxMultiply(MATRIX_SIZE, &grid, local_a, local_b, local_c);
    //mpi barrier for grid 
	MPI_Barrier(grid.grid_comm);
    if (grid.rank == 0)
    {
        finish = clock();
        clock_t result_time = finish - start;
        std::cout << "Time: " << double(result_time)/CLOCKS_PER_SEC << std::endl;
    }

    int *result_buff = new int[MATRIX_SIZE * MATRIX_SIZE];
    int *local_buff = new int[block_size * block_size];
    PackMatrix(local_buff, local_c, block_size);
    
     //mpi gather  
    MPI_Gather(local_buff, block_size * block_size, MPI_INT, result_buff, block_size * block_size, MPI_INT, 0, grid.grid_comm);
    MPI_Barrier(grid.grid_comm);
    if (grid.rank == 0) {
        int *data = new int[MATRIX_SIZE * MATRIX_SIZE];
        int k = 0;
        for (int bi = 0; bi < grid.dim; bi++)
            for (int bj = 0; bj < grid.dim; bj++)
                for (int i = bi * block_size; i < bi * block_size + block_size; i++)
                    for (int j = bj * block_size; j < bj * block_size + block_size; j++) {
                        data[i * MATRIX_SIZE + j] = result_buff[k];
                        k++;
                    }

        std::cout << "\nFox result check:\n";

        //it will call funciton to print packet matrix
        PrintPackedMatrix(data, MATRIX_SIZE);
        int ** res = TestSingleThread(first_matrix,second_matrix,N1);
        for(int i=0;i<MATRIX_SIZE;i++)
            for(int j=0;j<MATRIX_SIZE;j++)
            {
                if(res[i][j]!=data[i*MATRIX_SIZE+j])
                {
                    std::cout << "ERROR!" <<std::endl;
                    exit(1);
                }
            }
        std::cout << "OK" << std::endl;
    }
    //mpi will finalized here
    MPI_Finalize();
    exit(0);
    //end of program
}		
