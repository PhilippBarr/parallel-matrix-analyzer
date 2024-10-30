# MPI Parallel Diagonally Dominant Matrix Calculator


# Overview
This program uses MPI (Message Passing Interface) to distribute and process a 2D 
matrix among multiple processes.

The program performs the following tasks:
- Distributes the matrix among multiple processes.
- Checks if the matrix is diagonally dominant.
- Finds the maximum element of the diagonal of the matrix.
- Constructs a new matrix `B` based on the maximum element `m` (explained below).
- Finds the minimum element and its location from the matrix `B`

# Main Functions

`int* makeArray(int n)`
Allocates memory for an array of size `n` and returns a pointer to the array.

**Parameters:**
- `n`: Number of elements to allocate memory for.



`int *findSendCounts(int sendHeight, int sendWidth, int n_process)`

Calculates and returns an array of integers that indicates how many elements each process will recieve
when the data is distributed across `n_process` processes.

**Parameters:**

- `sendHeight`: The number of rows in the 2D data that needs to be distributed.
- `sendWidth`: The number of columns in the 2D data.
- `n_process`: The total number of processes to distribute the data across.

### Example:

**Suppose:**

- `sendHeight` = 10 (10 rows)
- `sendWidth` = 4 (4 columns)
- `n_process` = 3 (3 processes)

The function will assign `4 rows` to the first process (1 extra row) and
assign `3 rows` to the remaining two processes.

The final `sendcounts` array will look like this :
`sendcounts = [16, 12, 12]`


`int *findDispls(int n_process, int *sendcounts)`

Calculates the displacements (starting positions) for each process, indicating
from which index each process will start receiving its portion of the data
when the data is scattered among multiple processes.

###  Parameters:
- `n_process`: The total number of processes.
- `sendcounts`: The array calculated in `findSendCounts` function.

### Example:

 **Suppose**:

- `n_process = 3` (3 processes)
- `sendcounts = [16, 12, 12]` (the number of data elements each process will receive)

The final displacement array will be:

`displs = [0, 16, 28]`

This means:
- `Process 0` will start from index `0`
- `Process 1` will start from index `16`
- `Process 2` will start from index `28`


`void shareArray(int *sendArray, int sendHeight, int sendWidth, int **recvArray, int *recvHeight, int *recvWidth, int root, MPI_Comm comm)
`
This function distributes a the 2D array across multiple processes.

- `sendArray`: A pointer to the 2D array (matrix) in 1D form that will be distributed across the processes. This is located on the root process.
- `sendHeight`: The number of rows of the original matrix sendArray (total height).
- `sendWidth`: The number of columns in the original matrix (width).
- `recvArray`: A pointer to a pointer that will hold the chunk of the matrix received by each process.
- `recvHeight`: A pointer to an integer that will store the number of rows in the chunk that each process will receive.
- `recvWidth`: A pointer to an integer that will store the number of columns (the width) of the chunk each process will receive (this is the same for all processes).
- `root`: The rank of the root process (the process holding the full matrix).
- `comm`: The MPI communicator.


Let’s say we have 2 processes (processes 0 and 1), and we want to split this matrix between them:

Process 0 (root):

The matrix is split into two parts:

`Process 0` gets the first two rows `{7, 1, 2, 3}` and `{0, 5, 2, 2}`.

`Process 1` gets the last two rows `{2, 2, 8, 3}` and `{6, 7, 1, 15}`.

`Process 0` (root) initializes `sendcounts = {8, 8}` and `displs = {0, 8}` (each process gets 8 elements).


**From the initial array:**

`Process 0` gets `7, 1, 2, 3, 0, 5, 2, 2` (which corresponds to the first two rows).

`Process 1` gets `2, 2, 8, 3, 6, 7, 1, 15` (which corresponds to the last two rows).






`int giveFlag(int *A, int grammes, int stiles, int pos)`
Checks if a chunk of the matrix is diagonally dominant

**Parameters**
- `A`: The matrix chunk.
- `grammes`: Number of rows in the chunk.
- `stiles`: Number of columns in the chunk.
- `pos`: Starting position for checking diagonal dominance.



`int isItDesp(int *pinakas, int grammes, int stiles, 
int sendHeight, int sendWidth, int root, MPI_Comm comm)`

Determines if the matrix is diagonally dominant.

###  Parameters:

- `pinakas`: The matrix chunk.
- `grammes`: Number of rows in the chunk.
- `stiles`: Number of columns in the chunk.
- `sendHeight`: Height of the full matrix.
- `sendWidth`: Width of the full matrix.
- `root`: Rank of the root process.
- `comm`: The MPI communicator.

`int giveMax(int *A, int grammes, int stiles, int pos)`

Returns the maximum diagonal element from a matrix

###  Parameters:
- `A`: The matrix.
- `grammes`: Number of rows.
- `stiles`: Number of columns.
- `pos`: Starting position for checking the diagonal.

`int getMax(int *pinakas, int grammes, int stiles, int sendHeight, int sendWidth,
int root, MPI_Comm comm)`

Finds the maximum element of an array for the current process

###  Parameters:
- `pinakas`: The matrix.
- `grammes`: Number of rows.
- `stiles`: Number of columns.
- `pos`: Starting position for checking the diagonal.
- `root`: Rank of the root process.
- `comm`: The MPI communicator


`int *partOfBFind(int *pinakas, int grammes, int stiles, int pos, int max) `

This function constructs a chunk of matrix B from a given chunk of matrix A. The chunk is
built using the formula:

- `Bij = m - |Aij| for i != j` where `m` is the maximum diagonal element of matrix `A`
- `Bii = m for diagonal elements i = j`

###  Parameters:
- `pinakas`: The chunl pf matrix `A`
- `grammes`: Number of rows in the chunk of `A`
- `stiles`: Number of columns in the chunk of `A`
- `pos`: The starting position of the diagonal element
- `max` The maximum diagonal element `m` from matrix `A`


`int *givePartOfB(int *pinakas, int grammes, int stiles, int sendHeight, int sendWidth, int max, MPI_Comm comm)`

Calculates a chunck of matrix `B` for the current process on the corresponding
matrix `A` and the maximum diagonal element `m`. It internally uses `partOfBFind()`

###  Parameters:
- `pinakas`: The chunk of matrix `A` (received by the process).
- `grammes`: Number of rows in the chunk of  `A`. 
- `stiles`: Number of columns in the chunk of  `A`.
- `sendHeight`: The total height of matrix  `A`.
- `sendWidth`: The total width of matrix  `A`.
- `max`: The maximum diagonal element `m` from matrix `A`.
- comm: The MPI communicator.

`int *giveBackB(int *C, int h, int w, MPI_Comm comm)`

Gathers all the chunks of matrix `B` from the processes and merges them into a complete
matrix on the `root` process.

 ### Parameters:

- `C`: The "chunk" (this is why we named it `C`) of matrix `B` for the current process.
- `h`: The number of rows in the chunk of matrix `B`.
- `w`: The number of columns in the chunk of matrix `B`.
- `comm`: The MPI communicator.

`int giveMinLoc(int *A, int size)`

Returns the position of the minimum element contained in a array.

### Parameters:
-  `A`: Integer array
- `size`: The total number of elements in the array.

Finds the minimum location

`MPI_Scatter(&sendcounts[0], 1, MPI_INT, &smallSize, 1, MPI_INT, root, comm);`

Used to distribute data from one process to all others.

###  Parameters:

- `&sendcounts[0]`: Address of the first element of the sendcounts array
- `1`: Number of elements to send from `sendcounts` to each process
- `MPI_INT`: Data type of the element to be sent
- `&smallSize`: Here the recieved value will be stored
- `1`: Number of elements to be recieved
- `MPI_INT`: Data type of the element to be recieved
- `root`: Rank of the root process 
- `comm`: Communicator

**Example:**

`sendcounts = [16, 12, 12]`

`Process 0` will recieve `sendcounts[0] = 16` and store it in `smallSize`

`Process 1` will recieve `sendcounts[1] = 12` and store it in `smallSize`

`Process 2` will recieve `sendcounts[2] = 12` and store it in `smallSize`

`MPI_Scatterv(sendArray, sendcounts, displs, MPI_INT, *recvArray, smallSize, MPI_INT, root, comm);`

Unlike `MPI_Scatter()`, it sends different amounts of data to each process

# Run
To run this program i am using Linux. You need to have `Open MPI` installed
on your machine.

### To run the program:

- Initialize the 2D array inside the `parallel-matrix-analyzer.c` file:

We will initialize it with those values

| 7     | 1     | 2  | 3  |
|-------|-------|----|----|
| **0** | **5** | **2**  | **2**  |
| **2**     | **2**     | **8**  | **3**  |
| **6**     | **7**     | **1**  | **15** |


### On terminal:
1) `mpicc -o parallel-matrix-analyzer parallel-matrix-analyzer.c`
2) `mpirun -np 2 ./parallel-matrix-analyzer`

- The number `2` on the second terminal command is the number of processes
you would like to run the program with.

# Results:


**Diagonally Dominant**

**Final Max**: 15

### Full Matrix B:
| 15     | 14 | 13 | 12 |
|--------|----|----|----|
| **15** | **15** | **13** | **13** |
| **13**     | **13** | **15** | **12** |
| **9**      |  **8** | **14** | **15** |

**Min**: 8

**Min Location**:  
- **i** = 3  
- **j** = 1



