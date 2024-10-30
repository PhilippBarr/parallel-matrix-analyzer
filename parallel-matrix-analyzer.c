#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//Creating array of n size
int *makeArray(int n)
{
	int *array, i, x;
	array = (int *) malloc(n * sizeof(int));
	return array;
}

//Calculate and return array of integers that indicate how many elements each process will recieve.
int *findSendCounts(int sendHeight, int sendWidth, int n_process)
{
	int i, size;
	int *sendcounts;
	sendcounts = makeArray(n_process);//Create 1D array same size as the number of total processes.
	size = sendHeight / n_process + 1; //We know that sendHeight % n_procces number of processes will take sendHeight/n_procces +1 possitions.
	for (i = 0; i < sendHeight % n_process; i++) {
		sendcounts[i] = size * sendWidth;// Multiplying it with sendWidth because it is 2D and inserting the rows inside.
	}
	size--;//We know that n_procces-sendHeight % n_procces processes will take sendHeight/n_procces this is why we do size--.
	for (i = sendHeight % n_process; i < n_process; i++) {
		sendcounts[i] = size * sendWidth;//size = sendHeight/n_process.
	}
	return sendcounts;//return sendcounts array.
}

//Calculates starting positions for each process.
int *findDispls(int n_process, int *sendcounts)
{
	int i;
	int *displs;//Will store the displacement (starting index) for each process.
	displs = makeArray(n_process);
	displs[0] = 0;//The first process starts receiving data from the begining (index 0).
	for (i = 1; i < n_process; i++) {
		displs[i] = displs[i - 1] + sendcounts[i - 1];//the next process starts from where the previous started + the size of the previous.
	}
	return displs;//return displs array.
}

//Splits the array into chunks across processes.
void shareArray(int *sendArray, int sendHeight, int sendWidth, int **recvArray, int *recvHeight, int *recvWidth, int root, MPI_Comm comm)
{
	int myRank, n_process;
	int *sendcounts;
	int *displs;
	int smallSize;
	MPI_Comm_rank(comm, &myRank);//Curent procces.
	MPI_Comm_size(comm, &n_process);//Procces summary.

    //If i am on the root process.
	if (myRank == root) {
		sendcounts = findSendCounts(sendHeight, sendWidth, n_process);// Calculate how many elements each process will recieve.
		displs = findDispls(n_process, sendcounts);//Calculates starting positions for each process.
	}
	MPI_Scatter(&sendcounts[0], 1, MPI_INT, &smallSize, 1, MPI_INT, root, comm);//Sending each process it's own sendcounts so we can malloc the corresponding array.
	*recvArray = (int *) malloc(sizeof(int) * smallSize);//malloc array with smallSize positions (this is the size that the current process have.
	MPI_Scatterv(sendArray, sendcounts, displs, MPI_INT, *recvArray, smallSize, MPI_INT, root, comm);//Each process has it's own array.
	MPI_Bcast(&sendWidth, 1, MPI_INT, root, comm);//We broadcast on every process the columns (they are fixed).
	*recvWidth = sendWidth;//The columns of the new small array are the same as the old ones (they are fixed).
	*recvHeight = smallSize / sendWidth;//The lines of the small array that we have on each process are equal with the smallSize / the columns.

}

//Gives us flag if the chunk that we gave is (1) or is not (0) diagonally dominant.
int giveFlag(int *A, int grammes, int stiles, int pos)
{
	int sum;//Will get all the values ​​of each row and then subtract from the point that is diagonally dominant.
	int flag = 1;
	int i, j;
	//The loop continues only if i < the lines and if the flag == 1, because if it becomes - it means that is not diagonally dominant.
	for (i = 0; i < grammes && flag; i++) {
		sum = 0;
		for (j = 0; j < stiles; j++) {
			sum += abs(A[i * stiles + j]); //Takes the sum of the absolute value of the elements of the row i we are on.
		}
		sum -= abs(A[i * stiles + pos]); //The absolute value of the diagonal element is subtracted from it.

		//If the absolute value of the element of the diagonal is greater than the sum of the absolute values ​​of the remaining elements.
		//Then the flag takes the value 1 and continues.
		flag = flag && abs(A[i * stiles + pos]) > sum;
		pos++;//To find the diagonal element of the next line.
	}
	return flag;
}

// |Aii| > Σ |Aij|  i=/j for each line.
int isItDesp(int *pinakas, int grammes, int stiles, int sendHeight, int sendWidth, int root, MPI_Comm comm)
{
	//If each chunk is diagonally dominant then returns 1.
	int i, j, myRank, n_process;
	int flag, size;
	int *sendcounts;
	int *displs;
	MPI_Comm_rank(comm, &myRank);//Getting current procces.
	MPI_Comm_size(comm, &n_process);//Getting the procces sum.
	sendcounts = findSendCounts(sendHeight, sendWidth, n_process);
	displs = findDispls(n_process, sendcounts);
	int pos;
	pos = displs[myRank] / stiles; //pos = the position Aii we want to see. pos+1 = the next position.
	flag = giveFlag(pinakas, grammes, stiles, pos);
	return flag;
}

//Prints the array we give.
void print(void *Array, int height, int width)
{
	int i, j;
	int (*ptr)[width] = Array;
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			printf("%d ", ptr[i][j]);
		}
		printf("\n");
	}
}

//Returns the maximum element of the diagonal of the array we gave it.
int giveMax(int *A, int grammes, int stiles, int pos)
{
	int i, j, max; // A[i*stiles+pos] = element of the diagonal.
	max = A[pos]; //The first position of the diagonal of the piece we are on is initialized as maximum.
	for (i = 0; i < grammes; i++, pos++) {
		for (j = 0; j < stiles; j++) {
			if (max < A[i * stiles + pos]) //In case any diagonal is bigger then max, it takes its value.
				max = A[i * stiles + pos];
		}
	}
	return max;//returning the max.
}

// Calls giveMax and returns a flag.
int getMax(int *pinakas, int grammes, int stiles, int sendHeight, int sendWidth, int root, MPI_Comm comm)
{
	int i, j, myRank, n_process;
	int flag, size;
	int *sendcounts;
	int *displs;
	MPI_Comm_rank(comm, &myRank);//Current procces.
	MPI_Comm_size(comm, &n_process);//Procces sum.
	sendcounts = findSendCounts(sendHeight, sendWidth, n_process);
	displs = findDispls(n_process, sendcounts);
	int pos;
	pos = displs[myRank] / stiles; //pos = the Aii position we want to see . pos+1 = the next position.
	flag = giveMax(pinakas, grammes, stiles, pos);
	return flag;
}

// Finds the chunk of B array.
int *partOfBFind(int *pinakas, int grammes, int stiles, int pos, int max)
{
	int i, j;
	int *returnArray;//Array that at the end contains chunk of the array B found in each process.
	returnArray = makeArray(grammes * stiles); //Memory allocation for the array that will contain the fragment.
	for (i = 0; i < grammes; i++, pos++) { // A[i*stiles+pos] = the absolute value of the diagonal element.
		for (j = 0; j < stiles; j++) {
			returnArray[i * stiles + j] = max - abs(pinakas[i * stiles + j]); //Bij=m-|Aij|
		}
		returnArray[i * stiles + pos] = max; //Bij=max when we are in the diagonal element, we put the max
	}
	return returnArray;
}

//Returns a part of B that has been created by the partOfBFind.
int *givePartOfB(int *pinakas, int grammes, int stiles, int sendHeight, int sendWidth, int max, MPI_Comm comm)
{
	int i, j, myRank, n_process;
	int *sendcounts;
	int *displs;
	int *newPinakas;
	MPI_Comm_rank(comm, &myRank);//Current proccess.
	MPI_Comm_size(comm, &n_process);//Process summary.
	sendcounts = findSendCounts(sendHeight, sendWidth, n_process);
	displs = findDispls(n_process, sendcounts);
	int pos;//The position of the diagonal element.
	pos = displs[myRank] / stiles;
	newPinakas = partOfBFind(pinakas, grammes, stiles, pos, max);
	return newPinakas;
}

//Eeturns each chunk of array B found in our processes and merges it into one in process 0 as requested in the exercise.
int *giveBackB(int *C, int h, int w, MPI_Comm comm)
{
	int *rcvh;//Buffer that will contain the size of lines each process will have.
	int *rcvw, i; //Buffer that will contain the size of the columns each process will have.
	int myRank, n_process;
	int *sendcounts;
	int *displs;
	int *recvB = 0;//Final completed array.
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);//Current process.
	MPI_Comm_size(MPI_COMM_WORLD, &n_process);//Process summary.
	if (!myRank) {
		rcvh = makeArray(n_process);
		rcvw = makeArray(n_process);
	}
	MPI_Gather(&h, 1, MPI_INT, rcvh, 1, MPI_INT, 0, MPI_COMM_WORLD); //Send the h of each process to process 0 where it stores it in the table rcvh h=rows of the table the process has.
	MPI_Gather(&w, 1, MPI_INT, rcvw, 1, MPI_INT, 0, MPI_COMM_WORLD); //Send the w of each process to process 0 where it stores it in the table rcvw w=columns of the table the process has.
	if (!myRank) {
		recvB = (int *) malloc(sizeof(int) * 4 * 4);// The large array that 0 will have when it receives all the small ones from the processes.
		sendcounts = (int *)malloc(sizeof(int) * n_process);// Shows how many elements each process will send to 0.
		for (i = 0; i < n_process; i++) {
			sendcounts[i] = rcvh[i] * rcvw[i];//The i process will send its own h*w size of data.
			displs = findDispls(n_process, sendcounts);// From where it will start to receive and store the data from the processes in the big table.
		}
	}
	MPI_Gatherv(C, h * w, MPI_INT, recvB, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD); //We send the chunks to the recvB array.
	return recvB;
}

// Give an array and its size and return the position where the smallest element is.
int giveMinLoc(int *A, int size)
{
	int i, j, minloc;
	minloc = 0;
	for (i = 0; i < size; i++) {
		if (A[minloc] > A[i])
			minloc = i;
	}
	return minloc;//return of the position where the smallest element is located.
}

//Returns array of 3 positions. position[0]=min position[1]=row where min is located position[2] column where min is located.
int *giveMinLocParallila(int *smallpinaka, int grammes, int stiles, int sendHeight, int sendWidth, int root, MPI_Comm comm)
{
	int myRank, n_process;
	int send[2], recv[2];
	int *final;
	int *sendcounts;
	int *displs;
	int disp, minloc;
	int A[2];
	MPI_Comm_rank(comm, &myRank);//Curent process.
	MPI_Comm_size(comm, &n_process);//Process summary.
	if (myRank == root) {
		sendcounts = findSendCounts(sendHeight, sendWidth, n_process); //How many items we will send.
		displs = findDispls(n_process, sendcounts); //From where we start sending them.
	}
	MPI_Scatter(&displs[0], 1, MPI_INT, &disp, 1, MPI_INT, root, comm);
	send[1] = giveMinLoc(smallpinaka, grammes * stiles); //Has the position of min in the small array.
	send[0] = smallpinaka[send[1]];//Has the value of min in the small array.
	//Now it has the position of min in the big table due to the fact that we added the displ of each process where it contains from where we start sending.
	send[1] += disp;
	/*  We have the array send which contains [0] = the min element of the process
	[1] the position of min in the big array the MPI_Reduce I write below compares send[0] from all processes
	and stores in the recv table the recv[0] value and the recv[1] position of min in the big table.
	*/
	MPI_Reduce(&send, &recv, 1, MPI_2INT, MPI_MINLOC, root, comm); // MPI_2INT means compare between two integers.
	if (myRank == root) {
		final = makeArray(3);// Create an array of 3 positions.
		final[0] = recv[0];// Ctore the value of min.
		final[1] = recv[1] / sendHeight; //Lines of min on big array.
		final[2] = recv[1] % sendHeight; // Columns of min to big.
	}
	return final;
}

//Prints if it is diagonally dominant and the final maximum value.
void printifitis(int final_flag, int finalmax)
{
	int isit;
	if (final_flag) {
		printf("DIAGONALLY DOMINANT\n");
		fflush(stdout);
		printf("FINAL MAX %d\n", finalmax);
		fflush(stdout);
		isit = 1;
	}
	else {
		isit = 0;
		printf("NOT DIAGONALLY DOMINANT\n");
		fflush(stdout);
	}
}

int main(int argc, char **argv)
{
    //4X4 Array
	int A[4][4] = {{7, 1, 2, 3}, {0, 5, 2, 2}, {2, 2, 8, 3}, {6, 7, 1, 15}};
	int final_flag; // This tells if the array is strictly diagonally dominant.
	int final_max;  //Max value from A array.
	int myRank;     //Current rank.
	int flag;       //flag for each process.
	int n_process;  //Sum of processes.
	int i;          //Counter.
	int *recvB = 0; //Getting back the B array as a single pointer of integers.
	MPI_Init(&argc, &argv); // MPI initialization.
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank); //Current process.
	MPI_Comm_size(MPI_COMM_WORLD, &n_process);  //Process summary (number of total processes).
	int *B, h, w, max;// B array, h (Height) w (Width).
	int *C; // Support array.
	shareArray(A, 4, 4, &B, &h, &w, 0, MPI_COMM_WORLD); // Sharing array on the processes.
	flag = isItDesp(B, h, w, 4, 4, 0, MPI_COMM_WORLD); // Is it diagonally dominant? (each process).
	max = getMax(B, h, w, 4, 4, 0, MPI_COMM_WORLD); // Finds the maximum of the process.
	MPI_Reduce(&flag, &final_flag, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD); //Logical and of all the flags, so 1 if it is a diagonal of the dominants 0 if it is not, it sends it to 0.
	MPI_Reduce(&max, &final_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); // Takes the max of each process and returns the largest to process 0.

	//If we are in the main process (ie at 0).
	//We have final_flag and final_max i.e. if they are dominant diagonals and the absolute value of the maximum element of the diagonal.
	if (!myRank) {
		 printifitis(final_flag, final_max);
	}
	MPI_Bcast(&final_flag, 1, MPI_INT, 0, MPI_COMM_WORLD); // We send whether or not it is diagonally dominant to all processes.
	//if final flag = 1.
	if (final_flag) {
		MPI_Bcast(&final_max, 1, MPI_INT, 0, MPI_COMM_WORLD); //We send final max to all processes to use in their functions.
		C = givePartOfB(B, h, w, 4, 4, final_max, MPI_COMM_WORLD); //At this point we have part of table B that we ask for.
		recvB = giveBackB(C, h, w, MPI_COMM_WORLD); //We have obtained B.
	}
	if (!myRank) {
		if (final_flag) {
			printf("THIS IS FULL MATRIX B\n");
			fflush(stdout);
			print(recvB, 4, 4);
			fflush(stdout);
		}
	}
	int *pinakas;
	if (final_flag);
	pinakas = giveMinLocParallila(C, h, w, 4, 4, 0, MPI_COMM_WORLD); //Returns an array that contains the min, the position of the min in the rows, the position of the min in the columns.
	if (!myRank) {
		if (final_flag) {
			printf("MIN=%d\nMINLOC_i=%d\nMINLOC_j=%d\n", pinakas[0], pinakas[1], pinakas[2]);
			fflush(stdout);
		}
	}
	MPI_Finalize();
	return 0;
}
