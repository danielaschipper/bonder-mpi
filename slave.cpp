#include "readwfn.h"
#include "fill.h"
#include "output.h"
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include "line.h"

double res,cutoff;
int makeCube;
std::string outputFile;
wfnData *inputFile = 0;
const int Xsize = 600, Ysize = 600;

bool edgeComp(edgepoint i, edgepoint j)
{
        return i.Z < j.Z;
}

void slaveFill(double *data)
{
	bool sucsess;
	analysisBatch* batch = new analysisBatch(*inputFile);
	(*batch).setUpBatch(data[0], data[1], data[2], res);
	(*batch).cull(data[0], data[1], data[2]);
	grid results = fill(0,0,0,0,Xsize,Ysize,cutoff,&sucsess,batch);

	int maxX=0, maxY=0, maxZ=0, minX=0, minY=0, minZ=0;
        int vol = 0;
        gridPoint currentPoint;
        for (int i = -Xsize/2; i < Xsize/2; i++)
        {
                for (int j = -Ysize/2; j < Ysize/2; j++)
                {
                        currentPoint = *getPoint(&results, i, j); 
                        int size;
                        edgepoint* edges = (*(currentPoint.edges)).dump(&size);
                        if (size == 0)
                                continue;
                        if (i > maxX)
                                maxX = i;
                        if (i < minX)
                                minX = i;

                        if (j > maxY)
                                maxY = j;
                        if (j < minY)
                                minY = j;
    
                        std::sort(edges, edges + size, edgeComp);
                        int lastZ = INT32_MAX;
                        for (int k = 0; k < size; k++)
                        {
                                //printf("%d %d %d\n", i, j, edges[k].Z);
                                if (edges[k].Z < minZ)
                                        minZ = edges[k].Z;
                                if (edges[k].Z > maxZ)
                                        maxZ = edges[k].Z;

                                switch (edges[k].LR)
                                {
                                case 0:
                                        vol++;
                                        break;
                                case 1:
                                        lastZ = edges[k].Z;
                                        break;
                                case 2:
                                        vol += edges[k].Z - lastZ + 1;
                                        lastZ = INT32_MAX;
                                }
                        }
                        delete [] edges;

                }
        }
	delete batch;
        minX -= 10;
        minY -= 10;
        minZ -= 10;

        maxX += 10;
        maxY += 10;
        maxZ += 10;
	char signal = 'g';
	MPI_Send(&signal, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	data[3] = minX * res + data[0];
	data[4] = minY * res + data[1];
	data[5] = minZ * res + data[2];
	data[6] = maxX * res + data[0];
	data[7] = maxY * res + data[1];
	data[8] = maxZ * res + data[2];
	MPI_Send(data, 9, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
}



void slaveOuput(double* data)
{
	analysisBatch* batch = new analysisBatch(*inputFile);
	(*batch).setUpBatch(data[0], data[1], data[2], res);
	(*batch).cull(data[0], data[1], data[2]);
	outputCube(data[3], data[4], data[5],data[6], data[7], data[8], res, outputFile, *inputFile ,cutoff,batch,makeCube);
	char signal = 'd';
	MPI_Send(&signal, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	delete batch;
}

void slave(int argc, char *argv[])
{
	std::cout << "starting" << std::endl;
	//init
	try
	{
		inputFile = readFile(argv[2]);

	}
	catch (const std::invalid_argument& ia)
	{
		return;
	}
	res = std::stod(argv[3]);
	cutoff = std::stod(argv[4]);
	outputFile =  argv[5];
	if (argc == 6)
		makeCube = 1;
	else
		makeCube = std::stoi(argv[6]);


	std::cout << "waiting" << std::endl;
	for (;;)
	{
		int signal;
		MPI_Recv(&signal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//std::cout << "recived" << std::endl;
		if (signal ==3)
		{
			MPI_Finalize();
			return;
		}
		double data[9];
		MPI_Recv(data, 9, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (signal == 0)
			slaveFill(data);
		else if (signal == 1)
			slaveOuput(data);
		else
			drawLines((int)(data[0]),res,cutoff,inputFile);

	}
}



