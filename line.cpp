#include <queue>
#include "iface.h"
#include <stdio.h>
#include <mpi.h>

double* drawline(int a, int b, double res, double cutoff,  wfnData* inputFile)
{
    
        double lowX = (*inputFile).x[a];
        double lowY = (*inputFile).y[a];
        double lowZ = (*inputFile).z[a];

        double highX = (*inputFile).x[b];
        double highY = (*inputFile).y[b];
        double highZ = (*inputFile).z[b];

        if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
        {
                return 0;
        }

        printf("testing line between %d and %d\n",a,b);

        analysisBatch* batch = new analysisBatch(*inputFile);
        double jumpScaler = res * 5/ ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

        double dx = (highX - lowX) * jumpScaler;
        double dy = (highY - lowY) * jumpScaler;
        double dz = (highZ - lowZ) * jumpScaler;
        int reps = 1 / jumpScaler;
        batch->cullLine(lowX + reps * dx/2, lowY + dy * reps/2, lowZ + dz * reps/2);

        int flips = 1;
        for (int i = 0; i < reps; i++)
        {
                int k = reps / 2 + flips * i / 2;
                flips *= -1; 
                double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
                if (mesured <= cutoff)
                {

                        double* vals = new double[3];
                        vals[0] = lowX + k*dx;
			vals[1] = lowY + k*dy;
			vals[2] = lowZ + k*dz;
			printf("found point\n");
   			delete batch; 
			return vals;
                        
                }
        }


        delete batch;
	return 0;
}


void drawLines(int a,double res, double cutoff, wfnData* inputFile)
{
	std::queue<double*> *points = new std::queue<double*>;
	for(int i = 0; i < 10; i++)
	{
		if (a+i >= inputFile->nuc)
			break;
		for (int b= 0; b < a+i; b++)
		{
			double* center = drawline(a+i,b,res,cutoff,inputFile);
			if (center)
			{
				points->push(center);
			}
		}
	}
	char signal = 'c';
	MPI_Send(&signal, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	int numOfCent = points->size();
	MPI_Send(&numOfCent, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	double* toSend = new double[numOfCent*3];
	for (int i = 0; i < numOfCent; i++)
	{
		double* cent = points->front();
		toSend[3*i] = cent[0];
		toSend[3*i+1] = cent[1];
		toSend[3*i+2] = cent[2];
		points->pop();
		delete cent;
	}
	MPI_Send(toSend, 3*numOfCent, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
	delete toSend;
	delete points;



}
