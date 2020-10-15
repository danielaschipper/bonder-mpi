#include "stdafx.h"
#include <algorithm>
#include <vector>
#include <queue>
#include "analize.h"
#include <mutex>
#include <list>
#include <iostream>
#include <mutex>
#include <mpi.h>
#include <chrono>
#include <thread>



std::mutex mutexCube,mutexCenters,mutexJob,mutexRank;


struct grid
{
        double minX,minY,minZ,maxX,maxY,maxZ;
};





struct center
{
	double x, y, z;
	center(double X, double Y, double Z)
	{
		x = X;
		y = Y;
		z = Z;
	}
	center() {}
	double disi(center cent)
	{
		return (x - cent.x) * (x - cent.x) + (y - cent.y) * (y - cent.y) + (z - cent.z) * (z - cent.z);
	}
	inline bool operator==(const center lhs);
};

inline bool center::operator==(const center lhs) { return (lhs.x == x) && (lhs.y == y) && (lhs.z == z); }

struct cube
{
	double maxX, maxY, maxZ, minX, minY, minZ;
	cube(double NX,double NY,double NZ, double XX, double XY, double XZ)
	{
		maxX = XX;
		maxY = XY;
		maxZ = XZ;
		minX = NX;
		minY = NY;
		minZ = NZ;
	}
	cube(){}
};

//type =0 for fill type =1 for output
struct job
{
	center cent;
	cube cub;
	int type;
};

std::vector<cube> *spacesUsed = new std::vector<cube>;
std::list<center> *centers = new std::list<center>;
std::queue<int> *openRanks = new std::queue<int>;
std::queue<job> *jobQueue = new std::queue<job>;
unsigned int totalRanks;


void initRanks(int rank)
{
	for(int i = 1; i < rank; i++)
		(*openRanks).push(i);
	
	totalRanks = (*openRanks).size();
}


//moniters for this rank being freed
void* moniter(void* args)
{
	int rank;
	for (;;)
	{
		MPI_Status status;
		char signal = 'l';
		MPI_Recv(&signal, 1, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		rank = status.MPI_SOURCE;
		//std::cout << "recived" << std::endl;
		if (signal == 'g')
		{
			job newjob;
			double incomingData[9];
			//get raw data
			MPI_Recv(incomingData, 9, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			newjob.cent = center(incomingData[0],incomingData[1],incomingData[2]);
			newjob.cub = cube(incomingData[3],incomingData[4],incomingData[5],incomingData[6],incomingData[7],incomingData[8]);
			newjob.type = 1;
			
			mutexCenters.lock();
			(*centers).remove(newjob.cent);
			mutexCenters.unlock();



			bool dump = false;
			mutexCube.lock();
			for (size_t i = 0; i < (*spacesUsed).size(); i++)
			{
				cube var = (*spacesUsed)[i];
				if (newjob.cent.x < var.maxX && newjob.cent.x > var.minX && newjob.cent.y< var.maxY && newjob.cent.y > var.minY && newjob.cent.z < var.maxZ && newjob.cent.z > var.minZ)
				{
					dump = true;
					break;
				}
			}

			if (dump)
			{

				mutexCube.unlock();
				mutexRank.lock();
				(*openRanks).push(rank);
				mutexRank.unlock();

				continue;
			}
			else
			{
				(*spacesUsed).push_back(newjob.cub);
				mutexCube.unlock();
			}
			mutexJob.lock();
			(*jobQueue).push(newjob);
			mutexJob.unlock();
			std::cout << "jobs: " << (*jobQueue).size() << " ranks: " << (*openRanks).size() << " max ranks: " << totalRanks << std::endl;

		}
		else if (signal == 'c')
		{
			int newCents = 0;
			MPI_Recv(&newCents, 1, MPI_INT, rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			double* incomingData = new double[3*newCents];

			MPI_Recv(incomingData, 3*newCents, MPI_DOUBLE, rank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i = 0; i < newCents; i++)
			{
				center thisPoint(incomingData[i*3], incomingData[i*3+1], incomingData[i*3+2]);
				job newJob;
				newJob.type = 0;
				newJob.cent = thisPoint;
				mutexJob.lock();
				(*jobQueue).push(newJob);
				mutexJob.unlock();
			}
			delete[] incomingData;
			std::cout << "jobs: " << (*jobQueue).size() << " ranks: " << (*openRanks).size() << " max ranks: " << totalRanks << std::endl;
		}
		else
		{
			std::cout << "jobs: " << (*jobQueue).size() << " ranks: " << (*openRanks).size() << " max ranks: " << totalRanks << std::endl;
		}

		mutexRank.lock();
		(*openRanks).push(rank);
		mutexRank.unlock();
	}
	return 0;
}

void launchJob(job toLaunch)
{
	mutexRank.lock();
	int rank = (*openRanks).front();
	(*openRanks).pop();
	mutexRank.unlock();
	double data[9];
	data[0] = toLaunch.cent.x;
	data[1] = toLaunch.cent.y;
	data[2] = toLaunch.cent.z;
	data[3] = toLaunch.cub.minX;
	data[4] = toLaunch.cub.minY;
	data[5] = toLaunch.cub.minZ;
	data[6] = toLaunch.cub.maxX;
	data[7] = toLaunch.cub.maxY;
	data[8] = toLaunch.cub.maxZ;
	MPI_Send(&(toLaunch.type), 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
	MPI_Send(data, 9, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);

}

void* assignJobs(void* nah)
{
	std::cout << "assigning" << std::endl;
	for(;;)
	{
		mutexRank.lock();
		if (((*openRanks).size() != 0) && (!(*jobQueue).empty()))
		{

			mutexRank.unlock();
			mutexJob.lock();
			job nextJob = (*jobQueue).front();
			(*jobQueue).pop();
			mutexJob.unlock();
			bool dump =false;
			bool loop;
			if (nextJob.type == 0)//if the next job is a fill job
			{

				mutexCube.lock();
				for (size_t i = 0; i < (*spacesUsed).size(); i++)
				{
					cube var = (*spacesUsed)[i];
					if (nextJob.cent.x < var.maxX && nextJob.cent.x > var.minX && nextJob.cent.y< var.maxY && nextJob.cent.y > var.minY && nextJob.cent.z < var.maxZ && nextJob.cent.z > var.minZ)
					{
						dump = true;
						break;
					}
				}

				mutexCube.unlock();
				if (dump)
					continue;


				mutexCenters.lock();
				bool found = false;
				for (std::list<center>::iterator iti = ((*centers).begin()); iti != ((*centers).end()); ++iti)
				{
					if ((*iti).disi(nextJob.cent) < 100)
					{
						found = true;
						break;
					}
				}
				if (!found)
				{
					(*centers).push_back(nextJob.cent);
					loop = false;
					mutexCenters.unlock();
				}
				else
				{
					loop = true;
					mutexCenters.unlock();
				}
				if (loop)
				{
					mutexJob.lock();
					(*jobQueue).push(nextJob);
					mutexJob.unlock();
					std::this_thread::sleep_for(std::chrono::milliseconds(10));

				}
				else
				{
					launchJob(nextJob);
				}


			}
			else //if next job is a cube or lineset
			{
				launchJob(nextJob);
			}

		}
		else if ((*jobQueue).empty()&&(*openRanks).size() != 0)
		{

			std::this_thread::sleep_for(std::chrono::milliseconds(2000));
			if ((*jobQueue).empty()&&(*openRanks).size() != 0)
			{
				std::cout << "culling rank"  << std::endl;
				int rank = (*openRanks).front();
				(*openRanks).pop();
				int termsig = -1;
				totalRanks -= 1;
				MPI_Send(&termsig, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);

				if(totalRanks == 0)
				{
					//MPI_Finalize();
                        		return 0;	
				}
			}
			mutexRank.unlock();
		}
		else
		{
			mutexRank.unlock();
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
	}

	return 0;
}

//create the line jobs
void createLine(int n)
{
	for (int i =0; i < (n + 1)/10 +1; i++)
	{
		center thisPoint(i * 10, 0, 0);
		job newJob;
		newJob.type = 2;
		newJob.cent = thisPoint;
		mutexJob.lock();
		(*jobQueue).push(newJob);
		mutexJob.unlock();

	}
}

void analysis::setUpAnalysisBatch(double x, double y, double z, double resalution, analysisBatch* batch)
{
	//x = res * ((int)(x / res));
	//y = res * ((int)(y / res));
	//z = res * ((int)(z / res));
	(*batch).setUpBatch(x, y, z, resalution);
	offsetx = x;
	offsety = y;
	offsetz = z;
	res = resalution;
}

void analysis::anilizePoint()
{
	//std::cout << makeCube << std::endl;
	center thisPoint(offsetx, offsety, offsetz);
	job newJob;
	newJob.type = 0;
	newJob.cent = thisPoint;
	mutexJob.lock();
	(*jobQueue).push(newJob);
	mutexJob.unlock();

}

analysis::analysis()
{
}

analysis::~analysis()
{
}

