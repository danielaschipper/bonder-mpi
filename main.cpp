#include "stdafx.h"
#include "analize.h"
#include "readwfn.h"
#include "iface.h"
#include <string>
#include <string.h>
#include <iostream>
#include <thread>
#include <fstream>
#include "slave.h"
#include <mpi.h>
#include <vector>
#include <unistd.h>


const int SIZE = 600;

wfnData* init(std::string file)
{
	wfnData *inputFile = readFile(file);
	return inputFile;
}

void drawline(int a, int b, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{
	
	double lowX = (*inputFile).x[a];
	double lowY = (*inputFile).y[a];
	double lowZ = (*inputFile).z[a];

	double highX = (*inputFile).x[b];
	double highY = (*inputFile).y[b];
	double highZ = (*inputFile).z[b];

	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}

	printf("testing line between %d and %d\n",a,b);

	analysisBatch* batch = new analysisBatch(*inputFile);
	double jumpScaler = res * 5/ ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	batch->cullLine(lowX + reps * dx/2, lowY + dy * reps/2, lowZ + dz * reps/2);

	//printf("%d\n",reps);
	int flips = 1;
	for (int i = 0; i < reps; i++)
	{
		int k = reps / 2 + flips * i / 2;
		flips *= -1;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured <= cutoff)
		{

			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);

			analize.anilizePoint();
			
			break;
		}
	}


	delete batch;
}

void drawtrig(int a, int b,int c, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{

	
	analysisBatch* batch = new analysisBatch(*inputFile);
	double highX = ((*batch).atomx(a) + (*batch).atomx(b))/2;
	double highY = ((*batch).atomy(a) + (*batch).atomy(b))/2;
	double highZ = ((*batch).atomz(a) + (*batch).atomz(b))/2;

	double lowX = (*batch).atomx(c);
	double lowY = (*batch).atomy(c);
	double lowZ = (*batch).atomz(c);

	printf("testing line between centre of %d and %d and %d\n",a,b,c);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;

	printf("%d\n",reps);
	for (int i = 0; i < reps; i++)
	{
		int k = i;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured <= cutoff)
		{

			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);

			analize.anilizePoint();
			break;
		}
	}


	delete batch;
}

void drawquad(int a, int b,int c,int d, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{

	analysisBatch* batch = new analysisBatch(*inputFile);
	double highX = ((*batch).atomx(a) + (*batch).atomx(b))/2;
	double highY = ((*batch).atomy(a) + (*batch).atomy(b))/2;
	double highZ = ((*batch).atomz(a) + (*batch).atomz(b))/2;

	double lowX = ((*batch).atomx(c)+(*batch).atomx(d))/2;
	double lowY = ((*batch).atomy(c)+(*batch).atomx(d))/2;
	double lowZ = ((*batch).atomz(c)+(*batch).atomx(d))/2;

	printf("testing line between centre of %d and %d and %d\n",a,b,c);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;

	printf("%d\n",reps);
	for (int i = 0; i < reps; i++)
	{
		int k = i;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured <= cutoff)
		{

			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);

			analize.anilizePoint();
			break;
		}
	}


	delete batch;
}
struct pdrawArgs
{
	int a; int b; double res; double cutoff; std::string outputfile; int size; wfnData* inputFile; int makeCube;
	pdrawArgs(int A, int B, double Res, double Cutoff, std::string Outputfile, int Size, wfnData* InputFile,int MakeCube)
	{
		makeCube = MakeCube;
		a = A;
		b = B;
		res = Res;
		cutoff = Cutoff;
		outputfile = Outputfile;
		size = Size;
		inputFile = InputFile;
	}
};

void pDrawline(void *input)
{
	pdrawArgs* data = (pdrawArgs*)input;
	drawline((*data).a, (*data).b, (*data).res, (*data).cutoff, (*data).outputfile, (*data).size, (*data).inputFile, (*data).makeCube);
	//pthread_exit(NULL);
	delete data;
}

void runAll(double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{

	int world_size;
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	//set up multithreading
	initRanks(world_size);
	pthread_t* threads = new pthread_t[2];
	pthread_create(&threads[0], NULL, moniter, 0);
	createLine(inputFile->nuc);

	assignJobs(0);
	std::cout << "done" << std::endl;


}

std::vector<std::string> readFileLines(const char* filename)
{
	std::vector<std::string> file;
	std::ifstream input(filename);
	std::string line;
	while (getline(input, line)){
		file.push_back(line);
	}
	return file;
}
void useInputFile(char* filename)
{

	std::fstream inputFileTest(filename);
	if(!inputFileTest)
	{
		std::cout << "input file not found" << std::endl;
		return;
	}
	std::vector<std::string> lines;
	lines = readFileLines(filename);
	int lineNum = lines.size();
	if(lineNum == 0)
	{
		std::cout << "the input file needs text" <<  std::endl;
		printf("bonder h for help\n");
		return;
	}

	if(lineNum == 1)
	{
		std::cout << "please select option in the input file and the name of the wfn file";
		return;
	}

	wfnData *inputFile = 0;
	try
	{
		inputFile = init(lines[1]);

	}
	catch (const std::invalid_argument& ia) 
	{
		std::cout << "error in parssing wavefunction data, if you have more than 100 atoms 'bonder fixwfn' must be run" << std::endl;
		return;
	}

	//letter file res cutoff output
	if (lines[0] == "a")
	{
		if (lineNum != 6)
		{
			std::cout << "error in parsing input file\nall bonds file format is:\na\nwfn file\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}


		try
		{
			runAll(std::stod(lines[2]), std::stod(lines[3]), lines[4], SIZE, inputFile, std::stoi(lines[5]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;
	}


}


int master(int argc, char *argv[])
{
	if (argc == 1)
	{
		printf("bonder h for help\n");
		return 0;
	}


	if (argv[1][0] == 'h')
	{
		printf("the first letter determins the wht the program will do \n p looks at a point and determins the volume around it\n l looks at a line between two atoms\n a looks for all interactions\n g prints out a grid\n h displays this message\n for more detail on an operation type bonder letter");
		return 0;
	}

	if (argv[1][0] == 'f')
	{
		useInputFile(argv[2]);
		return 0;
	}

	wfnData *inputFile = 0;
	if (argc != 2)
	{
		try
		{
			inputFile = init(argv[2]);

		}
		catch (const std::invalid_argument& ia) 
		{
			std::cout << "error in parssing wavefunction data, if you have more than 100 atoms 'bonder fixwfn' must be run" << std::endl;
			return 1;
		}
	}

	std::cout << "data read" << std::endl;
	//letter file x y z res cutoff
	if (argv[1][0] == 'p')
	{
		if (argc != 9 && argc != 10 )
		{
			printf("arguments are bonder p inputFile x y z res cutoff outputFile\n");
			return 0;
		}
		bool sucsess = false;
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		try
		{
			analize.setUpAnalysisBatch( std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]),batch);
			printf("%f \n", (*batch).RDG(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5])));
			if (argc == 10)
			{
				analize.anilizePoint();
			}
			else
			{
				analize.anilizePoint();
			}
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}

		if (sucsess)
		{
			printf("point given is in region\n");
		}
		else
		{
			printf("point given is not in region\n");
		}
		return 0;
	}

	//letter file 1 2 res cutoff
	if (argv[1][0] == 'l')
	{
		if (!(argc == 8 || argc == 9))
		{
			printf("arguments are bonder l inputfile atom1 atom2 res cutoff outputfile\n");
			return 0;
		}
		try
		{
			if (argc == 8)
				drawline(std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), argv[7], SIZE, inputFile,1);
			else
				drawline(std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), argv[7], SIZE, inputFile, std::stoi(argv[8]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		return 0;

	}

	//letter file 1 2 res cutoff
	if (argv[1][0] == 't')
	{
		if (!(argc == 8 || argc == 9))
		{
			printf("arguments are bonder t inputfile atom1 atom2 atom3 res cutoff outputfile\n");
			return 0;
		}
		try
		{
			if (argc == 9)
				drawtrig(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), argv[8], SIZE, inputFile,1);
			else
				drawtrig(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), argv[8], SIZE, inputFile, std::stoi(argv[9]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		return 0;

	}
	
	if (argv[1][0] == 'q')
	{
		if (!(argc == 9 || argc == 10))
		{
			printf("arguments are bonder t inputfile atom1 atom2 atom3 res cutoff outputfile\n");
			return 0;
		}
		try
		{
			if (argc == 9)
				drawquad(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]),std::stoi(argv[6]), std::stod(argv[7]), std::stod(argv[8]), argv[1], SIZE, inputFile,1);
			else
				drawquad(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]),std::stoi(argv[6]), std::stod(argv[7]), std::stod(argv[8]), argv[1], SIZE, inputFile, std::stoi(argv[9]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		return 0;
	}
	

	//letter file res cutoff output
	if (argv[1][0] == 'a')
	{
		if (argc != 6 && argc != 7)
		{
			printf("arguments are bonder a inputFile res cutoff outputFile\n");
			return 0;
		}

		try
		{
			if (argc == 6)
				runAll(std::stod(argv[3]), std::stod(argv[4]), argv[5], SIZE, inputFile,true);
			else
				runAll(std::stod(argv[3]), std::stod(argv[4]), argv[5], SIZE, inputFile, std::stoi(argv[6]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		return 0;
	}

	//letter file minx miny minz maxx maxy maxz res outputFile
	if (argv[1][0] == 'g')
	{
		if (argc != 11 )
		{
			printf("arguments are bonder g inputFile minx miny minz maxx maxy maxz res outputFile\n");
			return 0;
		}
		analysisBatch* batch = new analysisBatch(*inputFile);
		try
		{
				outputCube(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), argv[10], *inputFile, 1.0, batch, 1);
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		printf("done");
		return 0;
	}

	//bool sucsess;
	//anilizePoint(0, 0, 0, 0, 2000, 2000, 2501, &sucsess);
	//printf("%d", sucsess);
	printf("bonder h for help\n");
	return 0;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	std::cout << "MPI rank: " << world_rank << " PID: " << getpid() << std::endl;
	if (!world_rank)
		master(argc,argv);
	else
		slave(argc,argv);
}

#ifdef preargs
int main()
{
	char* argv[] = {"","a","input.xyz","0.02","0.3","output"};
	int argc = 6;
	main2(argc,argv);
}
#endif

