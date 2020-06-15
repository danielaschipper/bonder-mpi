#include "stdafx.h"
#include "analize.h"
#include "readwfn.h"
#include "iface.h"
#include <string>
#include <string.h>
#include <iostream>
#include <thread>
#include <fstream>
#include "minion.h"
#include <mpi.h>
#include <vector>
#include <unistd.h>
#include <argp.h>


const int SIZE = 600;

struct arguments
{
	char* type;
	char *inputFile="input.wfn";
	double res=0.02;     
	double cutoff=0.3;
	char *output="output"; 
	int cubesize=1;
	double x1=0,x2=1,y1=0,y2=1,z1=0,z2=1;
	int atom1=0,atom2=1,atom3=2,atom4=3;
	//char* dir=".";
	char* configFile="config";
	double dist = 1;
};


static struct argp_option options[] =
{
	{"input",'i',"FILENAME",0,"The input wavefunction file"},
	{"resolution",'r',"VOXELSIZE",0,"The length in amstrongs of a single voxel in the intergration grid"},
	{"cutoff",'c',"RDGCUTTOFF",0,"The maximium value for RDG"},
	{"output",'o',"OUTPUTPREFIX",0,"The file prefix for all output files"},
	{"cubesize",'q',"CUBESCALE",0,"Reduces the resolution for output cubefiles by the given factor, 0 will prvent cube files from being written"},
	{"x1",1,"XSTART",0,"Either the x cordanate to start for point mode, the x coordanate of the sphere or the x cordante for one corner in cube mode"},
	{"y1",2,"YSTART",0,"Either the y cordanate to start for point mode, the y coordanate of the sphere or the y cordante for one corner in cube mode"},
	{"z1",3,"ZSTART",0,"Either the z cordanate to start for point mode, the z coordanate of the sphere or the z cordante for one corner in cube mode"},
	{"x2",4,"XEND",0,"The x cordanate for the second corner in cube mode"},
	{"y2",4,"YEND",0,"The y cordanate for the second corner in cube mode"},
	{"z2",4,"ZEND",0,"The z cordanate for the second corner in cube mode"},
	{"atom1",'1',"ATOMNUMBER",0,"The first atom for use in line, trinagle or quad mode"},
	{"atom2",'2',"ATOMNUMBER",0,"The first atom for use in line, trinagle or quad mode"},
	{"atom3",'3',"ATOMNUMBER",0,"The first atom for use in  trinagle or quad mode"},
	{"atom4",'4',"ATOMNUMBER",0,"The first atom for use in quad mode"},
	//{"directory",'d',"FOLDER",0,"The folder to put output files in"},
	{"config",'f',"CONFIGFILE",0,"A file containing all of the options to run bonder, only in input mode"},
	{"radius",'l',"RADIUS",0,"In sphere mode, defines the radius of the region to look at"},
	{0}
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	arguments *argument = (arguments*)(state->input);

	switch (key)
	{
		case 'i':
			argument->inputFile = arg;
			break;
		case 'r':
			argument->res = std::stod(arg);
			break;
		case 'c':
			argument->cutoff = std::stod(arg);
			break;
		case 'l':
			argument->dist = std::stod(arg);
			break;
		case 'o':
			argument->output = arg;
			break;
		case 'q':
			argument->cubesize = std::stoi(arg);
                        break;
		case 1:
			argument->x1 = std::stod(arg);
                        break;
		case 2:
                        argument->y1 = std::stod(arg);
                        break;
		case 3:
                        argument->z1 = std::stod(arg);
                        break;
		case 4:
                        argument->x2 = std::stod(arg);
                        break;
		case 5:
                        argument->y2 = std::stod(arg);
                        break;
		case 6:
                        argument->y1 = std::stod(arg);
                        break;
		case '1':
                        argument->atom1 = std::stoi(arg);
                        break;
		case '2':
                        argument->atom2 = std::stoi(arg);
                        break;
		case '3':
                        argument->atom3 = std::stoi(arg);
                        break;
		case '4':
                        argument->atom4 = std::stoi(arg);
                        break;
		//case 'd':
                  //      argument->dir = arg;
                    //    break;
		case 'f':
			argument->configFile=arg;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 1)
			{
				argp_usage(state);
			}
			else
			{
				argument->type = arg;
			}
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 1)
			{
				argp_usage (state);
			}
			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static char args_doc[] = "type";
static char doc[] = "Bonder, a program designed to identify and map non covalent interactions";
static struct argp argp = {options, parse_opt, args_doc, doc};



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

//yes it is exacly the same, all work is done by the minions, different signiture is required to ensure that file parsing is done correctly
void runAllCentre(double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube,double x, double y, double z,double dist)
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
	struct arguments arguments;
	argp_parse (&argp, argc, argv, 0, 0, &arguments);


	if (arguments.type[0] == 'f')
	{
		useInputFile(arguments.configFile);
		return 0;
	}

	wfnData *inputFile = 0;
		try
		{
			inputFile = init(arguments.inputFile);

		}
		catch (const std::invalid_argument& ia) 
		{
			std::cout << "error in parssing wavefunction data, if you have more than 100 atoms 'bonder fixwfn' must be run" << std::endl;
			return 1;
		}

	std::cout << "data read" << std::endl;
	//letter file x y z res cutoff
	if (arguments.type[0] == 'p')
	{
		bool sucsess;
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		analize.setUpAnalysisBatch( arguments.x1 , arguments.y1, arguments.z1, arguments.res,batch);

		analize.anilizePoint();

		return 0;
	}

	//letter file 1 2 res cutoff
	if (arguments.type[0] == 'l')
	{
		drawline(arguments.atom1, arguments.atom2, arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		return 0;

	}

	//letter file 1 2 res cutoff
	if (arguments.type[0] == 't')
	{
				drawtrig(arguments.atom1, arguments.atom2,arguments.atom3, arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		return 0;

	}

	if (arguments.type[0] == 'q')
	{
		drawquad(arguments.atom1, arguments.atom2,arguments.atom3,arguments.atom4, arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		return 0;
	}


	//letter file res cutoff output
	if (arguments.type[0] == 'a')
	{
		runAll(arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		return 0;
	}

	//letter file minx miny minz maxx maxy maxz res outputFile
	if (arguments.type[0] == 'g')
	{
		analysisBatch* batch = new analysisBatch(*inputFile);
		outputCube(arguments.x1, arguments.y1, arguments.z1, arguments.x2, arguments.y2, arguments.z2, arguments.res, arguments.output, *inputFile, 1.0, batch, arguments.cubesize);
		printf("done");
		return 0;
	}

	//bool sucsess;
	//anilizePoint(0, 0, 0, 0, 2000, 2000, 2501, &sucsess);
	//printf("%d", sucsess);
//letter file res cutoff output centrex centery centerz dist
	if (argv[1][0] == 's')
	{

		runAllCentre(arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile,arguments.x1, arguments.y1, arguments.z1,arguments.cubesize,arguments.dist);
		return 0;
	}


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
		minion(argc,argv);
}

#ifdef preargs
int main()
{
	char* argv[] = {"","a","input.xyz","0.02","0.3","output"};
	int argc = 6;
	main2(argc,argv);
}
#endif

