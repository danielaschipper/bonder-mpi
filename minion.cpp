#include "readwfn.h"
#include "fill.h"
#include "output.h"
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include "line.h"
#include <argp.h>


double res,cutoff;
int makeCube;
std::string outputFile;
wfnData *inputFile = 0;
const int Xsize = 600, Ysize = 600;

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
	char* dir=".";
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
	{"directory",'d',"FOLDER",0,"The folder to put output files in"},
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
		case 'd':
                        argument->dir = arg;
                        break;
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

bool edgeComp(edgepoint i, edgepoint j)
{
        return i.Z < j.Z;
}

void minionFill(double *data)
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



void minionOuput(double* data)
{
	std::cout << "output" << std::endl;
	analysisBatch* batch = new analysisBatch(*inputFile);
	(*batch).setUpBatch(data[0], data[1], data[2], res);
	(*batch).cull(data[0], data[1], data[2]);
	outputCube(data[3], data[4], data[5],data[6], data[7], data[8], res, outputFile, *inputFile ,cutoff,batch,makeCube);
	char signal = 'd';
	MPI_Send(&signal, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	delete batch;
}

void minion(int argc, char *argv[])
{
	struct arguments arguments;
	argp_parse (&argp, argc, argv, 0, 0, &arguments);
	std::cout << "starting" << std::endl;
	//init
	try
	{
		inputFile = readFile(arguments.inputFile);

	}
	catch (const std::invalid_argument& ia)
	{
		return;
	}
	res = arguments.res;
	cutoff = arguments.cutoff;
	outputFile =  arguments.output;
	makeCube = arguments.cubesize;
	bool sphere = false; 
	double centerx,centery,centerz,dist;
	if (arguments.type[0] == 's')
	{
		sphere = true; 
		centerx=arguments.x1;centery=arguments.y1;centerz=arguments.z1;dist=arguments.dist;
		makeCube = 1;
	}


	std::cout << "waiting" << std::endl;
	for (;;)
	{
		int signal;
		MPI_Recv(&signal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//std::cout << "recived" << std::endl;
		if (signal ==-1)
		{
			MPI_Finalize();
			return;
		}
		double data[9];
		MPI_Recv(data, 9, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (signal == 0)
			minionFill(data);
		else if (signal == 1)
		{
			minionOuput(data);
		}
		else
		{
			if (!sphere)
				drawLines((int)(data[0]),res,cutoff,inputFile);
			else
				drawLinesCent((int)(data[0]),res,cutoff,inputFile,centerx,centery,centerz,dist);
		}

	}
}



