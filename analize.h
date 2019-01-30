#pragma once
#include "readwfn.h"
#include "iface.h"


//this class takes an input point and anayses it
class  analysis
{
public:
	 analysis();
	~ analysis();
	void anilizePoint();
	void setUpAnalysisBatch(double x, double y, double z, double resalution, analysisBatch* batch);
private:
	double offsetx, offsety, offsetz, res;
};

void outputCube(double minx, double miny, double minz, double maxx, double maxy, double maxz, double res, std::string file, wfnData inputData,double cutoff,analysisBatch* batch,int makeCube);
void* assignJobs(void* nah);
void initRanks(int rank);
void* moniter(void* Vrank);
void createLine(int n);
