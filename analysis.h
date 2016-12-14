#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "seven-op.h"


#define CUBICLATTICEA 1
#define CUBICLATTICEB 0

//lattice coordinate to absolute address
int address(int x, int y, int z, int division);
int latticex(int address, int division);
int latticey(int address, int division);
int latticez(int address, int division);
void output_coordinates( const Vector3d &box, const sWaters &waters);
void optimize_water_positions(sWaters& waters, const Vector3d& grid);
void address_on_the_lattice(const sWaters& waters, const Vector3d& grid, int division,
			    int mycubiclattice[], int residentA[], int residentB[], int myaddress[]);
void hydrogen_bond(int residentA[], int residentB[], Vector3d& box,
		   const sWaters& waters, int division, vector<int> hb[],vector<int> dhb[]);
void list_neighbors(int division, vector<int> neighborA[], vector<int> neighborB[]);
void estimate_optimized_network(const sWaters& waters,
				int mycubiclattice[],
				vector<int> neighborA[],
				vector<int> neighborB[],
				const vector<float>& op,
				int residentA[],
				int residentB[],
				int myaddress[],
				vector<int> assign[]);


#endif /*ANALYSIS_H*/
