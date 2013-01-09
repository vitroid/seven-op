#include "analysis.h"


//lattice coordinate to absolute address
int address(int x, int y, int z, int division)
{
  x = (x+division) % division;
  y = (y+division) % division;
  z = (z+division) % division;
  return (x*division+y)*division+z;
}


int latticex(int address, int division)
{
  return address / (division*division);
}


int latticey(int address, int division)
{
  return (address % (division*division)) / division;
}


int latticez(int address, int division)
{
  return address % division;
}


void output_coordinates(const Vector3d &box, const sWaters &waters)
{
  cout << "@BOX3" << endl;
  cout << box[0] << " "
       << box[1] << " "
       << box[2] << endl;
  cout << "@AR3A" << endl;
  cout << waters.size() << endl;
  for(int i=0; i<waters.size(); i++){
    cout << waters[i].com[0] << " "
	 << waters[i].com[1] << " "
	 << waters[i].com[2] << endl;
  }
}


void optimize_water_positions(sWaters& waters, const Vector3d& grid)
{
  float delta = 1.0;
  while( delta > 0.001 ){
    Vector3d origin(0.0,0.0,0.0);
    for(int i=0;i<waters.size(); i++){
      for(int dim=0;dim<3;dim++){
	//displacement from the lattice points
	//this might fail when displacement is as large as grid/4.
	//so we solve it iteratively.
	origin[dim] += waters[i].com[dim] - rint( waters[i].com[dim] / (grid[dim]/2) ) * (grid[dim]/2);
      }
    }
    delta = 0.0;
    for(int dim=0;dim<3;dim++){
      origin[dim] /= waters.size();
      delta += origin[dim]*origin[dim];
    }
    //cout << delta << " Drift" << endl;
    for(int i=0;i<waters.size(); i++){
      for(int dim=0;dim<3;dim++){
	waters[i].com[dim] -= origin[dim];
      }
    }
  }    
}


void address_on_the_lattice(const sWaters& waters, const Vector3d& grid, int division, int mycubiclattice[], int residentA[], int residentB[], int myaddress[])
{
  float dsum = 0.0;
  for(int i=0;i<waters.size(); i++){
    float dA=0,dB=0;
    for(int dim=0;dim<3;dim++){
      float d = waters[i].com[dim] - rint( waters[i].com[dim] / grid[dim] ) * grid[dim];
      dA += d*d;
    }
    for(int dim=0;dim<3;dim++){
      float d = waters[i].com[dim] - (rint( (waters[i].com[dim]-grid[dim]/2) / grid[dim] ) * grid[dim] + grid[dim]/2);
      dB += d*d;
    }
    dsum += min(dA,dB);
    if ( dA < dB ){
      mycubiclattice[i] = CUBICLATTICEA;
    }
    else{
      mycubiclattice[i] = CUBICLATTICEB;
    }
    if ( mycubiclattice[i] == CUBICLATTICEA ){
      int ix = rint( waters[i].com[0] / grid[0] );
      int iy = rint( waters[i].com[1] / grid[1] );
      int iz = rint( waters[i].com[2] / grid[2] );
      int ia = address(ix,iy,iz,division);
      if ( residentA[ia] >= 0 ){
	int j = residentA[ia];
	fprintf(stderr, "Two residents in a cell %d: %d,%d.\n", ia, i, residentA[ia]);
	fprintf(stderr, "%f %f %f\n", waters[i].com[0], waters[i].com[1], waters[i].com[2]);
	fprintf(stderr, "%f %f %f\n", waters[i].com[0], waters[i].com[1], waters[i].com[2]);
	fprintf(stderr, "%d %d %d\n", ix,iy,iz);
	exit(1);
      }
      residentA[ia] = i;
      myaddress[i] = ia;
    }
    else{
      int ix = rint( (waters[i].com[0]-grid[0]/2) / grid[0] );
      int iy = rint( (waters[i].com[1]-grid[1]/2) / grid[1] );
      int iz = rint( (waters[i].com[2]-grid[2]/2) / grid[2] );
      int ia = address(ix,iy,iz,division);
      if ( residentB[ia] >= 0 ){
	int j = residentB[ia];
	fprintf(stderr, "Two residents in a cell %d: %d,%d.\n", ia, i, residentB[ia]);
	exit(1);
      }
      residentB[ia] = i;
      myaddress[i] = ia;
    }
  }
  //cerr << dsum << " Total square displacements from the lattice points." << endl;
}




void hydrogen_bond(int residentA[], int residentB[], Vector3d& box,
		   const sWaters& waters, int division,
		   vector<int> hb[],
		   vector<int> dhb[])
{
  //lattice B
  for(int ix=0;ix<division;ix++){
    for(int iy=0;iy<division;iy++){
      for(int iz=0;iz<division;iz++){
	int Lb = address(ix,iy,iz,division);
	int i = residentB[Lb];
	if ( i >= 0 ){
	  for(int jx=0;jx<2;jx++){
	    for(int jy=0;jy<2;jy++){
	      for(int jz=0;jz<2;jz++){
		int kx = (ix+jx) % division;
		int ky = (iy+jy) % division;
		int kz = (iz+jz) % division;
		int La = address(kx,ky,kz,division);
		int j = residentA[La];
		if ( j >= 0 ){
		  Vector3d delta = Wrap(waters[i].com-waters[j].com,box);
		  //#molecule pair is close enough
		  //#find the shortest OH pair distance among four possible combinations
		  float di1 = ( delta + waters[i].rH1 ).norm();
		  float di2 = ( delta + waters[i].rH2 ).norm();
		  float dj1 = ( delta - waters[j].rH1 ).norm();
		  float dj2 = ( delta - waters[j].rH2 ).norm();
		  float di = min(di1,di2);
		  float dj = min(dj1,dj2);
		  float d  = min(di,dj);
		  if (d < 2.2) {
		    //#ª©; threshold for O-H adjacency
		    //#register as a hydrogen bond.
		    //#bond direction is ignored.
		    hb[i].push_back(j); //#hb[i] contains the labels of the HB partners
		    hb[j].push_back(i);
		    if ( di < dj ){
		      dhb[i].push_back(j);
		    }
		    else{
		      dhb[j].push_back(i);
		    }		      
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


void list_neighbors(int division, vector<int> neighborA[], vector<int> neighborB[])
{
  for(int ix=0;ix<division;ix++){
    for(int iy=0;iy<division;iy++){
      for(int iz=0;iz<division;iz++){
	//lattice A
	int ia = address(ix,iy,iz,division);
	int jx = (ix - 1 + division) % division;
	int jy = (iy - 1 + division) % division;
	int jz = (iz - 1 + division) % division;
	neighborA[ia].push_back(address(jx,jy,jz,division));
	neighborA[ia].push_back(address(ix,iy,jz,division));
	neighborA[ia].push_back(address(jx,iy,iz,division));
	neighborA[ia].push_back(address(ix,jy,iz,division));
	neighborA[ia].push_back(address(ix,jy,jz,division));
	neighborA[ia].push_back(address(jx,iy,jz,division));
	neighborA[ia].push_back(address(jx,jy,iz,division));
	neighborA[ia].push_back(address(ix,iy,iz,division));
	//lattice B
	jx = (ix + 1) % division;
	jy = (iy + 1) % division;
	jz = (iz + 1) % division;
	neighborB[ia].push_back(address(ix,jy,jz,division));
	neighborB[ia].push_back(address(jx,iy,jz,division));
	neighborB[ia].push_back(address(jx,jy,iz,division));
	neighborB[ia].push_back(address(ix,iy,iz,division));
	neighborB[ia].push_back(address(jx,jy,jz,division));
	neighborB[ia].push_back(address(ix,iy,jz,division));
	neighborB[ia].push_back(address(jx,iy,iz,division));
	neighborB[ia].push_back(address(ix,jy,iz,division));
      }
    }
  }
}



void estimate_optimized_network(const sWaters& waters,
				int mycubiclattice[],
				vector<int> neighborA[],
				vector<int> neighborB[],
				const vector<float>& op,
				int residentA[],
				int residentB[],
				int myaddress[],
				vector<int> assign[])
{
  // 0==undetermined
  // 1==regular bond (directing the correct neighbor)
  // 2==regular vacancy
  for(int i=0;i<waters.size(); i++){
    assign[i].resize(waters.size());
  }
  for(int iop=4; iop>0; iop--){
    float fop = iop / 4.0;
    for(int i=0;i<waters.size(); i++){
      if ( op[i] == fop ){
	//op is positive: bond dirs are (+--, --+, -+-, +++)
	//even if op is 4, it should not always obey the ice rule.
	if ( mycubiclattice[i] == CUBICLATTICEA ){
	  //i am on lattice A
	  int Li = myaddress[i];
	  for(int k=4;k<8;k++){
	    int Lj = neighborA[Li][k];
	    int j = residentB[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 1;//regular bond
		assign[j][i] = 1;//regular bond
	      }
	    }
	  }
	  for(int k=0;k<4;k++){
	    int Lj = neighborA[Li][k];
	    int j = residentB[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 2;//regular vacancy
		assign[j][i] = 2;//regular vacancy
	      }
	    }
	  }
	}
	else{
	  //i am on lattice B
	  int Li = myaddress[i];
	  for(int k=4;k<8;k++){
	    int Lj = neighborB[Li][k];
	    int j = residentA[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 1;//regular bond
		assign[j][i] = 1;//regular bond
	      }
	    }
	  }
	  for(int k=0;k<4;k++){
	    int Lj = neighborB[Li][k];
	    int j = residentA[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 2;//regular vacancy
		assign[j][i] = 2;//regular vacancy
	      }
	    }
	  }
	}
      }
      else if ( op[i] == -fop ){
	//op is negative
	//even if op is 4, it should not always obey the ice rule.
	if ( mycubiclattice[i] == CUBICLATTICEA ){
	  //i am on lattice A
	  int Li = myaddress[i];
	  for(int k=4;k<8;k++){
	    int Lj = neighborA[Li][k];
	    int j = residentB[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 2;//regular vacancy
		assign[j][i] = 2;//regular vacancy
	      }
	    }
	  }
	  for(int k=0;k<4;k++){
	    int Lj = neighborA[Li][k];
	    int j = residentB[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 1;//regular bond
		assign[j][i] = 1;//regular bond
	      }
	    }
	  }
	}
	else{
	  //i am on lattice B
	  int Li = myaddress[i];
	  for(int k=4;k<8;k++){
	    int Lj = neighborB[Li][k];
	    int j = residentA[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 2;//regular vacancy
		assign[j][i] = 2;//regular vacancy
	      }
	    }
	  }
	  for(int k=0;k<4;k++){
	    int Lj = neighborB[Li][k];
	    int j = residentA[Lj];
	    if ( j >= 0 ){
	      //if undetermined
	      if ( assign[i][j] == 0 ){
		assign[i][j] = 1;//regular bond
		assign[j][i] = 1;//regular bond
	      }
	    }
	  }
	}
      }
    }
  }
}
