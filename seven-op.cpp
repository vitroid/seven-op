
/*#read water configurations in Gromacs format,
#calculate the pairwise products of local order parameter for ice 7,
#and output them with pair distances.

#2012-10-19 tested with ice VII crystal structure.
#2012-12-27 Renamed (formerly seven-op-gro3.cpp)
*/


#include "seven-op.h"


int oneframe(int mode = 0, int debug=0){
  //#crystal axes
  Vector3d ex(1.0,0.0,0.0);
  Vector3d ey(0.0,1.0,0.0);
  Vector3d ez(0.0,0.0,1.0);

  //#stage 1: load the structure
  pair<sWaters, Vector3d> pWB;
  pWB = Configure_gro(stdin);
  sWaters waters = pWB.first;
  Vector3d box  = pWB.second;
  
  if ( mode == 4 ){
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
    exit(0);
  }
  int lattice = latticeSize(waters.size()); //assume subic box
  Vector3d grid;
  for(int dim=0;dim<3;dim++){
    grid[dim] = box[dim] / lattice;
  }
  if (waters.size() ==0){
    return 0;
  }
  if (debug){
    //cout << waters;
  }
  //stage 1.5 optimize the origin of the lattice
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
  //#stage 2: determine the double lattices
  int mygroup[waters.size()];
  int latticeA[lattice*lattice*lattice];
  int latticeB[lattice*lattice*lattice];
  for(int i=0;i<lattice*lattice*lattice; i++){
    latticeA[i] = -1;
    latticeB[i] = -1;
  }
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
    mygroup[i] = (dA < dB);
    if ( mygroup[i] ){
      int ix = rint( waters[i].com[0] / grid[0] );
      ix = (ix + lattice) % lattice;
      int iy = rint( waters[i].com[1] / grid[1] );
      iy = (iy + lattice) % lattice;
      int iz = rint( waters[i].com[2] / grid[2] );
      iz = (iz + lattice) % lattice;
      int j = (ix*lattice+iy)*lattice+iz;
      if ( latticeA[j] >= 0 ){
	fprintf(stderr, "Two residents in a cell %d: %d,%d.\n", j, i, latticeA[j]);
	exit(1);
      }
      latticeA[j] = i;
    }
    else{
      int ix = rint( (waters[i].com[0]-grid[0]/2) / grid[0] );
      ix = (ix + lattice) % lattice;
      int iy = rint( (waters[i].com[1]-grid[1]/2) / grid[1] );
      iy = (iy + lattice) % lattice;
      int iz = rint( (waters[i].com[2]-grid[2]/2) / grid[2] );
      iz = (iz + lattice) % lattice;
      int j = (ix*lattice+iy)*lattice+iz;
      if ( latticeB[j] >= 0 ){
	fprintf(stderr, "Two residents in a cell %d: %d,%d.\n", j, i, latticeB[j]);
	exit(1);
      }
      latticeB[j] = i;
    }
  }
  //cerr << dsum << " Total square displacements from the lattice points." << endl;
  //#check the number of groups (should be 2)
  //#but you may have to loosen the condition.
  map<int, int> groups;
  for(int i=0;i<waters.size(); i++){
    int myg = mygroup[i];
    if (groups.find(myg)==groups.end()){
      groups[myg] = 0;
    }
    groups[myg] += 1;
  }
  typedef map<int, int>::const_iterator CI;
  for(CI p = groups.begin(); p != groups.end(); p++ ){
    cerr << p->second <<" ";
  }
  cerr << endl;
  if ( (groups.size() != 2) || (groups[0] != waters.size()/2) ){
    exit(1);
  }
  else{
    cerr << "No problem. Go ahead!" << endl;
  }
  //#stage 3: determine the HB networks
  //molecules must be on the lattice

  if ( mode == 3 ){
    cout << "@NGPH" << endl;
    cout << waters.size() << endl;
  }
  vector<int> hb[waters.size()];
  //lattice B
  for(int ix=0;ix<lattice;ix++){
    for(int iy=0;iy<lattice;iy++){
      for(int iz=0;iz<lattice;iz++){
	int Lb = (ix*lattice+iy)*lattice+iz;
	int i = latticeB[Lb];
	if ( i >= 0 ){
	  //lattice A
	  for(int jx=0;jx<2;jx++){
	    for(int jy=0;jy<2;jy++){
	      for(int jz=0;jz<2;jz++){
		int kx = (ix+jx) % lattice;
		int ky = (iy+jy) % lattice;
		int kz = (iz+jz) % lattice;
		int La = (kx*lattice+ky)*lattice+kz;
		int j = latticeA[La];
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
		    //#Å; threshold for O-H adjacency
		    //#register as a hydrogen bond.
		    //#bond direction is ignored.
		    if ( mode == 3 ){
		      if ( di < dj ){
			cout << i << " " << j << endl;
		      }
		      else{
			cout << j << " " << i << endl;
		      }
		    }
		    hb[i].push_back(j); //#hb[i] contains the labels of the HB partners
		    hb[j].push_back(i);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if ( mode == 3 ){
    cout << "-1 -1" << endl;
    exit(0);
  }
  if (debug){
    for(int i=0;i<hb[0].size();i++){
      cout << hb[0][i] << " ";
    }
    cout << endl;
  }
  //#stage 4: determine the order parameters and output it with distance
  cout << "@VMRK" << endl;
  cout << waters.size() << endl;
  vector<float> op;
  for(int i=0;i<waters.size(); i++){
    op.push_back(OrderParameter2(i,hb[i],waters,box));
    cout << op[i] << endl;
  }

    int cluster[waters.size()];
    for(int i=0;i<waters.size(); i++){
      //#-1 means an isolated node.
      cluster[i] = -1;
    }
    //lattice A
    for(int ix=0;ix<lattice;ix++){
      for(int iy=0;iy<lattice;iy++){
	for(int iz=0;iz<lattice;iz++){
	  int j = (ix*lattice+iy)*lattice+iz;
	  int res = latticeA[j];
	  if ( res >= 0 ){
	    int ixi = (ix+1) % lattice;
	    int jx = (ixi*lattice+iy)*lattice+iz;
	    int resx = latticeA[jx];
	    if ( (resx >= 0 ) && ( op[res] * op[resx] > 0 ) ){
	      //same sign
	      BindNodes(res,resx,cluster);
	    }
	    int iyi = (iy+1) % lattice;
	    int jy = (ix*lattice+iyi)*lattice+iz;
	    int resy = latticeA[jy];
	    if ( (resy >= 0 ) && ( op[res] * op[resy] > 0 ) ){
	      //same sign
	      BindNodes(res,resy,cluster);
	    }
	    int izi = (iz+1) % lattice;
	    int jz = (ix*lattice+iy)*lattice+izi;
	    int resz = latticeA[jz];
	    if ( (resz >= 0 ) && ( op[res] * op[resz] > 0 ) ){
	      //same sign
	      BindNodes(res,resz,cluster);
	    }
	  }
	  res = latticeB[j];
	  if ( res >= 0 ){
	    int ixi = (ix+1) % lattice;
	    int jx = (ixi*lattice+iy)*lattice+iz;
	    int resx = latticeB[jx];
	    if ( (resx >= 0 ) && ( op[res] * op[resx] > 0 ) ){
	      //same sign
	      BindNodes(res,resx,cluster);
	    }
	    int iyi = (iy+1) % lattice;
	    int jy = (ix*lattice+iyi)*lattice+iz;
	    int resy = latticeB[jy];
	    if ( (resy >= 0 ) && ( op[res] * op[resy] > 0 ) ){
	      //same sign
	      BindNodes(res,resy,cluster);
	    }
	    int izi = (iz+1) % lattice;
	    int jz = (ix*lattice+iy)*lattice+izi;
	    int resz = latticeB[jz];
	    if ( (resz >= 0 ) && ( op[res] * op[resz] > 0 ) ){
	      //same sign
	      BindNodes(res,resz,cluster);
	    }
	  }
	}
      }
    }
    if ( mode == 0 ){
      exit(0);
    }
  if ( mode >= 1 ) {
    cout << "@GRUP" << endl;
    cout << waters.size() << endl;
    for(int i=0;i<waters.size(); i++){
      //#-1 means an isolated node.
      cout << MyGroup(i,cluster) << endl;
    }
  }
  if ( mode >= 2 ){
    cout << "@OPSD" << endl;
    cout << waters.size() << endl;
    //lattice A
    for(int ix=0;ix<lattice;ix++){
      for(int iy=0;iy<lattice;iy++){
	for(int iz=0;iz<lattice;iz++){
	  int Li = (ix*lattice+iy)*lattice+iz;
	  int i = latticeA[Li];
	  if ( i >= 0 ){
	    float opi = op[i];
	    //lattice A
	    for(int jx=0;jx<lattice;jx++){
	      for(int jy=0;jy<lattice;jy++){
		for(int jz=0;jz<lattice;jz++){
		  int Lj = (jx*lattice+jy)*lattice+jz;
		  int j = latticeA[Lj];
		  if ( j > i ){
		    Vector3d delta = Wrap(waters[i].com-waters[j].com,box);
		    float opj = op[j];
		    //#i and j are on the same sublattice
		    double dist  = delta.norm();
		    if (dist < box(0)/2){
		      cout << i << " " << j << " " << dist << " " << opi*opj << endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    cout << "-1 -1 0 0\n";
  }
}





int main(int argc, char* argv[])
{
  int mode = 1;
  //mode 0 : raw OP
  //mode 1 : +cluster size
  //mode 2 : +spatial correlation
  //mode 3 : network topology only
  //mode 4 : CoM coordinate in @AR3A.
  if ( argc == 2 ){
    if ( strcmp(argv[1], "-0") == 0 ){
      mode = 0;
    }
    if ( strcmp(argv[1], "-r") == 0 ){
      mode = 2;
    }
    else if ( strcmp(argv[1], "-n") == 0 ){
      mode = 3;
    }
    else if ( strcmp(argv[1], "-c") == 0 ){
      mode = 4;
    }
  }
  oneframe(mode);
}
