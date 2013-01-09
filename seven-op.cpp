
/*#read water configurations in Gromacs format,
#calculate the pairwise products of local order parameter for ice 7,
#and output them with pair distances.

#2012-10-19 tested with ice VII crystal structure.
#2012-12-27 Renamed (formerly seven-op-gro3.cpp)
*/


#include "seven-op.h"
#include "analysis.h"




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
    output_coordinates(box, waters);
    exit(0);
  }
  int division = latticeSize(waters.size()); //assume subic box
  //cerr << division << " division" << endl;
  Vector3d grid;
  for(int dim=0;dim<3;dim++){
    grid[dim] = box[dim] / division;
  }
  if (waters.size() ==0){
    return 0;
  }
  if (debug){
    //cout << waters;
  }
  //stage 1.5 optimize the origin of the lattice
  optimize_water_positions(waters, grid);
  //#stage 2: determine the double lattices
  int mycubiclattice[waters.size()];
  int residentA[division*division*division];
  int residentB[division*division*division];
  int myaddress[waters.size()];
  for(int i=0;i<division*division*division; i++){
    residentA[i] = -1;
    residentB[i] = -1;
  }
  address_on_the_lattice( waters, grid, division, mycubiclattice, residentA, residentB, myaddress );

  //#check the number of groups (should be 2)
  //#but you may have to loosen the condition.
  map<int, int> groups;
  for(int i=0;i<waters.size(); i++){
    int myg = mycubiclattice[i];
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

  vector<int> uhb[waters.size()];//undirected
  vector<int> dhb[waters.size()];//directed
  hydrogen_bond(residentA, residentB, box, waters, division, uhb, dhb);
  if ( mode == 3 ){
    cout << "@NGPH" << endl;
    cout << waters.size() << endl;
    for(int i=0;i<waters.size(); i++){
      for(int k=0; k<dhb[i].size(); k++){
	int j = dhb[i][k];
	cout << i << " " << j << endl;
      }
    }
    cout << "-1 -1" << endl;
    exit(0);
  }
  if (debug){
    for(int i=0;i<dhb[0].size();i++){
      cout << dhb[0][i] << " ";
    }
    cout << endl;
  }
  //#stage 4: determine the order parameters and output it with distance
  cout << "@VMRK" << endl;
  cout << waters.size() << endl;
  vector<float> op;
  for(int i=0;i<waters.size(); i++){
    op.push_back(OrderParameter2(i,uhb[i],waters,box));
    cout << op[i] << endl;
  }
  //////topological quench (structure optimization)
  //8 neighbor lattice point list
  //first four are negative op directions, and the rest are positive op directions.
  vector<int> neighborA[division*division*division];
  vector<int> neighborB[division*division*division];
  list_neighbors(division, neighborA, neighborB);
  vector<int> assign[waters.size()];  //optimized network
  estimate_optimized_network(waters,mycubiclattice,neighborA,
			     neighborB,op,residentA,residentB,myaddress,assign);
  if ( mode == 5 ){
    cout << "@NGPH" << endl;
    cout << waters.size() << endl;
    //output only the "regular" bonds
    for(int i=0;i<waters.size(); i++){
      for(int k=0; k<dhb[i].size(); k++){
	int j = dhb[i][k];
	if ( assign[i][j] == 1 ){
	  cout << i << " " << j << endl;
	}
      }
    }
    cout << "-1 -1" << endl;
    exit(0);
  }
	
  ////////cluster analysis
  int cluster[waters.size()];
  for(int i=0;i<waters.size(); i++){
    //#-1 means an isolated node.
    cluster[i] = -1;
  }
  //lattice A
  for(int ix=0;ix<division;ix++){
    for(int iy=0;iy<division;iy++){
      for(int iz=0;iz<division;iz++){
	int j = address(ix,iy,iz,division);
	int res = residentA[j];
	if ( res >= 0 ){
	  int ixi = (ix+1) % division;
	  int jx = address(ixi,iy,iz,division);
	  int resx = residentA[jx];
	  if ( (resx >= 0 ) && ( op[res] * op[resx] > 0 ) ){
	    //same sign
	    BindNodes(res,resx,cluster);
	  }
	  int iyi = (iy+1) % division;
	  int jy = address(ix,iyi,iz,division);
	  int resy = residentA[jy];
	  if ( (resy >= 0 ) && ( op[res] * op[resy] > 0 ) ){
	    //same sign
	    BindNodes(res,resy,cluster);
	  }
	  int izi = (iz+1) % division;
	  int jz = address(ix,iy,izi,division);
	  int resz = residentA[jz];
	  if ( (resz >= 0 ) && ( op[res] * op[resz] > 0 ) ){
	    //same sign
	    BindNodes(res,resz,cluster);
	  }
	}
	res = residentB[j];
	if ( res >= 0 ){
	  int ixi = (ix+1) % division;
	  int jx = address(ixi,iy,iz,division);
	  int resx = residentB[jx];
	  if ( (resx >= 0 ) && ( op[res] * op[resx] > 0 ) ){
	    //same sign
	    BindNodes(res,resx,cluster);
	  }
	  int iyi = (iy+1) % division;
	  int jy = address(ix,iyi,iz,division);
	  int resy = residentB[jy];
	  if ( (resy >= 0 ) && ( op[res] * op[resy] > 0 ) ){
	    //same sign
	    BindNodes(res,resy,cluster);
	  }
	  int izi = (iz+1) % division;
	  int jz = address(ix,iy,izi,division);
	  int resz = residentB[jz];
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
  //spatial correlation
  if ( mode >= 2 ){
    cout << "@OPSD" << endl;
    cout << waters.size() << endl;
    //lattice A
    for(int ix=0;ix<division;ix++){
      for(int iy=0;iy<division;iy++){
	for(int iz=0;iz<division;iz++){
	  int Li = address(ix,iy,iz,division);
	  int i = residentA[Li];
	  if ( i >= 0 ){
	    float opi = op[i];
	    //lattice A
	    for(int jx=0;jx<division;jx++){
	      for(int jy=0;jy<division;jy++){
		for(int jz=0;jz<division;jz++){
		  int Lj = address(jx,jy,jz,division);
		  int j = residentA[Lj];
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
    if ( strcmp(argv[1], "-0") == 0 ){      // raw order parameter only
      mode = 0;
    }
    if ( strcmp(argv[1], "-r") == 0 ){      // +clustering & spatial corr
      mode = 2;
    }
    else if ( strcmp(argv[1], "-n") == 0 ){  // network only
      mode = 3;
    }
    else if ( strcmp(argv[1], "-c") == 0 ){ // coordinate only
      mode = 4;
    }
    else if ( strcmp(argv[1], "-q") == 0 ){ // topological quench
      mode = 5;
    }
  }
  oneframe(mode);
}
