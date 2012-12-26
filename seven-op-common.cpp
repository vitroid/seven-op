//############################################################################ Quaternion
#include "seven-op.h"

Matrix3d quat2rotmat(Vector4d q)
{
  double a,b,c,d;
  a = q[0];
  b = q[1];
  c = q[2];
  d = q[3];
  Matrix3d rot;
  rot(0,0) = (a*a+b*b-(c*c+d*d));
  rot(0,1) = -2.0*(a*d+b*c);
  rot(0,2) = 2.0*(b*d-a*c);
  rot(1,0) = 2.0*(a*d-b*c);
  rot(1,1) = a*a+c*c-(b*b+d*d);
  rot(1,2) = -2.0*(a*b+c*d);
  rot(2,0) = 2.0*(a*c+b*d);
  rot(2,1) = 2.0*(a*b-c*d);
  rot(2,2) = a*a+d*d-(b*b+c*c);
  return rot;
}


//############################################################### Periodic boundary utility
//#relative position vector at the periodic boundary condition
Vector3d Wrap( Vector3d vector, Vector3d& box )
{
  for(int dim=0; dim<3; dim++){
    vector(dim) -= rint( vector(dim) / box(dim) ) * box(dim);
  }
  return vector;
}




double OrderParameter2( int site, vector<int>& neighbors, sWaters& waters, Vector3d& box )
{
  Vector3d& com = waters[site].com;
  double op = 0.0;
  for(int n=0; n<neighbors.size(); n++){
    Vector3d& ncom = waters[neighbors[n]].com;
    Vector3d dir = Wrap(ncom - com, box);
    int parity = 1;
    for(int dim=0;dim<3;dim++){
      if (dir(dim) < 0){
	parity *= -1;
      }
    }
    op += parity;
  }
  return op / 4.0;
}


float min(vector<float> v)
{
  float x = v[0];
  for(int i=1; i<v.size(); i++){
    if ( v[i] < x ){
      x = v[i];
    }
  }
  return x;
}




pair<sWaters, Vector3d>
Configure_gro(FILE* in)
{
  char line[1024];
  sWaters waters;
  Vector3d box;
  //ignore the first line
  char* c = fgets(line,sizeof(line),in);
  c = fgets(line,sizeof(line),in);
  int nsite;
  sscanf(line,"%d", &nsite);
  for(int i=0;i<nsite;i++){
    sWater water;
    c = fgets(line,sizeof(line),in);
    char molnum[6], molname[6], elemname[6];
    int sitenum;
    float x,y,z;
    sscanf(line,"%5c%5c%5c%5d%8f%8f%8f", molnum, molname, elemname, &sitenum, &x,&y,&z);
    elemname[5] = '\0';
    //cout << elemname << endl; // waters.size();
    
    Vector3d xyz(x*10,y*10,z*10); //nm to AA
    if ( (0 == strcmp(elemname,"   OW")) || (0 == strcmp(elemname,"  OW ")) ){
      water.com = xyz;
    }
    else if (0 == strcmp(elemname,"  HW1")){
      water.rH1 = xyz;
    }
    else if (0 == strcmp(elemname,"  HW2")){
      water.rH2 = xyz;
      waters.push_back(water);
    }
  }
  c = fgets(line,sizeof(line),in);
  float bx,by,bz;
  sscanf(line,"%f %f %f", &bx,&by,&bz);
  box(0) = bx * 10; //nm to AA
  box(1) = by * 10;
  box(2) = bz * 10;
  //#convert hydrogen positions from absolute to relative
  for(int i=0; i<waters.size(); i++){
    sWater& water = waters[i];
    water.rH1 = Wrap(water.rH1-water.com,box);
    water.rH2 = Wrap(water.rH2-water.com,box);
  }
  return make_pair(waters, box);
}



int latticeSize(int nmol)
{
  float xmol = nmol;
  xmol /= 2;
  xmol = floor(exp((1./3.)*log(xmol))+0.5);
  if ( xmol*xmol*xmol*2 != nmol ){
    printf("Invalid lattice size: %d\n", nmol);
    exit(2);
  }
  return (int)xmol;
}
    



//############################################################################ Clustering
//#The cluster connectivity is saved in the "group dictionary" as a tree structure.
//#The value of the group dictionary, group[i], indicates:
//#group[i]>=0 ==> group[i] is the parent node for node i
//#group[i]<0  ==> i is the root node; -group[i] is number of nodes in the same tree.

//#x: node label
//#group[]: group dictionary
//#return value: root node of the group.



int MyGroup(int x,int group[])
{
  while (group[x] >= 0){
    x = group[x];
  }
  return x;
}


//#register the nodes x and y as connected
void BindNodes(int x, int y, int group[])
{
  int xg = MyGroup(x,group);
  int yg = MyGroup(y,group);
  if (xg != yg){
    //#merge the trees
    group[yg] += group[xg];
    group[xg] =  yg;
  }
}
