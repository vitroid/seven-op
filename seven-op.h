#ifndef SEVEN_OP_H
#define SEVEN_OP_H
#include <iostream>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>
#include <utility>   //for pair
#include <map>
using namespace Eigen;
using namespace std;
using Eigen::MatrixXd;

Matrix3d quat2rotmat(Vector4d q);
Vector3d Wrap( Vector3d vector, Vector3d& box );
typedef struct {
  Vector3d com;
  Vector3d rH1;
  Vector3d rH2;
} sWater;

typedef vector<sWater> sWaters;
double OrderParameter2( int site, vector<int>& neighbors, sWaters& waters, Vector3d& box );
float min(vector<float> v);
pair<sWaters, Vector3d> Configure_gro(FILE* in);
int latticeSize(int nmol);
int MyGroup(int x,int group[]);
void BindNodes(int x, int y, int group[]);

#endif 
