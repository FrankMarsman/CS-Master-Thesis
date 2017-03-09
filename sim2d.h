#ifndef SIM2D_H
#define SIM2D_H

#include "auxvar2d.h"

#include <Eigen>
#include <vector>
#include <QImage>
#include <iostream>
#include <Sparse>
#include <QString>
#include <SparseCholesky>
#include <omp.h>

using namespace Eigen;
using namespace std;

class Sim2D
{
public:
  typedef unsigned int uint;

  uint m; // number of vertices
  uint iterationsPerStep; // number of local/global iterations per step
  uint simStep; // how many steps computed

  double meshMass; // total mass of all vertices

  double timeInSim; // how many seconds have elapsed

  double h; // time step
  double springForceConstant; // [N/m] F_spring = -Cu (C is force constant)

  // stats:
  double qSolveTime; // [ms] time to solve lMatrix * q = b
  double pSolveTime; // [ms] time to compute p_i's
  double rMatrixMakeTime; // [ms] time to construct rMatrix
  double minSpringLen; // [m] smallest length of any spring ever

  // stats for during NextStep
  vector <double> STEP_disp; // total displacement after each iteration
  vector <double> STEP_totDiffE; // energy loss t.o.v. previous NextStep

  VectorXd q; // position vector in R^2m
  VectorXd v; // velocity vector in R^2m
  VectorXd F; // external forces

  SparseMatrix < double > M; // mass matrix of size 2m * 2m
  SparseMatrix < double > Minv; // inverse of mass matrix
  SparseMatrix < double > lMatrix; // constant system matrix

  // solvers:
  ConjugateGradient < SparseMatrix < double >, Lower|Upper > cg; // solver for lMatrix
  SimplicialLLT < SparseMatrix < double > > cholenskySolver;


  uint numEdges; // number of edges
  MatrixXi E; // numEdges * 2 matrix containing list of edges

  vector < SpringPotential2D > pVec; // list of auxiliary variables
  vector < double > restLen; // rest lengths of all springs

  // variables for drawing and manipulating
  vector < uint > selectedVertices; // list of selected vertices
  vector < uint > lockedVertices; // list of vertices that cannot move

  double imgCenterX, imgCenterY; // position of center of QImage
  double imgViewSize; // how large is area inside image

  // position functions
  double MinX( );
  double MinY( );
  double MaxX( );
  double MaxY( );

  double GetCenterOfMassY( ); // returns y-pos of COM
  double GetCenterOfMassX( ); // returns x-pos of COM

  // velocity functions
  double GetV(int k); // returns velocity of vertex n
  double MaxV( ); // max velocity
  double MinV( ); // min velocity

  // Energy functions:
  double GetKineticEnergy( ); // returns total kinetic energy of all vertices
  double GetSpringPotentialEnergy( ); // returns total potential energy of all potentials
  double GetGravPotEnergy( ); // gravitational potential energy w.r.t. gravRefHeight
  double GetTotEnergy( ); // sum of E_kin, E_spring and E_grav

  QString GetInfoString( ); // string with info about sim

  void InitMesh_Square(double _meshSize, bool dSprings); // arranges vertices in square grid
  void InitRestLen( );
  void InitPs( ); // initializes pVec
  void ComputePs( ); // computes auxilliary variables (updates pVec elements)
  void InitLMatrix( );

  double gravAcc; // gravitational acceleration
  double gravRefHeight; // measure Gravitational energy from this height

  void SetGravity(double g); // sets F with gravity and updates gravAcc
  void SetSpringForceConstant(double C); // updates springForceConstant (also in pVec)

  void SetSelectedVertices(double x1, double x2, double y1, double y2);
  void AddLockedVertices( );
  void ResetLockedVertices(bool forcesToo = true); // sets forces and velocity of locked vertices to 0

  void InitImgParams( ); // initializes imgCenterX, imgCenterY and imgViewSize

  // abuse functions:
  void AddNoise(double strength = 0.025); // displaces selected vertices randomly
  void AddVNoise(double strength = 0.025); // adds random velocity to selected vertices
  void AddVelocity(double vx, double vy); // adds given velocity to selected vertices
  void SqueezeX(double factor = 0.9); // squeezes mesh in x-direction by factor
  void SqueezeY(double factor = 0.9); // squeezes mesh in y-direction by factor

  void NextStep(bool useLLT = false, bool doMeasures = false); // compute next q

  QImage ToQImage(uint SIZE, bool useAA = true, bool drawEdges = true,
                  bool drawSelectedVertices = false, bool drawLockedVertices = false,
                  uint DRAW_MODE = 0);

  Sim2D(uint _m, double _meshMass = 1, double _meshSize = 1, bool dSprings = true);
  ~Sim2D( );
};

#endif // SIM2D_H
