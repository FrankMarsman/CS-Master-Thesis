#ifndef AUXVAR2D_H
#define AUXVAR2D_H

// this class represents auxiliary variables p_i's
// that are needed for the projective dynamics simulation
// it contains the 'desired' position of a number of nodes
// that follows from a certain potential

#include <Eigen>
#include <vector>
#include <QImage>
#include <iostream>
#include <Sparse>

using namespace Eigen;
using namespace std;

class AuxVar2D {
public:
  typedef unsigned int uint;

  double w; // w_i, weight of this potential

  virtual void AddToPosVec(VectorXd & q);

  AuxVar2D(double _w);
  ~AuxVar2D( );
}; // AuxVar2D


// represents potential describing two vertices (in 2D) that
// are connected by a spring (Hooke's law) of rest length r0
class SpringPotential2D : public AuxVar2D {
public:
  typedef unsigned int uint;

  uint v1, v2; // indices of connected vertices
  uint ix1, iy1, ix2, iy2; // indices of vertex positions in pos vector

  double r0; // rest length of spring

  // 'desired' positions of vertices
  double x1, y1, x2, y2;
  double dx, dy; // desired dx and dy

  virtual void AddToPosVec(VectorXd & q);

  double GetPotE(double _x1, double _y1, double _x2, double _y2);

  SpringPotential2D(uint _v1, uint _v2, double _r0, double _w = 1);
}; // SpringPotential2D



















#endif // AUXVAR2D_H















