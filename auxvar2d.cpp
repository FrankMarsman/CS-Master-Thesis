#include "auxvar2d.h"

#include <QDebug>
#include <cmath>

void AuxVar2D::AddToPosVec(VectorXd &q) {
  qDebug( ) << "AuxVar2D. AddToPosVec of size" << q.rows( );
  return;
} // AddToPosVec

AuxVar2D::AuxVar2D(double _w) {
  this->w = _w;
} // AuxVar2D

AuxVar2D::~AuxVar2D( ) {

} // ~AuxVar2D


//
//
// SPRING POTENTIAL
//
//


// adds 'desired' positions of vertices v1 and v2 to q,
// weighted by w. function does not compute positions!
// make sure x1, y1, x2, y2 are correct before calling!
void SpringPotential2D::AddToPosVec(VectorXd &q) {
  q(ix1) += w * x1;
  q(iy1) += w * y1;
  q(ix2) += w * x2;
  q(iy2) += w * y2;
} // AddToPosVec

// returns potential energy of spring system (E = 0.5 * C u^2), C = w
double SpringPotential2D::GetPotE(double _x1, double _y1, double _x2, double _y2) {
  double dx = _x1 - _x2;
  double dy = _y1 - _y2;
  double dist = sqrt(dx*dx + dy*dy);
  double u = dist - r0;
  return 0.25 * w * u * u;
} // GetPotE

SpringPotential2D::SpringPotential2D(uint _v1, uint _v2, double _r0, double _w) : AuxVar2D(_w) {
  this->v1 = _v1;
  this->v2 = _v2;
  this->r0 = _r0;

  // so that we can add positions faster
  this->ix1 = v1 * 2;
  this->iy1 = v1 * 2 + 1;
  this->ix2 = v2 * 2;
  this->iy2 = v2 * 2 + 1;
} // SpringPotential2D
























