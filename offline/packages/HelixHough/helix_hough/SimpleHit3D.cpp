#include "SimpleHit3D.h"

#include <iostream>

using namespace std;

SimpleHit3D::SimpleHit3D(float x, float dxx,
			 float y, float dyy,
			 float z, float dzz,
			 unsigned int id, int lyr) 
  : _id(id),
    _x(x),
    _y(y),
    _z(z),
    dx(dxx),    
    dy(dyy),
    dz(dzz),
    layer(lyr),
    _err()
{
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,0.0);
    }
  }
}			 

void SimpleHit3D::print(std::ostream& out) const {

  out << "SimpleHit3D: "
      << "id: " << get_id() << " layer: " << layer << " "
      << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ") "
      << "(dx,dy,dz) = (" << dx << "," << dy << "," << dz << ")" << endl;

  out << "       ( ";
  out << get_error(0,0) << " , ";
  out << get_error(0,1) << " , ";
  out << get_error(0,2) << " )" << endl; 
  out << "err  = ( ";
  out << get_error(1,0) << " , ";
  out << get_error(1,1) << " , ";
  out << get_error(1,2) << " )" << endl;
  out << "       ( ";
  out << get_error(2,0) << " , ";
  out << get_error(2,1) << " , ";
  out << get_error(2,2) << " )" << endl;

  return;
}

void SimpleHit3D::set_error(unsigned int i, unsigned int j, float value) {
  _err[covar_index(i,j)] = value;
  return;
}

float SimpleHit3D::get_error(unsigned int i, unsigned int j) const {
  return _err[covar_index(i,j)];
}

unsigned int SimpleHit3D::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}
