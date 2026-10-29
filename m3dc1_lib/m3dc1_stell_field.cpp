#include "m3dc1_stell_field.h"

#include <iostream>

m3dc1_stell_field::m3dc1_stell_field(m3dc1_mesh* mesh_in, m3dc1_coord_map* map_in)
  : m3dc1_3d_field(mesh_in)
{
  map = map_in;
}

bool m3dc1_stell_field::extent(double* X0, double* X1,
			 double* Phi0, double* Phi1,
			 double* Z0, double* Z1) const
{
  return map->extent(X0, X1, Phi0, Phi1, Z0, Z1);
}


bool m3dc1_stell_field::eval(const double R, const double Phi, const double Z, 
			     const m3dc1_field::m3dc1_get_op op, double* val, 
			     int* element)
{
  double R0[OP_NUM], Z0[OP_NUM];
  double x, phi, y;
  int refine=20;
  
  // find logical coordinates associated with R, Phi, Z
  if(!map->find_coordinates(R,Phi,Z,&x,&phi,&y,element,refine))
    return false;

  // evaluate field on logical mesh
  m3dc1_get_op logical_op = op;
  const bool second = ((op & GET_DDVAL) != 0) || ((op & GET_PPVAL) != 0)
    || (((op & GET_DVAL) != 0) && ((op & GET_PVAL) != 0));
  if(second) logical_op = GET_ALL;
  if(!m3dc1_3d_field::eval(x,phi,y,logical_op,val,element))
    return false;

  if(op==GET_VAL)
    return true;

  // translate derivatives (if requested)
  // This follows Zhou et al 2021 Nucl. Fusion 61 086015

  if(!map->eval_map_deriv(x,phi,y,R0,Z0,element)) {
    std::cerr << "Error evaluating map derivatives" << std::endl;
    return false;
  }

  // derivatives of v, R0, and Z0 here are wrt to logical coordinates
  double v[OP_NUM];
  for(int i=0; i<OP_NUM; i++) {
    v[i] = val[i];
    if(i != OP_1) val[i] = 0;
  }

  double A = R0[OP_DZ]*Z0[OP_DP] - R0[OP_DP]*Z0[OP_DZ];
  double B = R0[OP_DP]*Z0[OP_DR] - R0[OP_DR]*Z0[OP_DP];
  double D = R0[OP_DR]*Z0[OP_DZ] - R0[OP_DZ]*Z0[OP_DR];

  // First derivatives
  val[OP_DR] = ( Z0[OP_DZ]*v[OP_DR] - Z0[OP_DR]*v[OP_DZ])/D;
  val[OP_DZ] = (-R0[OP_DZ]*v[OP_DR] + R0[OP_DR]*v[OP_DZ])/D;
  val[OP_DP] = v[OP_DP] + (A*v[OP_DR] + B*v[OP_DZ])/D;

  // Hessian chain rule: H_q(v) = J^T H_X(v) J
  //                         + v_R H_q(R) + v_Z H_q(Z).
  // Columns of K are logical derivatives with respect to (R,Z,Phi).
  // First derivatives above retain the original arithmetic and conventions.
  const double K[3][3] = {{ Z0[OP_DZ]/D, -R0[OP_DZ]/D, A/D },
                         {-Z0[OP_DR]/D,  R0[OP_DR]/D, B/D },
                         { 0., 0., 1. }};
  const int index[3][3] = {{OP_DRR,OP_DRZ,OP_DRP},
                           {OP_DRZ,OP_DZZ,OP_DZP},
                           {OP_DRP,OP_DZP,OP_DPP}};
  double M[3][3], H[3][3];
  for(int i=0;i<3;++i) for(int j=0;j<3;++j)
    M[i][j] = v[index[i][j]] - val[OP_DR]*R0[index[i][j]]
                               - val[OP_DZ]*Z0[index[i][j]];
  for(int a=0;a<3;++a) for(int b=0;b<3;++b) {
    H[a][b] = 0.;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
      H[a][b] += K[i][a]*M[i][j]*K[j][b];
  }
  val[OP_DRR]=H[0][0]; val[OP_DRZ]=H[0][1]; val[OP_DZZ]=H[1][1];
  val[OP_DRP]=H[0][2]; val[OP_DZP]=H[1][2]; val[OP_DPP]=H[2][2];
  // Higher mixed derivatives are outside this correction's tested contract.

  return true;
}
