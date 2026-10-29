/*
 * Analytic regression for the mapped scalar derivative contract.
 * These are local mathematical charts, not plasma equilibria/coupling examples.
 * No HDF data, tracing, fitting, flux-surface reduction or physical closure.
 *
 * Derivation: let X(q)=(R(q),Z(q),phi), J_ai=dX_a/dq_i and K=J^{-1}.
 * For v(q)=F(X(q)), g_q=J^T g_X, and
 *   H_q(v)=J^T H_X(F) J + F_R H_q(R) + F_Z H_q(Z).
 * Thus H_X=K^T [H_q(v)-F_R H_q(R)-F_Z H_q(Z)] K.
 * Coordinate functions R and Z have exactly zero physical Hessians. Also
 * F=R^2+2RZ+3Z^2+phi R+2phi Z+phi^2 has the constant independent Hessian
 *   [[2,2,1],[2,6,2],[1,2,2]] in (R,Z,phi) order.
 * For R=3+x+alpha phi y, Z=y and v=x=R-3-alpha phi Z,
 *   d^2 v/(dZ dphi)=-alpha, with all other Hessian entries zero.
 *
 * Bound=4096*epsilon*max(1,abs(expected)) is a double-precision software
 * oracle tolerance on these well-conditioned fixtures. It is NOT a magnetic,
 * geometry, tangency or flux-acceptance threshold.
 */
#include "m3dc1_stell_field.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using Power=std::array<int,3>;
using Poly=std::map<Power,double>;
static Poly constant(double x){return {{{0,0,0},x}};}
static Poly monomial(int a,int b,int c,double x=1.){return {{{a,b,c},x}};}
static Poly add(Poly a,const Poly& b,double scale=1.){
  for(const auto& t:b)a[t.first]+=scale*t.second;return a;
}
static Poly mul(const Poly& a,const Poly& b){
  Poly c;
  for(const auto& x:a)for(const auto& y:b){
    Power p{{x.first[0]+y.first[0],x.first[1]+y.first[1],x.first[2]+y.first[2]}};
    c[p]+=x.second*y.second;
  }
  return c;
}
static double eval(const Poly& a,const std::array<double,3>& q){
  double x=0.;for(const auto& t:a)x+=t.second*std::pow(q[0],t.first[0])*std::pow(q[1],t.first[1])*std::pow(q[2],t.first[2]);return x;
}
static void fill(m3dc1_field& field,const Poly& p){
  std::fill(field.data,field.data+80,0.);
  for(const auto& t:p){
    if(t.second==0.)continue;
    bool found=false;
    for(int j=0;j<20;++j)if(m3dc1_field::mi[j]==t.first[0] && m3dc1_field::ni[j]==t.first[1]){
      if(t.first[2]<0 || t.first[2]>=4)throw std::runtime_error("test polynomial outside native toroidal basis");
      field.data[20*t.first[2]+j]=t.second;found=true;break;
    }
    if(!found)throw std::runtime_error("test polynomial outside native poloidal basis");
  }
}
static const int hessian_ops[6]={m3dc1_field::OP_DRR,m3dc1_field::OP_DRZ,m3dc1_field::OP_DZZ,
  m3dc1_field::OP_DRP,m3dc1_field::OP_DZP,m3dc1_field::OP_DPP};
struct Observation{std::string fixture,field;int point,request,op;double actual,expected,bound;bool pass;};
static std::vector<Observation> observations;
static int failures=0;
static void check(const std::string& fixture,const std::string& field,int point,int request,int op,double actual,double expected){
  double bound=4096*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(expected));
  bool pass=std::isfinite(actual) && std::abs(actual-expected)<=bound;
  observations.push_back({fixture,field,point,request,op,actual,expected,bound,pass});
  if(!pass)++failures;
}
static std::array<double,18> sample(m3dc1_stell_field& field,const std::array<double,3>& x,int request){
  std::array<double,18> out;out.fill(std::numeric_limits<double>::quiet_NaN());int hint=0;
  if(!field.eval(x[0],x[2],x[1],static_cast<m3dc1_field::m3dc1_get_op>(request),out.data(),&hint) || hint!=0)
    throw std::runtime_error("fixture inverse failed or changed cell");
  return out;
}
static void fixture(const std::string& name,bool curved,bool shear){
  m3dc1_3d_mesh mesh(1);
  mesh.a[0]=1.;mesh.b[0]=1.;mesh.c[0]=1.;mesh.co[0]=1.;mesh.sn[0]=0.;
  mesh.x[0]=0.;mesh.z[0]=0.;mesh.bound[0]=0;mesh.region[0]=0;
  mesh.d[0]=1.;mesh.phi[0]=0.;mesh.nplanes=1;mesh.nperiods=1;mesh.period=2*std::acos(-1.);mesh.toroidal=true;
  mesh.nneighbors=new int[1]{0};mesh.neighbor=new int*[1];mesh.neighbor[0]=new int[5]{};
  const Poly x=monomial(1,0,0),y=monomial(0,1,0),p=monomial(0,0,1);
  Poly R=add(constant(3.),x),Z=y;
  if(curved){R=add(R,mul(x,y),.12);R=add(R,mul(p,y),.18);Z=add(Z,mul(x,x),.07);Z=add(Z,mul(p,x),.09);}
  const double alpha=.2;
  if(shear)R=add(R,mul(p,y),alpha);
  Poly F=mul(R,R);F=add(F,mul(R,Z),2.);F=add(F,mul(Z,Z),3.);
  F=add(F,mul(p,R));F=add(F,mul(p,Z),2.);F=add(F,mul(p,p));
  m3dc1_3d_field rawR(&mesh),rawZ(&mesh);fill(rawR,R);fill(rawZ,Z);
  m3dc1_coord_map map(&mesh);if(!map.load(&rawR,&rawZ))throw std::runtime_error("fixture map construction failed");
  m3dc1_stell_field rfield(&mesh,&map),zfield(&mesh,&map),ffield(&mesh,&map),xfield(&mesh,&map);
  fill(rfield,R);fill(zfield,Z);fill(ffield,F);fill(xfield,x);
  const std::array<std::array<double,3>,3> points{{{{-.15,.15,.2}},{{.2,.25,.4}},{{-.25,.2,.6}}}};
  for(int i=0;i<3;++i){
    const auto& q=points[i];std::array<double,3> X{{eval(R,q),eval(Z,q),q[2]}};
    auto rv=sample(rfield,X,31),zv=sample(zfield,X,31);
    for(int field=0;field<2;++field){
      const auto& a=field==0?rv:zv;std::string label=field==0?"R":"Z";
      check(name,label,i,31,1,a[1],field==0?1.:0.);check(name,label,i,31,2,a[2],field==1?1.:0.);check(name,label,i,31,6,a[6],0.);
      for(int op:hessian_ops)check(name,label,i,31,op,a[op],0.);
    }
    const double r0=rv[0],z0=zv[0],phi=q[2];
    const double expected_value=r0*r0+2*r0*z0+3*z0*z0+phi*r0+2*phi*z0+phi*phi;
    const double expected_h[6]={2.,2.,6.,1.,2.,2.};
    auto all=sample(ffield,X,31);
    check(name,"quadratic",i,31,0,all[0],expected_value);
    check(name,"quadratic",i,31,1,all[1],2*r0+2*z0+phi);
    check(name,"quadratic",i,31,2,all[2],2*r0+6*z0+2*phi);
    check(name,"quadratic",i,31,6,all[6],r0+2*z0+2*phi);
    for(int j=0;j<6;++j)check(name,"quadratic",i,31,hessian_ops[j],all[hessian_ops[j]],expected_h[j]);
    for(int request:{4,10,16}){
      auto partial=sample(ffield,X,request);
      for(int j=0;j<6;++j){
        bool requested=(request==4 && j<3)||(request==10 && j>=3 && j<5)||(request==16 && j==5);
        if(requested)check(name,"quadratic",i,request,hessian_ops[j],partial[hessian_ops[j]],expected_h[j]);
      }
    }
    if(shear){
      auto a=sample(xfield,X,31);
      check(name,"logical_x",i,31,0,a[0],r0-3-alpha*phi*z0);
      check(name,"logical_x",i,31,1,a[1],1.);check(name,"logical_x",i,31,2,a[2],-alpha*phi);check(name,"logical_x",i,31,6,a[6],-alpha*z0);
      for(int j=0;j<6;++j)check(name,"logical_x",i,31,hessian_ops[j],a[hessian_ops[j]],j==4?-alpha:0.);
    }
  }
}
int main(int argc,char** argv){
  try{
    fixture("affine",false,false);fixture("curved",true,false);fixture("toroidal_shear",false,true);
    std::ostream* out=&std::cout;std::ofstream file;
    if(argc==2){file.open(argv[1]);if(!file)throw std::runtime_error("cannot open result file");out=&file;}
    *out<<std::setprecision(17)<<"fixture,field,point,request,operator,actual,expected,bound,pass\n";
    for(const auto& a:observations)*out<<a.fixture<<','<<a.field<<','<<a.point<<','<<a.request<<','<<a.op<<','<<a.actual<<','<<a.expected<<','<<a.bound<<','<<a.pass<<'\n';
    std::cerr<<"CHECKS="<<observations.size()<<" FAILURES="<<failures<<'\n';
    return failures?1:0;
  }catch(const std::exception& e){std::cerr<<"SETUP_ERROR: "<<e.what()<<'\n';return 2;}
}
