/*
 * Copyright 2019 Martins Brics <email>
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef INTEGRATOR1_H
#define INTEGRATOR1_H
#include <vector>
#include <complex>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>
#include <ctime>
#include <functional>
#include <sstream>
#include <chrono>
#include <sys/time.h>
#include <boost/numeric/odeint.hpp>
#include <sys/time.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>


typedef std::complex<double> dcmplx;
typedef std::vector<double> dvec;
typedef std::vector<size_t> ivec;


extern "C" void dgtsv_(const long *Np, const long *NRHSp, double *DL,
	     double *D, double *DU, double *B, const long *LDBp,
		     long *INFOp);

using namespace std;



/**
 * @todo write docs
 */
class Integrator1
{
  size_t Ns;
  size_t Np;
  size_t Nv;
  double alfa;
  double mu;//=p+1;
  double Cm;//=100;
  double Cr=-1.0;
  double H1;//=4.0;
  double om;//=50*16
  double Torq;
  double lam=-1.0;
  double theta=62.0*M_PI/180;
  Eigen::Vector3d hvec={1.0,0.0,0.0};
  Eigen::VectorXd J_Jt_u,J_Jt_l,J_Jt_d,J_DD_Jt_u,J_DD_Jt_l,J_DD_Jt_d, E_small_vec, rhs, ones, x_temp, K1, K2, K3, K4, rl, rll, FT;
  Eigen::MatrixXd A, E, E_small, P, Jt_J_inv, lhs;
  Eigen::VectorXd FA,Ftwist, A_r, J_r, temp_r;
  Eigen::SparseMatrix<double>A_sparse, DD_sparse, Jt, J,lhs_sparse, lhs_sparse1, lhs_sparse2 ;
  Eigen::SparseMatrix<double,Eigen::RowMajor>A_sparse_row, FA_mat, D1_sparse, D2_sparse;
  void fill_Fa(const Eigen::VectorXd &x, double t);
  void fill_Ftwist(const Eigen::VectorXd &x);
  void fill_Fr(const Eigen::VectorXd &x, double t);
  void fill_Fa_mat(const Eigen::VectorXd &x, double t);
  void fill_J(const Eigen::VectorXd &x);
  void fill_J_DD(const Eigen::VectorXd &x);
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_sparse, solver_sparse1, solver_sparse2;
  void make_vector_product(double *a, double *b, double *c);
  //VSIMEX
  double alpha0=0., alpha1=0., alpha2=0., alpha3=0., alpha4=1., beta0=0., beta1=0., beta2=0., beta3=0.,omega1=0.,omega2=0.,omega3=0.;
  
  
  
  long dgtsv(long Nx1, long NRHS, double *DL, double *D, double *DU, double *imag_t_Bx11, long LDB);
  double dist(const Eigen::VectorXd &x , int i);

  
public:
  Integrator1(size_t Ns1, double Cm1, double H11, double om1, double torq1);
  void do_CN(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_CN_propper(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_Euler(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void fill_lhs_Euler_impl_expl(const double dt);
  void fill_lhs_Euler_impl_expl_DD(const double dt);
  void fill_lhs1_Euler_impl_expl(const double dt);
  void fill_lhs2_Euler_impl_expl(const double dt);
  void do_Euler_impl_expl_DD(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_Euler_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_Euler1_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_midpoint(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_rk4_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);
  void do_IMEX_BDF2(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,  Eigen::VectorXd &new_x ,  Eigen::VectorXd &K1 , const double  t , const double dt);
  Eigen::VectorXd fill_K1_IMEX_BDF2(const Eigen::VectorXd &x,const double dt);
  Eigen::VectorXd fill_K1_IMEX_BDF3(const Eigen::VectorXd &x,const double dt, const double t=0);
  Eigen::VectorXd fill_K1_IMEX_BDF3_DD(const Eigen::VectorXd &x,const double dt, const double t=0);
  void do_IMEX_BDF3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt);
  void fill_lhs2_IMEX_BDF3(const double dt);
  void fill_lhs2_IMEX_TVB3(const double dt);
  void do_IMEX_TVB3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt);
  void do_IMEX_BDF5(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,const Eigen::VectorXd &x_old_old3 ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3,Eigen::VectorXd &K4 , const double  t , const double dt);
  void fill_lhs2_IMEX_BDF5(const double dt);
  void fill_lhs2_IMEX_BDF4(const double dt);
  void do_IMEX_BDF4(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3 , const double  t , const double dt);
  void do_IMEX_TVB5(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,const Eigen::VectorXd &x_old_old3 ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3,Eigen::VectorXd &K4 , const double  t , const double dt);
  void fill_lhs2_IMEX_TVB5(const double dt);
  void do_VSIMEX_BDF2(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1 , const double  t , const double dt);
  void fill_mat_VSIMEX_BDF2(const double dt, const double dt1);
  void fill_lhs1_IMEX_BDF2(const double dt, const double dt1);
  void fill_mat_VSIMEX_BDF4( const double dt, const double dt1, const double dt2, const double dt3);
  void fill_lhs1_VSIMEX_BDF4(const double dt, const double dt1, const double dt2, const double dt3);
  void do_VSIMEX_BDF4(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3 , const double  t , const double dt);
  void fill_lhs2_IMEX_BDF2(const double dt);
  void fill_mat_VSIMEX_BDF3( const double dt, const double dt1, const double dt2);
  void fill_lhs2_VSIMEX_BDF3(const double dt, const double dt1, const double dt2);
  void fill_lhs2_VSIMEX_BDF3_DD(double dt);
  void do_VSIMEX_BDF3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt);
   void do_VSIMEX_BDF3_DD(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt);
  void fill_lhs1_IMEX_ARS_232(const double dt);
  void fill_lhs2_IMEX_ARS_232(const double dt);
  void do_IMEX_ARS_232(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);

};

#endif // INTEGRATOR1_H
