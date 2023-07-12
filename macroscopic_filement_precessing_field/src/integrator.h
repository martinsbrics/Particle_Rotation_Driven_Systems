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

#ifndef INTEGRATOR_H
#define INTEGRATOR_H
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
class Integrator
{
  size_t Ns;
  size_t Np;
  size_t Nv;
  double alfa;
  double mu;//=p+1;
  double Cm;//=100;
  double H1;//=4.0;
  double om;//=50*16
  Eigen::VectorXd J_Jt_u,J_Jt_l,J_Jt_d, E_small_vec, rhs, ones;
  Eigen::MatrixXd A, E, E_small, P, Jt_J_inv, lhs;
  Eigen::VectorXd FA, A_r, J_r, temp_r, x_temp;
  Eigen::SparseMatrix<double> A_sparse ,Jt, J, lhs_sparse;
  void fill_Fa(const Eigen::VectorXd &x, double t);
  void fill_J(const Eigen::VectorXd &x);
  double dist(const  Eigen::VectorXd &x , int i);
  void renorm(const  Eigen::VectorXd &x1, Eigen::VectorXd &x2);
  
  long dgtsv(long Nx1, long NRHS, double *DL, double *D, double *DU, double *imag_t_Bx11, long LDB);

  
public:
  Integrator(size_t Ns1, double Cm1, double H11, double om1);
  void operator() (const Eigen::VectorXd &x ,  Eigen::VectorXd &dxdt , const double  t );
  \

};

#endif // INTEGRATOR_H
