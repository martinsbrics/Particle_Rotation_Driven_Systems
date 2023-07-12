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

#ifndef SIMULATION_H
#define SIMULATION_H

#include <boost/numeric/odeint.hpp>
#include "eigen_algebra.hpp"
#include "integrator.h"
#include "integrator1.h"
using namespace boost::numeric::odeint;

namespace boost {
namespace numeric {
namespace odeint {

template<typename B,int S1,int S2,int O, int M1, int M2>
struct algebra_dispatcher< Eigen::Matrix<B,S1,S2,O,M1,M2> >
{
    typedef vector_space_algebra algebra_type;
};

}}}



/**
 * @todo write docs
 */
class simulation
{
  size_t Ns;
  size_t Np;
  size_t Nv;
  double alfa;
  double mu;//=p+1;
  double Cm;//=100;
  double H1;//=4.0;
  double om;//=50*16
  double torq;
  // Integrator class instance
  Integrator *sys;Integrator1 *sys1;
  // time
  double t;
  //timestep
  double dt;
  // max to propagate
  double t_max;
  // time for analysis
  double extra_t;
  // iter per time step
  int iter_per_step;
//   void print_data(double step=0.1, double tmax=20);
//   void print_data_sigle_file( double step=0.1, double tmax=20);
//   void print_data(string &str);
//   //void print_data_(int i);
  void print_pars();
  std::chrono::time_point<std::chrono::system_clock> now;
  std::chrono::time_point<std::chrono::system_clock> foo;
  std::chrono::duration<double> diff1;
  string folder;
  double dist(const  Eigen::VectorXd &x , int i);
  vector<Eigen::MatrixXd> data; vector<Eigen::VectorXd> data1;
  dvec vec_t;
  void renorm(const  Eigen::VectorXd &x1, Eigen::VectorXd &x2);
  void print_data( double step, double tmax1);
  void print_data1( double step, double tmax1);
  dvec length;
  void calc_dist(const  Eigen::VectorXd &x);
  void anal_data();
  bool accept_step(double err=1.e-10);
  bool accept_step_beginning();
  double max_rel_error;
  int npart;
  void fill_0mega3();
//   double phi6(dvec & xlocal);
//   void do_delaunator( dvec &x);
//   vector <ivec> neighbours;
//   ivec neighbour_count;
//   dcmplx fi6(size_t point, dvec &xlocal);
//   double phi6_local( dvec & xlocal);
//   void print_data_w_phi6( double step=0.1, double tmax1=20);
//   void print_angle( double step=0.1, double tmax1=20);
  
public:
  simulation(size_t Ns1, double Cm1, double H11, double om1, double torq1, Eigen::VectorXd x_);
  simulation(size_t Ns1, double Cm1, double H11, double om1, double torq1);
  ~simulation();
   void init();
   Eigen::VectorXd x, x_temp, Omega3;
//   void rand_init();
//   void read(const char* fname);
  void propagate(double dt,  double tmax, double abs_tol = 1.0e-13,double rel_tol = 1.0e-13);
  void propagate_CN(double dt1, double tmax1);
  void propagate_Euler(double dt1, double tmax1);
  void propagate_CN_proper(double dt1, double tmax1);
  void propagate_Euler_impl_expl(double dt1, double tmax1);
  void propagate_Euler1_impl_expl(double dt1, double tmax1);
  void propagate_midp_impl_expl(double dt1, double tmax1);
  void propagate_rk4_impl_expl(double dt1, double tmax1);
  void propagate_IMEX_BDF2(double dt1, double tmax1);
  void propagate_IMEX_BDF3(double dt1, double tmax1);
  void propagate_IMEX_TVB3(double dt1, double tmax1);
  void propagate_IMEX_BDF4(double dt1, double tmax1);
  void propagate_IMEX_BDF5(double dt1, double tmax1);
  void propagate_IMEX_TVB5(double dt1, double tmax1);
  void propagate_IMEX_TVB5_adaptive(double dt1, double tmax1);
  void propagate_IMEX_BDF5_adaptive(double dt1, double tmax1);
  void propagate_IMEX_BDF4_adaptive(double dt1, double tmax1);
  void propagate_IMEX_BDF3_adaptive(double dt1, double tmax1);
  void propagate_VSIMEX_BDF2(double dt1, double tmax1, double tol=1e-6);
  void propagate_VSIMEX_BDF3(double dt1, double tmax1, double tol=1e-7);
  void propagate_VSIMEX_BDF3_DD(double dt1, double tmax1, double tol=1e-7);
  void propagate_VSIMEX_BDF4(double dt1, double tmax1);
  void propagate_VSIMEX_BDF4_adaptive(double dt1, double tmax1);
  void propagate_ARS_232(double dt1, double tmax1);
  string anal_data1();
  Eigen::VectorXd get_x(){return x;};
//   void propagate1(double dt, double tmax, double tmax1, double abs_tol = 1.0e-13,double rel_tol = 1.0e-13);
//   void print_phi6( double step, double tmax1);
};


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T>
void sort_indexes(const vector<T> &v,vector<size_t> &idx) {

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
}



#endif // SIMULATION_H
