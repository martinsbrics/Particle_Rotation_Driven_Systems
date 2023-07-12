// SPDX-FileCopyrightText: 2023 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#ifndef SIMULATION_H
#define SIMULATION_H

#include <boost/numeric/odeint.hpp>
#include "eigen_algebra.hpp"
#include "integrator.h"
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
  size_t Np;
  size_t Ns;
  size_t Nv;
  double omega_B;
  double beta;
  bool B_off;
  Integrator *sys;
  // time
  double t;
  //timestep
  double dt;
  // max to propagate
  double t_max;
  std::chrono::time_point<std::chrono::system_clock> now;
  std::chrono::time_point<std::chrono::system_clock> foo;
  std::chrono::duration<double> diff1;
  string folder;
  double dist(const  Eigen::VectorXd &x , int i);
  vector<Eigen::MatrixXd> data;
  dvec vec_t;
  void renorm(const  Eigen::VectorXd &x1, Eigen::VectorXd &x2);
  void print_data( double step, double tmax1);
  dvec length;
  void calc_dist(const  Eigen::VectorXd &x);
  void anal_data();
  bool accept_step(double err=1.e-10);
  bool accept_step_beginning();
  double max_rel_error;
  void print_pars();
public:
  simulation(int N_, double omega_B_, double beta_,  bool B_off, Eigen::VectorXd x_);
  ~simulation();
   void init();
   Eigen::VectorXd x, x_temp;
   void propagate_VSIMEX_BDF3_DD(double dt1, double tmax1, double tol=1e-7);



};

#endif // SIMULATION_H
