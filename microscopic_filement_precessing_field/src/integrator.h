/*
 * Copyright 2021 <copyright holder> <email>
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
#include <boost/numeric/odeint.hpp>
#include <sys/time.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <numeric>      // std::iota
#include <algorithm>
#include <iostream>
using namespace std;


// std::sort, std::stable_sort

typedef std::vector<double> dvec;
using namespace boost::numeric::odeint;



template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}




/**
 * @todo write docs
 */
class Integrator
{
private:
  int N; //number of magnetic particles
  int Ns; //numbe of segments
  int Nv; // number of variable
  double beta=62.0/180*M_PI;
  double a;
  double spring_l;
  double spring_k = 10.0;
  double k_bend = 2.6;

  Eigen::Vector3d gravitation_force={0.0,0.0,-1.0};
  Eigen::Vector3d reaction_force_dir={0.0,0.0,1.0};
  double bottom_pozition=/*-0.472624*/-M_SQRT1_2 /*+0.1505*/;
  double sigma_prefactor=std::pow(2.0,-1.0/6.0);
  std::vector<Eigen::Vector3d> pos; //position of centers
  Eigen::Vector3d magnetisation; //magnetization vectors
  std::vector<Eigen::Matrix3d> rot; //rotation matrix for amm cubes
  std::vector<Eigen::Quaterniond> quat;
  std::vector<Eigen::Vector3d> mag_Force; //vector of magnetic forces
  std::vector<Eigen::Vector3d> mag_Torque; //vector of magnetic torques
  std::vector<Eigen::Vector3d> ext_Force;
  std::vector<Eigen::Vector3d> ext_Torque; //vector of external Torques
  std::vector<Eigen::Vector3d> steric_Force; //steric forces
  std::vector<Eigen::Vector3d> steric_Torque; //steric torque
  double B_mag; //magnituede af external field
  double omega_B; // angular velocity of external field
  //Eigen::Vector3d B;
  Eigen::MatrixXd Omega_quat=Eigen::MatrixXd::Zero(4, 4);
  void fill_Omega_quat(Eigen::Vector3d & omega_vec);
  double mom_angle;
  void fill_support_vectors(double t);
  void fill_mag_foreces();
  void fill_external_foreces_and_torques(const double t);
  void fill_external_foreces_and_torques_no_stiffnes(const double t);
  std::vector<size_t> sort_indexes(const dvec &v);
  void fill_steric_foreces();
  void fill_steric_torques();
  //VSIMEX
  Eigen::MatrixXd A;
  Eigen::SparseMatrix<double>A_sparse, lhs_sparse2,lhs_sparse;
  Eigen::SparseLU<Eigen::SparseMatrix<double> >  solver_sparse2, solver_sparse;
  double alpha0=0., alpha1=0., alpha2=0., alpha3=0., alpha4=1., beta0=0., beta1=0., beta2=0., beta3=0.,omega1=0.,omega2=0.,omega3=0.;


    /**
     * @todo write docs
     */
public:
    Integrator(int N_, double omega_B_, double beta_,  bool B_off);
    Integrator(const Integrator& old);
    void operator()(const Eigen::VectorXd & x, Eigen::VectorXd & dxdt, const double t);
    Eigen::Vector3d get_HD_force(int i){return mag_Force[i];};
    Eigen::Vector3d get_steric_force(int i){return steric_Force[i];};
    Eigen::Vector3d get_B(double t){  Eigen::Vector3d B={cos(beta), sin(beta)*sin(2*M_PI*omega_B*t), sin(beta)*cos(2*M_PI*omega_B*t)};return B;};
    Eigen::Vector3d get_external_torque(){return ext_Torque[0];};
    double get_omega(){return omega_B;};
    void set_omega( double a){ omega_B=a;};
    void append_data( Eigen::VectorXd &x, std::vector <dvec> &data,double t);
    Eigen::VectorXd  x_temp, K1, K2, K3, K4, rhs;
    void fill_mat_VSIMEX_BDF3( const double dt, const double dt1, const double dt2);
    void fill_lhs2_VSIMEX_BDF3(const double dt, const double dt1, const double dt2);
    void do_VSIMEX_BDF3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,  Eigen::VectorXd &new_x, const double  t , const double dt);
    void fill_K1_IMEX_BDF3(const Eigen::VectorXd &x,const double dt, const double t);
    void fill_lhs_Euler_impl_expl(const double dt);
    void do_Euler_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt);

};

#endif // INTEGRATOR_Hvoid
