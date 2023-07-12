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
  int N_virtual=92/*8*/;// alternative is 44,20,8
  double a;
  double b;
  double alpha;
  double komega_kv;
  Eigen::Vector3d gravitation_force={0.0,0.0,-1.0};
  Eigen::Vector3d reaction_force_dir={0.0,0.0,1.0};
  double bottom_pozition=/*-0.472624*/-M_SQRT1_2 /*+0.1505*/;
  double sigma_prefactor=std::pow(2.0,-1.0/6.0);
  std::vector<Eigen::Vector3d> pos; //position of centers
  std::vector<Eigen::Vector3d> angle; //rotation angles
  std::vector<Eigen::Vector3d> HD_Force_and_Torque; //HD force and Torque
  std::vector<Eigen::Vector3d> init_magnetisation; //initial magnetization vectors
  std::vector<Eigen::Vector3d> magnetisation; //magnetization vectors
  std::vector<Eigen::Matrix3d> rot; //rotation matrix for amm cubes
  std::vector<Eigen::Quaterniond> quat;
  std::vector<Eigen::Vector3d> mag_Force; //vector of magnetic forces
  std::vector<Eigen::Vector3d> mag_Torque; //vector of magnetic torques
  std::vector<Eigen::Vector3d> ext_Torque; //vector of external Torques
  std::vector<Eigen::Vector3d> steric_Force; //steric forces
  std::vector<Eigen::Vector3d> steric_Torque; //steric torque
  double B_mag; //magnituede af external field
  double omega_B; // angular velocity of external field
  Eigen::Vector3d B;
  double phi0;
  double gravitation_perfactor=1.5;
  std::vector<Eigen::Vector3d> B_at_partivcle_pos;
  std::vector<std::vector<Eigen::Vector3d> > virtual_sites;
  std::vector<std::vector<Eigen::Vector3d> > new_virtual_sites;
  void update_virtual_pos();
  std::vector<double>  virtual_sites_radius;
  double radius_big_particle=0.0;
  Eigen::MatrixXd Omega_quat=Eigen::MatrixXd::Zero(4, 4);
  void fill_Omega_quat(Eigen::Vector3d & omega_vec);
  double mom_angle;
  dvec dist_inside;
  void fill_support_vectors();
  void fill_mag_foreces_and_torques();
  void fill_external_foreces_and_torques(const double t);
  void calc_B_field();
  void fill_virtual_sites();
  std::vector<size_t> sort_indexes(const dvec &v);
  void fill_steric_foreces();
  void fill_steric_foreces_print();
  void fill_steric_torques();
  void fill_steric_torques_cout();
    /**
     * @todo write docs
     */
public:
    Integrator(int N_, double omega_B_, Eigen::Vector3d& init_magnetisation_, double alpha_, double komega_kv_, double lag,double gravitation_perfactor1, double angle, int N_virtual1, bool B_off);
    void operator()(const Eigen::VectorXd & x, Eigen::VectorXd & dxdt, const double t);
    void debug(const Eigen::VectorXd & x, const double t);
    Eigen::Vector3d get_HD_force(int i){return mag_Force[i];};
    Eigen::Vector3d get_steric_force(int i){return steric_Force[i];};
    double get_Bx(double t){return B_mag*std::cos(t*omega_B+phi0);};
    void print_mathematica(const Eigen::VectorXd & x,std::string folder);
    void test_do_something(){cout<<"Number of virtual particles "<<N_virtual<<" "<<pos[0][0]<<" "<<pos[0][1]<<endl;};
    Eigen::Vector3d get_init_mom (size_t i){return init_magnetisation[i];};

   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;

   //std::cout <<"here3"<<std::endl;

//    std::cout <<"here4"<<std::endl;
    Eigen::Vector3d get_external_torque(){return ext_Torque[0];};
    double get_omega(){return omega_B;};
    void set_omega( double a){ omega_B=a;};
    void append_data_virtual( Eigen::VectorXd &x, std::vector <dvec> &data,double t);
    void append_data( Eigen::VectorXd &x, std::vector <dvec> &data,double t);
    void append_data_quat( Eigen::VectorXd &x, std::vector <dvec> &data, double t);
};

#endif // INTEGRATOR_Hvoid
