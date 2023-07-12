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
#include <vector>
#include <chrono>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

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
  double a;
  double b;
  double alpha;
  double komega_kv;
  std::vector<Eigen::Vector3d> pos; //position of centers
  std::vector<Eigen::Vector3d> angle; //rotation angles
  std::vector<Eigen::Vector3d> HD_Force_and_Torque; //HD force and Torque
  std::vector<Eigen::Vector3d> init_magnetisation; //initial magnetization vectors
  std::vector<Eigen::Vector3d> magnetisation; //magnetization vectors
  std::vector<Eigen::Matrix3d> rot; //rotation matrix for amm cubes
  std::vector<Eigen::Vector3d> mag_Force; //vector of magnetic forces
  std::vector<Eigen::Vector3d> mag_Torque; //vector of magnetic torques
  std::vector<Eigen::Vector3d> ext_Torque; //vector of external Torques
  std::vector<Eigen::Vector3d> steric_Force; //steric forces
  std::vector<Eigen::Vector3d> steric_Force1; //steric forces
  std::vector<Eigen::Vector3d> steric_Torque; //vector of steric torques
  std::vector<Eigen::Vector3d> steric_Torque1; //vector of steric torques
  std::vector<Eigen::Vector3d> perp; // perpendicul to top side
  std::vector<Eigen::Vector3d> perp1; // perpendicul to top side
  std::vector<Eigen::Vector3d> TR; // position of the top right corner
  std::vector<Eigen::Vector3d> TL; // position of the top left corner
  std::vector<Eigen::Vector3d> BR; // position of the bottom right corner
  std::vector<Eigen::Vector3d> BL; // position of the bottom left corner
  Eigen::Vector3d i_perp;
  Eigen::Vector3d i_perp1;
  Eigen::Vector3d i_TR;
  Eigen::Vector3d i_TL;
  Eigen::Vector3d i_BR;
  Eigen::Vector3d i_BL;
  double B_mag; //magnituede af external field
  double omega_B; // angular velocity of external field
  Eigen::Vector3d B;
  double phi0;
  double mom_angle;
  dvec dist_inside;
  std::vector<Eigen::Vector3d> steric_force_dir;
  dvec torque_dist;
  
  void fill_support_vectors();
  void fill_mag_foreces_and_torques();
  //double ScalProd(Eigen::Vector3d  a, Eigen::Vector3d b){double res=a[0]*b[0]+a[1]*b[1]; return res;};
  //double crossProd(Eigen::Vector3d  a, Eigen::Vector3d  b){return ( a[0]*b[1] - a[1]*b[0]);};
  void fill_external_foreces_and_torques(const double t);
  double is_point_inside_below(Eigen::Vector3d &point, int i);
  double is_point_inside_above(Eigen::Vector3d &point, int i);
  double is_point_inside(Eigen::Vector3d &point,Eigen::Vector3d &pointp,Eigen::Vector3d &pointp1, int i, Eigen::Vector3d &force_dir, double &torque_distance);
  void fill_steric_foreces_wp();
  double is_point_inside_wp(Eigen::Vector3d &point,Eigen::Vector3d &pointp,Eigen::Vector3d &pointp1, int i, Eigen::Vector3d &force_dir);
  double is_point_inside_TB(Eigen::Vector3d &point, int i, Eigen::Vector3d &force_dir);
  void fill_steric_foreces_TB();
  void fill_steric_torques();
  void fill_steric_torques1();
  void fill_steric_foreces_and_torques();
  void fill_steric_foreces_and_torques_old();
  void fill_steric_foreces_and_torques_print();
  std::vector<size_t> sort_indexes(const dvec &v);
  void fill_steric_foreces_and_torques_new();
public:
    /**
     * @todo write docs
     */
    Integrator(int N_, double a_, std::vector<Eigen::Vector3d> &pos_, std::vector<Eigen::Vector3d> angle_, double B_mag_, double omega_B_, Eigen::Vector3d& init_magnetisation_);
    Integrator(int N_, double omega_B_, Eigen::Vector3d& init_magnetisation_, double alpha_, double komega_kv_, double lag, bool B_off);
    //Integrator():N(2){std::cout<<"This newer should be used!!!!!!!"<<std::endl;abort();};
    void operator()(const Eigen::VectorXd & x, Eigen::VectorXd & dxdt, const double t);
    void fill_steric_foreces();
    Eigen::Vector3d get_HD_force(int i){return mag_Force[i];};
    Eigen::Vector3d get_steric_force(int i){return steric_Force[i];};
    double get_Bx(double t){return B_mag*std::cos(t*omega_B+phi0);};
    Eigen::Vector3d get_Mag_torque(Eigen::VectorXd &x){  for(size_t i=0; i<N; i++)
  {
    pos[i][0]=x[12*i];
    pos[i][1]=x[12*i+1];
    pos[i][2]=x[12*i+2];
    angle[i][0]=x[12*i+3];
    angle[i][1]=x[12*i+4];
    angle[i][2]=x[12*i+5];
  }
   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
   fill_support_vectors();
   //std::cout <<"here3"<<std::endl;
   fill_mag_foreces_and_torques();
//    std::cout <<"here4"<<std::endl;
     fill_steric_foreces_and_torques_new(); return mag_Torque[0];};
    Eigen::Vector3d get_Steric_torque(){return steric_Torque[0];};
    Eigen::Vector3d get_external_torque(){return ext_Torque[0];};
    double get_By(double t){return B_mag*std::sin(t*omega_B+phi0);};
    double get_phi(double t){return t*omega_B+phi0;};
    double get_phi0(){return phi0;};
    double get_mom_angle(){return mom_angle;};
    void print_forces(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t);
    void print_mag_foreces_and_torques();
    double get_omega(){return omega_B;}

};

#endif // INTEGRATOR_H
