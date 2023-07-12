/*
 * Copyright 2018 <copyright holder> <email>
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

#ifndef CLAC_H
#define CLAC_H

using namespace std;
#include <iostream>
#include <fstream>
#include <sstream> 

#include "optim.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>


//#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>
#include <Eigen/Geometry>

/**
 * @todo write docs
 */
class calc
{
private:
    int n;double B_mag;
public:
    /**
     * Default constructor
     */
    calc(int n1, double B_mag1, double phi);
    static double ackley_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_shift3(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_smart(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_smart_smart(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_smart_smart_second(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_with_grad1(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_with_grad(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
    static double fn_gsl(const gsl_vector * x, void * opt_data);
    static double fn_gsl_new(const gsl_vector * x, void * opt_data);
    static double fn_smart_smart_gsl(const gsl_vector * x, void* opt_data);
    static double fn_smart_smart_gsl_second(const gsl_vector * x, void* opt_data);
    void set_val(const arma::vec& vals_inp);
    void set_val_gsl(const gsl_vector * vals_inp);
    void set_val_shift3(const arma::vec& vals_inp);
    void fill_vec_x(const arma::vec& vals_inp);
    void set_val_smart(const arma::vec& vals_inp);
    void set_val_smart_smart(const arma::vec& vals_inp);
    void set_val_smart_smart_second(const arma::vec& vals_inp);
    void set_val_w_grad(const arma::vec& vals_inp);
    double calc_en();
    double calc_en_w_grad(arma::vec* grad_out);
    double dipole_energy(int i, int j);
    double dipole_energy_discret_cube(int i, int j,int ncube, int i1, int i2, int i3, int j1, int j2, int j3);
    double dipole_energy_w_grad(int i, int j, double & grad_ix, double & grad_iy, double & grad_jx, double & grad_jy, double & grad_iphi, double & grad_jphi);
    double dipole_energy_w_grad(int i, int j, double & grad_jx, double & grad_jy, double & grad_jphi);
    string get_data();
    double get_offset();
    double get_maxoffset();
    double get_angle();
    double get_rangle();
    void minimize_energy();
    void minimize_energy_first();
    void minimize_energy_second();
    double calc_en_ncube(int ncube);
    double orintation();
    void shift(double a);
    double xx;
    void minimize_energy_pso();
    void minimize_energy_pso_smart();
    void minimize_energy_de_smart();
    void improve_energy_first();
    void improve_energy_first_();
    void improve_energy_first_cout();
    double calc_en_w_grad1(arma::vec* grad_out);
    void minimize_energy_de_smart_smart(bool first=true);
    void minimize_energy_smart_smart_gsl( bool first=true);
    void improve_energy_smart_smart_gsl( bool first=true);
    void improve_energy_second();
    void improve_energy_second_();
    void improve_energy_second_cout();
    void minimize_energy_shit3(int kink1=1);
    void setB(double B_mag1);
    void minimize_energy_gsl(bool first1=true);
    void minimize_energy_gsl_new();
    void set_val_smart_smart_gsl(const gsl_vector * vals_inp);
    void improve_energy_gsl();
    double get_Bangle();
    void set_val_smart_smart_second_gsl(const gsl_vector * vals_inp);
    double get_beta();
    void set_val_gsl_new(const gsl_vector * vals_inp);
    void set_val_gsl_new_kink(const gsl_vector * vals_inp);
    static double fn_gsl_new_kink(const gsl_vector * x, void * opt_data);
    void minimize_energy_gsl_new_kink(int kink1=1);
    double Heviside(double x);
    void minimize_energy_gsl_new_kink_smart();
    double get_beta1();
    double get_beta2();
    void set_val_gsl_smart(const gsl_vector * vals_inp);
    void minimize_energy_gsl_smart(bool first1);
    static double fn_gsl_smart(const gsl_vector * x, void * opt_data);
    void set_val_gsl_smart_kink(const gsl_vector * vals_inp);
    static double fn_gsl_smart_kink(const gsl_vector * x, void * opt_data);
    void minimize_energy_gsl_smart_kink(int kink1);
    void set_val_gsl_superball(const gsl_vector * vals_inp);
    static double fn_gsl_superball(const gsl_vector * x, void * opt_data);
    void minimize_energy_gsl_superball( double q1);
    double q_par;
    vector <double> get_betas();
    double par_c_div_b;
    void set_val_gsl_cube(const gsl_vector * vals_inp);
    static double fn_gsl_cube(const gsl_vector * x, void * opt_data);
    void minimize_energy_gsl_cube( double q1);
    
    //
double theta, phi,pfi;
int rot_point, add_case; bool first; bool do_not_calc;
vector <Eigen::Vector3d> m,Rotoation_point, PRotoation_point; Eigen::Vector3d Xold,Yold,Zold,X,Y,Z,M, Centr;
vector <Eigen::Vector3d> m_grad;
vector<Eigen::Vector3d> r;
vector<Eigen::Matrix3d> rot;
Eigen::Matrix3d rotyz;
vector<Eigen::Matrix3d> rot_grad;
Eigen::Vector3d B;
Eigen::Vector3d r10;
Eigen::Vector3d tempr;
Eigen::Vector3d B_grad1;
Eigen::Vector3d B_grad2;
double percantage_sio2;
int kink;
vector < double> vec_x, vec_smart_smart_x;



};

#endif // CLAC_H
