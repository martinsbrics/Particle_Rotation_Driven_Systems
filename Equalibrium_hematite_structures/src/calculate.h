/*
 * Copyright 2018 Martins Brics <email>
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

#ifndef CALCULATE_H
#define CALCULATE_H

using namespace std;
#include <iostream>
#include <fstream>

#include "optim.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>

//#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>
#include <Eigen/Geometry>



/**
 * @todo write docs
 */
class calculate
{
private:
    static double aa;
public:
 calculate(double B_mag1);
 void print();
 static double ackley_fn1(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
  static double ackley_fn2(const arma::vec& vals_inp, void* opt_data);
 static double fn_with_grad(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);
 void angle_scan();
 
};

#endif // CALCULATE_H
