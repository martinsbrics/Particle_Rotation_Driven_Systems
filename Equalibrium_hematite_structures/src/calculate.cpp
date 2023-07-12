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

#include "calculate.h"

//Eigen::Vector3d m0(-0.5,0,0.5);
//Eigen::Vector3d m0(1,1,1);
Eigen::Vector3d m0(1,1,M_SQRT2*tan(47.0*M_PI/180));
Eigen::Vector3d m1;
Eigen::Vector3d Bx(1,0,0);
Eigen::Vector3d B=Bx;
Eigen::Matrix3d rot=Eigen::Matrix3d::Identity();
Eigen::Matrix3d rotBx=Eigen::Matrix3d::Identity();
Eigen::Matrix3d rotBy=Eigen::Matrix3d::Identity();
Eigen::Matrix3d rotBz=Eigen::Matrix3d::Identity();
Eigen::Matrix3d rot_grad=Eigen::Matrix3d::Zero();
Eigen::Vector3d r0;
Eigen::Vector3d r1;
Eigen::Vector3d r10;
double B_mag;


double calculate::ackley_fn1(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    r1[0]=vals_inp(0);
    r1[1]=vals_inp(1);
    rot(0,0)=cos(vals_inp(2));
    rot(0,1)=-sin(vals_inp(2));
    rot(1,0)=sin(vals_inp(2));
    rot(1,1)=cos(vals_inp(2));
    rotBz(0,0)=cos(vals_inp(5));
    rotBz(0,1)=-sin(vals_inp(5));
    rotBz(1,0)=sin(vals_inp(5));
    rotBz(1,1)=cos(vals_inp(5));
    rotBx(1,1)=cos(vals_inp(3));
    rotBx(1,2)=-sin(vals_inp(3));
    rotBx(2,1)=sin(vals_inp(3));
    rotBx(2,2)=cos(vals_inp(3));
    rotBy(0,0)=cos(vals_inp(4));
    rotBy(0,2)=sin(vals_inp(4));
    rotBy(2,0)=-sin(vals_inp(4));
    rotBy(2,2)=cos(vals_inp(4));
    B=(rotBx*rotBy*rotBz)*Bx;  
    m1=rot*m0;
    r10=r1-r0;
    double dr10=r10.norm();
    double dr10_2=dr10*dr10;
    double dr10_3=dr10_2*dr10;
    double dr10_5=dr10_3*dr10_2;
//     cout <<"Here"<<endl;
//     cout <<r1<<" " << dr10<<endl;;
//     abort();
    
    double obj_val = (m0.dot(rot*m0)/dr10_3-3*m0.dot(r10)*(rot*m0).dot(r10)/dr10_5);
    obj_val-=B_mag*(B.dot(m0)+B.dot(m1));

    //

    return obj_val;
}

double calculate::ackley_fn2(const arma::vec& vals_inp, void* opt_data)
{
    r1[0]=vals_inp(0);
    r1[1]=vals_inp(1);
    rot(0,0)=cos(vals_inp(2));
    rot(0,1)=-sin(vals_inp(2));
    rot(1,0)=sin(vals_inp(2));
    rot(1,1)=cos(vals_inp(2));
    m1=rot*m0;
    r10=r1-r0;
    double dr10=r10.norm();
    double dr10_2=dr10*dr10;
    double dr10_3=dr10_2*dr10;
    double dr10_5=dr10_3*dr10_2;
//     cout <<"Here"<<endl;
//     cout <<r1<<" " << dr10<<endl;;
//     abort();
    
    double obj_val = m0.dot(m1)/dr10_3-3*m0.dot(r10)*(m1).dot(r10)/dr10_5;

    //

    return obj_val;
}


double calculate::fn_with_grad(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    r1[0]=vals_inp(0);
    r1[1]=vals_inp(1);
    rot(0,0)=cos(vals_inp(2));
    rot(0,1)=-sin(vals_inp(2));
    rot(1,0)=sin(vals_inp(2));
    rot(1,1)=cos(vals_inp(2));
    rot_grad(0,0)=-sin(vals_inp(2));
    rot_grad(0,1)=-cos(vals_inp(2));
    rot_grad(1,0)=cos(vals_inp(2));
    rot_grad(1,1)=-sin(vals_inp(2));
    r10=r1-r0;
    double dr10=r10.norm();
    double dr10_2=dr10*dr10;
    double dr10_3=dr10_2*dr10;
    double dr10_5=dr10_3*dr10_2;
    double dr10_7=dr10_5*dr10_2;
    double obj_val = m0.dot(rot*m0)/dr10_3-3*m0.dot(r10)*(rot*m0).dot(r10)/dr10_5;
    
 
    //
    if (grad_out) {
        (*grad_out)(0) = -3*m0.dot(rot*m0)*r10(0)/dr10_5 +15*r10(0)*m0.dot(r10)*(rot*m0).dot(r10)/dr10_7
        -3*m0(0)*(rot*m0).dot(r10)/dr10_5-3*m0.dot(r10)*(rot*m0)(0)/dr10_5;
        (*grad_out)(1) = -3*m0.dot(rot*m0)*r10(1)/dr10_5 +15*r10(1)*m0.dot(r10)*(rot*m0).dot(r10)/dr10_7
        -3*m0(1)*(rot*m0).dot(r10)/dr10_5-3*m0.dot(r10)*(rot*m0)(1)/dr10_5;
        (*grad_out)(2)= m0.dot(rot_grad*m0)/dr10_3-3*m0.dot(r10)*(rot_grad*m0).dot(r10)/dr10_5;
    }
    //
    cout<<r1.transpose()<<" "<<obj_val<<" "<<(*grad_out)(0)<<" "<<(*grad_out)(1)<<" "<<(*grad_out)(2)<<endl;
    return obj_val;
}


calculate::calculate(double B_mag1)
{   
  cout.precision(16);
  double phi0=atan(1.0/sqrt(2.))*180/M_PI;
  m0=Eigen::Vector3d (1,1,M_SQRT2*tan((12+phi0)*M_PI/180));
    m0.normalize();
    r0=Eigen::Vector3d(0,0,0);
    //r1=Eigen::Vector3d(0.5,0.5,1);
    r1=Eigen::Vector3d(0,0,1);
    B_mag=B_mag1;
    cout <<r0<<endl;
    //a=23.718282L;
    
    cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(6,1) + 1.0; // (2,2)
    for(int i=0; i<2; i++)
    {
        x(i)=r1[i];
    }
    for(int i=2; i<6; i++)
    {
        x(i)=0.0;
    }
    cout <<x<<endl;

    //

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    //bool success = optim::de(x,ackley_fn1,nullptr);
    //bool success = optim::pso(x,ackley_fn1,nullptr);
    bool success =optim::bfgs(x,fn_with_grad,nullptr);

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    x(2)=fmod(x(2), 2*pi);
    arma::cout << "\n1de: solution to Ackley test:\n" << x << arma::endl;
    std::cout <<ackley_fn1(x, nullptr, nullptr) <<endl;
    std::cout <<r0.transpose()<<" "<<r1.transpose()<<" "<<m0.transpose()<<" "<<m1.transpose()<<" "<<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;
    
//     start = std::chrono::system_clock::now();
//      for(int i=0; i<2; i++)
//     {
//         x(i)=r1[i];
//     }
//     //x(2)=M_PI;
//     x(2)=0;
//     //bool success_2 = optim::cg(x,fn_with_grad,nullptr);
//     //bool success_2 = optim::bfgs(x,fn_with_grad,nullptr);
//     //bool success_2 = optim::gd(x,fn_with_grad,nullptr);
//     bool success_2 = optim::pso(x,ackley_fn1,nullptr);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end-start;
//     if (success_2) {
//         std::cout << "cg: Booth test completed successfully\n" 
//         << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     } else {
//         std::cout << "cg: Booth test completed unsuccessfully." << std::endl;
//     }
//  
//     arma::cout << "cg: solution to Booth test:\n" << x << arma::endl;
//     
//     cout <<m0.transpose()<<" "<<(rot*m0).transpose()<<endl;
// //     x(0)=-0.5;x(1)=0;x(2)=0;
// //     for(int i=0; i<200; i++)
// //     {
// //         cout << x(0)<<" "<< ackley_fn1(x, nullptr, nullptr) <<endl;
// //         x(0)+=0.005;
// //     }
//     
//  
    
}
void calculate::print()
{
    cout <<" I am here"<<endl;
}

void calculate::angle_scan()
{
//     double phi0=atan(1/sqrt(2))*180/M_PI;
    double phi0=0;
    ofstream file; file.open("energy_angle");
    arma::vec x = arma::ones(3,1) + 1.0; // (2,2)
    for (int i=0; i<90; i++)
    {
        m0=Eigen::Vector3d (1,1,M_SQRT2*tan((phi0+i)*M_PI/180));
        m0.normalize();
        for(int j=0; j<2; j++)
        {
            x(j)=0.5;
        }
        //x(2)=M_PI;
        x(2)=0;
        bool success_2 = optim::pso(x,ackley_fn1,nullptr);
        double ener=ackley_fn1(x, nullptr, nullptr);
        cout <<"For angle "<< i<<" E="<< ener<<endl;
        file<<i<<" "<<ener<<" "<<x(0)<<" "<<x(1)<<" "<<x(2)<<" "<<m0.transpose()<<" "<<(rot*m0).transpose()<<endl;
    }
    m0=Eigen::Vector3d (0,0,1); 
    for(int i=0; i<2; i++)
    {
        x(i)=0.5;
    }
    //x(2)=M_PI;
    x(2)=0;
    bool success_2 = optim::pso(x,ackley_fn1,nullptr);
    if(!success_2)
    {
        std::cout << "de: Ackley test completed unsuccessfully." << std::endl;
    }
    double ener=ackley_fn1(x, nullptr, nullptr);
    cout <<"For angle "<< 90<<" E="<< ener;
    file<<90<<" "<<ener<<" "<<x(0)<<" "<<x(1)<<" "<<x(2)<<" "<<m0.transpose()<<" "<<(rot*m0).transpose()<<endl;
    for (int i=1; i<90; i++)
    {
        m0=Eigen::Vector3d (-1,-1,M_SQRT2*tan((90+i)*M_PI/180));
        m0.normalize();
        for(int j=0; j<2; j++)
        {
            x(j)=0.5;
        }
        //x(2)=M_PI;
        x(2)=0;
        bool success_2 = optim::pso(x,ackley_fn1,nullptr);
        double ener=ackley_fn1(x, nullptr, nullptr);
        cout <<"For angle "<< 90+i<<" E="<< ener<<endl;;
        file<<90+i<<" "<<ener<<" "<<x(0)<<" "<<x(1)<<" "<<x(2)<<" "<<m0.transpose()<<" "<<(rot*m0).transpose()<<endl;
    }
    m0=Eigen::Vector3d (1,0,1);
    m0.normalize();
    for(int j=0; j<2; j++)
    {
        x(j)=0.5;
    }
    //x(2)=M_PI;
    x(2)=0;
    success_2 = optim::pso(x,ackley_fn1,nullptr);
    ener=ackley_fn1(x, nullptr, nullptr);
    cout <<"For angle "<< x<<" E="<< ener<<endl;;
    file<<x<<" "<<ener<<" "<<x(0)<<" "<<x(1)<<" "<<x(2)<<" "<<m0.transpose()<<" "<<(rot*m0).transpose()<<endl;
    file.close();
}

