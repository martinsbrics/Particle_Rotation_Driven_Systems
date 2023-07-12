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

#include "integrator.h"

Integrator::Integrator(int N_, double omega_B_, double beta_,  bool B_off):N(N_),Ns(N_-1), Nv(3*N_),
a(1.0),
spring_l(1.05*a),
pos(N),
mag_Force(N),
mag_Torque(N),
ext_Force(N),
ext_Torque(N),
steric_Force(N),
steric_Torque(N),
B_mag(B_off==true?1.0:0.0),
omega_B(omega_B_),
beta(beta_),
K1(Nv)
{
  std::cout<<"initalizeing"<<std::endl;
  
  magnetisation={cos(beta), sin(beta)*sin(2*M_PI*omega_B*0.0), sin(beta)*cos(2*M_PI*omega_B*0.0)};
  double kbend_hetero;
  double alfa=1.0/(3*M_PI*8.9e-04);
  A= Eigen::MatrixXd::Zero(Nv, Nv);
  for(size_t i=0; i<3; i++)
  {
    kbend_hetero=(0.6 + 0.4*(1.0*0/Ns));
    A(i,i)=-alfa*kbend_hetero;
    A(i,i+3)+=2*alfa*kbend_hetero;
    A(i,i+6)=-alfa*kbend_hetero;
    kbend_hetero=1.0;
    A(Nv-1-i,Nv-1-i)=-alfa*kbend_hetero;
    A(Nv-1-i,Nv-4-i)=2*alfa*kbend_hetero;
    A(Nv-1-i,Nv-7-i)=-alfa*kbend_hetero;
  }

  for(size_t i=3; i<6; i++)
  {
    kbend_hetero=(0.6 + 0.4*(1.0*1/Ns));
    A(i,i-3)=2*alfa*kbend_hetero;
    A(i,i)=-5*alfa*kbend_hetero;
    A(i,i+3)=4*alfa*kbend_hetero;
    A(i,i+6)=-alfa*kbend_hetero;
    kbend_hetero=(0.6 + 0.4*(1.0*(Ns-1)/Ns));
    A(Nv-1-i,Nv-1-i+3)=2*alfa*kbend_hetero;
    A(Nv-1-i,Nv-1-i)=-5*alfa*kbend_hetero;
    A(Nv-1-i,Nv-4-i)=4*alfa*kbend_hetero;
    A(Nv-1-i,Nv-7-i)=-alfa*kbend_hetero;
  }

  for(int i=2; i<Ns-1; i++)
  {
    kbend_hetero=(0.6 + 0.4*(1.0*i/Ns));
    A(3*i,3*(i+2))=-alfa*kbend_hetero;
    A(3*i+1,3*(i+2)+1)=-alfa*kbend_hetero;
    A(3*i+2,3*(i+2)+2)=-alfa*kbend_hetero;
    A(3*i,3*(i-2))=-alfa*kbend_hetero;
    A(3*i+1,3*(i-2)+1)=-alfa*kbend_hetero;
    A(3*i+2,3*(i-2)+2)=-alfa*kbend_hetero;
    A(3*i,3*(i+1))=4*alfa*kbend_hetero;
    A(3*i+1,3*(i+1)+1)=4*alfa*kbend_hetero;
    A(3*i+2,3*(i+1)+2)=4*alfa*kbend_hetero;
    A(3*i,3*(i-1))=4*alfa*kbend_hetero;
    A(3*i+1,3*(i-1)+1)=4*alfa*kbend_hetero;
    A(3*i+2,3*(i-1)+2)=4*alfa*kbend_hetero;
    A(3*i,3*(i))=-6*alfa*kbend_hetero;
    A(3*i+1,3*(i)+1)=-6*alfa*kbend_hetero;
    A(3*i+2,3*(i)+2)=-6*alfa*kbend_hetero;
  }
  A_sparse=A.sparseView();
  A_sparse.makeCompressed();

}

Integrator::Integrator(const Integrator& old)
{
  N=old.N;
  Ns=old.Ns;
  Nv=old.Nv;
  beta=old.beta;
  a=old.a;
  spring_l=old.spring_l;
  spring_k=old.spring_k;
  k_bend=old.k_bend;
  pos=old.pos;
  mag_Force=old.mag_Force;
  mag_Torque=old.mag_Torque;
  ext_Force=old.ext_Force;
  ext_Torque=old.ext_Torque;
  steric_Force=old.steric_Force;
  steric_Torque=old.steric_Torque;
  B_mag=old.B_mag;
  omega_B=old.omega_B;
  mom_angle=old.mom_angle;
}


void Integrator::fill_Omega_quat(Eigen::Vector3d& omega_vec)
{
  //the matrix has to be scew symmetric thus it is assumed that diagonals are always zero
//   Omega_quat(0,1)=0.5*omega_vec[2];
//   Omega_quat(0,2)=-0.5*omega_vec[0];
//   Omega_quat(0,3)=-0.5*omega_vec[1];
//   Omega_quat(1,0)=-0.5*omega_vec[2];
//   Omega_quat(1,2)=-0.5*omega_vec[1];
//   Omega_quat(1,3)=0.5*omega_vec[0];
//   Omega_quat(2,0)=0.5*omega_vec[0];
//   Omega_quat(2,1)=0.5*omega_vec[1];
//   Omega_quat(2,3)=0.5*omega_vec[2];
//   Omega_quat(3,0)=0.5*omega_vec[1];
//   Omega_quat(3,1)=-0.5*omega_vec[0];
//   Omega_quat(3,2)=-0.5*omega_vec[2];
  
  Omega_quat(0,1)=-0.5*omega_vec[2];
  Omega_quat(0,2)=0.5*omega_vec[1];
  Omega_quat(0,3)=0.5*omega_vec[0];
  Omega_quat(1,0)=0.5*omega_vec[2];
  Omega_quat(1,2)=-0.5*omega_vec[0];
  Omega_quat(1,3)=0.5*omega_vec[1];
  Omega_quat(2,0)=-0.5*omega_vec[1];
  Omega_quat(2,1)=0.5*omega_vec[0];
  Omega_quat(2,3)=0.5*omega_vec[2];
  Omega_quat(3,0)=-0.5*omega_vec[0];
  Omega_quat(3,1)=-0.5*omega_vec[1];
  Omega_quat(3,2)=-0.5*omega_vec[2];


//   Omega_quat(0,1)=0.5*omega_vec[2];
//   Omega_quat(0,2)=-0.5*omega_vec[1];
//   Omega_quat(0,3)=0.5*omega_vec[0];
//   Omega_quat(1,0)=-0.5*omega_vec[2];
//   Omega_quat(1,2)=0.5*omega_vec[0];
//   Omega_quat(1,3)=0.5*omega_vec[1];
//   Omega_quat(2,0)=0.5*omega_vec[1];
//   Omega_quat(2,1)=-0.5*omega_vec[0];
//   Omega_quat(2,3)=0.5*omega_vec[2];
//   Omega_quat(3,0)=-0.5*omega_vec[0];
//   Omega_quat(3,1)=-0.5*omega_vec[1];
//   Omega_quat(3,2)=-0.5*omega_vec[2];

}


//creepy flow
void Integrator::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t)
{
  dxdt=Eigen::VectorXd::Zero(3*N);// do not remoove this line important
  //std::cout<<"here t="<<t<<std::endl;
  //remove this
 // cout<<"here4"<<endl;
  double mu=3*M_PI*8.9e-04;
  double tens_force_prefactor=1.0/mu;
  double mag_force_prefactor=(t>0?-2.6/mu:0.0);
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(3*i);
  }
  //std::cout<<angle[0][0]<<" "<<angle[0][1]<<" "<<angle[0][2]<<std::endl<<std::endl;
  
   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
   fill_support_vectors(t);
   fill_mag_foreces();
   fill_external_foreces_and_torques(t);
//    for(size_t i=0; i<N; i++){
//     cout<<i<<" "<<mag_Force[i].transpose()<<endl;
//   }
   fill_steric_foreces();
//  fill_steric_torques();
//  print_forces(x,dxdt,t);abort();
  Eigen::Vector3d Force; Eigen::Vector3d Torque;
  for(size_t i=0; i<N; i++)
  {
    Force=(ext_Force[i]+steric_Force[i])*tens_force_prefactor+mag_force_prefactor*mag_Force[i];
    //Force=mag_Force[i]*dipole_force_prefactor+/*steric_force_prefactor**/(steric_Force[i])+gravitation_force*gravitation_perfactor;
    //Torque=mag_Torque[i]*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i]+/*steric_torque_prefactor**/steric_Torque[i];
    //fill_Omega_quat(Torque);
    //std::cout<<Torque.transpose()<<endl;

//     std::cout <<"i="<<i<<" "<< dxdt.size()<<" "<<x.size()<<std::endl;
    dxdt.segment<3>(3*i)=Force;

  }
    //std::cout <<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
  //std::cout <<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
}

void Integrator::fill_K1_IMEX_BDF3(const Eigen::VectorXd &x,const double dt, const double t)
{
  double mu=3*M_PI*8.9e-04;
  double tens_force_prefactor=1.0/mu*dt;
  double mag_force_prefactor=(t>0?-2.6/mu*dt:0.0);
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(3*i);
  }
  //std::cout<<angle[0][0]<<" "<<angle[0][1]<<" "<<angle[0][2]<<std::endl<<std::endl;

   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
   fill_support_vectors(t);
   fill_mag_foreces();
   fill_external_foreces_and_torques_no_stiffnes(t);
//    for(size_t i=0; i<N; i++){
//     cout<<i<<" "<<mag_Force[i].transpose()<<endl;
//   }
   fill_steric_foreces();
//  fill_steric_torques();
//  print_forces(x,dxdt,t);abort();
  Eigen::Vector3d Force; Eigen::Vector3d Torque;
  for(size_t i=0; i<N; i++)
  {
    Force=(ext_Force[i]+steric_Force[i])*tens_force_prefactor+mag_force_prefactor*mag_Force[i];
    //Force=mag_Force[i]*dipole_force_prefactor+/*steric_force_prefactor**/(steric_Force[i])+gravitation_force*gravitation_perfactor;
    //Torque=mag_Torque[i]*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i]+/*steric_torque_prefactor**/steric_Torque[i];
    //fill_Omega_quat(Torque);
    //std::cout<<Torque.transpose()<<endl;

//     std::cout <<"i="<<i<<" "<< dxdt.size()<<" "<<x.size()<<std::endl;
    K1.segment<3>(3*i)=Force;
  }
}

void Integrator::do_VSIMEX_BDF3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,  Eigen::VectorXd &new_x, const double  t , const double dt)
{
  int info;
  rhs=+(beta1/alpha3)*K1+(beta0/alpha3)*K2;
  K3=K2;
  K2=K1;
  fill_K1_IMEX_BDF3(x,dt,t);
  rhs+=-(alpha2/alpha3)*x-(alpha1/alpha3)*x_old-(alpha0/alpha3)*x_old_old+(beta2/alpha3)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}


void Integrator::do_Euler_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
//  cout <<"t="<<t<< " "<<dt<<endl;
  fill_K1_IMEX_BDF3(x,dt,t);
  rhs=K1+x;
  new_x= solver_sparse.solve(rhs);
  //rhs=(mu*dt)*(A_r-Jt*temp_r);
  //new_x= solver_sparse.solve(rhs)+x;
  //x_temp+=x;
//   cout<<"check norm!!!!!"<<endl;
//   for(size_t i=0; i<Ns; i++)
//   {
//     cout <<i<<" "<<dist(x_temp, i)<<endl;;
//   }

}



void Integrator::fill_support_vectors(double t)
{
  magnetisation={cos(beta), sin(beta)*sin(2*M_PI*omega_B*t), sin(beta)*cos(2*M_PI*omega_B*t)};
}

void Integrator::fill_mag_foreces()
{
  for(size_t i=0; i<N; i++){
    mag_Force[i]=Eigen::Vector3d::Zero();
    //mag_Torque[i]=Eigen::Vector3d::Zero();
  }
  
  for(size_t i=0; i<N; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d rad=pos[j]-pos[i];
      double radius=rad.norm();
      double rad2=radius*radius; //rad^2
      double rad5=rad2*rad2;//rad^4
      rad5*=radius;//rad^5
//       Eigen::Vector3d tau_ij=3.*magnetisation[i].dot(rad)*magnetisation[j].cross(rad)+rad2*magnetisation[i].cross(magnetisation[j]);
//       Eigen::Vector3d tau_ji=3.*magnetisation[j].dot(rad)*magnetisation[i].cross(rad)+rad2*magnetisation[j].cross(magnetisation[i]);
//       mag_Torque[i]+=tau_ji/rad5;
//       mag_Torque[j]+=tau_ij/rad5;
      double sc_pr=rad.dot(magnetisation);
      Eigen::Vector3d F_ij=magnetisation*(2.0*sc_pr)- (5.0/rad2*sc_pr*sc_pr-1) *rad;
      F_ij/=rad5;
      mag_Force[i]-=F_ij;
      mag_Force[j]+=F_ij; 
      //_radius_vec[j+_number_of_particles*i]=_particle_pos[j]-_particle_pos[i];
    }
  }
}


void Integrator::fill_external_foreces_and_torques(const double t){
  size_t i;
  for(i=0; i<N; i++){
    ext_Force[i]=Eigen::Vector3d::Zero();
    //ext_Torque[i]=Eigen::Vector3d::Zero();
  }
  //spring force
  Eigen::Vector3d rad;
  double radius;
  Eigen::Vector3d F_ij;
  for(i=0; i<N-1; i++){
    rad=pos[i+1]-pos[i];
    radius=rad.norm();
    F_ij= (spring_k * (1.0 - spring_l/radius)) * rad;
    ext_Force[i] += F_ij;
    ext_Force[i+1] -= F_ij;
  }
  //bending force
  double k_bend_hetero;
  i=0;
  k_bend_hetero=0.6*k_bend;
  ext_Force[i]-= k_bend_hetero*(-2*pos[1] + pos[0] + pos[2]);
  i++;
  k_bend_hetero =k_bend*(0.6 + 0.4*i/(N-1));
  ext_Force[i]-=k_bend_hetero*(-2*pos[0] + 5*pos[1] - 4*pos[2] + pos[3]);
  for(i=2; i<N-2; i++){
    k_bend_hetero =k_bend*(0.6 + 0.4*i/(N-1));
    ext_Force[i]-=k_bend_hetero*(pos[i-2] - 4*pos[i-1] + 6*pos[i] - 4*pos[i+1] + pos[i+2]);
  }
  i=N-2;
  k_bend_hetero =k_bend*(0.6 + 0.4*i/(N-1));
  ext_Force[i]-=k_bend_hetero*(-2*pos[N-1] + 5*pos[N-2] - 4*pos[N-3] + pos[N-4]);
  i=N-1;
  ext_Force[i]-=k_bend*(-2*pos[N-2] + pos[N-1] + pos[N-3]);  
  
}


void Integrator::fill_external_foreces_and_torques_no_stiffnes(const double t){
  size_t i;
  for(i=0; i<N; i++){
    ext_Force[i]=Eigen::Vector3d::Zero();
    //ext_Torque[i]=Eigen::Vector3d::Zero();
  }
  //spring force
  Eigen::Vector3d rad;
  double radius;
  Eigen::Vector3d F_ij;
  for(i=0; i<N-1; i++){
    rad=pos[i+1]-pos[i];
    radius=rad.norm();
    F_ij= (spring_k * (1.0 - spring_l/radius)) * rad;
    ext_Force[i] += F_ij;
    ext_Force[i+1] -= F_ij;
  }

}

std::vector<size_t> Integrator::sort_indexes(const dvec &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

void Integrator::fill_steric_foreces()
{

  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
    //steric_Torque[i]=Eigen::Vector3d::Zero();
  }
  size_t pos_it;
  double temp;
  Eigen::Vector3d vec_virt_dx=Eigen::Vector3d::Zero();
  Eigen::Vector3d reaction_force_dir=Eigen::Vector3d::Zero();
  double virt_dx=0.;
  Eigen::Vector3d vec_dx=Eigen::Vector3d::Zero();
  Eigen::Vector3d force_dir;
  double dx=0.0;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      vec_dx=pos[j]-pos[i];
      dx=vec_dx.norm();
      if(dx<a)
      {
//         double temp=a*sigma_prefactor/dx;
//         double coef2=temp*temp;
//         double coef4=coef2*coef2;
//         double coef8=coef4*coef4;
//         double coef13=coef8*coef4*temp;
//         double coef7=coef4*coef2*temp;
//         double force_magn=12.0* coef13-6.0*coef7;
//         steric_Force[i]-=force_dir*force_magn;
//         steric_Force[j]+=force_dir*force_magn;
        force_dir=vec_dx/dx;
        double coef=1.0/(1.0-2*dx/a);
        double coef2=coef*coef;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*coef;
        steric_Force[i]-=coef13*force_dir;
        steric_Force[j]+=coef13*force_dir;
      }
    }
  }
}


void Integrator::append_data( Eigen::VectorXd &x, std::vector <dvec> &data, double t)
{
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(3*i);
  }
  fill_support_vectors(t);
  data[0].push_back(t);
  size_t dsize=3;
  for(size_t i=0; i<N; i++)
  {

    data[1+dsize*i].push_back(pos[i][0]);
    data[2+dsize*i].push_back(pos[i][1]);
    data[3+dsize*i].push_back(pos[i][2]);
  }
//   Eigen::Vector3d temp_mag=magnetisation[0];
//   temp_mag[2]=0.0;temp_mag.normalize();
//   B[0]=B_mag*std::cos(t*omega_B+phi0);
//   B[2]=0;
//   B[1]=B_mag*std::sin(t*omega_B+phi0);
//   double alpha, beta;
//   double temp_angle=magnetisation[0].dot(temp_mag);
//   if(temp_angle>1.0)
//   {
//     temp_angle=1.0;
//   }
//   if(magnetisation[0][2]>0){
//     alpha=acos(temp_angle);
//   }
//   else{
//     alpha=-acos(temp_angle);
//   }
//   data[dsize*N+1].push_back(alpha);
//   alpha=(temp_mag.cross(B))[2];
//   if(alpha>0){
//     beta=acos(B.dot(temp_mag));
//   }
//   else{
//      beta=-acos(B.dot(temp_mag));
//   }
//   data[dsize*N+2].push_back(beta);
}

void Integrator::fill_mat_VSIMEX_BDF3( const double dt, const double dt1, const double dt2)
{
  omega2=dt/dt1;
  omega1=dt1/dt2;
  double A1=1+omega1*(1+omega2);
  alpha0=-(omega1*omega1*omega1)*(omega2*omega2)*(1+omega2)/((1+omega1)*A1);
  alpha1=(omega2*omega2)*(omega1+1.0/(1+omega2));
  alpha2=-1-omega2-omega1*omega2*(1+omega2)/(1+omega1);
  alpha3=1+omega2/(1+omega2)+omega1*omega2/A1;
  beta0=omega1*omega1*omega2*(1+omega2)/(1+omega1);
  beta1=-omega2*A1;
  beta2=(1+omega2)*A1/(1+omega1);
}

void Integrator::fill_lhs2_VSIMEX_BDF3(const double dt, const double dt1, const double dt2)
{
  fill_mat_VSIMEX_BDF3(dt, dt1, dt2);
  //double temp=1-M_SQRT1_2;
  double temp=(1.0/alpha3);
  lhs_sparse2=(-dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator::fill_lhs_Euler_impl_expl(const double dt)
{
  //cout<<"fill_lhs_Euler_impl_expl: dt=" <<dt<<endl;
  lhs_sparse=(-dt)*A_sparse;
  lhs_sparse.diagonal().array() += 1;
  solver_sparse.analyzePattern(lhs_sparse);
  solver_sparse.factorize(lhs_sparse);
}
