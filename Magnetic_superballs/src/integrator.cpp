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

Integrator::Integrator(int N_, double omega_B_, Eigen::Vector3d& init_magnetisation_, double alpha_, double komega_kv_, double lag, double gravitation_perfactor1, double angle1, int N_virtual1,  bool B_off):N(N_),
N_virtual(N_virtual1),
a(1.0),
b(M_SQRT2*a),
alpha(alpha_),
komega_kv(komega_kv_),
pos(N),
angle(N),
HD_Force_and_Torque(N),
init_magnetisation(N),
magnetisation(N),
rot(N),
quat(N),
mag_Force(N),
mag_Torque(N),
ext_Torque(N),
steric_Force(N),
steric_Torque(N),
B_mag(B_off==true?1.0:0.0),
omega_B(omega_B_),
B(init_magnetisation_),
phi0(lag),
gravitation_perfactor(gravitation_perfactor1),
B_at_partivcle_pos(N)
{
  std::cout<<"initalizeing"<<std::endl;
  
  //fill_support_vectors();
  phi0=lag;
  fill_virtual_sites();
  double tz = -0.125 * 2. * M_PI;
  double ty = 0.1312899719849809 * 2. * M_PI;
  double tx = 0.25 * 2. * M_PI;
  Eigen::Matrix3d m= Eigen::MatrixXd::Identity(3, 3);;
  m = Eigen::AngleAxisd(tx, Eigen::Vector3d::UnitX())
    * Eigen::AngleAxisd(ty, Eigen::Vector3d::UnitY())
    * Eigen::AngleAxisd(tz, Eigen::Vector3d::UnitZ());
  Eigen::Vector3d new_vec=m*Eigen::Vector3d::UnitZ();
  for(size_t j=0; j<N; j++){
    for(size_t i=0; i<N_virtual; i++)
    {
      //cout<<virtual_sites[j][i].transpose()<<" "<<(m *virtual_sites[j][i]).transpose()<<endl;
      virtual_sites[j][i]= m *virtual_sites[j][i];
    }
  }
  Eigen::Vector3d test=virtual_sites[0][6]-virtual_sites[0][7]; test.normalize();
  cout<<test.transpose()<<endl;
  cout<<new_vec.transpose()<<endl;
  //perturbation momentum out of plane
  ty =angle1 * 2. * M_PI;
  m =Eigen::AngleAxisd(ty, new_vec);
  for(size_t j=0; j<N; j++){
    for(size_t i=0; i<N_virtual; i++)
    {
      cout<<virtual_sites[j][i].transpose()<<" "<<(m *virtual_sites[j][i]).transpose()<<endl;
      virtual_sites[j][i]= m *virtual_sites[j][i];
    }
  }
  for(size_t i=0; i<N; i++)
  {
    init_magnetisation[i]=m*init_magnetisation_;
  }
 // cout<<init_magnetisation[0].transpose()<<endl;
 // cout<<virtual_sites[0][7][2]-virtual_sites_radius[7]<<endl;
  cout<<"Angle "<<angle1<<endl;;
  //abort();
  new_virtual_sites=virtual_sites;
  magnetisation=init_magnetisation;
  
  Eigen::Vector4d kv={0,0,1,0};
  Eigen::Quaterniond q(kv);

  std::cout << "This quaternion consists of a scalar " << q.w() << " and a vector " << std::endl << q.vec() << std::endl;



  q.normalize();

  std::cout << "To represent rotation, we need to normalize it such that its length is " << q.norm() << std::endl;
  Eigen::Quaterniond a = {0,1,0,0};
  //a.normalize();
  cout<<a.coeffs().transpose()<<" "<<a.w()<<endl;
  Eigen::Vector3d aa={0,0,1};
  fill_Omega_quat(aa);
  cout<<Omega_quat*a.coeffs()<<endl;
  //abort();
}

void Integrator::fill_virtual_sites()
{
  if(N_virtual>0){
    if(N_virtual>92){
      N_virtual=92;
    }
    radius_big_particle=0.5*a;
    double x_coner=1.0/std::pow(3,0.25)/3.0*2*radius_big_particle;
    double radius_coners=std::pow(3,0.25)/6.0*2*radius_big_particle;
 //   double radius_coners=0.01*2*radius_big_particle;
 //   double x_coner=0.5*2*radius_big_particle-radius_coners;
    double radius_midpoints=std::pow(2,0.25)/6.0*2*radius_big_particle;
    double radius_red=0.19782*2*radius_big_particle;
    double radius_red1=0.257166*2*radius_big_particle;;
    double radius_red2=0.237841*2*radius_big_particle;;
    double x_midpoints=1.0/std::pow(2,0.25)/3.0*2*radius_big_particle;
    double x_red=0.27953*2*radius_big_particle;
    double y_red=0.153645*2*radius_big_particle;
    double x_red1=0.15215*2*radius_big_particle;
    double y_red1=0.240383*2*radius_big_particle;
    double x_red2=0.178519*2*radius_big_particle;
    double y_red2=0.2598*2*radius_big_particle;
//     cout<<x_coner<<" "<<radius_midpoints<<" "<<x_midpoints<<" "<<radius_coners<<endl;
//     abort();

    std::vector <Eigen::Vector3d> Points; Points.resize(92);
    virtual_sites_radius.resize(92);
    size_t i=0;
    //midpoints
    for(i=0;i<8;i++){
      virtual_sites_radius[i]=radius_coners;
    }
    //corners
    for(i=8;i<8+12;i++){
      virtual_sites_radius[i]=radius_midpoints;
    }
    //red
    for(i=20;i<20+24;i++){
      virtual_sites_radius[i]=radius_red;
    }
    //red1
    for(i=44;i<44+24;i++){
      virtual_sites_radius[i]=radius_red1;
    }
    //red2
    for(i=68;i<68+24;i++){
      virtual_sites_radius[i]=radius_red2;
    }
    Points[0]={x_coner,x_coner,x_coner};
    Points[1]={x_coner,x_coner,-x_coner};
    Points[2]={x_coner,-x_coner,-x_coner};
    Points[3]={x_coner,-x_coner,x_coner};
    Points[4]={-x_coner,-x_coner,-x_coner};
    Points[5]={-x_coner,-x_coner,x_coner};
    Points[6]={-x_coner,x_coner,x_coner};
    Points[7]={-x_coner,x_coner,-x_coner};
    Points[8]={-x_midpoints,0.0,-x_midpoints};
    Points[9] ={-x_midpoints,0.0,x_midpoints};
    Points[10] ={x_midpoints,0.0,x_midpoints};
    Points[11] ={x_midpoints,0.0,-x_midpoints};
    Points[12] ={0.0,x_midpoints,x_midpoints};
    Points[13] ={0.0,x_midpoints,-x_midpoints};
    Points[14] ={0.0,-x_midpoints,-x_midpoints};
    Points[15] ={0.0,-x_midpoints,x_midpoints};
    Points[16] ={-x_midpoints,x_midpoints,0.0};
    Points[17] ={x_midpoints,x_midpoints,0.0};
    Points[18] ={x_midpoints,-x_midpoints,0.0};
    Points[19] ={-x_midpoints,-x_midpoints,0.0};
    Points[20] ={x_red, x_red, y_red};
    Points[21] ={x_red, x_red, -y_red};
    Points[22] ={x_red, -x_red, y_red};
    Points[23] ={-x_red, x_red, y_red};
    Points[24] ={-x_red, -x_red, y_red};
    Points[25] ={x_red, -x_red, -y_red};
    Points[26] ={-x_red, x_red, -y_red};
    Points[27] ={-x_red, -x_red, -y_red};
    Points[28] ={x_red, y_red, x_red};
    Points[29] ={x_red, y_red, -x_red};
    Points[30] ={x_red, -y_red, x_red};
    Points[31] ={-x_red, y_red, x_red};
    Points[32] ={-x_red, -y_red, x_red};
    Points[33] ={x_red, -y_red, -x_red};
    Points[34] ={-x_red, y_red, -x_red};
    Points[35] ={-x_red, -y_red, -x_red};
    Points[36] ={y_red, x_red,x_red};
    Points[37] ={y_red, x_red, -x_red};
    Points[38] ={y_red, -x_red,x_red};
    Points[39] ={-y_red, x_red,x_red};
    Points[40] ={-y_red, -x_red,x_red};
    Points[41] ={y_red, -x_red, -x_red};
    Points[42] ={-y_red, x_red, -x_red};
    Points[43] ={-y_red, -x_red, -x_red};
    Points[44] ={x_red1, x_red1, y_red1};
    Points[45] ={x_red1, x_red1, -y_red1};
    Points[46] ={x_red1, -x_red1, y_red1};
    Points[47] ={-x_red1, x_red1, y_red1};
    Points[48] ={-x_red1, -x_red1, y_red1};
    Points[49] ={x_red1, -x_red1, -y_red1};
    Points[50] ={-x_red1, x_red1, -y_red1};
    Points[51] ={-x_red1, -x_red1, -y_red1};
    Points[52] ={x_red1, y_red1, x_red1};
    Points[53] ={x_red1, y_red1, -x_red1};
    Points[54] ={x_red1, -y_red1, x_red1};
    Points[55] ={-x_red1, y_red1, x_red1};
    Points[56] ={-x_red1, -y_red1, x_red1};
    Points[57] ={x_red1, -y_red1, -x_red1};
    Points[58] ={-x_red1, y_red1, -x_red1};
    Points[59] ={-x_red1, -y_red1, -x_red1};
    Points[60] ={y_red1, x_red1,x_red1};
    Points[61] ={y_red1, x_red1, -x_red1};
    Points[62] ={y_red1, -x_red1,x_red1};
    Points[63] ={-y_red1, x_red1,x_red1};
    Points[64] ={-y_red1, -x_red1,x_red1};
    Points[65] ={y_red1, -x_red1, -x_red1};
    Points[66] ={-y_red1, x_red1, -x_red1};
    Points[67] ={-y_red1, -x_red1, -x_red1};
    Points[68] ={x_red2, y_red2, 0.0};
    Points[69] ={y_red2, x_red2, 0.0};
    Points[70] ={x_red2, 0.0, y_red2};
    Points[71] ={y_red2, 0.0, x_red2};
    Points[72] ={0.0, x_red2, y_red2};
    Points[73] ={0.0, y_red2, x_red2};
    Points[74] ={-x_red2, -y_red2, 0.0};
    Points[75] ={-y_red2, -x_red2, 0.0};
    Points[76] ={-x_red2, 0.0, -y_red2};
    Points[77] ={-y_red2, 0.0, -x_red2};
    Points[78] ={0.0, -x_red2, -y_red2};
    Points[79] ={0.0, -y_red2, -x_red2};
    Points[80] ={-x_red2, y_red2, 0.0};
    Points[81] ={-y_red2, x_red2, 0.0};
    Points[82] ={-x_red2, 0.0, y_red2};
    Points[83] ={-y_red2, 0.0, x_red2};
    Points[84] ={0.0, -x_red2, y_red2};
    Points[85] ={0.0, -y_red2, x_red2};
    Points[86] ={x_red2, -y_red2, 0.0};
    Points[87] ={y_red2, -x_red2, 0.0};
    Points[88] ={x_red2, 0.0, -y_red2};
    Points[89] ={y_red2, 0.0, -x_red2};
    Points[90] ={0.0, x_red2, -y_red2};
    Points[91] ={0.0, y_red2, -x_red2};
    if(N_virtual<92)
    {
      virtual_sites_radius.resize(N_virtual);
      Points.resize(N_virtual);
    }
    for(i=0;i<N; i++){
      virtual_sites.push_back(Points);
    }
    
    
  }//end  if Nvirtual>0
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


void Integrator::debug(const Eigen::VectorXd& x, const double t)
{
  
  cout<<"!!!!!!!!!!!! Start debug"<<endl;
  double dipole_torque_prefactor=alpha/3.0;
  double dipole_force_prefactor=alpha/3.0*komega_kv;
  double steric_force_prefactor=dipole_force_prefactor;
  double steric_torque_prefactor=dipole_torque_prefactor/*4*/;
  double ext_torque_prefactor=1.0;
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(7*i);
    quat[i].coeffs()=x.segment<4>(7*i+3);
  }
  //std::cout<<angle[0][0]<<" "<<angle[0][1]<<" "<<angle[0][2]<<std::endl<<std::endl;
  
   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
   fill_support_vectors();
   fill_external_foreces_and_torques(t);
 //  fill_mag_foreces_and_torques();
//    std::cout <<"here4"<<std::endl;
     //fill_steric_foreces_and_torques();
   fill_steric_foreces_print();
   fill_steric_torques_cout();
//   //std::cout <<"here6"<<std::endl;
//   print_forces(x,dxdt,t);abort();
   cout.precision(10);
   cout<<"Forces"<<" "<<mag_Force[0].transpose()<<" "<<steric_Force[0].transpose()<<" "<<gravitation_force.transpose()<<endl;
   cout<<"Torques"<<" "<<mag_Torque[0].transpose()<<" "<<ext_Torque[0].transpose()<<" "<<steric_Torque[0].transpose()<<endl;
   cout<<"B="<<B.transpose()<<" "<<magnetisation[0].transpose()<<endl;;
   for(size_t i=0; i<N; i++)
  {
    cout<<pos[i].transpose()<<endl;
  }
  if(N>1)
  {
    cout<<(pos[1]-pos[0]).norm()<<endl;
  }
  Eigen::Vector3d Force; Eigen::Vector3d Torque;
  for(size_t i=0; i<N; i++)
  {
    Force=mag_Force[i]*dipole_force_prefactor+/*steric_force_prefactor**/(steric_Force[i])+gravitation_force*gravitation_perfactor;
    Torque=mag_Torque[i]*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i]+/*steric_torque_prefactor**/steric_Torque[i];
    fill_Omega_quat(Torque);
    //std::cout<<Torque.transpose()<<endl;

//     std::cout <<"i="<<i<<" "<< dxdt.size()<<" "<<x.size()<<std::endl;
  }
    //std::cout <<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
  //std::cout <<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
}



//creepy flow
void Integrator::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t)
{
  dxdt=Eigen::VectorXd::Zero(7*N);// do not remoove this line important
  //std::cout<<"here t="<<t<<std::endl;
  //remove this
 // cout<<"here4"<<endl;
  double dipole_torque_prefactor=alpha/3.0;
  double dipole_force_prefactor=alpha/3.0*komega_kv;
  double steric_force_prefactor=dipole_force_prefactor;
  double steric_torque_prefactor=dipole_torque_prefactor/*4*/;
  double ext_torque_prefactor=1.0;
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(7*i);
    quat[i].coeffs()=x.segment<4>(7*i+3);
  }
  //std::cout<<angle[0][0]<<" "<<angle[0][1]<<" "<<angle[0][2]<<std::endl<<std::endl;
  
   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
   fill_support_vectors();
   fill_external_foreces_and_torques(t);
   fill_mag_foreces_and_torques();
     //fill_steric_foreces_and_torques();
   fill_steric_foreces();
   fill_steric_torques();
//   print_forces(x,dxdt,t);abort();
  Eigen::Vector3d Force; Eigen::Vector3d Torque;
  for(size_t i=0; i<N; i++)
  {
    Force=mag_Force[i]*dipole_force_prefactor+/*steric_force_prefactor**/(steric_Force[i])+gravitation_force*gravitation_perfactor;
    Torque=mag_Torque[i]*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i]+/*steric_torque_prefactor**/steric_Torque[i];
    fill_Omega_quat(Torque);
    //std::cout<<Torque.transpose()<<endl;

//     std::cout <<"i="<<i<<" "<< dxdt.size()<<" "<<x.size()<<std::endl;
    dxdt.segment<3>(7*i)=Force;

    dxdt.segment<4>(7*i+3)=Omega_quat*quat[i].coeffs();

  }
    //std::cout <<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
  //std::cout <<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
}

void Integrator::fill_support_vectors()
{
  for(size_t i=0; i<N; i++)
  {
    //rot[i] = Eigen::AngleAxisd(angle[i][0], Eigen::Vector3d::UnitX())
    //* Eigen::AngleAxisd(angle[i][1], Eigen::Vector3d::UnitY())
    //* Eigen::AngleAxisd(angle[i][2], Eigen::Vector3d::UnitZ());

    magnetisation[i]=quat[i]*init_magnetisation[i];
  }
}

void Integrator::fill_mag_foreces_and_torques()
{
  for(size_t i=0; i<N; i++){
    mag_Force[i]=Eigen::Vector3d::Zero();
    mag_Torque[i]=Eigen::Vector3d::Zero();
  }
  
  for(size_t i=0; i<N; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d rad=pos[j]-pos[i];
      double radius=rad.norm();
      double rad2=radius*radius; //rad^2
      double rad5=rad2*rad2;//rad^4
      rad5*=radius;//rad^5
      Eigen::Vector3d tau_ij=3.*magnetisation[i].dot(rad)*magnetisation[j].cross(rad)+rad2*magnetisation[i].cross(magnetisation[j]);
      Eigen::Vector3d tau_ji=3.*magnetisation[j].dot(rad)*magnetisation[i].cross(rad)+rad2*magnetisation[j].cross(magnetisation[i]);
      mag_Torque[i]+=tau_ji/rad5;
      mag_Torque[j]+=tau_ij/rad5;
      Eigen::Vector3d F_ij=rad*magnetisation[i].dot(magnetisation[j]) +magnetisation[i]*rad.dot(magnetisation[j])+ magnetisation[j]*magnetisation[i].dot(rad)- (5.0/rad2*rad.dot(magnetisation[i])*rad.dot(magnetisation[j])) *rad;
      F_ij/=rad5;
      mag_Force[i]-=F_ij;
      mag_Force[j]+=F_ij; 
      //_radius_vec[j+_number_of_particles*i]=_particle_pos[j]-_particle_pos[i];
    }
  }
}

void Integrator::calc_B_field()
{
  for(size_t i=0; i<N; i++){
    B_at_partivcle_pos[i]=Eigen::Vector3d::Zero();
  }
  for(size_t i=0; i<N; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d rad=pos[j]-pos[i];
      double radius=rad.norm();
      double rad2=radius*radius; //rad^2
      double rad5=rad2*rad2;//rad^4
      rad5*=radius;//rad^5
      B_at_partivcle_pos[i]+= (3.0/rad5*magnetisation[j].dot(rad)) * rad- (1.0/rad2/radius)*magnetisation[j];
      B_at_partivcle_pos[j]+= (3.0/rad5*magnetisation[i].dot(rad)) * rad- (1.0/rad2/radius)*magnetisation[i];
    }
  }
}


void Integrator::fill_external_foreces_and_torques(const double t)
{
  B[0]=B_mag*std::cos(t*omega_B+phi0);
  B[2]=0;
  B[1]=B_mag*std::sin(t*omega_B+phi0);
//   calc_B_field();
//   for(size_t i=0; i<N; i++){
//    magnetisation[i]=/*magnetisation[i]*0.5+0.5**/(0.5*(B/*+B_at_partivcle_pos[i]*/));
//    //std::cout<<"i="<<i<<" t="<<t<<" "<<B_at_partivcle_pos[i][0]<<" "<<B_at_partivcle_pos[i][1]<<" "<<B[0]<<" "<<B[1]<<" "<<magnetisation[i][0]<<" "<<magnetisation[i][1] <<std::endl;
//   }
  for(size_t i=0; i<N; i++){
    ext_Torque[i]=magnetisation[i].cross(B);
  }
  //std::cout<<ext_Torque[0].transpose()<<" "<<t<<" "<<magnetisation[0].transpose()<<" "<<B.transpose()<<" "<<+phi0<<std::endl; abort();
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



void Integrator::fill_steric_torques()
{
  if(N_virtual>=8)
  {
    update_virtual_pos();
    for(size_t i=0; i<N; i++)
    {
      if(pos[i][2]-bottom_pozition<radius_big_particle){
        double temp=radius_big_particle*sigma_prefactor/(pos[i][2]-bottom_pozition);
        double coef2=temp*temp;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*temp;
        double coef7=coef4*coef2*temp;
        steric_Force[i]+=reaction_force_dir*(12.0* coef13-6.0*coef7);
      }
      for(size_t j=0; j<N_virtual; j++)
      {
        if(new_virtual_sites[i][j][2]+pos[i][2]-bottom_pozition<virtual_sites_radius[j]){
          double temp=virtual_sites_radius[j]*sigma_prefactor/(new_virtual_sites[i][j][2]+pos[i][2]-bottom_pozition);
          double coef2=temp*temp;
          double coef4=coef2*coef2;
          double coef8=coef4*coef4;
          double coef13=coef8*coef4*temp;
          double coef7=coef4*coef2*temp;
          steric_Torque[i]+=new_virtual_sites[i][j].cross(reaction_force_dir)*(12.0* coef13-6.0*coef7);
          steric_Force[i]+=reaction_force_dir*(12.0* coef13-6.0*coef7);
        }
      }
    }
  }
}

void Integrator::fill_steric_torques_cout()
{
  for(size_t i=0; i<N; i++){
    steric_Torque[i]=Eigen::Vector3d::Zero();
  }
  if(N_virtual>=8)
  {
    update_virtual_pos();
    for(size_t i=0; i<N; i++)
    {      
      if(pos[i][2]-bottom_pozition<radius_big_particle){
        double temp=radius_big_particle*sigma_prefactor/(pos[i][2]-bottom_pozition);
        double coef2=temp*temp;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*temp;
        double coef7=coef4*coef2*temp;
        steric_Force[i]+=reaction_force_dir*(12.0* coef13-6.0*coef7);
        cout<<"Steric big:"<<i<<" "<<(reaction_force_dir*(12.0* coef13-6.0*coef7)).transpose()<<endl;;
      }
      cout<<"Steric big:"<<pos[i][2]-bottom_pozition<<" "<<radius_big_particle<<endl;
      for(size_t j=0; j<8; j++)
      {
        if(new_virtual_sites[i][j][2]+pos[i][2]-bottom_pozition<virtual_sites_radius[j]){
          double temp=virtual_sites_radius[j]*sigma_prefactor/(new_virtual_sites[i][j][2]+pos[i][2]-bottom_pozition);
          double coef2=temp*temp;
          double coef4=coef2*coef2;
          double coef8=coef4*coef4;
          double coef13=coef8*coef4*temp;
          double coef7=coef4*coef2*temp;
          cout<<"Steric:"<<j<<" "<<(new_virtual_sites[i][j].cross(reaction_force_dir)*(12.0* coef13-6.0*coef7)).transpose()<<
          " "<<(reaction_force_dir*(12.0* coef13-6.0*coef7)).transpose()<<" A "<<new_virtual_sites[i][j][2]+pos[i][2]<<" "<<bottom_pozition+virtual_sites_radius[j]<<endl;;
          steric_Torque[i]+=new_virtual_sites[i][j].cross(reaction_force_dir)*(12.0* coef13-6.0*coef7);
          steric_Force[i]+=reaction_force_dir*(12.0* coef13-6.0*coef7);
        }
      }
    }
  }
}


void Integrator::fill_steric_foreces()
{

  update_virtual_pos();
  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
    steric_Torque[i]=Eigen::Vector3d::Zero();
  }
  size_t pos_it;
  double temp;
  Eigen::Vector3d vec_virt_dx=Eigen::Vector3d::Zero();
  Eigen::Vector3d reaction_force_dir=Eigen::Vector3d::Zero();
  double virt_dx=0.;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d vec_dx=Eigen::Vector3d::Zero();
      vec_dx=pos[j]-pos[i];

      double dx=vec_dx.norm();
      Eigen::Vector3d force_dir=vec_dx/dx;
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
        double coef=1.0/(1.0-2*dx/a);
        double coef2=coef*coef;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*coef;
        steric_Force[i]=-coef13*force_dir;
        steric_Force[j]=coef13*force_dir;
      }
      //if(N_virtual>=8)
      {
         for(size_t k=0; k<N_virtual; k++){
          vec_virt_dx=vec_dx-new_virtual_sites[i][k];
          virt_dx=vec_virt_dx.norm();
          if(virt_dx<virtual_sites_radius[k]+radius_big_particle){
            reaction_force_dir=vec_virt_dx/virt_dx;
            double temp=(virtual_sites_radius[k]+radius_big_particle)*sigma_prefactor/(virt_dx);
            double coef2=temp*temp;
            double coef4=coef2*coef2;
            double coef8=coef4*coef4;
            double coef13=coef8*coef4*temp;
            double coef7=coef4*coef2*temp;
            double force_magn=12.0* coef13-6.0*coef7;
            steric_Torque[i]-=new_virtual_sites[i][k].cross(reaction_force_dir)*force_magn;
            steric_Force[i]-=reaction_force_dir*force_magn;
            steric_Force[j]+=reaction_force_dir*force_magn;
         }
         vec_virt_dx=vec_dx+new_virtual_sites[j][k];
         virt_dx=vec_virt_dx.norm();
          if(virt_dx<virtual_sites_radius[k]+radius_big_particle){
            reaction_force_dir=vec_virt_dx/virt_dx;
            double temp=(virtual_sites_radius[k]+radius_big_particle)*sigma_prefactor/(virt_dx);
            double coef2=temp*temp;
            double coef4=coef2*coef2;
            double coef8=coef4*coef4;
            double coef13=coef8*coef4*temp;
            double coef7=coef4*coef2*temp;
            double force_magn=12.0* coef13-6.0*coef7;
            steric_Torque[j]+=new_virtual_sites[i][k].cross(reaction_force_dir)*force_magn;
            steric_Force[i]-=reaction_force_dir*force_magn;
            steric_Force[j]+=reaction_force_dir*force_magn;
        }
          for(size_t l=0; l<N_virtual; l++)
          {
            vec_virt_dx=vec_dx+new_virtual_sites[j][l]-new_virtual_sites[i][k];
            virt_dx=vec_virt_dx.norm();
            if(virt_dx<virtual_sites_radius[k]+virtual_sites_radius[l]){
              reaction_force_dir=vec_virt_dx/virt_dx;
              double temp=(virtual_sites_radius[k]+virtual_sites_radius[l])*sigma_prefactor/(virt_dx);
              double coef2=temp*temp;
              double coef4=coef2*coef2;
              double coef8=coef4*coef4;
              double coef13=coef8*coef4*temp;
              double coef7=coef4*coef2*temp;
              double force_magn=12.0* coef13-6.0*coef7;
              steric_Torque[i]-=new_virtual_sites[i][k].cross(reaction_force_dir)*force_magn;
              steric_Force[i]-=reaction_force_dir*force_magn;
              steric_Force[j]+=reaction_force_dir*force_magn;
              steric_Torque[j]+=new_virtual_sites[j][l].cross(reaction_force_dir)*force_magn;
            }
          }
        }
      }
    }
  }
}

void Integrator::fill_steric_foreces_print()
{
  cout<<"start fill_steric_foreces_print"<<endl;
  update_virtual_pos();
  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
    steric_Torque[i]=Eigen::Vector3d::Zero();
  }
  size_t pos_it;
  double temp;
  Eigen::Vector3d vec_virt_dx=Eigen::Vector3d::Zero();
  Eigen::Vector3d reaction_force_dir=Eigen::Vector3d::Zero();
  double virt_dx=0.;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d vec_dx=Eigen::Vector3d::Zero();
      vec_dx=pos[j]-pos[i];

      double dx=vec_dx.norm();
      Eigen::Vector3d force_dir=vec_dx/dx;
      if(dx<a)
      {
        double temp=a*sigma_prefactor/dx;
        double coef2=temp*temp;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*temp;
        double coef7=coef4*coef2*temp;
        double force_magn=12.0* coef13-6.0*coef7;
        steric_Force[i]-=force_dir*force_magn;
        steric_Force[j]+=force_dir*force_magn;
//         double coef=1.0/(1.0-2*dx/a);
//         double coef2=coef*coef;
//         double coef4=coef2*coef2;
//         double coef8=coef4*coef4;
//         double coef13=coef8*coef4*coef;
        //steric_Force[i]=coef13*force_dir;
        //steric_Force[j]-=coef13*force_dir;
      }
      //if(N_virtual>=8)
      {
         for(size_t k=0; k<N_virtual; k++){
            vec_virt_dx=vec_dx-new_virtual_sites[i][k];
            virt_dx=vec_virt_dx.norm();
            if(virt_dx<virtual_sites_radius[k]+radius_big_particle){
              reaction_force_dir=vec_virt_dx/virt_dx;
              double temp=(virtual_sites_radius[k]+radius_big_particle)*sigma_prefactor/(virt_dx);
              double coef2=temp*temp;
              double coef4=coef2*coef2;
              double coef8=coef4*coef4;
              double coef13=coef8*coef4*temp;
              double coef7=coef4*coef2*temp;
              double force_magn=12.0* coef13-6.0*coef7;
              steric_Torque[i]-=new_virtual_sites[i][k].cross(reaction_force_dir)*force_magn;
              steric_Force[i]-=reaction_force_dir*force_magn;
              steric_Force[j]+=reaction_force_dir*force_magn;
          }
          vec_virt_dx=vec_dx+new_virtual_sites[j][k];
          virt_dx=vec_virt_dx.norm();
            if(virt_dx<virtual_sites_radius[k]+radius_big_particle){
              reaction_force_dir=vec_virt_dx/virt_dx;
              double temp=(virtual_sites_radius[k]+radius_big_particle)*sigma_prefactor/(virt_dx);
              double coef2=temp*temp;
              double coef4=coef2*coef2;
              double coef8=coef4*coef4;
              double coef13=coef8*coef4*temp;
              double coef7=coef4*coef2*temp;
              double force_magn=12.0* coef13-6.0*coef7;
              steric_Torque[j]+=new_virtual_sites[i][k].cross(reaction_force_dir)*force_magn;
              steric_Force[i]-=reaction_force_dir*force_magn;
              steric_Force[j]+=reaction_force_dir*force_magn;
          }
          for(size_t l=0; l<N_virtual; l++)
          {
            vec_virt_dx=vec_dx+new_virtual_sites[j][l]-new_virtual_sites[i][k];
            virt_dx=vec_virt_dx.norm();
            if(virt_dx<virtual_sites_radius[k]+virtual_sites_radius[l]){
              reaction_force_dir=vec_virt_dx/virt_dx;
              cout<<k<<" "<<l<<" "<<virt_dx<<" "<<virtual_sites_radius[k]<<" "<<virtual_sites_radius[k]<<endl;
              cout<<new_virtual_sites[i][k].transpose()<<" "<<new_virtual_sites[j][l].transpose()<<" "<<reaction_force_dir.transpose()<<endl;
              cout<<(new_virtual_sites[i][k]+pos[i]).transpose()<<" "<<(new_virtual_sites[j][l]+pos[j]).transpose()<<" "<<reaction_force_dir.transpose()<<endl;
              double temp=(virtual_sites_radius[k]+virtual_sites_radius[l])*sigma_prefactor/(virt_dx);
              double coef2=temp*temp;
              double coef4=coef2*coef2;
              double coef8=coef4*coef4;
              double coef13=coef8*coef4*temp;
              double coef7=coef4*coef2*temp;
              double force_magn=12.0* coef13-6.0*coef7;
              steric_Torque[i]-=new_virtual_sites[i][k].cross(reaction_force_dir)*force_magn;
              steric_Force[i]-=reaction_force_dir*force_magn;
              steric_Force[j]+=reaction_force_dir*force_magn;
              steric_Torque[j]+=new_virtual_sites[j][l].cross(reaction_force_dir)*force_magn;
            }
          }
        }
      }
    }
  }
  cout<<"The setric Force="<<steric_Force[0].transpose()<<endl;
}

void Integrator::update_virtual_pos()
{
  for(size_t j=0; j<N; j++){
    for(size_t i=0; i<N_virtual; i++)
    {
      new_virtual_sites[j][i]= quat[j] *virtual_sites[j][i];
    }
  }
}

void Integrator::print_mathematica(const Eigen::VectorXd & x, std::string folder)
{
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(7*i);
    quat[i].coeffs()=x.segment<4>(7*i+3);
  }
  update_virtual_pos();
  FILE *file;
  file=NULL;
  ostringstream os;
  os.str("");
  os<<folder<<".txt";
  file=fopen(os.str().c_str(),"w");
  size_t i,j;
  if(file!=NULL)
  {
    for(i=0; i<N; i++){
      fprintf(file,"part0%ld={{Blue,Ball[{%.15f, %.15f, %.15f}, %.15f]}",i,pos[i][0],pos[i][1],pos[i][2],0.5*a);
      
      if(N_virtual>=8){
        for(j=0;j<8; j++){
          fprintf(file,", {Orange,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
        }
        if(N_virtual>=20){
          for(j=8;j<20; j++){
            fprintf(file,", {White,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
          }
          if(N_virtual>=44){
            for(j=20;j<44; j++){
              fprintf(file,", {Magenta,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
            }
            if(N_virtual>=68){
              for(j=44;j<68; j++){
                fprintf(file,", {Red,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
              }
              if(N_virtual>=92){
                for(j=68;j<92; j++){
                  fprintf(file,", {Green,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
                }
              }//92
              else
              {
                for(j=68;j<N_virtual; j++){
                  fprintf(file,", {Green,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
                }
              }
            }//68
            else{
              for(j=44;j<N_virtual; j++){
                fprintf(file,", {Red,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
              }
            }
          }//44
          else{
            for(j=20;j<N_virtual; j++){
              fprintf(file,", {Magenta,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
            }
          }
        }
        else{
          for(j=8;j<N_virtual; j++){
            fprintf(file,", {White,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
          }
        }
      }
      else{
        for(j=0;j<N_virtual; j++){
          fprintf(file,", {Orange,Ball[{%.15f, %.15f, %.15f}, %.15f]}",new_virtual_sites[i][j][0],new_virtual_sites[i][j][1],new_virtual_sites[i][j][2], virtual_sites_radius[j]);
        }
      }
      fprintf(file,",{Green, Arrow[{{%.15f, %.15f, %.15f},am*{%.15f, %.15f, %.15f }+{%.15f, %.15f, %.15f }}]}};\n",pos[i][0],pos[i][1], pos[i][2], magnetisation[i][0],magnetisation[i][1],magnetisation[i][2],pos[i][0],pos[i][1], pos[i][2]);
    }
    fclose(file);
  }
}
void Integrator::append_data_virtual( Eigen::VectorXd &x, std::vector <dvec> &data, double t)
{
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(7*i);
    quat[i].coeffs()=x.segment<4>(7*i+3);
  }
  update_virtual_pos();
  fill_support_vectors();
  data[0].push_back(t);
  size_t dsize=6+N_virtual*4;
  for(size_t i=0; i<N; i++)
  {
    data[1+dsize*i].push_back(pos[i][0]);
    data[2+dsize*i].push_back(pos[i][1]);
    data[3+dsize*i].push_back(pos[i][2]);
    data[4+dsize*i].push_back(magnetisation[i][0]);
    data[5+dsize*i].push_back(magnetisation[i][1]);
    data[6+dsize*i].push_back(magnetisation[i][2]);
    for(size_t j=0; j<N_virtual; j++){
       data[7+dsize*i+4*j].push_back(new_virtual_sites[i][j][0]+pos[i][0]);
       data[8+dsize*i+4*j].push_back(new_virtual_sites[i][j][1]+pos[i][1]);
       data[9+dsize*i+4*j].push_back(new_virtual_sites[i][j][2]+pos[i][2]);
       data[10+dsize*i+4*j].push_back(virtual_sites_radius[j]);
    }
  }
}

void Integrator::append_data( Eigen::VectorXd &x, std::vector <dvec> &data, double t)
{
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(7*i);
    quat[i].coeffs()=x.segment<4>(7*i+3);
  }
  fill_support_vectors();
  data[0].push_back(t);
  size_t dsize=6;
  for(size_t i=0; i<N; i++)
  {

    data[1+dsize*i].push_back(pos[i][0]);
    data[2+dsize*i].push_back(pos[i][1]);
    data[3+dsize*i].push_back(pos[i][2]);
    data[4+dsize*i].push_back(magnetisation[i][0]);
    data[5+dsize*i].push_back(magnetisation[i][1]);
    data[6+dsize*i].push_back(magnetisation[i][2]);
  }
  Eigen::Vector3d temp_mag=magnetisation[0];
  temp_mag[2]=0.0;temp_mag.normalize();
  B[0]=B_mag*std::cos(t*omega_B+phi0);
  B[2]=0;
  B[1]=B_mag*std::sin(t*omega_B+phi0);
  double alpha, beta;
  double temp_angle=magnetisation[0].dot(temp_mag);
  if(temp_angle>1.0)
  {
    temp_angle=1.0;
  }
  if(magnetisation[0][2]>0){
    alpha=acos(temp_angle);
  }
  else{
    alpha=-acos(temp_angle);
  }
  data[dsize*N+1].push_back(alpha);
  alpha=(temp_mag.cross(B))[2];
  if(alpha>0){
    beta=acos(B.dot(temp_mag));
  }
  else{
     beta=-acos(B.dot(temp_mag));
  }
  data[dsize*N+2].push_back(beta);
}

void Integrator::append_data_quat( Eigen::VectorXd &x, std::vector <dvec> &data, double t)
{
  for(size_t i=0; i<N; i++)
  {
    pos[i] = x.segment<3>(7*i);
    quat[i].coeffs()=x.segment<4>(7*i+3);
  }
  fill_support_vectors();
  data[0].push_back(t);
//   size_t dsize=7;
//   for(size_t i=0; i<x.size(); i++)
//   {
//
//     data[1+i].push_back(x[i]);
//   }
  B[0]=B_mag*std::cos(t*omega_B+phi0);
  B[2]=0;
  B[1]=B_mag*std::sin(t*omega_B+phi0);
  size_t dsize=7+6;
  for(size_t i=0; i<N; i++)
  {
    data[1+dsize*i].push_back(pos[i][0]);
    data[2+dsize*i].push_back(pos[i][1]);
    data[3+dsize*i].push_back(pos[i][2]);
    data[4+dsize*i].push_back(quat[i].x());
    data[5+dsize*i].push_back(quat[i].y());
    data[6+dsize*i].push_back(quat[i].z());
    data[7+dsize*i].push_back(quat[i].w());
    data[8+dsize*i].push_back(magnetisation[i][0]);
    data[9+dsize*i].push_back(magnetisation[i][1]);
    data[10+dsize*i].push_back(magnetisation[i][2]);
    data[11+dsize*i].push_back(B[0]);
    data[12+dsize*i].push_back(B[1]);
    data[13+dsize*i].push_back(B[2]);
  }
  Eigen::Vector3d temp_mag=magnetisation[0];
  temp_mag[2]=0.0;temp_mag.normalize();
  double alpha, beta;
  double temp_angle=magnetisation[0].dot(temp_mag);
  if(temp_angle>1.0)
  {
    temp_angle=1.0;
  }
  if(magnetisation[0][2]>0){
    alpha=acos(temp_angle);
  }
  else{
    alpha=-acos(temp_angle);
  }
  data[dsize*N+1].push_back(alpha);
  alpha=(temp_mag.cross(B))[2];
  if(alpha>0){
    beta=acos(B.dot(temp_mag));
  }
  else{
     beta=-acos(B.dot(temp_mag));
  }
  data[dsize*N+2].push_back(beta);
}
