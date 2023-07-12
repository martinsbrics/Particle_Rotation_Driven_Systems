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

Integrator::Integrator(int N_, double a_, std::vector<Eigen::Vector3d> &pos_, std::vector<Eigen::Vector3d> angle_, double B_mag_, double omega_B_,Eigen::Vector3d& init_magnetisation_ ):N(N_), a(a_), b(M_SQRT2*a_), pos(pos_), angle(angle_), HD_Force_and_Torque(N), init_magnetisation(N), magnetisation(N), rot(N), mag_Force(N), mag_Torque(N), ext_Torque(N), steric_Force(N), steric_Force1(N), steric_Torque(N),  perp(N),  perp1(N), TR(N), TL(N), BR(N), BL(N), i_perp(0.0,1.0,0.0), i_perp1(1.0,0.0,0.0), i_TR(0.5*a,0.5*a,0.0), i_TL(-0.5*a,0.5*a,0.0), i_BR(0.5*a,-0.5*a,0.0), i_BL(-0.5*a,-0.5*a,0.0),B_mag(B_mag_), omega_B(omega_B_), B(init_magnetisation_), dist_inside(8),steric_force_dir(8),torque_dist(8)
{
  std::cout<<"This newer should be used!!!!!!!"<<std::endl;abort();
  for(size_t i=0; i<N; i++)
  {
    init_magnetisation[i]=init_magnetisation_;
  }
  fill_support_vectors();
  phi0=std::atan2(init_magnetisation_[1],init_magnetisation_[0]);
  mom_angle=std::atan2(init_magnetisation_[1],init_magnetisation_[0]);
}

Integrator::Integrator(int N_, double omega_B_, Eigen::Vector3d& init_magnetisation_, double alpha_, double komega_kv_, double lag , bool B_off):N(N_),
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
mag_Force(N),
mag_Torque(N),
ext_Torque(N),
steric_Force(N),
steric_Force1(N),
steric_Torque(N),
steric_Torque1(N),
perp(N),
perp1(N),
TR(N),
TL(N),
BR(N),
BL(N),
i_perp(0.0,1.0,0.0),
i_perp1(1.0,0.0,0.0),
i_TR(0.5*a,0.5*a,0.0),
i_TL(-0.5*a,0.5*a,0.0),
i_BR(0.5*a,-0.5*a,0.0),
i_BL(-0.5*a,-0.5*a,0.0),
B_mag(B_off==true?1.0:0.0),
omega_B(omega_B_),
B(init_magnetisation_),
dist_inside(8),
steric_force_dir(8),
torque_dist(8)
{
  std::cout<<"initalizeing"<<std::endl;
  for(size_t i=0; i<N; i++)
  {
    init_magnetisation[i]=init_magnetisation_;
  }
  //fill_support_vectors();
  phi0=std::atan2(init_magnetisation_[1],init_magnetisation_[0])+lag;
  mom_angle=std::atan2(init_magnetisation_[1],init_magnetisation_[0]);
  std::cout<<"initalizeing here"<<std::endl;
}


// void Integrator::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t)
// {
//   dxdt=Eigen::VectorXd::Zero(6*N);
//  // std::cout<<"here"<<std::endl;
//   //remove this
//   for(size_t i=0; i<N; i++)
//   {
//     HD_Force_and_Torque[i][0]=-x[6*i+3];
//     HD_Force_and_Torque[i][1]=-x[6*i+4];
//     HD_Force_and_Torque[i][2]=-x[6*i+5];
//   }
// 
//   double inv_mass=6.0;
//   double inv_momentum_inertia=1.0;
//   double dipole_force_prefactor=inv_mass*3;
//   double dipole_torque_prefactor=0.0;
//   double steric_force_prefactor=inv_mass;
//   double steric_torque_prefactor=inv_momentum_inertia;
//   double ext_torque_prefactor=inv_momentum_inertia;
//   for(size_t i=0; i<N; i++)
//   {
//     pos[i][0]=x[6*i];
//     pos[i][1]=x[6*i+1];
//     angle[i]=x[6*i+2];
//   }
//    //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
//    fill_support_vectors();
//    //std::cout <<"here3"<<std::endl;
//    fill_mag_foreces_and_torques();
// //    std::cout <<"here4"<<std::endl;
// //    fill_steric_foreces_and_torques();
//    fill_steric_foreces();
// //   fill_steric_foreces_TB();
// // //   fill_steric_foreces_and_torques();
// //    if((steric_Force[0]-steric_Force1[0]).norm()>1.0e-13){
// //      std::cout<<"wrong "<<steric_Force1[0].transpose()<<" "<<steric_Force1[1].transpose()<<std::endl;
// //      std::cout<<"correct "<<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
// //      std::cout<<"mag "<<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
// //      std::cout<<"mag "<<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
// //      std::cout<<t<<std::endl;
// //      fill_steric_foreces_wp();
// //      //abort();
// //    }
// //   std::cout <<"here5"<<std::endl;
//    fill_external_foreces_and_torques(t);
// //   //std::cout <<"here6"<<std::endl;
//   for(size_t i=0; i<N; i++)
//   {
// //     std::cout <<"i="<<i<<" "<< dxdt.size()<<" "<<x.size()<<std::endl;
//     dxdt[6*i]=x[6*i+3];
//     dxdt[6*i+1]=x[6*i+4];
//     dxdt[6*i+2]=x[6*i+5];
//     dxdt[6*i+3]=HD_Force_and_Torque[i][0]*inv_mass+mag_Force[i][0]*dipole_force_prefactor+steric_Force[i][0]*steric_force_prefactor;
//     dxdt[6*i+4]=HD_Force_and_Torque[i][1]*inv_mass+mag_Force[i][1]*dipole_force_prefactor+steric_Force[i][1]*steric_force_prefactor;
//     dxdt[6*i+5]=HD_Force_and_Torque[i][2]*inv_momentum_inertia + mag_Torque[i]*dipole_torque_prefactor /*+ steric_Torque[i]*steric_torque_prefactor*/ + ext_torque_prefactor*ext_Torque[i];
//   }
//   //std::cout <<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
//   //std::cout <<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
// }

//creepy flow
void Integrator::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t)
{
  dxdt=Eigen::VectorXd::Zero(12*N);
  //std::cout<<"here"<<std::endl;
  //remove this

  double dipole_torque_prefactor=alpha/3.0;
  double dipole_force_prefactor=alpha/3.0*komega_kv;
  double steric_force_prefactor=dipole_force_prefactor;
  double steric_torque_prefactor=dipole_torque_prefactor/*4*/;
  double ext_torque_prefactor=1.0;
  for(size_t i=0; i<N; i++)
  {
    pos[i][0]=x[12*i];
    pos[i][1]=x[12*i+1];
    pos[i][2]=x[12*i+2];
    angle[i][0]=x[12*i+3];
    angle[i][1]=x[12*i+4];
    angle[i][2]=x[12*i+5];
  }
  std::cout<<angle[0][0]<<" "<<angle[0][1]<<" "<<angle[0][2]<<std::endl<<std::endl;
  
   //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
   fill_support_vectors();
   //std::cout <<"here3"<<std::endl;
   fill_mag_foreces_and_torques();
//    std::cout <<"here4"<<std::endl;
     //fill_steric_foreces_and_torques();
     fill_steric_foreces_and_torques_new();
     fill_steric_torques1();
//   fill_steric_foreces();
//   fill_steric_torques();
//   fill_steric_foreces_TB();
// //   fill_steric_foreces_and_torques();
//    if((steric_Force[0]-steric_Force1[0]).norm()>1.0e-13){
//      std::cout<<"wrong "<<steric_Force1[0].transpose()<<" "<<steric_Force1[1].transpose()<<std::endl;
//      std::cout<<"correct "<<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
//      std::cout<<"mag "<<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
//      std::cout<<"mag "<<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
//      std::cout<<t<<std::endl;
//      fill_steric_foreces_wp();
//      //abort();
//    }
//   std::cout <<"here5"<<std::endl;
   fill_external_foreces_and_torques(t);
//   //std::cout <<"here6"<<std::endl;
//   print_forces(x,dxdt,t);abort(); 
  for(size_t i=0; i<N; i++)
  {
//     std::cout <<"i="<<i<<" "<< dxdt.size()<<" "<<x.size()<<std::endl;
    dxdt[12*i]=mag_Force[i][0]*dipole_force_prefactor+steric_Force[i][0]*steric_force_prefactor;
    dxdt[12*i+1]=mag_Force[i][1]*dipole_force_prefactor+steric_Force[i][1]*steric_force_prefactor;
    dxdt[12*i+2]=mag_Force[i][2]*dipole_force_prefactor+steric_Force[i][2]*steric_force_prefactor;
    dxdt[12*i+3]=(mag_Torque[i][0]+steric_Torque[i][0]+steric_Torque1[i][0])*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i][0];;
    dxdt[12*i+4]=(mag_Torque[i][1]+steric_Torque[i][1]+steric_Torque1[i][1])*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i][1];
    dxdt[12*i+5]=(mag_Torque[i][2]+steric_Torque[i][2]+steric_Torque1[i][2])*dipole_torque_prefactor + ext_torque_prefactor*ext_Torque[i][2];
    dxdt[12*i+6]=0;
    dxdt[12*i+7]=0;
    dxdt[12*i+8]=0;
    dxdt[12*i+9]=0;
    dxdt[12*i+10]=0;
    dxdt[12*i+11]=0;
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
    rot[i] =Eigen::AngleAxisd(angle[i][2], Eigen::Vector3d::UnitZ()) 
    * Eigen::AngleAxisd(angle[i][1], Eigen::Vector3d::UnitY())
    *Eigen::AngleAxisd(angle[i][0], Eigen::Vector3d::UnitX()); 
    
    std::vector<Eigen::Matrix3d> rot1; rot1.resize(N);
    
    double yaw=angle[i][2];     //Yaw   angle (radians)
    double pitch=angle[i][1];   //Pitch angle (radians)
    double roll=angle[i][0];    //Roll  angle (radians)
    double su = sin(roll);
    double cu = cos(roll);
    double sv = sin(pitch);
    double cv = cos(pitch);
    double sw = sin(yaw);
    double cw = cos(yaw);
    
    rot1[i](0, 0) = cv*cw;
    rot1[i](0, 1) = su*sv*cw - cu*sw;
    rot1[i](0, 2) = su*sw + cu*sv*cw;
    rot1[i](1, 0) = cv*sw;
    rot1[i](1, 1) = cu*cw + su*sv*sw;
    rot1[i](1, 2) = cu*sv*sw - su*cw;
    rot1[i](2, 0) = -sv;
    rot1[i](2, 1) = su*cv;
    rot1[i](2, 2) = cu*cv;
    
    std::vector<Eigen::Matrix3d> rot2; rot2.resize(N);
    
    std::cout<<angle[0][0]<<" "<<angle[0][1]<<" "<<angle[0][2]<<std::endl<<std::endl;
    std::cout<<rot[0]<<std::endl<<std::endl;
    std::cout<<rot1[0]<<std::endl<<std::endl;
    abort();
    magnetisation[i]=rot[i]*init_magnetisation[i];
    perp[i]=rot[i]*i_perp;
    perp1[i]=rot[i]*i_perp1;
    TR[i]=rot[i]*i_TR+pos[i];
    TL[i]=rot[i]*i_TL+pos[i];
    BR[i]=rot[i]*i_BR+pos[i];
    BL[i]=rot[i]*i_BL+pos[i];
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

void Integrator::print_mag_foreces_and_torques()
{
  
  std::vector<Eigen::Vector3d> imag_Force;imag_Force.resize(N*N);
  for(size_t i=0; i<N*N; i++){
    imag_Force[i]=Eigen::Vector3d::Zero();
  }
  
  for(size_t i=0; i<N; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d rad=pos[j]-pos[i];
      double radius=rad.norm();
      double rad2=radius*radius; //rad^2
      double rad5=rad2*rad2;//rad^4
      rad5*=radius;//rad^
      Eigen::Vector3d F_ij=rad*magnetisation[i].dot(magnetisation[j]) +magnetisation[i]*rad.dot(magnetisation[j])+ magnetisation[j]*magnetisation[i].dot(rad)- (5.0/rad2*rad.dot(magnetisation[i])*rad.dot(magnetisation[j])) *rad;
      F_ij/=rad5;
      imag_Force[i*N+j]=-F_ij;
      imag_Force[j*N+i]=F_ij; 
      //_radius_vec[j+_number_of_particles*i]=_particle_pos[j]-_particle_pos[i];
    }
  }
  for(size_t i=0; i<N*N; i++){
    std::cout<< i<<" "<<imag_Force[i].transpose()<<std::endl;
  }
}


void Integrator::fill_external_foreces_and_torques(const double t)
{
  B[0]=B_mag*std::cos(t*omega_B+phi0);
  B[1]=B_mag*std::sin(t*omega_B+phi0);
  B[2]=0;
  for(size_t i=0; i<N; i++){
    ext_Torque[i]=magnetisation[i].cross(B);
  }
  //std::cout<<ext_Torque[0].transpose()<<" "<<t<<" "<<magnetisation[0].transpose()<<" "<<B.transpose()<<" "<<+phi0<<std::endl; abort();
}


double Integrator::is_point_inside_below(Eigen::Vector3d &point, int i){
  Eigen::Vector3d AB, AD, AM;
  AB=BR[i]-BL[i];
  AD=TL[i]-BL[i];
  AM=point-BL[i];
  double v1=AM.dot(AB);
  double v2=AM.dot(AD);
  if((v1>-1.0e-13&&v1<a*a+1.0e-13)&&(v2>-1.0e-13&&v2<a*a+1.0e-13))
  {
    return fabs((BL[i][0]-BR[i][0])*(BR[i][1]-point[1])-(BR[i][0]-point[0])*(BL[i][1]-BR[i][1]))/a;
  }
  else
  {
    return -1.0;
  }
};
double Integrator::is_point_inside_above(Eigen::Vector3d &point, int i){
  Eigen::Vector3d AB, AD, AM;
  AB=BR[i]-BL[i];
  AD=TL[i]-BL[i];
  AM=point-BL[i];
  double v1=AM.dot(AB);
  double v2=AM.dot(AD);
  if((v1>-1.0e-13&&v1<a*a+1.0e-13)&&(v2>-1.0e-13&&v2<a*a+1.0e-13))
  {
    return fabs((TL[i][0]-TR[i][0])*(TR[i][1]-point[1])-(TR[i][0]-point[0])*(TL[i][1]-TR[i][1]))/a;
  }
  else
  {
    return -1.0;
  }
};

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

double Integrator::is_point_inside(Eigen::Vector3d &point,Eigen::Vector3d &pointp,Eigen::Vector3d &pointp1, int i, Eigen::Vector3d &force_dir, double &torque_distance){
  Eigen::Vector3d AB, AD, AM;
  AB=BR[i]-BL[i];
  AD=TL[i]-BL[i];
  AM=point-BL[i];
  dvec dist(4);
  std::vector<size_t> indexes;
  std::vector<Eigen::Vector3d> f_dir(4);
  double v1=AM.dot(AB);
  double v2=AM.dot(AD);
  int pos_it;
  //std::cout <<" v1, v2="<<v1<<" "<<v2<<" "<<point.transpose()<<" "<<AM.transpose()<<" "<<AB.transpose()<<" "<<AD.transpose()<<std::endl;
  if((v1>-1.0e-13&&v1<a*a+1.0e-13)&&(v2>-1.0e-13&&v2<a*a+1.0e-13)){
    dist[0]=fabs((TL[i][0]-TR[i][0])*(TR[i][1]-point[1])-(TR[i][0]-point[0])*(TL[i][1]-TR[i][1]))/a;
    f_dir[0]=perp[i];
    dist[1]=fabs((BR[i][0]-TR[i][0])*(TR[i][1]-point[1])-(TR[i][0]-point[0])*(BR[i][1]-TR[i][1]))/a;
    f_dir[1]=perp1[i];
    dist[2]=fabs((BR[i][0]-BL[i][0])*(BL[i][1]-point[1])-(BL[i][0]-point[0])*(BR[i][1]-BL[i][1]))/a;
    f_dir[2]=-perp[i];
    dist[3]=fabs((TL[i][0]-BL[i][0])*(BL[i][1]-point[1])-(BL[i][0]-point[0])*(TL[i][1]-BL[i][1]))/a;
    f_dir[3]=-perp1[i];
    indexes=sort_indexes(dist);
    if(dist[indexes[1]]>0.02){
      if(indexes[0]%2==0){
        force_dir=f_dir[indexes[0]];
        if(f_dir[1].dot(pointp1-point)<0){
          torque_distance=dist[1];
        }
        else{
          torque_distance=dist[3];
        }
        return dist[indexes[0]];
      }
      else{
        force_dir=f_dir[indexes[0]];
        if(f_dir[0].dot(pointp-point)<0){
          torque_distance=dist[0];
        }
        else{
          torque_distance=dist[2];
        }
        return dist[indexes[0]];
      }
    }
    else{
      if(indexes[0]%2==0){
        force_dir=f_dir[indexes[0]];
        if(f_dir[1].dot(pointp1-point)<0){
          torque_distance=dist[1];
        }
        else{
          torque_distance=dist[3];
        }
        return dist[indexes[0]];
      }
      else{
        force_dir=f_dir[indexes[1]];
        if(f_dir[1].dot(pointp1-point)<0){
          torque_distance=dist[1];
        }
        else{
          torque_distance=dist[3];
        }
        return dist[indexes[1]];
      }
    }
  }
//     pos_it=min_element(dist.begin(), dist.end())-dist.begin();
//     if(pos_it%2==0){
//       if(f_dir[pos_it].dot(pointp-point)>0){
//         force_dir=f_dir[pos_it];
//         if(f_dir[1].dot(pointp1-point)<0){
//           torque_distance=dist[1];
//         }
//         else{
//           torque_distance=dist[3];
//         }
//         return dist[pos_it];
//       }
//       else{
//         dist[pos_it]=2;
//         pos_it=min_element(dist.begin(), dist.end())-dist.begin();
//         force_dir=f_dir[pos_it];
//         if(f_dir[0].dot(pointp-point)<0){
//           torque_distance=dist[0];
//         }
//         else{
//           torque_distance=dist[2];
//         }
//         return dist[pos_it];
//       }
//     }
//     else{
//       if(f_dir[pos_it].dot(pointp1-point)>0){
//         force_dir=f_dir[pos_it];
//         if(f_dir[0].dot(pointp-point)<0){
//           torque_distance=dist[0];
//         }
//         else{
//           torque_distance=dist[2];
//         }
//         return dist[pos_it];
//       }
//       else{
//         dist[pos_it]=2;
//         pos_it=min_element(dist.begin(), dist.end())-dist.begin();
//         force_dir=f_dir[pos_it];
//         if(f_dir[1].dot(pointp1-point)<0){
//           torque_distance=dist[1];
//         }
//         else{
//           torque_distance=dist[3];
//         }
//         return dist[pos_it];
//       }
//     }
//   }
  else{
    return -1.0;
  }
};

double Integrator::is_point_inside_TB(Eigen::Vector3d &point, int i, Eigen::Vector3d &force_dir){
  Eigen::Vector3d AB, AD, AM;
  AB=BR[i]-BL[i];
  AD=TL[i]-BL[i];
  AM=point-BL[i];
  dvec dist(2);
  std::vector<Eigen::Vector3d> f_dir(2);
  double v1=AM.dot(AB);
  double v2=AM.dot(AD);
  int pos_it;
  //std::cout <<" v1, v2="<<v1<<" "<<v2<<" "<<point.transpose()<<" "<<AM.transpose()<<" "<<AB.transpose()<<" "<<AD.transpose()<<std::endl;
  if((v1>-1.0e-13&&v1<a*a+1.0e-13)&&(v2>-1.0e-13&&v2<a*a+1.0e-13))
  {
    dist[0]=fabs((TL[i][0]-TR[i][0])*(TR[i][1]-point[1])-(TR[i][0]-point[0])*(TL[i][1]-TR[i][1]))/a;
    f_dir[0]=perp[i];
    dist[1]=fabs((BR[i][0]-BL[i][0])*(BL[i][1]-point[1])-(BL[i][0]-point[0])*(BR[i][1]-BL[i][1]))/a;
    f_dir[1]=-perp[i];
    pos_it=min_element(dist.begin(), dist.end())-dist.begin();
    force_dir=f_dir[pos_it];
    return dist[pos_it];
  }
  else
  {
    return -1.0;
  }
};


double Integrator::is_point_inside_wp(Eigen::Vector3d &point,Eigen::Vector3d &pointp,Eigen::Vector3d &pointp1, int i, Eigen::Vector3d &force_dir){
  Eigen::Vector3d AB, AD, AM;
  AB=BR[i]-BL[i];
  AD=TL[i]-BL[i];
  AM=point-BL[i];
  dvec dist(4);
  std::vector<Eigen::Vector3d> f_dir(4);
  double v1=AM.dot(AB);
  double v2=AM.dot(AD);
  int pos_it;
  std::cout <<" v1, v2="<<v1<<" "<<v2<<" "<<point.transpose()<<" "<<AM.transpose()<<" "<<AB.transpose()<<" "<<AD.transpose()<<std::endl;
  if((v1>-1.0e-13&&v1<a*a+1.0e-13)&&(v2>-1.0e-13&&v2<a*a+1.0e-13))
  {
    dist[0]=fabs((TL[i][0]-TR[i][0])*(TR[i][1]-point[1])-(TR[i][0]-point[0])*(TL[i][1]-TR[i][1]))/a;
    f_dir[0]=perp[i];
    dist[1]=fabs((BR[i][0]-TR[i][0])*(TR[i][1]-point[1])-(TR[i][0]-point[0])*(BR[i][1]-TR[i][1]))/a;
    f_dir[1]=perp1[i];
    dist[2]=fabs((BR[i][0]-BL[i][0])*(BL[i][1]-point[1])-(BL[i][0]-point[0])*(BR[i][1]-BL[i][1]))/a;
    f_dir[2]=-perp[i];
    dist[3]=fabs((TL[i][0]-BL[i][0])*(BL[i][1]-point[1])-(BL[i][0]-point[0])*(TL[i][1]-BL[i][1]))/a;
    f_dir[3]=-perp1[i];
    pos_it=min_element(dist.begin(), dist.end())-dist.begin();
    if(pos_it%2==0){
      if(f_dir[pos_it].dot(pointp-point)>0){
        force_dir=f_dir[pos_it];
        return dist[pos_it];
      }
      else{
        std::cout <<dist[0]<<" "<<dist[1]<<" "<<dist[2]<<" "<<dist[3]<<std::endl;
        dist[pos_it]=2;
        pos_it=min_element(dist.begin(), dist.end())-dist.begin();
        force_dir=f_dir[pos_it];
        std::cout <<dist[0]<<" "<<dist[1]<<" "<<dist[2]<<" "<<dist[3]<<std::endl;
        return dist[pos_it];
      }
    }
    else{
      if(f_dir[pos_it].dot(pointp1-point)>0){
        force_dir=f_dir[pos_it];
        return dist[pos_it];
      }
      else{
        std::cout <<dist[0]<<" "<<dist[1]<<" "<<dist[2]<<" "<<dist[3]<<std::endl;
        dist[pos_it]=2;
        pos_it=min_element(dist.begin(), dist.end())-dist.begin();
        force_dir=f_dir[pos_it];
        std::cout <<dist[0]<<" "<<dist[1]<<" "<<dist[2]<<" "<<dist[3]<<std::endl;
        return dist[pos_it];
      }
    }
  }
  else
  {
    return -1.0;
  }
};

void Integrator::fill_steric_foreces_and_torques_old()
{
  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
    steric_Torque[i]=Eigen::Vector3d::Zero();
  }
  for(size_t i=0; i<N-1; i++){
    size_t j=i+1;
    double dx11=is_point_inside_below(TL[i],j); double dx12=is_point_inside_below(TR[i],j);
    double dx1=std::max(dx11,dx12);
    double dx21=is_point_inside_above(BL[j],i); double dx22=is_point_inside_above(BR[j],i);
    double dx2=std::max(dx21,dx22);
   // cout <<"dx1="<<dx1<<" dx2="<<dx2<<endl;
    if(std::max(dx1,dx2)>0)
    {
      //cout <<"moving centrs!!!!"<<endl;
      if(dx1>dx2){
        double coef=1.0/(1.0-2*dx1/a);
        double coef2=coef*coef;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*coef;
        steric_Force[i]-=coef13*perp[j];
        steric_Force[j]+=coef13*perp[j];   
        if(fmod(angle[i][2]-angle[j][2], 2*M_PI)>0.0001){
//        if(false){
          if(dx11>dx12){
            steric_Torque[i]-=i_TL.cross(perp[j])*coef13;
            steric_Torque[j]+=(TL[i]-pos[j]).cross(perp[j])*coef13;
            
          }
          else{
            steric_Torque[i]-=i_TR.cross(perp[j])*coef13;
            steric_Torque[j]+=(TR[i]-pos[j]).cross(perp[j])*coef13;
          }
            
        }
      }
      else{
        double coef=1.0/(1.0-2*dx2/a);
        double coef2=coef*coef;
        double coef4=coef2*coef2;
        double coef8=coef4*coef4;
        double coef13=coef8*coef4*coef;
        steric_Force[i]-=coef13*perp[i];
        steric_Force[j]+=coef13*perp[i];
        if(fmod(angle[i][2]-angle[j][2], 2*M_PI)>0.0001){
//        if(false){
          if(dx21>dx22){
            steric_Torque[i]-=(BL[j]-pos[i]).cross(perp[i])*coef13;
            steric_Torque[j]+=i_BL.cross(perp[i])*coef13;
            
          }
          else{
            steric_Torque[i]-=(BR[j]-pos[i]).cross(perp[i])*coef13;
            steric_Torque[j]+=i_BR.cross(perp[i])*coef13;
          }
        }
      }
    }
    else
    {
      //do 
    }
  }
}

void Integrator::fill_steric_foreces()
{
/*  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
  }
  Eigen::Vector3d force_dir=Eigen::Vector3d::Zero();
  size_t pos_it;
  double dx=-1.0;
  double temp;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      if(ScalProd(pos[i]-pos[j],pos[i]-pos[j])<2*a*a){
        dist_inside[0]=is_point_inside(TR[j],BR[j],TL[j],i,steric_force_dir[0],temp);
        dist_inside[1]=is_point_inside(TL[j],BL[j],TR[j],i,steric_force_dir[1],temp);
        dist_inside[2]=is_point_inside(BL[j],TL[j],BR[j],i,steric_force_dir[2],temp);
        dist_inside[3]=is_point_inside(BR[j],TR[j],BL[j],i,steric_force_dir[3],temp);
        dist_inside[4]=is_point_inside(TR[i],BR[i],TL[i],j,steric_force_dir[4],temp);
        dist_inside[5]=is_point_inside(TL[i],BL[i],TR[i],j,steric_force_dir[5],temp);
        dist_inside[6]=is_point_inside(BL[i],TL[i],BR[i],j,steric_force_dir[6],temp);
        dist_inside[7]=is_point_inside(BR[i],TR[i],BL[i],j,steric_force_dir[7],temp);
        pos_it=max_element(dist_inside.begin(), dist_inside.end())- dist_inside.begin();
        if(pos_it>3)
        {
          force_dir=-steric_force_dir[pos_it];
        }
        else
        {
          force_dir=steric_force_dir[pos_it];
        }
        dx=dist_inside[pos_it];
//         std::cout <<pos_it<<" "<<dx<<" "<<force_dir.transpose()<<" "<<*max_element(dist_inside.begin(), dist_inside.end())<<std::endl;
//         for(int i=0; i<8; i++)
//         {
//           std::cout <<dist_inside[i]<<std::endl;
//         }
//         abort();
        if(dx>0){
          double coef=1.0/(1.0-2*dx/a);
          double coef2=coef*coef;
          double coef4=coef2*coef2;
          double coef8=coef4*coef4;
          double coef13=coef8*coef4*coef;
          steric_Force[i]-=coef13*force_dir;
          steric_Force[j]+=coef13*force_dir;  
        }
        else{
          //DO
        }
      }
    }
  }
     */  
}

void Integrator::fill_steric_foreces_and_torques()
{
/*  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
    steric_Torque[i]=0.0;
  }
  Eigen::Vector3d force_dir=Eigen::Vector3d::Zero();
  size_t pos_it;
  double dx=-1.0;
  double torque_dir=1.0, temp_torque;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      if(ScalProd(pos[i]-pos[j],pos[i]-pos[j])<2*a*a){
        dist_inside[0]=is_point_inside(TR[j],BR[j],TL[j],i,steric_force_dir[0], torque_dist[0]);
        dist_inside[1]=is_point_inside(TL[j],BL[j],TR[j],i,steric_force_dir[1], torque_dist[1]);
        dist_inside[2]=is_point_inside(BL[j],TL[j],BR[j],i,steric_force_dir[2], torque_dist[2]);
        dist_inside[3]=is_point_inside(BR[j],TR[j],BL[j],i,steric_force_dir[3], torque_dist[3]);
        dist_inside[4]=is_point_inside(TR[i],BR[i],TL[i],j,steric_force_dir[4], torque_dist[4]);
        dist_inside[5]=is_point_inside(TL[i],BL[i],TR[i],j,steric_force_dir[5], torque_dist[5]);
        dist_inside[6]=is_point_inside(BL[i],TL[i],BR[i],j,steric_force_dir[6], torque_dist[6]);
        dist_inside[7]=is_point_inside(BR[i],TR[i],BL[i],j,steric_force_dir[7], torque_dist[7]);
        pos_it=max_element(dist_inside.begin(), dist_inside.end())- dist_inside.begin();
        if(pos_it>3){
          force_dir=-steric_force_dir[pos_it];
        }
        else{
          force_dir=steric_force_dir[pos_it];
        }
        if(pos_it%2==0){
          torque_dir=-1.0;        
        }
        else{
          torque_dir=1.0;
        }
        dx=dist_inside[pos_it];
//         std::cout <<pos_it<<" "<<dx<<" "<<force_dir.transpose()<<" "<<*max_element(dist_inside.begin(), dist_inside.end())<<std::endl;
//         for(int i=0; i<8; i++)
//         {
//           std::cout <<dist_inside[i]<<std::endl;
//         }
//         abort();
        if(dx>0){
          double coef=1.0/(1.0-4*dx/a);
          double coef2=coef*coef;
          double coef4=coef2*coef2;
          double coef8=coef4*coef4;
          double coef13=coef8*coef4*coef;
          //double coef7=coef4*coef2*coef;
          coef13-=1;
          steric_Force[i]-=coef13*force_dir;
          steric_Force[j]+=coef13*force_dir;
          //steric_Torque[i]-=crossProd(pos[j]-pos[i],force_dir)/ScalProd(pos[j]-pos[i],force_dir)*coef13*1.5*0.984809;
          //steric_Torque[j]-=crossProd(pos[j]-pos[i],force_dir)/ScalProd(pos[j]-pos[i],force_dir)*coef13*1.5*0.984809;
//          steric_Torque[i]-=((pos[j][0]-pos[i][0])*std::cos(angle[i])+(pos[j][1]-pos[i][1])*std::sin(angle[i]))*coef13*1.5;
//          steric_Torque[j]-=((pos[j][0]-pos[i][0])*std::cos(angle[i])+(pos[j][1]-pos[i][1])*std::sin(angle[i]))*coef13*1.5;
          temp_torque=torque_dir*torque_dist[pos_it]*coef13*1.5;
          steric_Torque[i]+=temp_torque;
          steric_Torque[j]+=temp_torque;
        }
        else{
          //DO
        }
      }
    }
  }
    */   
}


void Integrator::fill_steric_foreces_and_torques_new()
{
  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
    steric_Torque[i]=Eigen::Vector3d::Zero();
  }
  Eigen::Vector3d force_dir=Eigen::Vector3d::Zero();
  size_t pos_it;
  double dx=-1.0;
  double torque_dir=1.0; Eigen::Vector3d temp_torque;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      Eigen::Vector3d rij=pos[j]-pos[i];
      if(rij.squaredNorm()<3*a*a){
        double yrij= perp[i].dot(rij);
        double xrij= perp1[i].dot(rij);
        Eigen::Vector3d vyrij= rij.cross(perp[i]);
        Eigen::Vector3d vxrij= rij.cross(perp1[i]);
        double ay=a-fabs(yrij);
        double ax=b-fabs(xrij);
        //std::cout <<"From the new: "<<ax<<" "<<ay<<std::endl;
        //double mina,maxa;
        if(ax>ay){
          if(ay>0){
            double coef=1.0/(1.0-4*ay/a);
            double coef2=coef*coef;
            double coef4=coef2*coef2;
            double coef8=coef4*coef4;
            double coef13=coef8*coef4*coef;
            coef13-=1;
            force_dir=(signum(yrij)*coef13)*perp[i];
            steric_Force[i]-=force_dir;
            steric_Force[j]+=force_dir;
            temp_torque=vyrij*signum(yrij)*coef13*1.5;
            steric_Torque[i]-=temp_torque;
            steric_Torque[j]-=temp_torque;
          }
          else{
            //do something;
          }
        }
        else
        {
          if(ax>0){
            if(ay<0.02){
              double coef=1.0/(1.0-4*ay/a);
              double coef2=coef*coef;
              double coef4=coef2*coef2;
              double coef8=coef4*coef4;
              double coef13=coef8*coef4*coef;
              coef13-=1;
              force_dir=(signum(yrij)*coef13)*perp[i];
              steric_Force[i]-=force_dir;
              steric_Force[j]+=force_dir;
              temp_torque=vyrij*coef13*1.5*signum(yrij);
              steric_Torque[i]-=temp_torque;
              steric_Torque[j]-=temp_torque;
            }
            else{
              //std::cout <<"Am I here "<<ax<<" "<<ay<<std::endl;
              double coef=1.0/(1.0-4*ax/a);
              double coef2=coef*coef;
              double coef4=coef2*coef2;
              double coef8=coef4*coef4;
              double coef13=coef8*coef4*coef;
              coef13-=1;
              force_dir=(signum(xrij)*coef13)*perp1[i];
              steric_Force[i]-=force_dir;
              steric_Force[j]+=force_dir;
              temp_torque=vxrij*coef13*1.5*signum(xrij);
              steric_Torque[i]-=temp_torque;
              steric_Torque[j]-=temp_torque;
            }
          }
          else{
            //do something;
          }
        }
      }
    }
  }
       
}


void Integrator::fill_steric_foreces_and_torques_print()
{
//   std::cout<<"Printing fill_steric_foreces_and_torques_prin()"<<std::endl;
//   for(size_t i=0; i<N; i++){
//     steric_Force[i]=Eigen::Vector3d::Zero();
//     steric_Torque[i]=0.0;
//   }
//   Eigen::Vector3d force_dir=Eigen::Vector3d::Zero();
//   size_t pos_it;
//   double dx=-1.0;
//   double torque_dir=1.0, temp_torque;
//   for(size_t i=0; i<N-1; i++){
//     for(size_t j=i+1; j<N; j++){
//       if(ScalProd(pos[i]-pos[j],pos[i]-pos[j])<2*a*a){
//         dist_inside[0]=is_point_inside(TR[j],BR[j],TL[j],i,steric_force_dir[0], torque_dist[0]);
//         dist_inside[1]=is_point_inside(TL[j],BL[j],TR[j],i,steric_force_dir[1], torque_dist[1]);
//         dist_inside[2]=is_point_inside(BL[j],TL[j],BR[j],i,steric_force_dir[2], torque_dist[2]);
//         dist_inside[3]=is_point_inside(BR[j],TR[j],BL[j],i,steric_force_dir[3], torque_dist[3]);
//         dist_inside[4]=is_point_inside(TR[i],BR[i],TL[i],j,steric_force_dir[4], torque_dist[4]);
//         dist_inside[5]=is_point_inside(TL[i],BL[i],TR[i],j,steric_force_dir[5], torque_dist[5]);
//         dist_inside[6]=is_point_inside(BL[i],TL[i],BR[i],j,steric_force_dir[6], torque_dist[6]);
//         dist_inside[7]=is_point_inside(BR[i],TR[i],BL[i],j,steric_force_dir[7], torque_dist[7]);
//         pos_it=max_element(dist_inside.begin(), dist_inside.end())- dist_inside.begin();
//         if(pos_it>3){
//           force_dir=-steric_force_dir[pos_it];
//         }
//         else{
//           force_dir=steric_force_dir[pos_it];
//         }
//         if(pos_it%2==0){
//           torque_dir=-1.0;        
//         }
//         else{
//           torque_dir=1.0;
//         }         
//         dx=dist_inside[pos_it];
// //         std::cout <<pos_it<<" "<<dx<<" "<<force_dir.transpose()<<" "<<*max_element(dist_inside.begin(), dist_inside.end())<<std::endl;
// //         for(int i=0; i<8; i++)
// //         {
// //           std::cout <<dist_inside[i]<<std::endl;
// //         }
// //         abort();
//         if(dx>0){
//           std::cout <<"From the old dx=: "<<dx<<std::endl;
//           double coef=1.0/(1.0-4*dx/a);
//           double coef2=coef*coef;
//           double coef4=coef2*coef2;
//           double coef8=coef4*coef4;
//           double coef13=coef8*coef4*coef;
//           coef13-=1;
//           steric_Force[i]-=coef13*force_dir;
//           steric_Force[j]+=coef13*force_dir;
//           //steric_Torque[i]-=crossProd(pos[j]-pos[i],force_dir)/ScalProd(pos[j]-pos[i],force_dir)*coef13*1.5*0.984809;
//           //steric_Torque[j]-=crossProd(pos[j]-pos[i],force_dir)/ScalProd(pos[j]-pos[i],force_dir)*coef13*1.5*0.984809;
// //          steric_Torque[i]-=((pos[j][0]-pos[i][0])*std::cos(angle[i])+(pos[j][1]-pos[i][1])*std::sin(angle[i]))*coef13*1.5;
// //          steric_Torque[j]-=((pos[j][0]-pos[i][0])*std::cos(angle[i])+(pos[j][1]-pos[i][1])*std::sin(angle[i]))*coef13*1.5;
//           temp_torque=torque_dir*torque_dist[pos_it]*coef13*1.5;
//           steric_Torque[i]+=temp_torque;
//           steric_Torque[j]+=temp_torque;
//           std::cout<<i<<" "<<j<<" "<<torque_dir <<" "<<pos_it<<std::endl;
//         }
//         else{
//           //DO
//         }
//       }
//     }
//   }
       
}


void Integrator::fill_steric_torques1()
{
  Eigen::Vector3d avg_torque=Eigen::Vector3d::Zero();
  for(size_t i=0; i<N; i++){
    avg_torque+=mag_Torque[i]+steric_Torque[i];
  }
  avg_torque/=N;
  for(size_t i=0; i<N; i++){
    steric_Torque1[i]=avg_torque-mag_Torque[i]-steric_Torque[i];
  }
       
}

void Integrator::fill_steric_torques()
{
//   for(size_t i=0; i<N; i++){
//     steric_Torque[i]=0.0;
//   }
//   
//   for(size_t i=0; i<N; i++){
//     for(size_t j=i+1; j<N; j++){
//       Eigen::Vector3d rad=pos[j]-pos[i];
//       double radius=rad.norm();
//       double rad2=radius*radius; //rad^2
//       double rad5=rad2*rad2;//rad^4
//       rad5*=radius;//rad^5
//       double tau_ij=1.5*ScalProd(magnetisation[i],rad)*crossProd(magnetisation[j],rad)+1.5*ScalProd(magnetisation[j],rad)*crossProd(magnetisation[i],rad);
//       steric_Torque[i]-=tau_ij/rad5;
//       steric_Torque[j]-=tau_ij/rad5;
//       //_radius_vec[j+_number_of_particles*i]=_particle_pos[j]-_particle_pos[i];
//     }
//   }
}


void Integrator::fill_steric_foreces_TB()
{
/*  for(size_t i=0; i<N; i++){
    steric_Force[i]=Eigen::Vector3d::Zero();
  }
  Eigen::Vector3d force_dir=Eigen::Vector3d::Zero();
  size_t pos_it;
  double dx=-1.0;
  for(size_t i=0; i<N-1; i++){
    for(size_t j=i+1; j<N; j++){
      if(ScalProd(pos[i]-pos[j],pos[i]-pos[j])<2*a*a){
        dist_inside[0]=is_point_inside_TB(TR[j],i,steric_force_dir[0]);
        dist_inside[1]=is_point_inside_TB(TL[j],i,steric_force_dir[1]);
        dist_inside[2]=is_point_inside_TB(BL[j],i,steric_force_dir[2]);
        dist_inside[3]=is_point_inside_TB(BR[j],i,steric_force_dir[3]);
        dist_inside[4]=is_point_inside_TB(TR[i],j,steric_force_dir[4]);
        dist_inside[5]=is_point_inside_TB(TL[i],j,steric_force_dir[5]);
        dist_inside[6]=is_point_inside_TB(BL[i],j,steric_force_dir[6]);
        dist_inside[7]=is_point_inside_TB(BR[i],j,steric_force_dir[7]);
        pos_it=max_element(dist_inside.begin(), dist_inside.end())- dist_inside.begin();
        if(pos_it>3)
        {
          force_dir=-steric_force_dir[pos_it];
        }
        else
        {
          force_dir=steric_force_dir[pos_it];
        }
        dx=dist_inside[pos_it];
//         std::cout <<pos_it<<" "<<dx<<" "<<force_dir.transpose()<<" "<<*max_element(dist_inside.begin(), dist_inside.end())<<std::endl;
//         for(int i=0; i<8; i++)
//         {
//           std::cout <<dist_inside[i]<<std::endl;
//         }
//         abort();
        if(dx>0){
          double coef=1.0/(1.0-2*dx/a);
          double coef2=coef*coef;
          double coef4=coef2*coef2;
          double coef8=coef4*coef4;
          double coef13=coef8*coef4*coef;
          steric_Force[i]-=coef13*force_dir;
          steric_Force[j]+=coef13*force_dir;  
        }
        else{
          //DO
        }
      }
    }
  }
    */   
}

void Integrator::fill_steric_foreces_wp()
{
//   for(size_t i=0; i<N; i++){
//     steric_Force1[i]=Eigen::Vector3d::Zero();
//   }
//   Eigen::Vector3d force_dir=Eigen::Vector3d::Zero();
//   size_t pos_it;
//   double dx=-1.0;
//   for(size_t i=0; i<N-1; i++){
//     for(size_t j=i+1; j<N; j++){
//       if(ScalProd(pos[i]-pos[j],pos[i]-pos[j])<2*a*a){
//         dist_inside[0]=is_point_inside_wp(TR[j],BR[j],TL[j],i,steric_force_dir[0]);
//         dist_inside[1]=is_point_inside_wp(TL[j],BL[j],TR[j],i,steric_force_dir[1]);
//         dist_inside[2]=is_point_inside_wp(BL[j],TL[j],BR[j],i,steric_force_dir[2]);
//         dist_inside[3]=is_point_inside_wp(BR[j],TR[j],BL[j],i,steric_force_dir[3]);
//         dist_inside[4]=is_point_inside_wp(TR[i],BR[i],TL[i],j,steric_force_dir[4]);
//         dist_inside[5]=is_point_inside_wp(TL[i],BL[i],TR[i],j,steric_force_dir[5]);
//         dist_inside[6]=is_point_inside_wp(BL[i],TL[i],BR[i],j,steric_force_dir[6]);
//         dist_inside[7]=is_point_inside_wp(BR[i],TR[i],BL[i],j,steric_force_dir[7]);
//         pos_it=max_element(dist_inside.begin(), dist_inside.end())- dist_inside.begin();
//         if(pos_it>3)
//         {
//           force_dir=-steric_force_dir[pos_it];
//         }
//         else
//         {
//           force_dir=steric_force_dir[pos_it];
//         }
//         dx=dist_inside[pos_it];
//         std::cout <<pos_it<<" "<<dx<<" "<<force_dir.transpose()<<" "<<*max_element(dist_inside.begin(), dist_inside.end())<<std::endl;
//         for(int i=0; i<8; i++)
//         {
//           std::cout <<dist_inside[i]<<" "<<steric_force_dir[i].transpose()<<std::endl;
//         }
//         std::cout <<perp[0].transpose()<<" "<<perp1[0].transpose()<<std::endl;
//         //abort();
//         if(dx>0){
//           double coef=1.0/(1.0-2*dx/a);
//           double coef2=coef*coef;
//           double coef4=coef2*coef2;
//           double coef8=coef4*coef4;
//           double coef13=coef8*coef4*coef;
//           steric_Force1[i]-=coef13*force_dir;
//           steric_Force1[j]+=coef13*force_dir;  
//         }
//         else{
//           //DO
//         }
//       }
//     }
//   }
       
}


void Integrator::print_forces(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t)
{
//   dxdt=Eigen::VectorXd::Zero(6*N);
//   std::vector<Eigen::Vector3d> HD_Force;HD_Force.resize(N);
//  // std::cout<<"here"<<std::endl;
//   //remove this
//   for(size_t i=0; i<N; i++)
//   {
//     HD_Force_and_Torque[i][0]=-x[6*i+3];
//     HD_Force[i][0]=-x[6*i+3];
//     HD_Force_and_Torque[i][1]=-x[6*i+4];
//     HD_Force[i][1]=-x[6*i+4];
//     HD_Force_and_Torque[i][2]=-x[6*i+5];
//   }
// //   double inv_mass=6.0;
// //   double inv_momentum_inertia=1.0;
// //   double dipole_force_prefactor=6.0*3;
// //   double dipole_torque_prefactor=1.0*3;
// //   double steric_force_prefactor=6.0;
// //   double steric_torque_prefactor=1.0;
// //   double ext_torque_prefactor=1.0;
//   for(size_t i=0; i<N; i++)
//   {
//     pos[i][0]=x[6*i];
//     pos[i][1]=x[6*i+1];
//     angle[i]=x[6*i+2];
//   }
//    //std::cout <<pos[0].transpose()<<" "<<pos[1].transpose()<<" "<<x.transpose()<<std::endl;
//    fill_support_vectors();
//    //std::cout <<"here3"<<std::endl;
//    fill_mag_foreces_and_torques();
// //    std::cout <<"here4"<<std::endl;
//     fill_steric_foreces_and_torques_print();
// //   fill_steric_foreces();
// //   fill_steric_torques();
// //   fill_steric_foreces_TB();
// // //   fill_steric_foreces_and_torques();
// //    if((steric_Force[0]-steric_Force1[0]).norm()>1.0e-13){
// //      std::cout<<"wrong "<<steric_Force1[0].transpose()<<" "<<steric_Force1[1].transpose()<<std::endl;
// //      std::cout<<"correct "<<steric_Force[0].transpose()<<" "<<steric_Force[1].transpose()<<std::endl;
// //      std::cout<<"mag "<<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
// //      std::cout<<"mag "<<mag_Force[0].transpose()<<" "<<mag_Force[1].transpose()<<std::endl;
// //      std::cout<<t<<std::endl;
// //      fill_steric_foreces_wp();
// //      //abort();
// //    }
// //   std::cout <<"here5"<<std::endl;
//    fill_external_foreces_and_torques(t);
//    for(size_t i=0; i<N; i++){
//      double cosa=std::cos(angle[i]);
//      double sina=std::sin(-angle[i]);
//      rot[i](0,0)=cosa;
//      rot[i](0,1)=-sina;
//      rot[i](1,0)=sina;
//      rot[i](1,1)=cosa;
// //      HD_Force[i]=rot[i]*HD_Force[i];
// //      mag_Force[i]=rot[i]*mag_Force[i];
// //      steric_Force[i]=rot[i]*steric_Force[i];
//      //pos[i]-=pos[0];
//      //pos[i]=rot[0]*pos[i];
//   }
//   std::cout<<"Forces: HD, Mag, Steric, po, angle"<<std::endl;
//   for(size_t i=0; i<N; i++){
//     std::cout<<(rot[i]*HD_Force[i]).transpose()<<"  "<<(rot[i]*mag_Force[i]).transpose()<<"  "<<(rot[i]*steric_Force[i]).transpose()<<"  "<<(rot[i]*(pos[i]-pos[i-1>0?i-1:0])).transpose()<<"  "<<-fmod(angle[i],2*M_PI)<<"  "<<angle[i]-get_phi(t)+phi0<<std::endl; 
//   }
//   std::cout<<"Torques:  Mag, Steric"<<std::endl;
//   for(size_t i=0; i<N; i++){
//     std::cout<<mag_Torque[i]<<"  "<<steric_Torque[i]<<" "<<mag_Torque[i]+steric_Torque[i]<<"  "<<ext_Torque[i]<<std::endl; 
//   }
//   fill_steric_foreces_and_torques_new();
//   std::cout<<"Forces: HD, Mag, Steric, po, angle"<<std::endl;
//   for(size_t i=0; i<N; i++){
//     std::cout<<(rot[i]*HD_Force[i]).transpose()<<"  "<<(rot[i]*mag_Force[i]).transpose()<<"  "<<(rot[i]*steric_Force[i]).transpose()<<"  "<<(rot[i]*(pos[i]-pos[i-1>0?i-1:0])).transpose()<<"  "<<-fmod(angle[i],2*M_PI)<<"  "<<angle[i]-get_phi(t)+phi0<<std::endl; 
//   }
//   std::cout<<"Torques:  Mag, Steric"<<std::endl;
//   for(size_t i=0; i<N; i++){
//     std::cout<<mag_Torque[i]<<"  "<<steric_Torque[i]<<" "<<mag_Torque[i]+steric_Torque[i]<<"  "<<ext_Torque[i]<<std::endl; 
//   }
//   
//   print_mag_foreces_and_torques(); 
}
