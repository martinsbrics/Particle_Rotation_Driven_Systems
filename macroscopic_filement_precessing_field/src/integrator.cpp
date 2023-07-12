/*
 * Copyright 2019 Martins Brics <email>
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


Integrator::Integrator(size_t Ns1, double Cm1, double H11, double om1):Ns(Ns1), Np(Ns1+1),Nv(3*(Np)), alfa(Np*Np*Np), mu(Np) ,  Cm(Cm1),H1{H11}, om{om1},J_Jt_u(Ns-1),J_Jt_l(Ns-1),J_Jt_d(Ns), FA(Nv),A_r(Nv),J_r(Ns), temp_r(Ns), Jt(Nv,Ns), J(Ns, Nv) 
{
  E= Eigen::MatrixXd::Identity(Nv, Nv);
  E_small= Eigen::MatrixXd::Identity(Ns, Ns);
  ones=Eigen::VectorXd::Ones(Nv);
  A= Eigen::MatrixXd::Zero(Nv, Nv);

  for(size_t i=0; i<3; i++)
  {
    A(i,i)=-alfa;
    A(Nv-1-i,Nv-1-i)=-alfa;
    A(i,i+3)+=2*alfa;
    A(Nv-1-i,Nv-4-i)=2*alfa;
    A(i,i+6)=-alfa;
    A(Nv-1-i,Nv-7-i)=-alfa;
  }

  for(size_t i=3; i<6; i++)
  {
    A(i,i-3)=2*alfa;
    A(Nv-1-i,Nv-1-i+3)=2*alfa;
    A(i,i)=-5*alfa;
    A(Nv-1-i,Nv-1-i)=-5*alfa;
    A(i,i+3)=4*alfa;
    A(Nv-1-i,Nv-4-i)=4*alfa;
    A(i,i+6)=-alfa;
    A(Nv-1-i,Nv-7-i)=-alfa;
  }

  for(int i=2; i<Ns-1; i++)
  {
    A(3*i,3*(i+2))=-alfa;
    A(3*i+1,3*(i+2)+1)=-alfa;
    A(3*i+2,3*(i+2)+2)=-alfa;
    A(3*i,3*(i-2))=-alfa;
    A(3*i+1,3*(i-2)+1)=-alfa;
    A(3*i+2,3*(i-2)+2)=-alfa;
    A(3*i,3*(i+1))=4*alfa;
    A(3*i+1,3*(i+1)+1)=4*alfa;
    A(3*i+2,3*(i+1)+2)=4*alfa;
    A(3*i,3*(i-1))=4*alfa;
    A(3*i+1,3*(i-1)+1)=4*alfa;
    A(3*i+2,3*(i-1)+2)=4*alfa;
    A(3*i,3*(i))=-6*alfa;
    A(3*i+1,3*(i)+1)=-6*alfa;
    A(3*i+2,3*(i)+2)=-6*alfa;
  }
  
  A_sparse=A.sparseView();
  A_sparse.makeCompressed();
  
  for(size_t i=0; i<Ns; i++)
  {
    J.insert(i, 3*i)=-1;
    J.insert(i, 3*i+1)=-2;
    J.insert(i, 3*i+2)=-3;
    J.insert(i, 3*(i+1))=1;
    J.insert(i, 3*(i+1)+1)=2;
    J.insert(i, 3*(i+1)+2)=3;
  }
  J.makeCompressed();
  Jt=J.transpose();  
  
};

void Integrator::fill_J( const Eigen::VectorXd & x)
{
  double *Jtt=Jt.valuePtr();
  size_t i=0;
  for(i=0; i<Ns; i++)
  {
    Jtt[6*i]=-2*(x[3*(i+1)]-x[3*i]);
    Jtt[6*i+1]=-2*(x[3*(i+1)+1]-x[3*i+1]);
    Jtt[6*i+2]=-2*(x[3*(i+1)+2]-x[3*i+2]);
    Jtt[6*i+3]=2*(x[3*(i+1)]-x[3*i]);
    Jtt[6*i+4]=2*(x[3*(i+1)+1]-x[3*i+1]);
    Jtt[6*i+5]=2*(x[3*(i+1)+2]-x[3*i+2]);
  }
  J=Jt.transpose();
  for(i=0; i<Ns-1; i++)
  {
    J_Jt_d[i]=8*((x[3*(i+1)]-x[3*i])*(x[3*(i+1)]-x[3*i])+(x[3*(i+1)+1]-x[3*i+1])*(x[3*(i+1)+1]-x[3*i+1])+(x[3*(i+1)+2]-x[3*i+2])*(x[3*(i+1)+2]-x[3*i+2]));
    J_Jt_u[i]=-4*((x[3*(i+1)]-x[3*i])*(x[3*(i+2)]-x[3*(i+1)])+(x[3*(i+1)+1]-x[3*i+1])*(x[3*(i+2)+1]-x[3*(i+1)+1])+(x[3*(i+1)+2]-x[3*i+2])*(x[3*(i+2)+2]-x[3*(i+1)+2]));
    J_Jt_l[i]=J_Jt_u[i];
  }
  i=Ns-1;
  J_Jt_d[i]=8*((x[3*(i+1)]-x[3*i])*(x[3*(i+1)]-x[3*i])+(x[3*(i+1)+1]-x[3*i+1])*(x[3*(i+1)+1]-x[3*i+1])+(x[3*(i+1)+2]-x[3*i+2])*(x[3*(i+1)+2]-x[3*i+2])); 
}


void Integrator::fill_Fa(const Eigen::VectorXd & x, double t)
{
  double res;double t1=t+M_PI/2/om;
  
  for(size_t i=1; i<Ns; i++)
  {
    res=(x(3*(i+1))+x(3*(i-1))-2*x(3*i))+H1*cos(om*t1)*(x(3*(i+1)+1)+x(3*(i-1)+1)-2*x(3*i+1))+H1*sin(om*t1)*(x(3*(i+1)+2)+x(3*(i-1)+2)-2*x(3*i+2));
    FA(3*i)=-Cm*mu*res;
    FA(3*i+1)=-Cm*mu*res*H1*cos(om*t1);
    FA(3*i+2)=-Cm*mu*res*H1*sin(om*t1);
  }
  res=(8*x(3)-7*x(0)-x(6))+H1*cos(om*t1)*(8*x(4)-7*x(1)-x(7))+H1*sin(om*t1)*(8*x(5)-7*x(2)-x(8));
  FA(0)=-Cm*res/6*mu;
  FA(1)=-Cm*res*H1*cos(om*t1)/6*mu;
  FA(2)=-Cm*res*H1*sin(om*t1)/6*mu;
  res=7*x(3*Ns)+x(3*(Ns-2))-8*x(3*(Ns-1))+H1*cos(om*t1)*(7*x(3*Ns+1)+x(3*(Ns-2)+1)-8*x(3*(Ns-1)+1))+H1*sin(om*t1)*(7*x(3*Ns+2)+x(3*(Ns-2)+2)-8*x(3*(Ns-1)+2));
  FA(3*(Ns))=Cm*res/6*mu;
  FA(3*(Ns)+1)=Cm*res/6*mu*H1*cos(om*t1);
  FA(3*(Ns)+2)=Cm*res/6*mu*H1*sin(om*t1);
}

long Integrator::dgtsv(long Nx1, long NRHS, double *DL, double *D, double *DU, double *imag_t_Bx11, long LDB)
{
    long info;
    dgtsv_(&Nx1, &NRHS, DL, D, DU, imag_t_Bx11, &LDB, &info);
    return info;
}


 void Integrator::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t)
{
  int info;
  renorm(x,x_temp);
  fill_Fa(x_temp, t);
  fill_J(x_temp);
  A_r=A_sparse*x_temp+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  dxdt=mu*(A_r-Jt*temp_r);
//  cout.precision(14);
//   fill_J(x);
//   Eigen::MatrixXd AE= Eigen::MatrixXd::Identity(Ns, Ns);
//   Eigen::VectorXd b(Eigen::Map<Eigen::VectorXd>(AE.data(), AE.cols()*AE.rows()));
//   info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &b[0], Ns);
//   cout<<dxdt<<endl<<endl;;
//   Eigen::MatrixXd B;
//   B=Eigen::Map<Eigen::MatrixXd> (b.data(), Ns,Ns);
//   //Eigen::Map<Eigen::MatrixXd> B(b.data(), Ns,Ns);
//   cout<<mu*(A_r-Jt*B*J*A_r)<<endl;
//   abort();
  
  
  
//    cout.precision(14);
//   cout<<"t="<<(t-M_PI/(2*om))<<endl;
//   cout<<dxdt<<endl<<endl;;
//   abort();
//   cout <<A_r<<endl<<endl;
//   cout <<Jt*temp_r<<endl<<endl;
//    cout <<FA<<endl<<endl;
//   for(int i=0; i<Np; i++)
//   {
//     dxdt[i]=0;
//     dxdt[i+Np]=0;
//     dxdt[i+2*Np]=0;
//     for(int j=0; j<i; j++)
//     {
//       lcosj=cos(t-x[j+2*Np]);
//       lsinj=sin(t-x[j+2*Np]);
//       b=sqrt((x[i]-x[j])*(x[i]-x[j])+(x[i+Np]-x[j+Np])*(x[i+Np]-x[j+Np]));
//       b3=b*b*b;
//       b5=b3*b*b;
//       b13=b5*b5*b3;
//       a=lcosj*(x[i]-x[j])+lsinj*(x[i+Np]-x[j+Np]);
//       c=lcosj*(x[i+Np]-x[j+Np])-lsinj*(x[i]-x[j]);
//       dxdt[i]+=-(x[i]-x[j])/(b3)+3.*a*a*(x[i]-x[j])/b5-(x[i]-x[j])/(b13);
//       dxdt[i+Np]+=-(x[i+Np]-x[j+Np])/(b3)+3.*a*a*(x[i+Np]-x[j+Np])/b5-(x[i+Np]-x[j+Np])/(b13);
//       dxdt[i+2*Np]+=-3.*a*c/b5;
//     }
//     for(int j=i+1; j<Np; j++)
//     {
//       lcosj=cos(t-x[j+2*Np]);
//       lsinj=sin(t-x[j+2*Np]);
//       b=sqrt((x[i]-x[j])*(x[i]-x[j])+(x[i+Np]-x[j+Np])*(x[i+Np]-x[j+Np]));
//       b3=b*b*b;
//       b5=b3*b*b;
//       b13=b5*b5*b3;
//       a=lcosj*(x[i]-x[j])+lsinj*(x[i+Np]-x[j+Np]);
//       c=lcosj*(x[i+Np]-x[j+Np])-lsinj*(x[i]-x[j]);
//       dxdt[i]+=-(x[i]-x[j])/(b3)+3.*a*a*(x[i]-x[j])/b5-(x[i]-x[j])/(b13);
//       dxdt[i+Np]+=-(x[i+Np]-x[j+Np])/(b3)+3.*a*a*(x[i+Np]-x[j+Np])/b5-(x[i+Np]-x[j+Np])/(b13);
//       dxdt[i+2*Np]+=-3.*a*c/b5;
//     }
//    // cout<<dxdt[0]<<" "<<dxdt[1]<<" "<<dxdt[2]<<" "<<dxdt[3]<<" "<<dxdt[4]<<" "<<dxdt[5]<<" "<<endl;
//     dxdt[i]*=lambda;
//     dxdt[i+Np]*=lambda;
//     dxdt[i+2*Np]*=lambda;
//     dxdt[i]+=cos(t-x[i+2*Np]);
//     dxdt[i+Np]+=sin(t-x[i+2*Np]);
//     dxdt[i+2*Np]+=1-gamma*sin(x[i+2*Np]);
//   }
}
  
double Integrator::dist(const  Eigen::VectorXd &x , int i)
{
  return sqrt((x(3*(i+1))-x(3*i))*(x(3*(i+1))-x(3*i))+(x(3*(i+1)+1)-x(3*i+1))*(x(3*(i+1)+1)-x(3*i+1))+(x(3*(i+1)+2)-x(3*i+2))*(x(3*(i+1)+2)-x(3*i+2)));
}

void Integrator::renorm(const  Eigen::VectorXd &x1, Eigen::VectorXd &x2)
{
  double att, renorm, alpha1, alpha2;
  x2=x1;
  for (size_t i=0; i<Ns; i++)
  {
    att=dist(x1,i);
    renorm=1.-1.0/(Ns)/att;
    alpha1=(Np-(i+1))*renorm/Np;
    alpha2=(i+1)*renorm/Np;
    for(size_t j=0; j<=i; j++)
    {
      x2(3*j)+=alpha1*(x1(3*(i+1))-x1(3*i));
      x2(3*j+1)+=alpha1*(x1(3*(i+1)+1)-x1(3*i+1));
      x2(3*j+2)+=alpha1*(x1(3*(i+1)+2)-x1(3*i+2));
    }
    for(size_t j=i+1; j<Np; j++)
    {
      x2(3*j)-=alpha2*(x1(3*(i+1))-x1(3*i));
      x2(3*j+1)-=alpha2*(x1(3*(i+1)+1)-x1(3*i+1));
      x2(3*j+2)-=alpha2*(x1(3*(i+1)+2)-x1(3*i+2));
    }
  }
}
