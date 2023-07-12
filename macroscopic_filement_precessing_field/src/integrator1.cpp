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

#include "integrator1.h"


Integrator1::Integrator1(size_t Ns1, double Cm1, double H11, double om1, double torq1):Ns(Ns1), Np(Ns1+1),Nv(3*(Np)), alfa(Np*Np*Np), mu(Np) ,  Cm(Cm1),H1{H11}, om{om1}, Torq(torq1),J_Jt_u(Ns-1),J_Jt_l(Ns-1),J_Jt_d(Ns),J_DD_Jt_u(Ns-1),J_DD_Jt_l(Ns-1),J_DD_Jt_d(Ns), rl(Nv), rll(Nv), FT(Nv), FA(Nv), Ftwist(Nv) /*, FR(Nv)*/,A_r(Nv),J_r(Ns), temp_r(Ns),DD_sparse(Nv,Nv),D1_sparse(Nv,Nv),D2_sparse(Nv,Nv), Jt(Nv,Ns), J(Ns, Nv), FA_mat(Nv, Nv) 
{
  cout<<"!!!! initilizing Cm="<<Cm<<" Troq="<<Torq<<endl;
  //abort();
  Cr=Torq;
  E= Eigen::MatrixXd::Identity(Nv, Nv);
  E_small= Eigen::MatrixXd::Identity(Ns, Ns);
  ones=Eigen::VectorXd::Ones(Nv);
  A= Eigen::MatrixXd::Zero(Nv, Nv);
  double kbend_hetero;

for(size_t i=0; i<3; i++)
  {
    kbend_hetero=1.0;
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
    kbend_hetero=1.0;
    A(i,i-3)=2*alfa*kbend_hetero;
    A(i,i)=-5*alfa*kbend_hetero;
    A(i,i+3)=4*alfa*kbend_hetero;
    A(i,i+6)=-alfa*kbend_hetero;
    kbend_hetero=1.0;
    A(Nv-1-i,Nv-1-i+3)=2*alfa*kbend_hetero;
    A(Nv-1-i,Nv-1-i)=-5*alfa*kbend_hetero;
    A(Nv-1-i,Nv-4-i)=4*alfa*kbend_hetero;
    A(Nv-1-i,Nv-7-i)=-alfa*kbend_hetero;
  }

  for(int i=2; i<Ns-1; i++)
  {
    kbend_hetero=1.0;
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
  
//   for(size_t i=0; i<3; i++)
//   {
//     kbend_hetero=(0.6 + 0.4*(1.0*0/Ns));
//     A(i,i)=-alfa*kbend_hetero;
//     A(i,i+3)+=2*alfa*kbend_hetero;
//     A(i,i+6)=-alfa*kbend_hetero;
//     kbend_hetero=1.0;
//     A(Nv-1-i,Nv-1-i)=-alfa*kbend_hetero;
//     A(Nv-1-i,Nv-4-i)=2*alfa*kbend_hetero;
//     A(Nv-1-i,Nv-7-i)=-alfa*kbend_hetero;
//   }
// 
//   for(size_t i=3; i<6; i++)
//   {
//     kbend_hetero=(0.6 + 0.4*(1.0*1/Ns));
//     A(i,i-3)=2*alfa*kbend_hetero;
//     A(i,i)=-5*alfa*kbend_hetero;
//     A(i,i+3)=4*alfa*kbend_hetero;
//     A(i,i+6)=-alfa*kbend_hetero;
//     kbend_hetero=(0.6 + 0.4*(1.0*(Ns-1)/Ns));
//     A(Nv-1-i,Nv-1-i+3)=2*alfa*kbend_hetero;
//     A(Nv-1-i,Nv-1-i)=-5*alfa*kbend_hetero;
//     A(Nv-1-i,Nv-4-i)=4*alfa*kbend_hetero;
//     A(Nv-1-i,Nv-7-i)=-alfa*kbend_hetero;
//   }
// 
//   for(int i=2; i<Ns-1; i++)
//   {
//     kbend_hetero=(0.6 + 0.4*(1.0*i/Ns));
//     A(3*i,3*(i+2))=-alfa*kbend_hetero;
//     A(3*i+1,3*(i+2)+1)=-alfa*kbend_hetero;
//     A(3*i+2,3*(i+2)+2)=-alfa*kbend_hetero;
//     A(3*i,3*(i-2))=-alfa*kbend_hetero;
//     A(3*i+1,3*(i-2)+1)=-alfa*kbend_hetero;
//     A(3*i+2,3*(i-2)+2)=-alfa*kbend_hetero;
//     A(3*i,3*(i+1))=4*alfa*kbend_hetero;
//     A(3*i+1,3*(i+1)+1)=4*alfa*kbend_hetero;
//     A(3*i+2,3*(i+1)+2)=4*alfa*kbend_hetero;
//     A(3*i,3*(i-1))=4*alfa*kbend_hetero;
//     A(3*i+1,3*(i-1)+1)=4*alfa*kbend_hetero;
//     A(3*i+2,3*(i-1)+2)=4*alfa*kbend_hetero;
//     A(3*i,3*(i))=-6*alfa*kbend_hetero;
//     A(3*i+1,3*(i)+1)=-6*alfa*kbend_hetero;
//     A(3*i+2,3*(i)+2)=-6*alfa*kbend_hetero;
//   }
  
  for(int i=1; i<Ns; i++)
  {
    A(3*i+1,3*(i+1)+2)+=Torq*mu;
    A(3*i+1,3*(i-1)+2)+=Torq*mu;
    A(3*i+1,3*(i)+2)-=2*Torq*mu;
    A(3*i+2,3*(i+1)+1)-=Torq*mu;
    A(3*i+2,3*(i-1)+1)-=Torq*mu;
    A(3*i+2,3*(i)+1)+=2*Torq*mu;
  }

  A(1,5)+=8*Torq/6*mu;
  A(1,8)-=Torq*mu/6;
  A(1,2)-=Torq*mu*7/6;
  A(2,4)-=Torq*8*mu/6;
  A(2,7)+=Torq*mu/6;
  A(2,1)+=7*Torq*mu/6;
  A(Nv-2,Nv-1)-=Torq*7*mu/6;
  A(Nv-2,Nv-7)-=Torq*mu/6;
  A(Nv-2,Nv-4)+=Torq*8*mu/6;
  A(Nv-1,Nv-2)+=Torq*7*mu/6;
  A(Nv-1,Nv-8)+=Torq*mu/6;
  A(Nv-1,Nv-5)-=Torq*8*mu/6;
  
//   Eigen::EigenSolver<Eigen::MatrixXd> es(A);
//   cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
//   abort();

  
  A_sparse=A.sparseView();
  A_sparse_row=A.sparseView();
  A_sparse.makeCompressed();
  A_sparse_row.makeCompressed();
  
  DD_sparse.insert(0,0)=2.0;
  DD_sparse.insert(0,1)=-1.0;
  DD_sparse.insert(0,2)=-1.0;
  DD_sparse.insert(1,0)=-1.0;
  DD_sparse.insert(1,1)=2.0;
  DD_sparse.insert(1,2)=-1.0;
  DD_sparse.insert(2,0)=-1.0;
  DD_sparse.insert(2,1)=-1.0;
  DD_sparse.insert(2,2)=2.0;
  for(size_t i=1; i<Ns; i++)
  {
    DD_sparse.insert(3*i+0,3*i+0)=2.0;
    DD_sparse.insert(3*i+0,3*i+1)=-1.0;
    DD_sparse.insert(3*i+0,3*i+2)=-1.0;
    DD_sparse.insert(3*i+1,3*i+0)=-1.0;
    DD_sparse.insert(3*i+1,3*i+1)=2.0;
    DD_sparse.insert(3*i+1,3*i+2)=-1.0;
    DD_sparse.insert(3*i+2,3*i+0)=-1.0;
    DD_sparse.insert(3*i+2,3*i+1)=-1.0;
    DD_sparse.insert(3*i+2,3*i+2)=2.0;
  }
  DD_sparse.insert(3*Ns+0,3*Ns+0)=2.0;
  DD_sparse.insert(3*Ns+0,3*Ns+1)=-1.0;
  DD_sparse.insert(3*Ns+0,3*Ns+2)=-1.0;
  DD_sparse.insert(3*Ns+1,3*Ns+0)=-1.0;
  DD_sparse.insert(3*Ns+1,3*Ns+1)=2.0;
  DD_sparse.insert(3*Ns+1,3*Ns+2)=-1.0;
  DD_sparse.insert(3*Ns+2,3*Ns+0)=-1.0;
  DD_sparse.insert(3*Ns+2,3*Ns+1)=-1.0;
  DD_sparse.insert(3*Ns+2,3*Ns+2)=2.0;
  DD_sparse.makeCompressed();
  
  D1_sparse.insert(0,0)=-1.5;
  D1_sparse.insert(0,3)=2.0;
  D1_sparse.insert(0,6)=-0.5;
  D1_sparse.insert(1,1)=-1.5;
  D1_sparse.insert(1,4)=2.0;
  D1_sparse.insert(1,7)=-0.5;
  D1_sparse.insert(2,2)=-1.5;
  D1_sparse.insert(2,5)=2.0;
  D1_sparse.insert(2,8)=-0.5;
  for(size_t i=1; i<Ns; i++){
  
    D1_sparse.insert(3*i,3*(i-1))=-1.0*0.5;
    D1_sparse.insert(3*i,3*(i+1))=1.0*0.5;
    D1_sparse.insert(3*i+1,3*(i-1)+1)=-1.0*0.5;
    D1_sparse.insert(3*i+1,3*(i+1)+1)=1.0*0.5;
    D1_sparse.insert(3*i+2,3*(i-1)+2)=-1.0*0.5;
    D1_sparse.insert(3*i+2,3*(i+1)+2)=1.0*0.5;
  }
  D1_sparse.insert(3*Ns,3*(Ns-2))=1.0*0.5;
  D1_sparse.insert(3*Ns,3*(Ns-1))=-4.0*0.5;
  D1_sparse.insert(3*Ns,3*(Ns))=1.5;
  D1_sparse.insert(3*Ns+1,3*(Ns-2)+1)=1.0*0.5;
  D1_sparse.insert(3*Ns+1,3*(Ns-1)+1)=-4.0*0.5;
  D1_sparse.insert(3*Ns+1,3*(Ns)+1)=1.5;
  D1_sparse.insert(3*Ns+2,3*(Ns-2)+2)=1.0*0.5;
  D1_sparse.insert(3*Ns+2,3*(Ns-1)+2)=-4.0*0.5;
  D1_sparse.insert(3*Ns+2,3*(Ns)+2)=3.0*0.5;
  D1_sparse.makeCompressed();
  
  for(size_t i=1; i<Ns; i++){
    D2_sparse.insert(3*i,3*(i-1))=1.0;
    D2_sparse.insert(3*i,3*(i))=-2.0;
    D2_sparse.insert(3*i,3*(i+1))=1.0;
    D2_sparse.insert(3*i+1,3*(i-1)+1)=1.0;
    D2_sparse.insert(3*i+1,3*(i)+1)=-2.0;
    D2_sparse.insert(3*i+1,3*(i+1)+1)=1.0;
    D2_sparse.insert(3*i+2,3*(i-1)+2)=1.0;
    D2_sparse.insert(3*i+2,3*(i)+2)=-2.0;
    D2_sparse.insert(3*i+2,3*(i+1)+2)=1.0;
  }
  D2_sparse.makeCompressed();

  
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
  
  for(size_t i=1; i<Ns; i++)
  {
    FA_mat.insert(3*i,3*(i-1))=1;
    FA_mat.insert(3*i,3*(i-1)+1)=1;
    FA_mat.insert(3*i,3*(i-1)+2)=1;
    FA_mat.insert(3*i,3*(i))=1;
    FA_mat.insert(3*i,3*(i)+1)=1;
    FA_mat.insert(3*i,3*(i)+2)=1;
    FA_mat.insert(3*i,3*(i+1))=1;
    FA_mat.insert(3*i,3*(i+1)+1)=1;
    FA_mat.insert(3*i,3*(i+1)+2)=1;
    
    FA_mat.insert(3*i+1,3*(i-1))=1;
    FA_mat.insert(3*i+1,3*(i-1)+1)=1;
    FA_mat.insert(3*i+1,3*(i-1)+2)=1;
    FA_mat.insert(3*i+1,3*(i))=1;
    FA_mat.insert(3*i+1,3*(i)+1)=1;
    FA_mat.insert(3*i+1,3*(i)+2)=1;
    FA_mat.insert(3*i+1,3*(i+1))=1;
    FA_mat.insert(3*i+1,3*(i+1)+1)=1;
    FA_mat.insert(3*i+1,3*(i+1)+2)=1;
    
    FA_mat.insert(3*i+2,3*(i-1))=1;
    FA_mat.insert(3*i+2,3*(i-1)+1)=1;
    FA_mat.insert(3*i+2,3*(i-1)+2)=1;
    FA_mat.insert(3*i+2,3*(i))=1;
    FA_mat.insert(3*i+2,3*(i)+1)=1;
    FA_mat.insert(3*i+2,3*(i)+2)=1;
    FA_mat.insert(3*i+2,3*(i+1))=1;
    FA_mat.insert(3*i+2,3*(i+1)+1)=1;
    FA_mat.insert(3*i+2,3*(i+1)+2)=1;
  }
  FA_mat.insert(0,0)=1;
  FA_mat.insert(0,1)=1;
  FA_mat.insert(0,2)=1;
  FA_mat.insert(0,3)=1;
  FA_mat.insert(0,4)=1;
  FA_mat.insert(0,5)=1;
  FA_mat.insert(0,6)=1;
  FA_mat.insert(0,7)=1;
  FA_mat.insert(0,8)=1;
  
  FA_mat.insert(1,0)=1;
  FA_mat.insert(1,1)=1;
  FA_mat.insert(1,2)=1;
  FA_mat.insert(1,3)=1;
  FA_mat.insert(1,4)=1;
  FA_mat.insert(1,5)=1;
  FA_mat.insert(1,6)=1;
  FA_mat.insert(1,7)=1;
  FA_mat.insert(1,8)=1;
  
  FA_mat.insert(2,0)=1;
  FA_mat.insert(2,1)=1;
  FA_mat.insert(2,2)=1;
  FA_mat.insert(2,3)=1;
  FA_mat.insert(2,4)=1;
  FA_mat.insert(2,5)=1;
  FA_mat.insert(2,6)=1;
  FA_mat.insert(2,7)=1;
  FA_mat.insert(2,8)=1;
  
  FA_mat.insert(3*Ns,3*(Ns-2))=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+1)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+2)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+3)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+4)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+5)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+6)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+7)=1;
  FA_mat.insert(3*Ns,3*(Ns-2)+8)=1;
  
  FA_mat.insert(3*Ns+1,3*(Ns-2))=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+1)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+2)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+3)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+4)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+5)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+6)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+7)=1;
  FA_mat.insert(3*Ns+1,3*(Ns-2)+8)=1;
  
  FA_mat.insert(3*Ns+2,3*(Ns-2))=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+1)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+2)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+3)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+4)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+5)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+6)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+7)=1;
  FA_mat.insert(3*Ns+2,3*(Ns-2)+8)=1;
  
  FA_mat.makeCompressed();

};

void Integrator1::fill_J( const Eigen::VectorXd & x)
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


void Integrator1::fill_J_DD( const Eigen::VectorXd & x)
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
  double *DD=DD_sparse.valuePtr();
  DD[0]=1-lam*(4.0*x[3*(1)]-x[3*(2)]-3.0*x[3*(0)])*(4.0*x[3*(1)]-x[3*(2)]-3.0*x[3*(0)])*0.25*mu*mu;
  DD[1]=-lam*(4.0*x[3*(1)]-x[3*(2)]-3.0*x[3*(0)])*(4.0*x[3*(1)+1]-x[3*(2)+1]-3.0*x[3*(0)+1])*0.25*mu*mu;
  DD[2]=-lam*(4.0*x[3*(1)]-x[3*(2)]-3.0*x[3*(0)])*(4.0*x[3*(1)+2]-x[3*(2)+2]-3.0*x[3*(0)+2])*0.25*mu*mu;
  DD[3]=-lam*(4.0*x[3*(1)]-x[3*(2)]-3.0*x[3*(0)])*(4.0*x[3*(1)+1]-x[3*(2)+1]-3.0*x[3*(0)+1])*0.25*mu*mu;
  DD[4]=1-lam*(4.0*x[3*(1)+1]-x[3*(2)+1]-3.0*x[3*(0)+1])*(4.0*x[3*(1)+1]-x[3*(2)+1]-3.0*x[3*(0)+1])*0.25*mu*mu;
  DD[5]=-lam*(4.0*x[3*(1)+1]-x[3*(2)+1]-3.0*x[3*(0)+1])*(4.0*x[3*(1)+2]-x[3*(2)+2]-3.0*x[3*(0)+2])*0.25*mu*mu;
  DD[6]=-lam*(4.0*x[3*(1)]-x[3*(2)]-3.0*x[3*(0)])*(4.0*x[3*(1)+2]-x[3*(2)+2]-3.0*x[3*(0)+2])*0.25*mu*mu;
  DD[7]=-lam*(4.0*x[3*(1)+1]-x[3*(2)+1]-3.0*x[3*(0)+1])*(4.0*x[3*(1)+2]-x[3*(2)+2]-3.0*x[3*(0)+2])*0.25*mu*mu;
  DD[8]=1-lam*(4.0*x[3*(1)+2]-x[3*(2)+2]-3.0*x[3*(0)+2])*(4.0*x[3*(1)+2]-x[3*(2)+2]-3.0*x[3*(0)+2])*0.25*mu*mu;
  for(size_t i=1; i<Ns; i++)
  {
    DD[9*i+0]=1-lam*(x[3*(i+1)]-x[3*(i-1)])*(x[3*(i+1)]-x[3*(i-1)])*0.25*mu*mu;
    DD[9*i+1]=-lam*(x[3*(i+1)]-x[3*(i-1)])*(x[3*(i+1)+1]-x[3*(i-1)+1])*0.25*mu*mu;
    DD[9*i+2]=-lam*(x[3*(i+1)]-x[3*(i-1)])*(x[3*(i+1)+2]-x[3*(i-1)+2])*0.25*mu*mu;
    DD[9*i+3]=-lam*(x[3*(i+1)]-x[3*(i-1)])*(x[3*(i+1)+1]-x[3*(i-1)+1])*0.25*mu*mu;
    DD[9*i+4]=1-lam*(x[3*(i+1)+1]-x[3*(i-1)+1])*(x[3*(i+1)+1]-x[3*(i-1)+1])*0.25*mu*mu;
    DD[9*i+5]=-lam*(x[3*(i+1)+1]-x[3*(i-1)+1])*(x[3*(i+1)+2]-x[3*(i-1)+2])*0.25*mu*mu;
    DD[9*i+6]=-lam*(x[3*(i+1)]-x[3*(i-1)])*(x[3*(i+1)+2]-x[3*(i-1)+2])*0.25*mu*mu;
    DD[9*i+7]=-lam*(x[3*(i+1)+2]-x[3*(i-1)+2])*(x[3*(i+1)+1]-x[3*(i-1)+1])*0.25*mu*mu;
    DD[9*i+8]=1-lam*(x[3*(i+1)+2]-x[3*(i-1)+2])*(x[3*(i+1)+2]-x[3*(i-1)+2])*0.25*mu*mu;
  }
  DD[9*Ns+0]=1-lam*(x[3*(Ns-2)]-4.0*x[3*(Ns-1)]+3.0*x[3*(Ns)])*(x[3*(Ns-2)]-4.0*x[3*(Ns-1)]+3.0*x[3*(Ns)])*0.25*mu*mu;
  DD[9*Ns+1]=-lam*(x[3*(Ns-2)]-4.0*x[3*(Ns-1)]+3.0*x[3*(Ns)])*(x[3*(Ns-2)+1]-4.0*x[3*(Ns-1)+1]+3.0*x[3*(Ns)+1])*0.25*mu*mu;
  DD[9*Ns+2]=-lam*(x[3*(Ns-2)]-4.0*x[3*(Ns-1)]+3.0*x[3*(Ns)])*(x[3*(Ns-2)+2]-4.0*x[3*(Ns-1)+2]+3.0*x[3*(Ns)+2])*0.25*mu*mu;
  DD[9*Ns+3]=-lam*(x[3*(Ns-2)]-4.0*x[3*(Ns-1)]+3.0*x[3*(Ns)])*(x[3*(Ns-2)+1]-4.0*x[3*(Ns-1)+1]+3.0*x[3*(Ns)+1])*0.25*mu*mu;
  DD[9*Ns+4]=1-lam*(x[3*(Ns-2)+1]-4.0*x[3*(Ns-1)+1]+3.0*x[3*(Ns)+1])*(x[3*(Ns-2)+1]-4.0*x[3*(Ns-1)+1]+3.0*x[3*(Ns)+1])*0.25*mu*mu;
  DD[9*Ns+5]=-lam*(x[3*(Ns-2)+1]-4.0*x[3*(Ns-1)+1]+3.0*x[3*(Ns)+1])*(x[3*(Ns-2)+2]-4.0*x[3*(Ns-1)+2]+3.0*x[3*(Ns)+2])*0.25*mu*mu;
  DD[9*Ns+6]=-lam*(x[3*(Ns-2)]-4.0*x[3*(Ns-1)]+3.0*x[3*(Ns)])*(x[3*(Ns-2)+2]-4.0*x[3*(Ns-1)+2]+3.0*x[3*(Ns)+2])*0.25*mu*mu;
  DD[9*Ns+7]=-lam*(x[3*(Ns-2)+1]-4.0*x[3*(Ns-1)+1]+3.0*x[3*(Ns)+1])*(x[3*(Ns-2)+2]-4.0*x[3*(Ns-1)+2]+3.0*x[3*(Ns)+2])*0.25*mu*mu;
  DD[9*Ns+8]=1-lam*(x[3*(Ns-2)+2]-4.0*x[3*(Ns-1)+2]+3.0*x[3*(Ns)+2])*(x[3*(Ns-2)+2]-4.0*x[3*(Ns-1)+2]+3.0*x[3*(Ns)+2])*0.25*mu*mu;
  
  Eigen::SparseMatrix<double>J_DD_Jt=J*DD_sparse*Jt;
  double *mat=J_DD_Jt.valuePtr();
 
  J_DD_Jt_d[0]=mat[0];
  J_DD_Jt_l[0]=mat[1];
  for(i=1; i<Ns-1; i++)
  {
    J_DD_Jt_u[i-1]=mat[3*i-1];
    J_DD_Jt_d[i]=mat[3*i];
    J_DD_Jt_l[i]=mat[3*i+1];
  }
  i=Ns-1;
  J_DD_Jt_u[i-1]=mat[3*Ns-4];
  J_DD_Jt_d[i]=mat[3*Ns-3];
  
  
//   for(i=0; i<Ns-1; i++)
//   {
//     J_Jt_d[i]=8*((x[3*(i+1)]-x[3*i])*(x[3*(i+1)]-x[3*i])+(x[3*(i+1)+1]-x[3*i+1])*(x[3*(i+1)+1]-x[3*i+1])+(x[3*(i+1)+2]-x[3*i+2])*(x[3*(i+1)+2]-x[3*i+2]));
//     J_Jt_u[i]=-4*((x[3*(i+1)]-x[3*i])*(x[3*(i+2)]-x[3*(i+1)])+(x[3*(i+1)+1]-x[3*i+1])*(x[3*(i+2)+1]-x[3*(i+1)+1])+(x[3*(i+1)+2]-x[3*i+2])*(x[3*(i+2)+2]-x[3*(i+1)+2]));
//     J_Jt_l[i]=J_Jt_u[i];
//   }
//   i=Ns-1;
//   J_Jt_d[i]=8*((x[3*(i+1)]-x[3*i])*(x[3*(i+1)]-x[3*i])+(x[3*(i+1)+1]-x[3*i+1])*(x[3*(i+1)+1]-x[3*i+1])+(x[3*(i+1)+2]-x[3*i+2])*(x[3*(i+1)+2]-x[3*i+2])); 
//   cout<<"!!!!"<<endl;
//   cout<<J_Jt_d.transpose()<<" "<<J_Jt_u.transpose()<<" "<<J_Jt_l.transpose()<<endl;
//   cout<<J_DD_Jt_d.transpose()<<" "<<J_DD_Jt_u.transpose()<<" "<<J_DD_Jt_l.transpose()<<endl;
  
//   cout<<DD_sparse<<endl<<endl<<J<<endl<<x.transpose()<<endl;
//   cout <<J_DD_Jt<<endl<<J_DD_Jt.nonZeros()<<endl;
//   abort();
}

 void Integrator1::fill_Fa_mat(const Eigen::VectorXd & x, double t)
{
  double t1=(t>0?t+M_PI/2/om :M_PI/2/om );
  double *FA_sparse=FA_mat.valuePtr();
  
  FA_sparse[0]=-Cm/6*mu*(-7);
  FA_sparse[1]=-Cm/6*mu*(-7)*H1*cos(om*t1);
  FA_sparse[2]=-Cm/6*mu*(-7)*H1*sin(om*t1);
  FA_sparse[3]=-Cm/6*mu*(8);
  FA_sparse[4]=-Cm/6*mu*(8)*H1*cos(om*t1);
  FA_sparse[5]=-Cm/6*mu*(8)*H1*sin(om*t1);
  FA_sparse[6]=-Cm/6*mu*(-1);
  FA_sparse[7]=-Cm/6*mu*(-1)*H1*cos(om*t1);
  FA_sparse[8]=-Cm/6*mu*(-1)*H1*sin(om*t1);
  
  FA_sparse[9]=-Cm*H1*cos(om*t1)/6*mu*(-7);
  FA_sparse[10]=-Cm*H1*cos(om*t1)/6*mu*(-7)*H1*cos(om*t1);
  FA_sparse[11]=-Cm*H1*cos(om*t1)/6*mu*(-7)*H1*sin(om*t1);
  FA_sparse[12]=-Cm*H1*cos(om*t1)/6*mu*(8);
  FA_sparse[13]=-Cm*H1*cos(om*t1)/6*mu*(8)*H1*cos(om*t1);
  FA_sparse[14]=-Cm*H1*cos(om*t1)/6*mu*(8)*H1*sin(om*t1);
  FA_sparse[15]=-Cm*H1*cos(om*t1)/6*mu*(-1);
  FA_sparse[16]=-Cm*H1*cos(om*t1)/6*mu*(-1)*H1*cos(om*t1);
  FA_sparse[17]=-Cm*H1*cos(om*t1)/6*mu*(-1)*H1*sin(om*t1);
  
  FA_sparse[18]=-Cm*H1*sin(om*t1)/6*mu*(-7);
  FA_sparse[19]=-Cm*H1*sin(om*t1)/6*mu*(-7)*H1*cos(om*t1);
  FA_sparse[20]=-Cm*H1*sin(om*t1)/6*mu*(-7)*H1*sin(om*t1);
  FA_sparse[21]=-Cm*H1*sin(om*t1)/6*mu*(8);
  FA_sparse[22]=-Cm*H1*sin(om*t1)/6*mu*(8)*H1*cos(om*t1);
  FA_sparse[23]=-Cm*H1*sin(om*t1)/6*mu*(8)*H1*sin(om*t1);
  FA_sparse[24]=-Cm*H1*sin(om*t1)/6*mu*(-1);
  FA_sparse[25]=-Cm*H1*sin(om*t1)/6*mu*(-1)*H1*cos(om*t1);
  FA_sparse[26]=-Cm*H1*sin(om*t1)/6*mu*(-1)*H1*sin(om*t1);
  for(size_t i=1; i<Ns; i++)
  {
    FA_sparse[27*i+0]=-Cm*mu*(1);
    FA_sparse[27*i+1]=-Cm*mu*(1)*H1*cos(om*t1);
    FA_sparse[27*i+2]=-Cm*mu*(1)*H1*sin(om*t1);
    FA_sparse[27*i+3]=-Cm*mu*(-2);
    FA_sparse[27*i+4]=-Cm*mu*(-2)*H1*cos(om*t1);
    FA_sparse[27*i+5]=-Cm*mu*(-2)*H1*sin(om*t1);
    FA_sparse[27*i+6]=-Cm*mu*(1);
    FA_sparse[27*i+7]=-Cm*mu*(1)*H1*cos(om*t1);
    FA_sparse[27*i+8]=-Cm*mu*(1)*H1*sin(om*t1);
    
    FA_sparse[27*i+9]=-Cm*H1*cos(om*t1)*mu*(1);
    FA_sparse[27*i+10]=-Cm*H1*cos(om*t1)*mu*(1)*H1*cos(om*t1);
    FA_sparse[27*i+11]=-Cm*H1*cos(om*t1)*mu*(1)*H1*sin(om*t1);
    FA_sparse[27*i+12]=-Cm*H1*cos(om*t1)*mu*(-2);
    FA_sparse[27*i+13]=-Cm*H1*cos(om*t1)*mu*(-2)*H1*cos(om*t1);
    FA_sparse[27*i+14]=-Cm*H1*cos(om*t1)*mu*(-2)*H1*sin(om*t1);
    FA_sparse[27*i+15]=-Cm*H1*cos(om*t1)*mu*(1);
    FA_sparse[27*i+16]=-Cm*H1*cos(om*t1)*mu*(1)*H1*cos(om*t1);
    FA_sparse[27*i+17]=-Cm*H1*cos(om*t1)*mu*(1)*H1*sin(om*t1);
    
    FA_sparse[27*i+18]=-Cm*H1*sin(om*t1)*mu*(1);
    FA_sparse[27*i+19]=-Cm*H1*sin(om*t1)*mu*(1)*H1*cos(om*t1);
    FA_sparse[27*i+20]=-Cm*H1*sin(om*t1)*mu*(1)*H1*sin(om*t1);
    FA_sparse[27*i+21]=-Cm*H1*sin(om*t1)*mu*(-2);
    FA_sparse[27*i+22]=-Cm*H1*sin(om*t1)*mu*(-2)*H1*cos(om*t1);
    FA_sparse[27*i+23]=-Cm*H1*sin(om*t1)*mu*(-2)*H1*sin(om*t1);
    FA_sparse[27*i+24]=-Cm*H1*sin(om*t1)*mu*(1);
    FA_sparse[27*i+25]=-Cm*H1*sin(om*t1)*mu*(1)*H1*cos(om*t1);
    FA_sparse[27*i+26]=-Cm*H1*sin(om*t1)*mu*(1)*H1*sin(om*t1);
  }
  FA_sparse[27*Ns+0]=Cm/6*mu*(1);
  FA_sparse[27*Ns+1]=Cm/6*mu*(1)*H1*cos(om*t1);
  FA_sparse[27*Ns+2]=Cm/6*mu*(1)*H1*sin(om*t1);
  FA_sparse[27*Ns+3]=Cm/6*mu*(-8);
  FA_sparse[27*Ns+4]=Cm/6*mu*(-8)*H1*cos(om*t1);
  FA_sparse[27*Ns+5]=Cm/6*mu*(-8)*H1*sin(om*t1);
  FA_sparse[27*Ns+6]=Cm/6*mu*(7);
  FA_sparse[27*Ns+7]=Cm/6*mu*(7)*H1*cos(om*t1);
  FA_sparse[27*Ns+8]=Cm/6*mu*(7)*H1*sin(om*t1);
  
  FA_sparse[27*Ns+9]=Cm*H1*cos(om*t1)/6*mu*(1);
  FA_sparse[27*Ns+10]=Cm*H1*cos(om*t1)/6*mu*(1)*H1*cos(om*t1);
  FA_sparse[27*Ns+11]=Cm*H1*cos(om*t1)/6*mu*(1)*H1*sin(om*t1);
  FA_sparse[27*Ns+12]=Cm*H1*cos(om*t1)/6*mu*(-8);
  FA_sparse[27*Ns+13]=Cm*H1*cos(om*t1)/6*mu*(-8)*H1*cos(om*t1);
  FA_sparse[27*Ns+14]=Cm*H1*cos(om*t1)/6*mu*(-8)*H1*sin(om*t1);
  FA_sparse[27*Ns+15]=Cm*H1*cos(om*t1)/6*mu*(7);
  FA_sparse[27*Ns+16]=Cm*H1*cos(om*t1)/6*mu*(7)*H1*cos(om*t1);
  FA_sparse[27*Ns+17]=Cm*H1*cos(om*t1)/6*mu*(7)*H1*sin(om*t1);
  
  FA_sparse[27*Ns+18]=Cm*H1*sin(om*t1)/6*mu*(1);
  FA_sparse[27*Ns+19]=Cm*H1*sin(om*t1)/6*mu*(1)*H1*cos(om*t1);
  FA_sparse[27*Ns+20]=Cm*H1*sin(om*t1)/6*mu*(1)*H1*sin(om*t1);
  FA_sparse[27*Ns+21]=Cm*H1*sin(om*t1)/6*mu*(-8);
  FA_sparse[27*Ns+22]=Cm*H1*sin(om*t1)/6*mu*(-8)*H1*cos(om*t1);
  FA_sparse[27*Ns+23]=Cm*H1*sin(om*t1)/6*mu*(-8)*H1*sin(om*t1);
  FA_sparse[27*Ns+24]=Cm*H1*sin(om*t1)/6*mu*(7);
  FA_sparse[27*Ns+25]=Cm*H1*sin(om*t1)/6*mu*(7)*H1*cos(om*t1);
  FA_sparse[27*Ns+26]=Cm*H1*sin(om*t1)/6*mu*(7)*H1*sin(om*t1);
  
  
  
}

//  void Integrator1::fill_Fa(const Eigen::VectorXd & x, double t)
// {
//   double res;double t1=t+M_PI/2/om;
//   for(size_t i=1; i<Ns; i++)
//   {
//     res=(x(3*(i+1))+x(3*(i-1))-2*x(3*i))+H1*cos(om*t1)*(x(3*(i+1)+1)+x(3*(i-1)+1)-2*x(3*i+1))+H1*sin(om*t1)*(x(3*(i+1)+2)+x(3*(i-1)+2)-2*x(3*i+2));
//     FA(3*i)=-Cm*mu*res;
//     FA(3*i+1)=-Cm*mu*res*H1*cos(om*t1);
//     FA(3*i+2)=-Cm*mu*res*H1*sin(om*t1);
//   }
//   res=(8*x(3)-7*x(0)-x(6))+H1*cos(om*t1)*(8*x(4)-7*x(1)-x(7))+H1*sin(om*t1)*(8*x(5)-7*x(2)-x(8));
//   FA(0)=-Cm*res/6*mu;
//   FA(1)=-Cm*res*H1*cos(om*t1)/6*mu;
//   FA(2)=-Cm*res*H1*sin(om*t1)/6*mu;
//   res=7*x(3*Ns)+x(3*(Ns-2))-8*x(3*(Ns-1))+H1*cos(om*t1)*(7*x(3*Ns+1)+x(3*(Ns-2)+1)-8*x(3*(Ns-1)+1))+H1*sin(om*t1)*(7*x(3*Ns+2)+x(3*(Ns-2)+2)-8*x(3*(Ns-1)+2));
//   FA(3*(Ns))=Cm*res/6*mu;
//   FA(3*(Ns)+1)=Cm*res/6*mu*H1*cos(om*t1);
//   FA(3*(Ns)+2)=Cm*res/6*mu*H1*sin(om*t1);
// }



 void Integrator1::fill_Fa(const Eigen::VectorXd & x, double t)
{
  double res;double t1=t+M_PI/2/om;
  //hvec={cos(theta),sin(theta)*cos(om*t1),sin(theta)*sin(om*t1)};
  hvec={1.0,0.0,0.0};
  for(size_t i=1; i<Ns; i++)
  {
    res=(x(3*(i+1))+x(3*(i-1))-2*x(3*i))*hvec[0]+hvec[1]*(x(3*(i+1)+1)+x(3*(i-1)+1)-2*x(3*i+1))+hvec[2]*(x(3*(i+1)+2)+x(3*(i-1)+2)-2*x(3*i+2));
    FA(3*i)=-Cm*mu*res*hvec[0];
    FA(3*i+1)=-Cm*mu*res*hvec[1];
    FA(3*i+2)=-Cm*mu*res*hvec[2];
  }
  res=(8*x(3)-7*x(0)-x(6))*hvec[0]+hvec[1]*(8*x(4)-7*x(1)-x(7))+hvec[2]*(8*x(5)-7*x(2)-x(8));
  FA(0)=-Cm*res/6*mu*hvec[0];
  FA(1)=-Cm*res*hvec[1]/6*mu;
  FA(2)=-Cm*res*hvec[2]/6*mu;
  res=(7*x(3*Ns)+x(3*(Ns-2))-8*x(3*(Ns-1)))*hvec[0]+hvec[1]*(7*x(3*Ns+1)+x(3*(Ns-2)+1)-8*x(3*(Ns-1)+1))+hvec[2]*(7*x(3*Ns+2)+x(3*(Ns-2)+2)-8*x(3*(Ns-1)+2));
  FA(3*(Ns))=Cm*res/6*mu*hvec[0];
  FA(3*(Ns)+1)=Cm*res/6*mu*hvec[1];
  FA(3*(Ns)+2)=Cm*res/6*mu*hvec[2];
}

void Integrator1::make_vector_product(double *a, double *b, double *c){
  //cout<<"!!!here"<<endl; 
  for(size_t i=0; i<Np; i++){
   
     c[3*i]=a[3*i+1]*b[3*i+2]-a[3*i+2]*b[3*i+1];
     c[3*i+1]=a[3*i+2]*b[3*i]-a[3*i]*b[3*i+2];
     c[3*i+2]=a[3*i]*b[3*i+1]-a[3*i+1]*b[3*i];
  }  
}
  
  
void Integrator1::fill_Ftwist(const Eigen::VectorXd & x){
//    Eigen::VectorXd x1=x;
//    x1[0]=-0.0004995892312681977;
//    x1[1]=-0.49958923126819754;
//    x1[2]=0.0011172324723727313;
//    x1[3]=-0.0002497986081704593;
//    x1[4]=-0.2497986081704593;
//    x1[5]=-0.009109312256335413;
//    x1[6]=-3.241714864066196e-9;
//    x1[7]=-3.241714864066195e-6;
//    x1[8]=0.0010007075057537314;
//    x1[9]=0.0002497813526934264;
//    x1[10]=0.24978135269342644;
//    x1[11]=0.01137346046445289;
//    x1[12]=0.0004996097284600945;
//    x1[13]=0.49960972846009444;
//    x1[14]=0.002114937495252403;
  double c1=x[0]-x[3*Ns];
  //double c1=x1[0]-x1[3*Ns];
  double c2=-0.5*(x[0]+x[3*Ns]);
  //double c2=-0.5*(x1[0]+x1[3*Ns]);
  double h=1.0/Ns;
  double mu3=mu*mu*mu;
  double tmp;
   rl=D1_sparse*x;
   rll=D2_sparse*x;
  //rl=D1_sparse*x1;
  //rll=D2_sparse*x1;
  make_vector_product(&rll[0],&rl[0], &FT[0]);
  //cout<<mu3*FT.transpose()<<endl<<endl;
  for(size_t i=0; i<Np; i++){
    tmp=Cr*mu3*(x[3*i]+c1*(-0.5+h*i)+c2);
    //tmp=Cr*mu3*(x1[3*i]+c1*(-0.5+h*i)+c2);
    //cout<<"i="<<i <<" om="<<Cr*(x1[3*i]+c1*(-0.5+h*i)+c2)/**(x[3*i])*/<<" "<<c1<<" "<<c2 <<endl;
    FT[3*i]*=tmp;
    FT[3*i+1]*=tmp;
    FT[3*i+2]*=tmp;
  }
  Ftwist=D1_sparse*FT;
  Ftwist[0]=0.0;Ftwist[1]=0.0;Ftwist[2]=0.0;Ftwist[3*Ns]=0.0;Ftwist[3*Ns+1]=0.0;Ftwist[3*Ns+2]=0.0;
//   cout<<x.transpose()<<endl<<endl;
//   cout<<Ftwist.transpose()<<endl<<endl;
//   cout<<FT.transpose()<<endl<<endl;
//   cout<<rl.transpose()<<endl<<endl;
//   cout<<rll.transpose()<<endl<<endl;
//   cout<<D1_sparse<<endl;
//   abort();
}




//  void Integrator1::fill_Fr(const Eigen::VectorXd & x, double t)
// {
//   for(size_t i=1; i<Ns; i++)
//   {
//     FR(3*i)=0.0;
//     FR(3*i+1)=Cr*mu*(x(3*(i+1)+2)+x(3*(i-1)+2)-2*x(3*i+2));
//     FR(3*i+2)=-Cr*mu*(x(3*(i+1)+1)+x(3*(i-1)+2)-2*x(3*i+1));
//   }
//   FR(0)=0.0;
//   FR(1)=Cr*mu*0.5*(4*x(5)-3*x(2)-x(8));
//   FR(2)=-Cr*mu*0.5*(4*x(4)-3*x(1)-x(7));;
//   FR(3*(Ns))=0.0;
//   FR(3*(Ns)+1)=-Cr*mu*0.5*(3*x(3*Ns+2)+x(3*(Ns-2)+2)-4*x(3*(Ns-1)+2));
//   FR(3*(Ns)+2)=Cr*mu*0.5*(3*x(3*Ns+1)+x(3*(Ns-2)+1)-4*x(3*(Ns-1)+1));
// }



long Integrator1::dgtsv(long Nx1, long NRHS, double *DL, double *D, double *DU, double *imag_t_Bx11, long LDB)
{
    long info;
    dgtsv_(&Nx1, &NRHS, DL, D, DU, imag_t_Bx11, &LDB, &info);
    return info;
}

void Integrator1::do_Euler(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  int info;
  fill_Fa(x, t+dt*0.5);
  fill_J(x);
  A_r=A_sparse*x+FA;
  E_small_vec=(Eigen::Map<Eigen::VectorXd>(E_small.data(), E_small.cols()*E_small.rows()));
  info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &E_small_vec[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  Jt_J_inv=Eigen::Map<Eigen::MatrixXd> (E_small_vec.data(), Ns,Ns);
  rhs=/*(mu*dt)*(A_r-Jt*Jt_J_inv*J*A_r)*/(mu*dt)*(FA-Jt*Jt_J_inv*J*FA)+x;
  lhs=(-Jt*Jt_J_inv*J*A_sparse);
  lhs+=E*A_sparse;
  lhs*=(-mu*dt);
  lhs.diagonal().array() += 1;
  lhs_sparse=lhs.sparseView(1, 1e-13);
  //x_temp=lhs.fullPivLu().solve(rhs);
  //cout<<(mu*0.5*dt)*(A_r-Jt*Jt_J_inv*J*A_r)<<endl<<endl;abort();
  solver_sparse.analyzePattern(lhs_sparse);
  solver_sparse.factorize(lhs_sparse);
  new_x= solver_sparse.solve(rhs);
  //x_temp+=x;
//   cout<<"check norm!!!!!"<<endl;
//   for(size_t i=0; i<Ns; i++)
//   {
//     cout <<i<<" "<<dist(x_temp, i)<<endl;;
//   }

}

void Integrator1::fill_lhs_Euler_impl_expl(const double dt)
{
  //cout<<"fill_lhs_Euler_impl_expl: dt=" <<dt<<endl;
  lhs_sparse=(-mu*dt)*A_sparse;
  lhs_sparse.diagonal().array() += 1;
  solver_sparse.analyzePattern(lhs_sparse);
  solver_sparse.factorize(lhs_sparse);
}

void Integrator1::fill_lhs1_Euler_impl_expl(const double dt)
{
  lhs_sparse1=(-mu*dt*0.5)*A_sparse;
  lhs_sparse1.diagonal().array() += 1;
  solver_sparse1.analyzePattern(lhs_sparse1);
  solver_sparse1.factorize(lhs_sparse1);
}

void Integrator1::fill_lhs2_Euler_impl_expl(const double dt)
{
  double temp=1-M_SQRT1_2;
  //double temp=2./3.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::fill_lhs2_IMEX_BDF2(const double dt)
{
  double temp=2./3.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::do_Euler_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
//  cout <<"t="<<t<< " "<<dt<<endl;
  int info;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
//   E_small_vec=(Eigen::Map<Eigen::VectorXd>(E_small.data(), E_small.cols()*E_small.rows()));
//   info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &E_small_vec[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   }
//   Jt_J_inv=Eigen::Map<Eigen::MatrixXd> (E_small_vec.data(), Ns,Ns);
//   rhs=/*(mu*dt)*(A_r-Jt*Jt_J_inv*J*A_r)*/(mu*dt)*(FA-Jt*Jt_J_inv*J*A_r)+x;
  temp_r=J*A_r;  
//    cout<<temp_r<<endl<<endl;abort();
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  //cout<<temp_r<<endl<<endl;abort();
  rhs=(mu*dt)*(FA-Jt*temp_r)+x;  
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

void Integrator1::do_Euler_impl_expl_DD(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
//  cout <<"t="<<t<< " "<<dt<<endl;
  int info;
  fill_Fa(x, t);
  //fill_Fr(x, t);
  fill_J_DD(x);
  fill_Ftwist(x);
  A_r=DD_sparse*(A_sparse*x+FA+Ftwist)/*+FR*/;
//   E_small_vec=(Eigen::Map<Eigen::VectorXd>(E_small.data(), E_small.cols()*E_small.rows()));
//   info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &E_small_vec[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   }
//   Jt_J_inv=Eigen::Map<Eigen::MatrixXd> (E_small_vec.data(), Ns,Ns);
//   rhs=/*(mu*dt)*(A_r-Jt*Jt_J_inv*J*A_r)*/(mu*dt)*(FA-Jt*Jt_J_inv*J*A_r)+x;
  temp_r=J*A_r;  
//    cout<<temp_r<<endl<<endl;abort();
  info = dgtsv(Ns, 1, &J_DD_Jt_l[0], &J_DD_Jt_d[0], &J_DD_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  //cout<<temp_r<<endl<<endl;abort();
  rhs=(mu*dt)*(A_r-A_sparse*x/*+FR*/-DD_sparse*Jt*temp_r)+x;  
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


void Integrator1::do_rk4_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  int info;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
//   E_small_vec=(Eigen::Map<Eigen::VectorXd>(E_small.data(), E_small.cols()*E_small.rows()));
//   info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &E_small_vec[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   }
//   Jt_J_inv=Eigen::Map<Eigen::MatrixXd> (E_small_vec.data(), Ns,Ns);
//   rhs=/*(mu*dt)*(A_r-Jt*Jt_J_inv*J*A_r)*/(mu*dt)*(FA-Jt*Jt_J_inv*J*A_r)+x;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  rhs=(mu*dt)*(A_r-Jt*temp_r);  
  K1= solver_sparse.solve(rhs);
  x_temp=x+0.5*K1;
  fill_Fa(x_temp, t+0.5*dt);
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
  rhs=(mu*dt)*(A_r-Jt*temp_r);  
  K2= solver_sparse.solve(rhs);
  x_temp=x+0.5*K2;
  fill_Fa(x_temp, t+0.5*dt);
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
  rhs=(mu*dt)*(A_r-Jt*temp_r);  
  K3= solver_sparse.solve(rhs);
  x_temp=x+K3;
  fill_Fa(x_temp, t+dt);
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
  rhs=(mu*dt)*(A_r-Jt*temp_r);
  K4= solver_sparse.solve(rhs);
  new_x=x+1.0/6.0*(K1+2*K2+2*K3+K4);
  //x_temp+=x;
//   cout<<"check norm!!!!!"<<endl;
//   for(size_t i=0; i<Ns; i++)
//   {
//     cout <<i<<" "<<dist(x_temp, i)<<endl;;
//   }

}

// void Integrator1::do_Euler1_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
// {
//   int info;
//   fill_Fa(x, t+dt);
//   fill_J(x);
//   A_r=A_sparse*x+FA;
//   temp_r=J*A_r;  
//   info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   }
//   rhs=(mu*dt)*(FA-Jt*temp_r)+x;  
//   x_temp= solver_sparse.solve(rhs);
//   x_temp=0.5*(x_temp+x) ;
//   fill_Fa(x_temp, t+dt);
//   A_r=A_sparse*x_temp+FA;
//   fill_J(x_temp);
//   temp_r=J*A_r;  
//   info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   } 
//   rhs=(mu*dt)*(FA-Jt*temp_r)+x;  
//   x_temp= solver_sparse.solve(rhs);
//   x_temp=0.5*(x_temp+x) ;
//   fill_Fa(x_temp, t+dt);
//   A_r=A_sparse*x_temp+FA;
//   fill_J(x_temp);
//   temp_r=J*A_r;  
//   info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   } 
//   rhs=(mu*dt)*(FA-Jt*temp_r)+x;  
//   x_temp= solver_sparse.solve(rhs);
//   x_temp=0.5*(x_temp+x) ;
//   fill_Fa(x_temp, t+dt);
//   A_r=A_sparse*x_temp+FA;
//   fill_J(x_temp);
//   temp_r=J*A_r;  
//   info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   } 
//   rhs=(mu*dt)*(FA-Jt*temp_r)+x;
//   new_x= solver_sparse.solve(rhs);
// }

void Integrator1::do_Euler1_impl_expl(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  int info;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  rhs=(mu*dt)*(FA-Jt*temp_r)+x;  
  x_temp= solver_sparse.solve(rhs);
  fill_Fa(x_temp, t+dt);
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
  new_x=(mu*dt)*(A_r-Jt*temp_r)+x;
  
}

// void Integrator1::do_midpoint(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
// {
//   int info;
//   fill_Fa(x, t);
//   fill_J(x);
//   A_r=A_sparse*x+FA;
//   temp_r=J*A_r;  
//   info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   }
//   rhs=(mu*0.5*dt)*(FA-Jt*temp_r)+x;  
//   x_temp= solver_sparse1.solve(rhs);
//   fill_Fa(x_temp, t+0.5*dt);
//   fill_J(x_temp);
//   A_r=A_sparse*x_temp+FA;
//   temp_r=J*A_r;  
//   info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
//   if (info !=0)
//   {
//     cout << info<<endl;
//     fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
//     abort();
//   }
//   rhs=(mu*dt)*(FA-Jt*temp_r)+x;  
//   new_x= solver_sparse.solve(rhs);
//   //new_x=(mu*dt)*(A_r-Jt*temp_r)+x;  
//   
// }


void Integrator1::do_midpoint(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  double gamma=1-M_SQRT1_2;
  double delta=1-1.0/(2*gamma);
  int info;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K3=(FA-Jt*temp_r);
  rhs=(gamma*dt*mu)*K3+x;  
  x_temp= solver_sparse2.solve(rhs);
  K1=A_sparse*x_temp;
  fill_Fa(x_temp, t+gamma*dt);
  fill_J(x_temp);
  A_r=K1+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K4=(FA-Jt*temp_r);
  rhs=((1-gamma)*dt*mu)*K1+(delta*dt*mu)*K3+((1-delta)*mu*dt)*K4+x;
  x_temp= solver_sparse2.solve(rhs);
  K2=A_sparse*x_temp;
  new_x=dt*mu*((1-gamma)*K1+gamma*K2+(1-delta)*K3+delta*K4)+x;  
  
}


void Integrator1::do_IMEX_ARS_232(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  double gamma=1-M_SQRT1_2;
  double delta=-2*M_SQRT2/3;
  int info;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K3=(FA-Jt*temp_r);
  rhs=(gamma*dt*mu)*K3+x;  
  x_temp= solver_sparse1.solve(rhs);
  K1=A_sparse*x_temp;
  fill_Fa(x_temp, t+gamma*dt);
  fill_J(x_temp);
  A_r=K1+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K4=(FA-Jt*temp_r);
  rhs=((1-gamma)*dt*mu)*K1+(delta*dt*mu)*K3+((1-delta)*mu*dt)*K4+x;
  x_temp= solver_sparse1.solve(rhs);
  fill_Fa(x_temp, t+dt);
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
  K2=(A_r-Jt*temp_r);
  new_x=dt*mu*((1-gamma)*(K1+K4)+gamma*K2)+x;  
}


void Integrator1::do_CN(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  int info;
  fill_Fa(x, t+dt*0.5);
  fill_J(x);
  A_r=A_sparse*x+FA;
  E_small_vec=(Eigen::Map<Eigen::VectorXd>(E_small.data(), E_small.cols()*E_small.rows()));
  info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &E_small_vec[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  Jt_J_inv=Eigen::Map<Eigen::MatrixXd> (E_small_vec.data(), Ns,Ns);
  rhs=(mu*dt)*(A_r-Jt*Jt_J_inv*J*A_r);
  lhs=(-Jt*Jt_J_inv*J*A_sparse);
  lhs+=E*A_sparse;
  lhs*=(-mu*0.5*dt);
  lhs.diagonal().array() += 1;
  lhs_sparse=lhs.sparseView(1, 1e-13);
  //x_temp=lhs.fullPivLu().solve(rhs);
  //cout<<(mu*0.5*dt)*(A_r-Jt*Jt_J_inv*J*A_r)<<endl<<endl;abort();
  solver_sparse.analyzePattern(lhs_sparse);
  solver_sparse.factorize(lhs_sparse);
  x_temp= solver_sparse.solve(rhs);
  x_temp+=x;
//   cout<<"check norm!!!!!"<<endl;
//   for(size_t i=0; i<Ns; i++)
//   {
//     cout <<i<<" "<<dist(x_temp, i)<<endl;;
//   }
   double att, renorm, alpha1, alpha2;
   new_x=x_temp;
   for (size_t i=0; i<Ns; i++)
   {
     att=dist(x_temp, i);
     renorm=1.-1.0/(Ns)/(att);
    // cout<<i<<" "<<renorm<<endl;
     alpha1=(Np-(i+1))*renorm/Np;
     alpha2=(i+1)*renorm/Np;
      for(size_t j=0; j<=i; j++)
      {
        new_x(3*j)+=alpha1*(x_temp(3*(i+1))-x_temp(3*i));
        new_x(3*j+1)+=alpha1*(x_temp(3*(i+1)+1)-x_temp(3*i+1));
        new_x(3*j+2)+=alpha1*(x_temp(3*(i+1)+2)-x_temp(3*i+2));
      }
      for(size_t j=i+1; j<Np; j++)
      {
        new_x(3*j)-=alpha2*(x_temp(3*(i+1))-x_temp(3*i));
        new_x(3*j+1)-=alpha2*(x_temp(3*(i+1)+1)-x_temp(3*i+1));
        new_x(3*j+2)-=alpha2*(x_temp(3*(i+1)+2)-x_temp(3*i+2));
      }
   }
  // cout<<(Eigen::Map<Eigen::MatrixXd> (new_x.data(), 3,Np)).transpose()<<endl<<endl;;
//    double temp;
//    double norm=0;
//    for(size_t i=0; i<Ns; i++)
//    {
//      temp=dist(new_x, i);
//      norm+=temp;
//      cout <<i<<" "<<temp<<endl;;
//    }
//    cout<<"norm="<<norm<<endl;
//   abort();

}

void Integrator1::do_CN_propper(const Eigen::VectorXd &x ,  Eigen::VectorXd &new_x , const double  t , const double dt)
{
  int info;
  fill_Fa_mat(x, t+dt*0.5);
//   fill_Fa(x, t+dt*0.5);
//   if((FA_mat*x-FA).norm()>1)
//   {
//     cout <<FA<<endl<<endl;
//     cout <<x<<endl<<endl;
//     cout <<FA_mat<<endl<<endl;
//     cout <<FA_mat*x-FA<<endl<<endl;
//     cout<<(FA_mat*x-FA).norm()<<endl;
//    double *FA_sparse=FA_mat.valuePtr();
//    FA_sparse[27+27-3]=11;
//   for (int k=0; k<FA_mat.outerSize(); ++k)
//   {
//     int ii=0;
//     for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(FA_mat,k); it; ++it)
//     {
//       cout<<it.value()<<" ";
//       ii++;
// //       it.row();   // row index
// //       it.col();   // col index (here it is equal to k)
// //       it.index(); // inner index, here it is equal to it.row()
//     }
//     cout<<" ii="<<ii<<endl;
//   }
//     
//     abort();
//   }
  fill_J(x);
  //FA_mat+=A_sparse_row;
  A_r=(FA_mat+A_sparse_row)*x;
  E_small_vec=(Eigen::Map<Eigen::VectorXd>(E_small.data(), E_small.cols()*E_small.rows()));
  info = dgtsv(Ns, Ns, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &E_small_vec[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  Jt_J_inv=Eigen::Map<Eigen::MatrixXd> (E_small_vec.data(), Ns,Ns);
  rhs=(mu*0.5*dt)*(A_r-Jt*Jt_J_inv*J*A_r)+x;
  lhs=(-Jt*Jt_J_inv*J*(FA_mat+A_sparse_row));
  lhs+=E*FA_mat;
   lhs*=(-mu*0.5*dt);
   lhs.diagonal().array() += 1;
   lhs_sparse=lhs.sparseView(1, 1e-13);
   solver_sparse.analyzePattern(lhs_sparse);
   solver_sparse.factorize(lhs_sparse);
   x_temp= solver_sparse.solve(rhs);
   double att, renorm, alpha1, alpha2;
   new_x=x_temp;
   for (size_t i=0; i<Ns; i++)
   {
     att=dist(x_temp, i);
     renorm=1.-1.0/(Ns)/att;
     //cout<<i<<" "<<renorm<<endl;
     alpha1=(Np-(i+1))*renorm/Np;
     alpha2=(i+1)*renorm/Np;
      for(size_t j=0; j<=i; j++)
      {
        new_x(3*j)+=alpha1*(x_temp(3*(i+1))-x_temp(3*i));
        new_x(3*j+1)+=alpha1*(x_temp(3*(i+1)+1)-x_temp(3*i+1));
        new_x(3*j+2)+=alpha1*(x_temp(3*(i+1)+2)-x_temp(3*i+2));
      }
      for(size_t j=i+1; j<Np; j++)
      {
        new_x(3*j)-=alpha2*(x_temp(3*(i+1))-x_temp(3*i));
        new_x(3*j+1)-=alpha2*(x_temp(3*(i+1)+1)-x_temp(3*i+1));
        new_x(3*j+2)-=alpha2*(x_temp(3*(i+1)+2)-x_temp(3*i+2));
      }
   }
//   cout <<(Eigen::Map<Eigen::MatrixXd> (new_x.data(), 3,Np)).transpose()<<endl<<endl;;
//    double temp;
//    double norm=0;
//    for(size_t i=0; i<Ns; i++)
//    {
//      temp=dist(new_x, i);
//      norm+=temp;
//      cout <<i<<" "<<temp<<endl;;
//    }
//    cout<<"norm="<<norm<<endl;

}

double Integrator1::dist(const  Eigen::VectorXd &x , int i)
{
  return sqrt((x(3*(i+1))-x(3*i))*(x(3*(i+1))-x(3*i))+(x(3*(i+1)+1)-x(3*i+1))*(x(3*(i+1)+1)-x(3*i+1))+(x(3*(i+1)+2)-x(3*i+2))*(x(3*(i+1)+2)-x(3*i+2)));
}

Eigen::VectorXd Integrator1::fill_K1_IMEX_BDF2(const Eigen::VectorXd &x,const double dt)
{
  int info;
  fill_Fa(x,0);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", 0., info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  return K1;
  
}


void Integrator1::fill_mat_VSIMEX_BDF2( const double dt, const double dt1)
{
  omega1=dt/dt1;
  alpha2=(1+2*omega1)/(1+omega1);
  alpha1=-(1+omega1);
  alpha0=omega1*omega1/(1+omega1);
  beta1=(1+omega1);
  beta0=-omega1;
}

void Integrator1::fill_mat_VSIMEX_BDF3( const double dt, const double dt1, const double dt2)
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

void Integrator1::fill_mat_VSIMEX_BDF4( const double dt, const double dt1, const double dt2, const double dt3)
{
  omega3=dt/dt1;
  omega2=dt1/dt2;
  omega1=dt2/dt3;
  double A1=1+omega1*(1+omega2);
  double A2=1+omega2*(1+omega3);
  double A3=1+omega1*A2;
  alpha0=(1+omega3)/(1+omega1)*(A2/A1)*(omega1*omega1*omega1*omega1)*(omega2*omega2*omega2)*(omega3*omega3)/A3;
  alpha1=-(omega2*omega2*omega2)*(omega3*omega3)*(1+omega3)/(1+omega2)*(A3/A2);
  alpha2=omega3*(omega3/(1+omega3)+omega2*omega3*(A3+omega1)/(1+omega1));
  alpha3=-1-omega3*(1+omega2*(1+omega3)/(1+omega2)*(1+omega1*A2/A1));
  alpha4=1+omega3/(1+omega3)+omega2*omega3/A2+omega1*omega2*omega3/A3;;
  beta0=-omega1*omega1*omega1*omega2*omega2*omega3*A2/A1*(1+omega3)/(1+omega1);
  beta1=omega2*omega2*omega3*A3*(1+omega3)/(1+omega2);
  beta2=-A2*A3*(omega3)/(1+omega1);
  beta3=omega2*(1+omega3)/(1+omega2)/A1*((1+omega3)*(A3+omega1)+(1+omega1)/omega2);
}

void Integrator1::fill_lhs1_IMEX_BDF2(const double dt, const double dt1)
{
  fill_mat_VSIMEX_BDF2(dt, dt1);
  //double temp=1-M_SQRT1_2;
  double temp=(1.0/alpha2);
  lhs_sparse1=(-mu*dt*temp)*A_sparse;
  lhs_sparse1.diagonal().array() += 1;
  solver_sparse1.analyzePattern(lhs_sparse1);
  solver_sparse1.factorize(lhs_sparse1);
}

void Integrator1::fill_lhs1_VSIMEX_BDF4(const double dt, const double dt1, const double dt2, const double dt3)
{
  fill_mat_VSIMEX_BDF4(dt, dt1, dt2, dt3);
  //double temp=1-M_SQRT1_2;
  double temp=(1.0/alpha4);
  lhs_sparse1=(-mu*dt*temp)*A_sparse;
  lhs_sparse1.diagonal().array() += 1;
  solver_sparse1.analyzePattern(lhs_sparse1);
  solver_sparse1.factorize(lhs_sparse1);
}

void Integrator1::fill_lhs1_IMEX_ARS_232(const double dt)
{
  
  double temp=1-M_SQRT1_2;
  lhs_sparse1=(-mu*dt*temp)*A_sparse;
  lhs_sparse1.diagonal().array() += 1;
  solver_sparse1.analyzePattern(lhs_sparse1);
  solver_sparse1.factorize(lhs_sparse1);
}

void Integrator1::fill_lhs2_IMEX_ARS_232(const double dt)
{
  
  double temp=-2*M_SQRT2/3;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::fill_lhs2_VSIMEX_BDF3(const double dt, const double dt1, const double dt2)
{
  fill_mat_VSIMEX_BDF3(dt, dt1, dt2);
  //double temp=1-M_SQRT1_2;
  double temp=(1.0/alpha3);
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}


void Integrator1::do_VSIMEX_BDF2(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1 , const double  t , const double dt)
{
  int info;
  rhs=(beta0/alpha2)*K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=(-alpha1/alpha2)*x+(-alpha0/alpha2)*x_old+(beta1/alpha2)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse1.solve(rhs);
}

void Integrator1::do_VSIMEX_BDF3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt)
{
  int info;
  rhs=+(beta1/alpha3)*K1+(beta0/alpha3)*K2;
  K3=K2;
  K2=K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=-(alpha2/alpha3)*x-(alpha1/alpha3)*x_old-(alpha0/alpha3)*x_old_old+(beta2/alpha3)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}

void Integrator1::do_VSIMEX_BDF3_DD(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt)
{
  int info;
  rhs=+(beta1/alpha3)*K1+(beta0/alpha3)*K2;
  K3=K2;
  K2=K1;
  fill_Fa(x, t);
  //fill_Fr(x, t);
  fill_J_DD(x);
  fill_Ftwist(x);
  A_r=DD_sparse*(A_sparse*x+FA+Ftwist)/*+FR*/;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_DD_Jt_l[0], &J_DD_Jt_d[0], &J_DD_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(A_r-A_sparse*x/*+FR*/-DD_sparse*Jt*temp_r);
  rhs+=-(alpha2/alpha3)*x-(alpha1/alpha3)*x_old-(alpha0/alpha3)*x_old_old+(beta2/alpha3)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}


void Integrator1::do_VSIMEX_BDF4(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3 , const double  t , const double dt)
{
  int info;
  rhs=+(beta2/alpha4)*K1+(beta1/alpha4)*K2+(beta0/alpha4)*K3;
  K3=K2;
  K2=K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=-(alpha3/alpha4)*x-(alpha2/alpha4)*x_old-(alpha1/alpha4)*x_old_old-(alpha0/alpha4)*x_old_old2+(beta3/alpha4)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse1.solve(rhs);
}


void Integrator1::do_IMEX_BDF2(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1 , const double  t , const double dt)
{
  int info;
  rhs=-(2/3.)*K1;
  fill_Fa(x, t+dt*0.5);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=4/3.*x-1/3.*x_old+(4/3.)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}

void Integrator1::do_IMEX_BDF3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt)
{
  int info;
  rhs=+(6.0/11.0)*K2-(18.0/11.0)*K1;
  K2=K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=(18./11.)*x-(9./11.)*x_old+(2./11.)*x_old_old+(18./11.)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}

void Integrator1::do_IMEX_BDF5(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,const Eigen::VectorXd &x_old_old3 ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3,Eigen::VectorXd &K4 , const double  t , const double dt)
{
  int info;
  rhs=-(600.0/137.0)*K1+(600.0/137.0)*K2-(300.0/137.0)*K3+(60.0/137.0)*K4 ;
  K4=K3;
  K3=K2;
  K2=K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=(300./137.)*x-(300./137.)*x_old+(200./137.)*x_old_old-(75./137.)*x_old_old2+(12./137.)*x_old_old3 +(300./137.)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}

void Integrator1::do_IMEX_BDF4(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3 , const double  t , const double dt)
{
  int info;
  rhs=-(72.0/25.0)*K1+(48.0/25.0)*K2-(12.0/25.0)*K3;
  K3=K2;
  K2=K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=(48./25.)*x-(36./25.)*x_old+(16./25.)*x_old_old-(3./25.)*x_old_old2+(48./25.)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}



Eigen::VectorXd Integrator1::fill_K1_IMEX_BDF3(const Eigen::VectorXd &x,const double dt, const double t)
{
  int info;
  fill_Fa(x,t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", 0., info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  return K1;
}

Eigen::VectorXd Integrator1::fill_K1_IMEX_BDF3_DD(const Eigen::VectorXd &x,const double dt, const double t)
{
  int info;
  fill_Fa(x,t);
 // fill_Fr(x,t);
  fill_J_DD(x);
  fill_Ftwist(x);
  A_r=DD_sparse*(A_sparse*x+FA+Ftwist)/*+FR*/;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_DD_Jt_l[0], &J_DD_Jt_d[0], &J_DD_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", 0., info);
    abort();
  }
  K1=(mu*dt)*(A_r-A_sparse*x/*+FR*/-DD_sparse*Jt*temp_r);
  return K1;
}

void Integrator1::fill_lhs2_IMEX_BDF3(const double dt)
{
  //double temp=1-M_SQRT1_2;
  double temp=6./11.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::fill_lhs2_IMEX_BDF4(const double dt)
{
  //double temp=1-M_SQRT1_2;
  double temp=12./25.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::fill_lhs2_IMEX_BDF5(const double dt)
{
  //double temp=1-M_SQRT1_2;
  double temp=60./137.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::fill_lhs2_IMEX_TVB3(const double dt)
{
  //double temp=1-M_SQRT1_2;
  double temp=1089./2048.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::fill_lhs2_IMEX_TVB5(const double dt)
{
  //double temp=1-M_SQRT1_2;
  double temp=4007./8192.;
  lhs_sparse2=(-mu*dt*temp)*A_sparse;
  lhs_sparse2.diagonal().array() += 1;
  solver_sparse2.analyzePattern(lhs_sparse2);
  solver_sparse2.factorize(lhs_sparse2);
}

void Integrator1::do_IMEX_TVB3(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2 , const double  t , const double dt)
{
  int info;
  rhs=+(8233./12288.)*K2-(1271.0/768.0)*K1+A_sparse*((-1139./12288.*mu*dt)*x+(-367./6144.*mu*dt)*x_old+(1699./12288.*mu*dt)*x_old_old);
  K2=K1;
  fill_Fa(x, t+dt*0.5);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=(3909./2048.)*x-(1367./1024.)*x_old+(873./2048.)*x_old_old+(18463./12288.)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}

void Integrator1::do_IMEX_TVB5(const Eigen::VectorXd &x ,const Eigen::VectorXd &x_old ,const Eigen::VectorXd &x_old_old,const Eigen::VectorXd &x_old_old2,const Eigen::VectorXd &x_old_old3 ,  Eigen::VectorXd &new_x,  Eigen::VectorXd &K1,Eigen::VectorXd &K2,Eigen::VectorXd &K3,Eigen::VectorXd &K4 , const double  t , const double dt)
{
  int info;
  rhs=+-(13656497.0/2949120.0)*K1+(1249949.0/245760.0)*K2-(7937687.0/2949120.0)*K3+(3387361.0/5898240.0)*K4+A_sparse*((-4118249./5898240.*mu*dt)*x+(768703./2949120.*mu*dt)*x_old+(47849./245760.*mu*dt)*x_old_old+(-725087./2949120.*mu*dt)*x_old_old2 +(502321./5898240.*mu*dt)*x_old_old3);
  K4=K3;
  K3=K2;
  K2=K1;
  fill_Fa(x, t);
  fill_J(x);
  A_r=A_sparse*x+FA;
  temp_r=J*A_r;  
  info = dgtsv(Ns, 1, &J_Jt_l[0], &J_Jt_d[0], &J_Jt_u[0], &temp_r[0], Ns);
  if (info !=0)
  {
    cout << info<<endl;
    fprintf(stderr, "at t= %f dgtsv failure with error %d\n", t, info);
    abort();
  }
  K1=(mu*dt)*(FA-Jt*temp_r);
  rhs+=(13553./4096.)*x-(38121./8192.)*x_old+(7315./2048.)*x_old_old-(6161./4096.)*x_old_old2+(2269./8192.)*x_old_old3 +(10306951./5898240.)*K1;
//  cout<<"A!!!!!!";
  new_x= solver_sparse2.solve(rhs);
}

