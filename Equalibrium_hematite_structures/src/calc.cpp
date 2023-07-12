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

#include "calc.h"


void calc::set_val(const arma::vec& vals_inp)
{
  for(int i=1; i<n; i++)
  {
    r[i][0]=vals_inp(2*i-2);
    r[i][1]=0;
    rot[i](0,0)=cos(2*M_PI*vals_inp(2*n+i-3));
    rot[i](0,1)=-sin(2*M_PI*vals_inp(2*n+i-3));
    rot[i](1,0)=sin(2*M_PI*vals_inp(2*n+i-3));
    rot[i](1,1)=cos(2*M_PI*vals_inp(2*n+i-3));
    m[i]=rot[i]*m[0];
  }
  B[0]=cos(2*M_PI*vals_inp(3*n-3)); B[1]=sin(2*M_PI*vals_inp(3*n-3))*cos(2*M_PI*vals_inp(3*n-2));
  B[2]=sin(2*M_PI*vals_inp(3*n-3))*sin(2*M_PI*vals_inp(3*n-2));
  r10=r[n-1]-r[0];
}


void calc::set_val_gsl_new(const gsl_vector * vals_inp)
{
  double int_part, temp;
  for(int i=1; i<n; i++)
  {
    r[i][0]=r[i-1][0]+ fabs(gsl_vector_get(vals_inp,(i-1)));
  }
  B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(n-1))); B[1]=0;
  B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(n-1)));
  //   B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(3*n-3))); B[1]=sin(2*M_PI*gsl_vector_get(vals_inp,(3*n-3)))*cos(2*M_PI*gsl_vector_get(vals_inp,(3*n-2)));
  //   B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(3*n-3)))*sin(2*M_PI*gsl_vector_get(vals_inp,(3*n-2)));
  r10=r[n-1]-r[0];
}

double calc::Heviside(double x)
{
  if(x>0)
  {
    return x;
  }
  else
  {
    return 0.0;
  }
}

void calc::set_val_gsl_new_kink(const gsl_vector * vals_inp)
{
//    for(int i=1; i< n-kink; i++)
//   {
//     r[i][0]=r[i-1][0]+(gsl_vector_get(vals_inp,(i-1)));
//     r[i][1]=0;
//     r[i][2]=r[i-1][2]+1;
//   }
  for(int i=1; i<n-kink; i++)
  {
    r[i][0]=r[i-1][0]+fabs(gsl_vector_get(vals_inp,(i-1)));
    r[i][1]=0;
    r[i][2]=r[i-1][2]+1;
  }
  for(int i=n-kink; i<n; i++)
  {
    r[i][0]=r[i-1][0]+1;
    r[i][1]=0;
    r[i][2]=r[i-1][2]+gsl_vector_get(vals_inp,(i-1));
  }
  B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(n-1))); B[1]=0;
  B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(n-1)));
  //   B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(3*n-3))); B[1]=sin(2*M_PI*gsl_vector_get(vals_inp,(3*n-3)))*cos(2*M_PI*gsl_vector_get(vals_inp,(3*n-2)));
  //   B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(3*n-3)))*sin(2*M_PI*gsl_vector_get(vals_inp,(3*n-2)));
  r10=r[n-1]-r[0];
}

 


void calc::set_val_gsl(const gsl_vector * vals_inp)
{
  do_not_calc=false;
  //theta=2*M_PI*(-0.125);//gsl_vector_get(vals_inp,(0));
  theta=2*M_PI*gsl_vector_get(vals_inp,(0));
  phi=2*M_PI*gsl_vector_get(vals_inp,(1))/*0.09795663771925478*/;
  //pfi=2*M_PI*(0.25);//gsl_vector_get(vals_inp,(2));
  pfi=2*M_PI*gsl_vector_get(vals_inp,(2));

  rotyz(0,0)=cos(phi)*cos(theta); 
  rotyz(0,1)=-cos(phi)*sin(theta); 
  rotyz(0,2)=sin(phi);
  rotyz(1,0)=cos(theta)*sin(pfi)*sin(phi) +cos(pfi)*sin(theta); 
  rotyz(1,1)=cos(pfi)*cos(theta)-sin(pfi)*sin(phi)*sin(theta);
  rotyz(1,2)=-cos(phi)*sin(pfi);
  rotyz(2,0)=-cos(pfi)*cos(theta)*sin(phi) +sin(pfi)*sin(theta); 
  rotyz(2,1)=cos(theta)*sin(pfi) +cos(pfi)*sin(phi)*sin(theta); 
  rotyz(2,2)=cos(pfi)*cos(phi);
  X=rotyz*Xold;Y=rotyz*Yold;Z=rotyz*Zold;
  rot_point=-1; double zmin=1e100;
  for(int i=0; i<8; i++)
  {
    PRotoation_point[i]=rotyz*(Rotoation_point[i]);
    if(PRotoation_point[i][2]<zmin)
    {
      rot_point=i;
      zmin=PRotoation_point[i][2];
    }
  }
  r[0]=-PRotoation_point[rot_point];
  m[0]=rotyz*M;
  double int_part, temp;
  //Eigen::Matrix3d rot_matrix;
  
  switch(add_case)
  {
    case 1:
      for(int i=1; i<n; i++)
      {
        temp=(X[2] + gsl_vector_get(vals_inp,(2+i))*Y[2])/Z[2];
        r[i]=r[i-1]+X + gsl_vector_get(vals_inp,(2+i))*Y - temp*Z;
        if(first)
        {
          m[i]=Eigen::AngleAxisd(M_PI, X)*m[i-1];
        }
        else
        {
          m[i]=m[i-1];
        }
      }
      break;
    case 2:
      for(int i=1; i<n; i++)
      {temp=(Y[2] + gsl_vector_get(vals_inp,(2+i))*X[2])/Z[2];
        r[i]=r[i-1]+Y + gsl_vector_get(vals_inp,(2+i))*X - temp*Z;
        if(first)
        {
          m[i]=Eigen::AngleAxisd(M_PI, Y)*m[i-1];
        }
        else
        {
          m[i]=m[i-1];
        }
      }
      break;
    case 3:
      for(int i=1; i<n; i++)
      {
        if(fabs(X[2])<1e-10)
        {
           do_not_calc=true;
        }
        temp=(Z[2] + gsl_vector_get(vals_inp,(2+i))*Y[2])/X[2];
        r[i]=r[i-1]+Z + gsl_vector_get(vals_inp,(2+i))*Y - temp*X;
        if(first)
        {
          m[i]=Eigen::AngleAxisd(M_PI, Z)*m[i-1];
        }
        else
        {
          m[i]=m[i-1];
        }
      }
      break;
      case 4:
      for(int i=1; i<n; i++)
      {
        if(fabs(Y[2])<1e-10)
        {
           do_not_calc=true;
        }
        temp=(Z[2] + gsl_vector_get(vals_inp,(2+i))*X[2])/Y[2];
        r[i]=r[i-1]+Z + gsl_vector_get(vals_inp,(2+i))*X - temp*Y;
        if(first)
        {
          m[i]=Eigen::AngleAxisd(M_PI, Z)*m[i-1];
        }
        else
        {
          m[i]=m[i-1];
        }
      }
      break;
  }
  
//   for(int i=1; i<n; i++)
//   {
//     //if((fmod(theta,0.5)>=0.375)||(fmod(theta,1)<=-1e-3&&fmod(theta,1)>=-0.125001))
//     if(true)
//     {
//       //if((fmod(phi,0.5)>=-1e-3&&fmod(phi,0.5)<=0.126)||(fmod(theta,1)<=-0.375))
//       if(true)
//       {
//         add_case=1;
//         temp=(X[2] + gsl_vector_get(vals_inp,(2+i))*Y[2])/Z[2];
//         r[i]=r[i-1]+X + gsl_vector_get(vals_inp,(2+i))*Y - temp*Z;
//         if(first)
//         {
//           m[i]=Eigen::AngleAxisd(M_PI, X)*m[i-1];
//         }
//         else
//         {
//           m[i]=m[i-1];
//         }
//       }
// //       else
// //       {
// //         add_case=2;
// //         temp=(Z[2] + gsl_vector_get(vals_inp,(2+i))*Y[2])/X[2];
// //         r[i]=r[i-1]+Z + gsl_vector_get(vals_inp,(2+i))*Y - temp*X;
// //         if(first)
// //         {
// //           m[i]=Eigen::AngleAxisd(M_PI, Z)*m[i-1];
// //         }
// //         else
// //         {
// //           m[i]=m[i-1];
// //         }
// //       }
//     }
// //     else
// //     {
// //       add_case=3;
// //       if((fmod(phi,0.5)>=0&&fmod(phi,0.5)<=0.125)||(fmod(theta,1)<=-0.375))
// //       {
// //         temp=(Y[2] + gsl_vector_get(vals_inp,(2+i))*X[2])/Z[2];
// //         r[i]=r[i-1]+Y + gsl_vector_get(vals_inp,(2+i))*X - temp*Z;
// //         if(first)
// //         {
// //           m[i]=Eigen::AngleAxisd(M_PI, Y)*m[i-1];
// //         }
// //         else
// //         {
// //           m[i]=m[i-1];
// //         }
// //       }
// //       else
// //       {
// //         add_case=4;
// //         temp=(Z[2] + gsl_vector_get(vals_inp,(2+i))*X[2])/Y[2];
// //         r[i]=r[i-1]+Z + gsl_vector_get(vals_inp,(2+i))*X - temp*Y;
// //         if(first)
// //         {
// //           m[i]= Eigen::AngleAxisd(M_PI, Z)*m[i-1];
// //         }
// //         else
// //         {
// //           m[i]=m[i-1];
// //         }
// //       }
// //     }
//   }
  r10=r[n-1]-r[0];
}

void calc::set_val_gsl_smart(const gsl_vector * vals_inp)
{
  do_not_calc=false;
  theta=2*M_PI*(-0.125);
  phi=2*M_PI*gsl_vector_get(vals_inp,(0))/*0.09795663771925478*/;
  pfi=2*M_PI*(-0.25);

  rotyz(0,0)=cos(phi)*cos(theta); 
  rotyz(0,1)=-cos(phi)*sin(theta); 
  rotyz(0,2)=sin(phi);
  rotyz(1,0)=cos(theta)*sin(pfi)*sin(phi) +cos(pfi)*sin(theta); 
  rotyz(1,1)=cos(pfi)*cos(theta)-sin(pfi)*sin(phi)*sin(theta);
  rotyz(1,2)=-cos(phi)*sin(pfi);
  rotyz(2,0)=-cos(pfi)*cos(theta)*sin(phi) +sin(pfi)*sin(theta); 
  rotyz(2,1)=cos(theta)*sin(pfi) +cos(pfi)*sin(phi)*sin(theta); 
  rotyz(2,2)=cos(pfi)*cos(phi);
  X=rotyz*Xold;Y=rotyz*Yold;Z=rotyz*Zold;
  rot_point=-1; double zmin=1e100;
  for(int i=0; i<8; i++)
  {
    PRotoation_point[i]=rotyz*(Rotoation_point[i]);
    if(PRotoation_point[i][2]<zmin)
    {
      rot_point=i;
      zmin=PRotoation_point[i][2];
    }
  }
  r[0]=-PRotoation_point[rot_point];
  m[0]=rotyz*M;
  double int_part, temp;
  //Eigen::Matrix3d rot_matrix;
  
  for(int i=1; i<n; i++)
  {
    if(fabs(X[2])<1e-10)
    {
        do_not_calc=true;
    }
    temp=(Z[2] + gsl_vector_get(vals_inp,(i))*Y[2])/X[2];
    r[i]=r[i-1]+Z + gsl_vector_get(vals_inp,(i))*Y - temp*X;
    if(first)
    {
      m[i]=Eigen::AngleAxisd(M_PI, Z)*m[i-1];
    }
    else
    {
      m[i]=m[i-1];
    }
  }
  r10=r[n-1]-r[0];
}


void calc::set_val_gsl_smart_kink(const gsl_vector * vals_inp)
{
  do_not_calc=false;
  theta=2*M_PI*(-0.125);
  phi=2*M_PI*gsl_vector_get(vals_inp,(0))/*0.09795663771925478*/;
  pfi=2*M_PI*(-0.25);

  rotyz(0,0)=cos(phi)*cos(theta); 
  rotyz(0,1)=-cos(phi)*sin(theta); 
  rotyz(0,2)=sin(phi);
  rotyz(1,0)=cos(theta)*sin(pfi)*sin(phi) +cos(pfi)*sin(theta); 
  rotyz(1,1)=cos(pfi)*cos(theta)-sin(pfi)*sin(phi)*sin(theta);
  rotyz(1,2)=-cos(phi)*sin(pfi);
  rotyz(2,0)=-cos(pfi)*cos(theta)*sin(phi) +sin(pfi)*sin(theta); 
  rotyz(2,1)=cos(theta)*sin(pfi) +cos(pfi)*sin(phi)*sin(theta); 
  rotyz(2,2)=cos(pfi)*cos(phi);
  X=rotyz*Xold;Y=rotyz*Yold;Z=rotyz*Zold;
  rot_point=-1; double zmin=1e100;
  for(int i=0; i<8; i++)
  {
    PRotoation_point[i]=rotyz*(Rotoation_point[i]);
    if(PRotoation_point[i][2]<zmin)
    {
      rot_point=i;
      zmin=PRotoation_point[i][2];
    }
  }
  r[0]=-PRotoation_point[rot_point];
  m[0]=rotyz*M;
  double int_part, temp;
  //Eigen::Matrix3d rot_matrix;
  
  for(int i=1; i<n-kink; i++)
  {
    if(fabs(X[2])<1e-10)
    {
        do_not_calc=true;
    }
    temp=(Z[2] + fabs(gsl_vector_get(vals_inp,(i)))*Y[2])/X[2];
    r[i]=r[i-1]+Z + fabs(gsl_vector_get(vals_inp,(i)))*Y - temp*X;
    m[i]=m[i-1];
  }
  for(int i=n-kink; i<n; i++)
  {
    if(fabs(X[2])<1e-10)
    {
        do_not_calc=true;
    }
    temp=(Y[2] + gsl_vector_get(vals_inp,(i))*Z[2])/X[2];
    r[i]=r[i-1]+Y + gsl_vector_get(vals_inp,(i))*Z - temp*X;
    m[i]=m[i-1];
  }
  
  r10=r[n-1]-r[0];
}


void calc::set_val_gsl_superball(const gsl_vector * vals_inp)
{
  //cout<<"Here"<<endl;
  do_not_calc=false;
  r[0]=Eigen::Vector3d(0,0,0);
  m[0]=M;
  B=M;
  double int_part, temp;
  //Eigen::Matrix3d rot_matrix;
  for(int i=1; i<n; i++)
  {
        

    double b=gsl_vector_get(vals_inp,(i-1));
    double c=par_c_div_b*b;
   if(1 - pow((b*b),q_par)-pow((c*c),q_par)>0)//if(fabs(gsl_vector_get(vals_inp,(i-1))0)
    {
      r[i]=r[i-1]+Eigen::Vector3d(b,c,pow(1 - pow((b*b),q_par)-pow((c*c),q_par),0.5/q_par));//, 2*b, 
      m[i]=m[0];
    }
    else
    {
      do_not_calc=true;
    }
  }
  
  
  r10=r[n-1]-r[0];
}

void calc::set_val_gsl_cube(const gsl_vector * vals_inp)
{
  //cout<<"Here"<<endl;
  do_not_calc=false;
  r[0]=Eigen::Vector3d(0,0,0);
  m[0]=M;
  B=M;
  double int_part, temp;
  //Eigen::Matrix3d rot_matrix;
  for(int i=1; i<n; i++)
  {
        

    double b=gsl_vector_get(vals_inp,(i-1));
    double c=par_c_div_b*b;
   if(true)//if(fabs(gsl_vector_get(vals_inp,(i-1))0)
    {
      r[i]=r[i-1]+Eigen::Vector3d(b,c,1);//, 2*b, 
      m[i]=m[0];
    }
    else
    {
      do_not_calc=true;
    }
  }
  
  
  r10=r[n-1]-r[0];
}


// void calc::set_val_smart_smart_gsl(const gsl_vector * vals_inp)
// {
//   if(n%2==0)
//   {
//     for(int i=1; i<n/2; i++)
//     {
//       r[i][0]=gsl_vector_get(vals_inp,(i-1));
//       r[i][1]=0;
//       r[n-1-i][0]=r[i][0];
//       r[n-1-i][1]=r[i][1];
//     }
//   }
//   else
//   {
//     int n2=(n-1)/2;
//     for(int j=1; j<=n2; j++)
//     {
//       r[j][0]=gsl_vector_get(vals_inp,(j-1));
//       r[j][1]=0;
//     }
//     for(int j=1; j<=n2; j++)
//     {
//       r[n2+j][0]=2*r[n2][0]-r[n2-j][0];
//       r[n2+j][1]=0;
//     }
//   }
//   if(B_mag>1e-6)
//   {
//     if(n%2==0)
//     {
//       B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(n/2-1))); B[1]=sin(2*M_PI*gsl_vector_get(vals_inp,(n/2-1)))*cos(2*M_PI*gsl_vector_get(vals_inp,(n/2)));
//       B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(n/2-1)))*sin(2*M_PI*gsl_vector_get(vals_inp,(n/2)));
//     }
//     else
//     {
//       B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2))); B[1]=sin(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2)))*cos(2*M_PI*gsl_vector_get(vals_inp,((n+1)/2)));
//       B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2)))*sin(2*M_PI*gsl_vector_get(vals_inp,((n+1)/2)));
//     }
//   }
//   else
//   {
//     B[0]=0; B[1]=0; B[2]=0;
//   }
//   r10=r[n-1]-r[0];
// }

//2D magnetic field
void calc::set_val_smart_smart_gsl(const gsl_vector * vals_inp)
{
  if(n%2==0)
  {
    for(int i=1; i<n/2; i++)
    {
      r[i][0]=gsl_vector_get(vals_inp,(i-1));
      r[i][1]=gsl_vector_get(vals_inp,(i-1))*par_c_div_b;
      r[n-1-i][0]=r[i][0];
      r[n-1-i][1]=r[i][1];
    }
  }
  else
  {
    int n2=(n-1)/2;
    for(int j=1; j<=n2; j++)
    {
      r[j][0]=gsl_vector_get(vals_inp,(j-1));
      r[j][1]=gsl_vector_get(vals_inp,(j-1))*par_c_div_b;
    }
    for(int j=1; j<=n2; j++)
    {
      r[n2+j][0]=2*r[n2][0]-r[n2-j][0];
      r[n2+j][1]=2*r[n2][1]-r[n2-j][1];
    }
  }
  if(B_mag>1e-6)
  {
    if(n%2==0)
    {
      B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(n/2-1))); B[1]=0;
      B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(n/2-1)));
    }
    else
    {
      B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2))); B[1]=0;
      B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2)));
    }
  }
  else
  {
    B[0]=0; B[1]=0; B[2]=0;
  }
  r10=r[n-1]-r[0];
}



//2D limit
void calc::set_val_smart_smart_second_gsl(const gsl_vector * vals_inp)
{
  if(n%2==0)
  {
    for(int i=1; i<=n/2; i++)
    {
      r[i][0]=gsl_vector_get(vals_inp,(i-1));
      r[i][1]=gsl_vector_get(vals_inp,(i-1))*par_c_div_b;
    }
    int n2=n/2;
    for(int j=1; j<n2; j++)
    {
      r[n2+j][0]=r[n2][0]+r[n2-1][0]-r[n2-j-1][0];
      r[n2+j][1]=r[n2][1]+r[n2-1][1]-r[n2-j-1][1];
    }
  }
  else
  {
    int n2=(n-1)/2;
    for(int j=1; j<=n2; j++)
    {
      r[j][0]=gsl_vector_get(vals_inp,(j-1));
      r[j][1]=gsl_vector_get(vals_inp,(j-1))*par_c_div_b;
    }
    for(int j=1; j<=n2; j++)
    {
      r[n2+j][0]=2*r[n2][0]-r[n2-j][0];
      r[n2+j][1]=2*r[n2][1]-r[n2-j][1];
    }
  }
  if(B_mag>1e-6)
  {
    if(n%2==0)
    {
      B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,(n/2))); B[1]=0;
      B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,(n/2)));
    }
    else
    {
      B[0]=cos(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2))); B[1]=0;
      B[2]=sin(2*M_PI*gsl_vector_get(vals_inp,((n-1)/2)));
    }
  }
  else
  {
    B[0]=0; B[1]=0; B[2]=0;
  }
  r10=r[n-1]-r[0];
}


void calc::set_val_shift3(const arma::vec& vals_inp)
{
  for(int i=1; i<n-kink; i++)
  {
    r[i][0]=r[i-1][0]+fabs(vals_inp(2*i-2));
    r[i][1]=0;
    r[i][2]=r[i-1][2]+1;
  }
  for(int i=n-kink; i<n; i++)
  {
    r[i][0]=r[i-1][0]+1;
    r[i][1]=0;
    r[i][2]=r[i-1][2]+vals_inp(2*i-1);
  }
  r10=r[n-1]-r[0];
}


void calc::set_val_smart(const arma::vec& vals_inp)
{
  for(int i=1; i<n; i++)
  {
    r[i][0]=vals_inp(i-1);
    r[i][1]=0;
    rot[i](0,0)=cos(2*M_PI*vals_inp(n+i-2));
    rot[i](0,1)=-sin(2*M_PI*vals_inp(n+i-2));
    rot[i](1,0)=sin(2*M_PI*vals_inp(n+i-2));
    rot[i](1,1)=cos(2*M_PI*vals_inp(n+i-2));
    m[i]=rot[i]*m[0];
  }
  if(B_mag>1e-6)
  {
    B[0]=cos(2*M_PI*vals_inp(2*n-2)); B[1]=sin(2*M_PI*vals_inp(2*n-2))*cos(2*M_PI*vals_inp(2*n-1));
    B[2]=sin(2*M_PI*vals_inp(2*n-2))*sin(2*M_PI*vals_inp(2*n-1));
  }
  else
  {
    B[0]=0; B[1]=0; B[2]=0;
  }
  r10=r[n-1]-r[0];
}


void calc::set_val_smart_smart_second(const arma::vec& vals_inp)
{
  for(int i=1; i<n; i++)
  {
    r[i][0]=i*vals_inp(0);
    r[i][1]=0;
  }
  if(B_mag>1e-6)
  {
    B[0]=cos(2*M_PI*vals_inp(1)); B[1]=sin(2*M_PI*vals_inp(1))*cos(2*M_PI*vals_inp(2));
    B[2]=sin(2*M_PI*vals_inp(1))*sin(2*M_PI*vals_inp(2));
  }
  else
  {
    B[0]=0; B[1]=0; B[2]=0;
  }
  r10=r[n-1]-r[0];
}

void calc::set_val_smart_smart(const arma::vec& vals_inp)
{
  if(n%2==0)
  {
    for(int i=1; i<n/2; i++)
    {
      r[i][0]=vals_inp(i-1);
      r[i][1]=0;
      r[n-1-i][0]=vals_inp(i-1);
      r[n-1-i][1]=0;
    }
  }
  else
  {
    int n2=(n-1)/2;
    for(int j=1; j<=n2; j++)
    {
      r[j][0]=vals_inp(j-1);
      r[j][1]=0;
    }
    for(int j=1; j<=n2; j++)
    {
      r[n2+j][0]=2*r[n2][0]-r[n2-j][0];
      r[n2+j][1]=0;
    }
  }
  if(B_mag>1e-6)
  {
    if(n%2==0)
    {
      B[0]=cos(2*M_PI*vals_inp(n/2-1)); B[1]=sin(2*M_PI*vals_inp(n/2-1))*cos(2*M_PI*vals_inp(n/2));
      B[2]=sin(2*M_PI*vals_inp(n/2-1))*sin(2*M_PI*vals_inp(n/2));
    }
    else
    {
      B[0]=cos(2*M_PI*vals_inp((n-1)/2)); B[1]=sin(2*M_PI*vals_inp((n-1)/2))*cos(2*M_PI*vals_inp((n+1)/2));
      B[2]=sin(2*M_PI*vals_inp((n-1)/2))*sin(2*M_PI*vals_inp((n+1)/2));
    }
  }
  else
  {
    B[0]=0; B[1]=0; B[2]=0;
  }
  r10=r[n-1]-r[0];
}

void calc::set_val_w_grad(const arma::vec& vals_inp)
{
  for(int i=1; i<n; i++)
  {
    r[i][0]=vals_inp(2*i-2);
    r[i][1]=0;
    rot[i](0,0)=cos(2*M_PI*vals_inp(2*n+i-3));
    rot[i](0,1)=-sin(2*M_PI*vals_inp(2*n+i-3));
    rot[i](1,0)=sin(2*M_PI*vals_inp(2*n+i-3));
    rot[i](1,1)=cos(2*M_PI*vals_inp(2*n+i-3));
    rot_grad[i](0,0)=-sin(2*M_PI*vals_inp(2));
    rot_grad[i](0,1)=-cos(2*M_PI*vals_inp(2));
    rot_grad[i](1,0)=cos(2*M_PI*vals_inp(2));
    rot_grad[i](1,1)=-sin(2*M_PI*vals_inp(2));
    m[i]=rot[i]*m[0];
    m_grad[i]=rot_grad[i]*m[0];
  }
  B[0]=cos(2*M_PI*vals_inp(3*n-3)); B[1]=sin(2*M_PI*vals_inp(3*n-3))*cos(2*M_PI*vals_inp(3*n-2));
  B[2]=sin(2*M_PI*vals_inp(3*n-3))*sin(2*M_PI*vals_inp(3*n-2));
  B_grad1[0]=-sin(2*M_PI*vals_inp(3*n-3)); B_grad1[1]=cos(2*M_PI*vals_inp(3*n-3))*cos(2*M_PI*vals_inp(3*n-2));
  B_grad1[2]=cos(2*M_PI*vals_inp(3*n-3))*sin(2*M_PI*vals_inp(3*n-2));
  B_grad2[0]=0; B_grad2[1]=-sin(2*M_PI*vals_inp(3*n-3))*sin(2*M_PI*vals_inp(3*n-2));
  B_grad2[2]=sin(2*M_PI*vals_inp(3*n-3))*cos(2*M_PI*vals_inp(3*n-2));
  r10=r[n-1]-r[0];
}


double calc::calc_en()
{
  double obj_val=0;
  double tempx=0;
  for(int i=0; i<n; i++)
  {
    for(int j=i+1; j<n; j++)
    {
      obj_val+=dipole_energy(i,j);
    }
    tempx+=B.dot(m[i]);
  }
  obj_val-=B_mag*tempx;
  return obj_val;

}

double calc::calc_en_ncube(int ncube)
{
  double obj_val=0;
  double tempx=0;
  for(int i=0; i<n; i++)
  {
    for(int j=i+1; j<n; j++)
    {
      for(int j1=0; j1<ncube; j1++)
      {
        for(int j2=0; j2<ncube; j2++)
        {
          for(int j3=0; j3<ncube; j3++)
          {
            for(int i1=0; i1<ncube; i1++)
            {
              for(int i2=0; i2<ncube; i2++)
              {
                for(int i3=0; i3<ncube; i3++)
                {
                  obj_val+=dipole_energy_discret_cube(i,j, ncube, i1,i2,i3,j1,j2,j3);
                }
              }
            }
          }
        }
      }
    }
    tempx+=B.dot(m[i]);
  }
  obj_val/=pow(ncube,6);
  obj_val-=B_mag*tempx;
  return obj_val;
}




double calc::dipole_energy(int i, int j)
{
  tempr=r[i]-r[j];
  double dr10=tempr.norm();
  double dr10_2=dr10*dr10;
  double dr10_3=dr10_2*dr10;
  double dr10_5=dr10_3*dr10_2;
  return (m[i].dot(m[j])/dr10_3-3*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_5);
}

double calc::dipole_energy_discret_cube(int i, int j,int ncube, int i1, int i2, int i3, int j1, int j2, int j3)
{
  tempr=r[i]-r[j];
  tempr[0]+=percantage_sio2*(j1-i1)/ncube;
  tempr[1]+=percantage_sio2*(j2-i2)/ncube;
  tempr[2]+=percantage_sio2*(j3-i3)/ncube;
  double dr10=tempr.norm();
  double dr10_2=dr10*dr10;
  double dr10_3=dr10_2*dr10;
  double dr10_5=dr10_3*dr10_2;
  return (m[i].dot(m[j])/dr10_3-3*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_5);
}




double calc::fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val(vals_inp);
  return p->calc_en();
}

double calc::fn_gsl(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl(x);
  if(p->do_not_calc)
  {
    return 1e100;
  }
  return p->calc_en();
}

double calc::fn_gsl_smart(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl_smart(x);
  if(p->do_not_calc)
  {
    return 1e100;
  }
  return p->calc_en();
}
double calc::fn_gsl_smart_kink(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl_smart_kink(x);
  if(p->do_not_calc)
  {
    return 1e100;
  }
  return p->calc_en();
}

double calc::fn_gsl_superball(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl_superball(x);
  if(p->do_not_calc)
  {
    return 1e100;
  }
  return p->calc_en();
}

double calc::fn_gsl_cube(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl_cube(x);
  if(p->do_not_calc)
  {
    return 1e100;
  }
  return p->calc_en();
}



double calc::fn_gsl_new(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl_new(x);
  return p->calc_en();
}

double calc::fn_gsl_new_kink(const gsl_vector * x, void * opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_gsl_new_kink(x);
  return p->calc_en();
}


double calc::fn_shift3(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_shift3(vals_inp);
  return p->calc_en();
}

double calc::fn_smart(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_smart(vals_inp);
  return p->calc_en();
}
double calc::fn_smart_smart(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_smart_smart(vals_inp);
  return p->calc_en();
}

double calc::fn_smart_smart_gsl(const gsl_vector * x, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_smart_smart_gsl(x);
  return p->calc_en();
}

double calc::fn_smart_smart_gsl_second(const gsl_vector * x, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_smart_smart_second_gsl(x);
  return p->calc_en();
}

double calc::fn_smart_smart_second(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_smart_smart_second(vals_inp);
  return p->calc_en();
}

double calc::fn_with_grad(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_w_grad(vals_inp);
  return p->calc_en_w_grad(grad_out);
}
double calc::fn_with_grad1(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* p = reinterpret_cast<calc*>(opt_data);
  p->set_val_w_grad(vals_inp);
  for(int i=0; i<3*p->n-1; i++)
  { 
    cout<< vals_inp(i)<<" ";
  }
  cout<<endl;
  return p->calc_en_w_grad1(grad_out);
}


double calc::calc_en_w_grad(arma::vec* grad_out)
{
  double obj_val=0;
  for(int i=0; i<3*n-1; i++)
  {
    (*grad_out)(i)=0.;
  }
  double tempx=m[0][0];
  double tempy=m[0][1];
  double tempz=m[0][2];
  for(int j=1; j<n; j++)
  {
    obj_val+=dipole_energy_w_grad(0,j,(*grad_out)[2*j-2],(*grad_out)[2*j-1] ,(*grad_out)[2*n+j-1]);
  }
  for(int i=1; i<n; i++)
  {
    for(int j=i+1; j<n; j++)
    {
      obj_val+=dipole_energy_w_grad(i,j,(*grad_out)[2*i-2],(*grad_out)[2*i-1],(*grad_out)[2*j-2],(*grad_out)[2*j-1],(*grad_out)[2*n+i-1], (*grad_out)[2*n+j-1]);
    }
    tempx+=m[i][0];
    tempy+=m[i][1];
    tempz+=m[i][2];
  }
  obj_val-=B_mag*(B(0)*tempx+B(1)*tempy+B(2)*tempz);
  (*grad_out)[3*n-3]=-B_mag*(B_grad1[0]*tempx+B_grad1[1]*tempy+B_grad1[2]*tempz);
  (*grad_out)[3*n-2]=-B_mag*(B_grad2[0]*tempx+B_grad2[1]*tempy+B_grad2[2]*tempz);
  
//  cout<<r[1].transpose()<<" "<< obj_val<<" "<<(*grad_out)(0)<<" "<<(*grad_out)(1)<<" "<<(*grad_out)(2)<<" "<<(*grad_out)(3)<<" "<<(*grad_out)(4)<<endl; /*abort();*/

  return obj_val;

}

double calc::calc_en_w_grad1(arma::vec* grad_out)
{
  double obj_val=0;
  for(int i=0; i<3*n-1; i++)
  {
    (*grad_out)(i)=0.;
  }
  double tempx=m[0][0];
  double tempy=m[0][1];
  double tempz=m[0][2];
  for(int j=1; j<n; j++)
  {
    obj_val+=dipole_energy_w_grad(0,j,(*grad_out)[2*j-2],(*grad_out)[2*j-1] ,(*grad_out)[2*n+j-1]);
  }
  for(int i=1; i<n; i++)
  {
    for(int j=i+1; j<n; j++)
    {
      obj_val+=dipole_energy_w_grad(i,j,(*grad_out)[2*i-2],(*grad_out)[2*i-1],(*grad_out)[2*j-2],(*grad_out)[2*j-1],(*grad_out)[2*n+i-1], (*grad_out)[2*n+j-1]);
    }
    tempx+=m[i][0];
    tempy+=m[i][1];
    tempz+=m[i][2];
  }
  obj_val-=B_mag*(B(0)*tempx+B(1)*tempy+B(2)*tempz);
  (*grad_out)[3*n-3]=-B_mag*(B_grad1[0]*tempx+B_grad1[1]*tempy+B_grad1[2]*tempz);
  (*grad_out)[3*n-2]=-B_mag*(B_grad2[0]*tempx+B_grad2[1]*tempy+B_grad2[2]*tempz);
  cout<<"Here!!!!!"<<endl;
  for(int i=0; i<3*n-1; i++)
  { 
    cout<< (*grad_out)(i)<<" ";
  }
  cout<<endl;
 // abort();
  
//  cout<<r[1].transpose()<<" "<< obj_val<<" "<<(*grad_out)(0)<<" "<<(*grad_out)(1)<<" "<<(*grad_out)(2)<<" "<<(*grad_out)(3)<<" "<<(*grad_out)(4)<<endl; /*abort();*/

  return obj_val;

}

// double calc::fn_with_grad1(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
// {
//     calc* p = reinterpret_cast<calc*>(opt_data);
//     p->set_val_w_grad(vals_inp);
//     double dr10=p->r10.norm();
//     double dr10_2=dr10*dr10;
//     double dr10_3=dr10_2*dr10;
//     double dr10_5=dr10_3*dr10_2;
//     double dr10_7=dr10_5*dr10_2;
//     
//     double obj_val = p->m[0].dot(p->m[1])/dr10_3-3*p->m[0].dot(p->r10)*(p->m[1]).dot(p->r10)/dr10_5 -p->B_mag*(p->B.dot(p->m[0])+p->B.dot(p->m[1]));
//     
//  
//     //
//     if (grad_out) {
//         (*grad_out)(0) = -3*p->m[0].dot(p->m[1])*p->r10(0)/dr10_5 +15*p->r10(0)*p->m[0].dot(p->r10)*(p->m[1]).dot(p->r10)/dr10_7
//         -3*p->m[0](0)*(p->m[1]).dot(p->r10)/dr10_5-3*p->m[0].dot(p->r10)*p->m[1](0)/dr10_5;
//         (*grad_out)(1) = -3*p->m[0].dot(p->m[1])*p->r10(1)/dr10_5 +15*p->r10(1)*p->m[0].dot(p->r10)*p->m[1].dot(p->r10)/dr10_7
//         -3*p->m[0](1)*p->m[1].dot(p->r10)/dr10_5-3*p->m[0].dot(p->r10)*p->m[1](1)/dr10_5;
//         p->B[0]=-sin(vals_inp(3*p->n-3)); p->B[1]=cos(vals_inp(3*p->n-3))*cos(vals_inp(3*p->n-2));
//         p->B[2]=cos(vals_inp(3*p->n-3))*sin(vals_inp(3*p->n-2));
//         (*grad_out)(3)=-p->B_mag*(p->B.dot(p->m[0])+p->B.dot(p->m[1]));
//         p->B[0]=0.0; p->B[1]=-sin(vals_inp(3*p->n-3))*sin(vals_inp(3*p->n-2));
//         p->B[2]=sin(vals_inp(3*p->n-3))*cos(vals_inp(3*p->n-2));
//         (*grad_out)(4)=-p->B_mag*(p->B.dot(p->m[0])+p->B.dot(p->m[1]));        
//     }
//     //
//     cout<<p->r[1].transpose()<<" "<< obj_val<<" "<<(*grad_out)(0)<<" "<<(*grad_out)(1)<<" "<<(*grad_out)(2)<<" "<<(*grad_out)(3)<<" "<<(*grad_out)(4)<<endl; /*abort();*/
//     return obj_val;
// }


double calc::dipole_energy_w_grad(int i, int j, double & grad_ix, double & grad_iy, double & grad_jx, double & grad_jy,double & grad_iphi,double & grad_jphi)
{
  tempr=r[i]-r[j];
  double dr10=tempr.norm();
  double dr10_2=dr10*dr10;
  double dr10_3=dr10_2*dr10;
  double dr10_5=dr10_3*dr10_2;
  double dr10_7=dr10_5*dr10_2;
  double tempx=-3*m[i].dot(m[j])*tempr(0)/dr10_5 +15*tempr(0)*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_7
        -3*m[i](0)*(m[j]).dot(tempr)/dr10_5-3*m[i].dot(tempr)*m[j](0)/dr10_5;
  double tempy=-3*m[i].dot(m[j])*tempr(1)/dr10_5 +15*tempr(1)*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_7
        -3*m[i](1)*(m[j]).dot(tempr)/dr10_5-3*m[i].dot(tempr)*m[j](1)/dr10_5;
  grad_iphi+=m_grad[i].dot(m[j])/dr10_3-3*m_grad[i].dot(tempr)*m[j].dot(tempr)/dr10_5;
  grad_jphi+=m[i].dot(m_grad[j])/dr10_3-3*m[i].dot(tempr)*m_grad[j].dot(tempr)/dr10_5;
  grad_ix +=tempx;
  grad_iy +=tempy;
  grad_jx -=tempx;
  grad_jy -=tempy;
  return (m[i].dot(m[j])/dr10_3-3*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_5);  
}

double calc::dipole_energy_w_grad(int i, int j, double & grad_jx, double & grad_jy, double & grad_jphi)
{
  tempr=r[i]-r[j];
  double dr10=tempr.norm();
  double dr10_2=dr10*dr10;
  double dr10_3=dr10_2*dr10;
  double dr10_5=dr10_3*dr10_2;
  double dr10_7=dr10_5*dr10_2;
  double tempx=-3*m[i].dot(m[j])*tempr(0)/dr10_5 +15*tempr(0)*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_7
        -3*m[i](0)*(m[j]).dot(tempr)/dr10_5-3*m[i].dot(tempr)*m[j](0)/dr10_5;
  double tempy=-3*m[i].dot(m[j])*tempr(1)/dr10_5 +15*tempr(1)*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_7
        -3*m[i](1)*(m[j]).dot(tempr)/dr10_5-3*m[i].dot(tempr)*m[j](1)/dr10_5;
  grad_jx -=tempx;
  grad_jy -=tempy;
  grad_jphi+=m[i].dot(m_grad[j])/dr10_3-3*m[i].dot(tempr)*m_grad[j].dot(tempr)/dr10_5;
  return (m[i].dot(m[j])/dr10_3-3*m[i].dot(tempr)*(m[j]).dot(tempr)/dr10_5);  
}

calc::calc(int n1, double B_mag1, double phi):n(n1), B_mag(B_mag1), m(n), m_grad(n), r(n), rot(n), rot_grad(n), percantage_sio2(/*0.9*/0.1), kink(1),vec_x(3*n-1)
{
  double phi0=atan(1.0/sqrt(2.))*180/M_PI;
  cout.precision(16);
  double s;
  cout.precision(16);
  if(phi>55||phi<-35)
  {
    cout <<"angles from -35 till 55 deg. allowed"<<endl;
    abort();
  }
  else
  {
    s=1.0/sqrt(2.)*cos((phi0+phi)*M_PI/180);
    if(phi>0)
    {
      M=Eigen::Vector3d (s,s,sqrt(1-2*s*s));
      par_c_div_b=1.0;
    }
    else
    {
      M=Eigen::Vector3d (s,sqrt(1-2*s*s),s);
      par_c_div_b=M[1]/M[2];
    }
  }
  m[0]=M;
  r[0]=Eigen::Vector3d(0,0,0);
  theta=0; phi=0; pfi=0;
  rotyz=Eigen::Matrix3d::Identity();
  B=Eigen::Vector3d (1,0,0);
  Xold =Eigen::Vector3d (1,0,0);
  Yold = Eigen::Vector3d (0,1,0); Zold = Eigen::Vector3d (0,0,1); Centr = Eigen::Vector3d (0.,0.,0.); Rotoation_point.resize(8);
  Rotoation_point[0]=Eigen::Vector3d (0.5,0.5,-0.5);Rotoation_point[1]=Eigen::Vector3d (0.5,-0.5,-0.5);Rotoation_point[2]=Eigen::Vector3d (-0.5,-0.5,-0.5);Rotoation_point[3]=Eigen::Vector3d (-0.5,0.5,-0.5);
  Rotoation_point[4]=Eigen::Vector3d (0.5,0.5,0.5);Rotoation_point[5]=Eigen::Vector3d (0.5,-0.5,0.5);Rotoation_point[6]=Eigen::Vector3d (-0.5,-0.5,0.5);Rotoation_point[7]=Eigen::Vector3d (-0.5,0.5,0.5);
  PRotoation_point.resize(8);
  add_case=1;
  for(int i=1; i<n; i++)
  {
    r[i]=Eigen::Vector3d(0.25*i,0.25*i,i);
  }
  q_par=2.;
}

void calc::improve_energy_first()
{
  improve_energy_first_cout();
  double a=calc_en();
  improve_energy_first_();
  double b=calc_en();
  while(a-b>1e-14)
  {
    a=b;
    improve_energy_first_();
    b=calc_en();
  }
  improve_energy_first_cout();
}

void calc::improve_energy_second()
{
  cout<< "!!!!!!!!Here1"<<endl;
  improve_energy_second_();
  double a=calc_en();
  improve_energy_second_();
  double b=calc_en();
  cout<< "!!!!!!!!Here2"<<endl;
  while(a-b>1e-14)
  {
    a=b;
    improve_energy_second_();
    b=calc_en();
  }
  cout<< "!!!!!!!!Here3"<<endl;
  improve_energy_second_cout();
}


void calc::improve_energy_first_()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
    for(int j=1; j<n; j++)
    {
//         x(2*j-2)=0.5*(r[j][0]+r[j][1]);
//         x(2*j-2+1)=0.5*(r[j][0]+r[j][1]);
        x(2*j-2)=(r[j][0]);
        x(2*j-2+1)=(r[j][1]);
    }
    for(int i=2*n-2; i<3*n-3;i++)
    {
      if(i%2==0)
      {
        x(i)=0.5;
      }
      else
      {
        x(i)=0;
      }
    }
    if(n%2==0)
    {
      x(3*n-3)=-0.25;x(3*n-2)=-0.25;
    }
    else
    {
      x(3*n-3)=acos(B[0])/(2*M_PI);x(3*n-2)=acos(B[1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    }
 
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
   
    //bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

//     if (success) {
//         std::cout << "1de: Ackley test completed successfully.\n"
//                   << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     } else {
//         std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
//     }
//     const double pi = arma::datum::pi;
//     x(2)=fmod(x(2), 2*pi);
//     cout << "\n1de: solution :\n" << x << endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    /*for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;   */ 
}

void calc::improve_energy_second_()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
    for(int j=1; j<2; j++)
    {
        x(2*j-2)=0.5*(r[j][0]+r[j][1]);
        x(2*j-2+1)=0.5*(r[j][0]+r[j][1]);
    }
    for(int j=2; j<n; j++)
    {
      x(2*j-2)=j*x(0);
        x(2*j-2+1)=j*x(0);
    }
    for(int i=2*n-2; i<3*n-3;i++)
    {
      x(i)=0;
    }
    if(n%2==0)
    {
      x(3*n-3)=acos(m[0][0])/(2*M_PI);x(3*n-2)=acos(m[0][1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    }
    else
    {
      x(3*n-3)=acos(B[0])/(2*M_PI);x(3*n-2)=acos(B[1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    }
    //x(3*n-3)=0.5*M_PI;x(3*n-2)=0.5*M_PI;
 
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
   
    //bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

//     if (success) {
//         std::cout << "1de: Ackley test completed successfully.\n"
//                   << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     } else {
//         std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
//     }
//     const double pi = arma::datum::pi;
//     x(2)=fmod(x(2), 2*pi);
//     cout << "\n1de: solution :\n" << x << endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    /*for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;   */ 
}

void calc::improve_energy_second_cout()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
    for(int j=1; j<n; j++)
    {
        x(2*j-2)=0.5*(r[j][0]+r[j][1]);
        x(2*j-2+1)=0.5*(r[j][0]+r[j][1]);
    }
    for(int i=2*n-2; i<3*n-3;i++)
    {
      x(i)=0;
    }
    if(n%2==0)
    {
      x(3*n-3)=acos(m[0][0])/(2*M_PI);x(3*n-2)=acos(m[0][1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    }
    else
    {
      x(3*n-3)=acos(B[0])/(2*M_PI);x(3*n-2)=acos(B[1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    }    //x(3*n-3)=0.5*M_PI;x(3*n-2)=0.5*M_PI;
 
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
   
    //bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    cout << "\n1de: solution :\n";
    for(int i=0; i<x.n_elem; i++)
    { 
      cout<< x(i)<<" "<<x(i)/x(0)<<endl;
    } cout<<endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=1; i<n; i++)
    {
        cout <<r[i][0]-r[i-1][0]<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::improve_energy_first_cout()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
    for(int j=1; j<n; j++)
    {
        x(2*j-2)=r[j][0];
        x(2*j-2+1)=r[j][1];
    }
    for(int i=2*n-2; i<3*n-3;i++)
    {
      if(i%2==0)
      {
        x(i)=0.5;
      }
      else
      {
        x(i)=0;
      }
    }
    if(n%2==0)
    {
      x(3*n-3)=-0.25;x(3*n-2)=-0.25;
    }
    else
    {
      x(3*n-3)=acos(B[0])/(2*M_PI);x(3*n-2)=acos(B[1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    }
//     if(n==4)
//     {
//       x(0)=0.5*(r[0][0]+r[1][0]);
//       x(1)=x(0);
//       x(2)=x(0);
//       x(3)=x(0);
//       x(4)=0;
//       x(5)=0;
//     }
 
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
   
    //bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    //set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    for(int i=2*n-2; i<3*n-1;i++)
    {  
      x(i)=fmod(x(i), 1);
    }
    cout << "\n1de: solution :\n";
    for(int i=0; i<x.n_elem; i++)
    { 
      cout<< x(i)<<endl;
    }
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
    for(int j=1; j<n; j++)
    {
      for(int i=0; i<2; i++)
      {
        x(2*j-2+i)=r[j][i];
      }
    }
    for(int i=2*n-2; i<3*n-1; i++)
    {
        x(i)=0;
    }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    optim::algo_settings_t settings_1;
 
  settings_1.pso_center_particle = false;
//     settings_1.vals_bound = true;
//     arma::vec x_1 = arma::zeros(3*n,1);
//     arma::vec x_2 = arma::zeros(3*n,1);
//     for(int i=0; i<n-1; i++)
//     {
//       x_1(2*i)=-0.5*(n-1);
//       x_2(2*i)=0.5*(n-1);
//       x_1(2*i+1)=-0.5*(n-1);
//       x_2(2*i+1)=0.5*(n-1);
//     }
//     for(int i=2*n-2; i<3*n-1; i++)
//     {
//       x_1(i)=-M_PI;
//       x_2(i)=M_PI;
//     }
//      settings_1.lower_bounds = x_1;
//      settings_1.upper_bounds = x_2;
    //cout<<settings_1.pso_n_pop<<" "<<settings_1.pso_n_gen<<endl; abort();
    //cout<<settings_1.de_n_pop<<" "<<settings_1.de_n_gen<<" "<<settings_1.de_check_freq<<" "<<settings_1.de_mutation_method<<endl; abort();
    settings_1.de_n_pop=500;
    settings_1.de_n_gen=10000;
    //settings_1.de_check_freq=1000;
    settings_1.de_mutation_method=1;
    
//     settings_1.pso_n_pop = 4000;
//     settings_1.pso_n_gen = 4000;
    
    //bool success = optim::pso(x,fn,this,settings_1);
    bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0;x(1)=0; x(2)=M_PI; x(3)=0.5*M_PI;x(4)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    //bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    for(int i=2*n-2; i<3*n-1;i++)
    {  
      x(i)=fmod(x(i), 2*pi);
    }
    cout << "\n1de: solution :\n" << x << endl;
    fill_vec_x(x);
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}




void calc::minimize_energy_smart_smart_gsl( bool first)
{
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  int xsize;
  double size, old_size=1.0e100;
  cout << "start to initialize"<<endl;;
  if (first)
  {
    if(abs(B_mag)<1e-6)
    {
      if(n%2==0)
      {
        xsize=(n/2-1);
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
      else
      {
        xsize=((n-1)/2);
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
    }
    else
    {
      if(n%2==0)
      {
        xsize=((n/2));
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
      else
      {
        xsize=((n+1)/2);
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
    }
    if(n%2==0)
    {
//       for(int j=1; j<n/2; j++)
//       {
//           gsl_vector_set(x,(j-1),r[j][0]);
//           r[n-j-1][0]=r[j][0];
//       }
      r[n-1][0]=r[0][0];r[n-1][1]=r[0][1];
    }
    else
    {
//       int n2=(n-1)/2;
//       for(int j=1; j<=n2; j++)
//       {
//           gsl_vector_set(x,(j-1),r[j][0]);
//       }
//       for(int j=1; j<=n2; j++)
//       {
//         r[n2+j][0]=2*r[n2][0]-r[n2-j][0];
//       }
    }
  }
  else
  {
    //cout<<"Aborting!!"<<endl;abort();
    if(abs(B_mag)<1e-6)
    {
      if(n%2==0)
      {
        xsize=((n/2));
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
      else
      {
        xsize=((n+1)/2);
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
    }
    else
    {
      if(n%2==0)
      {
        xsize=((n/2)+1);
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
      else
      {
        xsize=((n+1)/2);
        x = gsl_vector_alloc (xsize);gsl_vector_set_all (x, 0.0);
      }
    }
  }
  if (first)
  {
    for(int j=1; j<n; j++)
    {
      if(j%2==0)
      {
        m[j]=m[0];
      }
      else
      {
        m[j][0]=-m[0][0];m[j][1]=-m[0][1];m[j][2]=m[0][2];
      }
    }
  }
  else
  {
    for(int j=1; j<n; j++)
    {
      if(j%2==0)
      {
        m[j]=m[0];
      }
      else
      {
        m[j][0]=m[0][0];m[j][1]=-m[0][1];m[j][2]=m[0][2];
      }
    }
  }
//   for(int i=n-1; i<2*n-2; i++)
//   {
//       x(i)=0;
//   }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

      /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (xsize);
  gsl_vector_set_all (ss, 0.01);
  /* Initialize method and iterate */
    
    
    
    if(first)
    {
      minex_func.n = xsize;
      minex_func.f = fn_smart_smart_gsl;
      minex_func.params = this;
      s = gsl_multimin_fminimizer_alloc (T, xsize);
      gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
      do
      {
        iter++;
        // cout <<iter<<endl;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status)
        {
          break;
        }

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, tol);

        if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }
        if(iter%20000==0)
        {
          if(fabs(old_size-calc_en())<1e-16)
          {
            status =GSL_SUCCESS; 
            //cout<<old_size<<" "<<calc_en()<<endl;
          }
          else
          {
            //cout<<old_size<<" "<<calc_en()<<endl;
            old_size=calc_en();
          }
        }
    //     if(iter%10==0)
    //     {
    //       cout <<iter<<" f() ="<<s->fval<<" "<<size<<endl;
    // //          printf ("%5d f() = %7.10f size = %.10f\n", 
    // //                   iter,
    // //                   s->fval, size);
    //     }
      }
      while (status == GSL_CONTINUE && iter < n_gen);
    }
    else
    {
      minex_func.n = xsize;
      minex_func.f = fn_smart_smart_gsl_second;
      minex_func.params = this;
      s = gsl_multimin_fminimizer_alloc (T, xsize);
      gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
      do
      {
        iter++;
        // cout <<iter<<endl;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status)
        {
          break;
        }

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, tol);

        if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }
        if(iter%20000==0)
        {
          if(fabs(old_size-calc_en())<1e-16)
          {
            status =GSL_SUCCESS; 
            //cout<<old_size<<" "<<calc_en()<<endl;
          }
          else
          {
            //cout<<old_size<<" "<<calc_en()<<endl;
            old_size=calc_en();
          }
        }
    //     if(iter%10==0)
    //     {
    //       cout <<iter<<" f() ="<<s->fval<<" "<<size<<endl;
    // //          printf ("%5d f() = %7.10f size = %.10f\n", 
    // //                   iter,
    // //                   s->fval, size);
    //     }
      }
      while (status == GSL_CONTINUE && iter < n_gen);

    }
    
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    if(first)
    {
      set_val_smart_smart_gsl(s->x);
    }
    else
    {
      //
    }
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "gsl SMART minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1gsl smart: solution :\n" << endl;
    cout<<"momnets"<<endl;
  for(int i=0; i<n; i++)
  {
    cout<< m[i].transpose()<<endl;
  }
    std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<endl;
    }
    cout<<"Bfield"<<endl;
     cout <<B.transpose()<<" "<<get_beta()<<" "<<acos(get_beta())*180/M_PI<<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<r[i].transpose()<<endl;
//     }
//     cout <<endl;
//     for(int i=1; i<n; i++)
//     {
//         cout <<r[i][0]-r[i-1][0]<<endl;
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//     cout <<B.transpose()<<endl;
//     cout<<"M dot B" <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].dot(B)<<" ";
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//        cout << acos(m[i].dot(B))*180/M_PI<<" ";
//     }
//     cout <<endl;
//     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;   
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);    
}


void calc::improve_energy_smart_smart_gsl( bool first)
{
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  int xsize=vec_smart_smart_x.size();
  double size, old_size=1.0e100;
  cout << "start to initialize"<<endl;;
  x = gsl_vector_alloc (xsize);
  for(int i=0; i<xsize; i++)
  {
    gsl_vector_set(x,i,vec_smart_smart_x[i]);
  }
  cout<< "initial energy="<<calc_en()<<endl;
  set_val_smart_smart_gsl(x);


  cout<< "initial energy="<<calc_en()<<endl;
  
//   if (first)
//   {
//     for(int j=1; j<n; j++)
//     {
//       if(j%2==0)
//       {
//         m[j]=m[0];
//       }
//       else
//       {
//         m[j][0]=-m[0][0];m[j][1]=-m[0][1];m[j][2]=m[0][2];
//       }
//     }
//   }
//   else
//   {
//     for(int j=1; j<n; j++)
//     {
//       m[j]=m[0];
//     }
//   }
//   for(int i=n-1; i<2*n-2; i++)
//   {
//       x(i)=0;
//   }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

      /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (xsize);
  gsl_vector_set_all (ss, 0.01);
  /* Initialize method and iterate */
  minex_func.n = xsize;
  minex_func.f = fn_smart_smart_gsl;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, xsize);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    
    
    if(first)
    {
      do
      {
        iter++;
        // cout <<iter<<endl;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status)
        {
          break;
        }

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, tol);

        if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }
        if(iter%20000==0)
        {
          if(fabs(old_size-calc_en())<1e-16)
          {
            status =GSL_SUCCESS; 
            //cout<<old_size<<" "<<calc_en()<<endl;
          }
          else
          {
            //cout<<old_size<<" "<<calc_en()<<endl;
            old_size=calc_en();
          }
        }
    //     if(iter%10==0)
    //     {
    //       cout <<iter<<" f() ="<<s->fval<<" "<<size<<endl;
    // //          printf ("%5d f() = %7.10f size = %.10f\n", 
    // //                   iter,
    // //                   s->fval, size);
    //     }
      }
      while (status == GSL_CONTINUE && iter < n_gen);
    }
    else
    {
      //
    }
    
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    if(first)
    {
      set_val_smart_smart_gsl(s->x);
    }
    else
    {
      //
    }
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "De minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1de: solution :\n" << endl;
  for(int i=0; i<xsize; i++)
  {
    cout<< gsl_vector_get(s->x,i)<<endl;
  }
  cout << "\n1gsl: solution :\n" << endl;
  for(int i=0; i<n; i++)
  {
    cout<< r[i].transpose()<<endl;
  }
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);    
}

void calc::minimize_energy_gsl_new()
{
  int n_gen=1000000; double tol=1e-18;
  //cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2/*rand*/;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n);
  /* Set initial step sizes to 1 */
  /* Initialize method and iterate */
  minex_func.n = n;
  minex_func.f = fn_gsl_new;
  minex_func.params = this;
   ss = gsl_vector_alloc (n);
  s = gsl_multimin_fminimizer_alloc (T, n);
 
  int conf_count= pow(4,n-1);
  vector <double> conf_energy; conf_energy.resize(conf_count);
  vector <double> rotangle;  rotangle.resize(n-1);
  for(int ii=0; ii<conf_count; ii++)
  {
    iter = 0;
    old_size=1.0e100;
    gsl_vector_set_all (x, 0.0);
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    for(int i=1; i<n; i++)
    {
      r[i]=Eigen::Vector3d(0.25*i,0,i);
    }
     
    for(int j=1; j<n; j++)
    {
      gsl_vector_set(x,(j-1),r[j][0]-r[j-1][0]);
    }
    gsl_vector_set(x,(n-1),0.25);
    B[0]=cos(2*M_PI*gsl_vector_get(x,(n-1))); B[1]=0;
    B[2]=sin(2*M_PI*gsl_vector_get(x,(n-1)));
    int l_ii=ii;
    for(int j=0; j<n-1; j++)
    {
      rotangle[j]=l_ii%4;
//      cout<<ii<<" "<< j<<" "<<l_ii%4<<" "<<rotangle[j]<<endl;
      l_ii-=rotangle[j];
      l_ii/=4;
      rot[j+1](0,0)=cos(0.5*M_PI*rotangle[j]);
      rot[j+1](0,1)=-sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,0)=sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,1)=cos(0.5*M_PI*rotangle[j]);
      m[j+1]=rot[j+1]*m[0];
    }
//    cout<< ii<<" :";
//     for(int i=0; i<n-1; i++)
//     {
//       cout<<rotangle[i]<<" ";
//     }
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    
    set_val_gsl_new(s->x);
    conf_energy[ii]=calc_en();
//     cout<<" "<<calc_en()<<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<r[i].transpose()<<endl;
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//    cout <<B.transpose()<<endl;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout << "GSL minimization completed successfully.\n"
//                     << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     cout << "\n1gsl: solution :\n" << endl;
  }
  
  for(int i=0; i<conf_count; i++)
  {
    cout <<i<< " "<<conf_energy[i]<<" :";
    int l_ii=i;
    for(int j=0; j<n-1; j++)
    {
      rotangle[j]=l_ii%4;
      l_ii-=rotangle[j];
      l_ii/=4;
    }
    for(int j=0; j<n-1; j++)
    {
      cout <<rotangle[j]<<" ";
    }
    cout<<endl;
 }
  int minElementIndex = std::min_element(conf_energy.begin(),conf_energy.end()) - conf_energy.begin();
  cout<< "min enegy index="<<minElementIndex<<endl;
  
  int ii=minElementIndex;
  {
        iter = 0;
        old_size=1.0e100;
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    for(int i=1; i<n; i++)
    {
      r[i]=Eigen::Vector3d(0.25*i,0,i);
    }
     
    for(int j=1; j<n; j++)
    {
      gsl_vector_set(x,(j-1),r[j][0]);
    }
    gsl_vector_set(x,(n-1),0.25);
    B[0]=cos(2*M_PI*gsl_vector_get(x,(n-1))); B[1]=0;
    B[2]=sin(2*M_PI*gsl_vector_get(x,(n-1)));
    int l_ii=ii;
    for(int j=0; j<n-1; j++)
    {
      rotangle[j]=l_ii%4;
      l_ii-=rotangle[j];
      l_ii/=4;
      rot[j+1](0,0)=cos(0.5*M_PI*rotangle[j]);
      rot[j+1](0,1)=-sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,0)=sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,1)=cos(0.5*M_PI*rotangle[j]);
      m[j+1]=rot[j+1]*m[0];
    }
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    
      cout<<"!!!!!!!!!!!!done"<<endl;
    set_val_gsl_new(s->x);
    cout<<" Energy"<<calc_en()<<endl;;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "GSL minimization completed successfully.\n"
                    << "elapsed time: " << elapsed_seconds.count() << "s\n";
    cout << "\n1gsl: solution :\n" << endl;
    for(int i=0; i<n; i++)
  {
    cout<< r[i].transpose()<<endl;
  }
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<"   "<<r[i][0]-r[i-1][0]<<endl;
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<endl;
    }
    cout <<endl;
    cout <<B.transpose()<<" "<<get_Bangle()<<" "<<acos(get_Bangle())*180/M_PI<<endl;
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].dot(B)<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
       cout << acos(m[i].dot(B))*180/M_PI<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}


void calc::minimize_energy_gsl_new_kink(int kink1)
{
  kink=kink1;
  int n_gen=1000000; double tol=1e-18;
  //cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n);
  /* Set initial step sizes to 1 */
  /* Initialize method and iterate */
  minex_func.n = n;
  minex_func.f = fn_gsl_new_kink;
  minex_func.params = this;
   ss = gsl_vector_alloc (n);
  s = gsl_multimin_fminimizer_alloc (T, n);
 
  int conf_count= pow(4,n-1);
  vector <double> conf_energy; conf_energy.resize(conf_count);
  vector <double> rotangle;  rotangle.resize(n-1);
  for(int ii=0; ii<conf_count; ii++)
  {
    iter = 0;
    old_size=1.0e100;
    gsl_vector_set_all (x, 0.0);
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    for(int i=1; i<n-kink; i++)
    {
      r[i]=r[i-1]+Eigen::Vector3d(0.0,0,i);
    }
    for(int i=n-kink; i<n; i++)
    {
      r[i][0]=r[i-1][0]+1;
      r[i][1]=0;
      r[i][2]=r[i-1][2]+0.25;
    }
     
//     for(int j=1; j<n; j++)
//     {
//       gsl_vector_set(x,(j-1),r[j][0]-r[j-1][0]);
//     }
    gsl_vector_set(x,(n-1),0.25);
    B[0]=cos(2*M_PI*gsl_vector_get(x,(n-1))); B[1]=0;
    B[2]=sin(2*M_PI*gsl_vector_get(x,(n-1)));
    int l_ii=ii;
    for(int j=0; j<n-1; j++)
    {
      rotangle[j]=l_ii%4;
//      cout<<ii<<" "<< j<<" "<<l_ii%4<<" "<<rotangle[j]<<endl;
      l_ii-=rotangle[j];
      l_ii/=4;
      rot[j+1](0,0)=cos(0.5*M_PI*rotangle[j]);
      rot[j+1](0,1)=-sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,0)=sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,1)=cos(0.5*M_PI*rotangle[j]);
      m[j+1]=rot[j+1]*m[0];
    }
//     for(int j=1; j<n; j++)
//     {
//       if(j%2==0)
//       {
//         m[j]=m[0];
//       }
//       else
//       {
//         m[j][0]=m[0][0];m[j][1]=-m[0][1];m[j][2]=m[0][2];
//       }
//     }
//    cout<< ii<<" :";
//     for(int i=0; i<n-1; i++)
//     {
//       cout<<rotangle[i]<<" ";
//     }
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    
    set_val_gsl_new_kink(s->x);
    conf_energy[ii]=calc_en();
//     cout<<" "<<calc_en()<<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<r[i].transpose()<<endl;
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//    cout <<B.transpose()<<endl;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout << "GSL minimization completed successfully.\n"
//                     << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     cout << "\n1gsl: solution :\n" << endl;
  }
  
  for(int i=0; i<conf_count; i++)
  {
    cout <<i<< " "<<conf_energy[i]<<" :";
    int l_ii=i;
    for(int j=0; j<n-1; j++)
    {
      rotangle[j]=l_ii%4;
      l_ii-=rotangle[j];
      l_ii/=4;
    }
    for(int j=0; j<n-1; j++)
    {
      cout <<rotangle[j]<<" ";
    }
    cout<<endl;
 }
  int minElementIndex = std::min_element(conf_energy.begin(),conf_energy.end()) - conf_energy.begin();
  cout<< "min enegy index="<<minElementIndex<<endl;
  
  int ii=minElementIndex;
  {
        iter = 0;
        old_size=1.0e100;
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    for(int i=1; i<n; i++)
    {
      r[i]=Eigen::Vector3d(0.25*i,0,i);
    }
     
    for(int j=1; j<n; j++)
    {
      gsl_vector_set(x,(j-1),r[j][0]);
    }
    gsl_vector_set(x,(n-1),0.25);
    B[0]=cos(2*M_PI*gsl_vector_get(x,(n-1))); B[1]=0;
    B[2]=sin(2*M_PI*gsl_vector_get(x,(n-1)));
    int l_ii=ii;
    for(int j=0; j<n-1; j++)
    {
      rotangle[j]=l_ii%4;
      l_ii-=rotangle[j];
      l_ii/=4;
      rot[j+1](0,0)=cos(0.5*M_PI*rotangle[j]);
      rot[j+1](0,1)=-sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,0)=sin(0.5*M_PI*rotangle[j]);
      rot[j+1](1,1)=cos(0.5*M_PI*rotangle[j]);
      m[j+1]=rot[j+1]*m[0];
    }
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    
      cout<<"!!!!!!!!!!!!done"<<endl;
    set_val_gsl_new_kink(s->x);
    cout<<" Energy"<<calc_en()<<endl;;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "GSL minimization completed successfully.\n"
                    << "elapsed time: " << elapsed_seconds.count() << "s\n";
    cout << "\n1gsl: solution :\n" << endl;
    for(int i=0; i<n; i++)
  {
    cout<< r[i].transpose()<<endl;
  }
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<"   "<<r[i][0]-r[i-1][0]<<endl;
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<endl;
    }
    cout <<endl;
    cout <<B.transpose()<<" "<<get_Bangle()<<" "<<acos(get_Bangle())*180/M_PI<<endl;
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].dot(B)<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
       cout << acos(m[i].dot(B))*180/M_PI<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;
    cout <<acos(get_beta())*180/M_PI<<" "<<acos(get_beta1())*180/M_PI<<" "<<acos(get_beta2())*180/M_PI<<endl;
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}



void calc::minimize_energy_gsl_new_kink_smart()
{
  int n_gen=1000000; double tol=1e-18;
  //cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n);
  /* Set initial step sizes to 1 */
  /* Initialize method and iterate */
  minex_func.n = n;
  minex_func.f = fn_gsl_new_kink;
  minex_func.params = this;
   ss = gsl_vector_alloc (n);
  s = gsl_multimin_fminimizer_alloc (T, n);
  for(int j=1; j<n; j++)
  {
    if(j%2==0)
    {
      m[j]=m[0];
    }
    else
    {
      m[j][0]=m[0][0];m[j][1]=-m[0][1];m[j][2]=m[0][2];
    }
  }
 
  int conf_count= n-1;
  vector <double> conf_energy; conf_energy.resize(conf_count);
  vector <double> vbeta1, vbeta2; vbeta1.resize(conf_count); vbeta2.resize(conf_count);
  vector <double> rotangle;  rotangle.resize(n-1);
  for(int ii=0; ii<conf_count; ii++)
  {
    kink=ii;
    iter = 0;
    old_size=1.0e100;
    gsl_vector_set_all (x, 0.0);
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    for(int i=1; i<n-kink; i++)
    {
      r[i]=r[i-1]+Eigen::Vector3d(0.0,0,i);
    }
    for(int i=n-kink; i<n; i++)
    {
      r[i][0]=r[i-1][0]+1;
      r[i][1]=0;
      r[i][2]=r[i-1][2]+0.25;
    }
     
//     for(int j=1; j<n; j++)
//     {
//       gsl_vector_set(x,(j-1),r[j][0]-r[j-1][0]);
//     }
    gsl_vector_set(x,(n-1),0.25);
    B[0]=cos(2*M_PI*gsl_vector_get(x,(n-1))); B[1]=0;
    B[2]=sin(2*M_PI*gsl_vector_get(x,(n-1)));
    int l_ii=ii;

//    cout<< ii<<" :";
//     for(int i=0; i<n-1; i++)
//     {
//       cout<<rotangle[i]<<" ";
//     }
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    
    set_val_gsl_new_kink(s->x);
    conf_energy[ii]=calc_en();
    vbeta1[ii]=acos(get_beta1())*180/M_PI;
    vbeta2[ii]=acos(get_beta2())*180/M_PI;
//     cout<<" "<<calc_en()<<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<r[i].transpose()<<endl;
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//    cout <<B.transpose()<<endl;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout << "GSL minimization completed successfully.\n"
//                     << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     cout << "\n1gsl: solution :\n" << endl;
  }
  
  for(int i=0; i<conf_count; i++)
  {
    cout <<i<< " "<<conf_energy[i]<<" "<<vbeta1[i]<<" "<<vbeta2[i];
    cout<<endl;
 }
  int minElementIndex = std::min_element(conf_energy.begin(),conf_energy.end()) - conf_energy.begin();
  cout<< "min enegy index="<<minElementIndex<<endl;
  
  kink=minElementIndex;
  {
        iter = 0;
        old_size=1.0e100;
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    for(int i=1; i<n; i++)
    {
      r[i]=Eigen::Vector3d(0.25*i,0,i);
    }
     
    for(int j=1; j<n; j++)
    {
      gsl_vector_set(x,(j-1),r[j][0]);
    }
    gsl_vector_set(x,(n-1),0.25);
    B[0]=cos(2*M_PI*gsl_vector_get(x,(n-1))); B[1]=0;
    B[2]=sin(2*M_PI*gsl_vector_get(x,(n-1)));
    
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    
      cout<<"!!!!!!!!!!!!done"<<endl;
    set_val_gsl_new_kink(s->x);
    cout<<" Energy"<<calc_en()<<endl;;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "GSL minimization completed successfully.\n"
                    << "elapsed time: " << elapsed_seconds.count() << "s\n";
    cout << "\n1gsl: solution :\n" << endl;
    for(int i=0; i<n; i++)
  {
    cout<< r[i].transpose()<<endl;
  }
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<"   "<<r[i][0]-r[i-1][0]<<endl;
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<endl;
    }
    cout <<endl;
    cout <<B.transpose()<<" "<<get_Bangle()<<" "<<acos(get_Bangle())*180/M_PI<<endl;
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].dot(B)<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
       cout << acos(m[i].dot(B))*180/M_PI<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;
    cout <<acos(get_beta())*180/M_PI<<" "<<acos(get_beta1())*180/M_PI<<" "<<acos(get_beta2())*180/M_PI<<endl;
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

void calc::minimize_energy_gsl(bool first1)
{
  cout<<"!!!!!!!!!minimize_energy_gsl"<<endl;
  first=first1;
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (2+n);
  
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (2+n);
  /* Initialize method and iterate */
  minex_func.n =2+n;
  minex_func.f = fn_gsl;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, 2+n);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  int conf_count= 4;
  vector <double> conf_energy; conf_energy.resize(conf_count);
  for(int ii=0; ii<conf_count; ii++)
  {
    cout <<"!!!!!testnig case="<<ii+1<<endl;
    add_case=ii+1;
    iter = 0;
    old_size=1.0e100;
    gsl_vector_set_all (x, 0.0);
    gsl_vector_set_all (ss, 0.01);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
   
    do
    {
      iter++;
      // cout <<iter<<endl;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status)
      {
        break;
      }

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, tol);

      if (status == GSL_SUCCESS)
      {
        printf ("converged to minimum at\n");
      }
      if(iter%20000==0)
      {
        if(fabs(old_size-calc_en())<1e-16)
        {
          status =GSL_SUCCESS; 
        }
        else
        {
          old_size=calc_en();
        }
      }
    }
    while (status == GSL_CONTINUE && iter < n_gen);
    set_val_gsl(s->x);
    conf_energy[ii]=calc_en();
  }
  cout<<"configuration energies!"<<endl;
  for(int i=0; i<conf_count; i++)
  {
    cout <<i<< " "<<conf_energy[i]<<" :";
    cout<<endl;
  }
  int minElementIndex = std::min_element(conf_energy.begin(),conf_energy.end()) - conf_energy.begin();
  cout<< "min enegy index="<<minElementIndex<<endl;
  
  
  add_case=minElementIndex+1;
  iter = 0;
  old_size=1.0e100;
  gsl_vector_set_all (x, 0.0);
  gsl_vector_set_all (ss, 0.01);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  do
  {
    iter++;
    // cout <<iter<<endl;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status)
    {
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
    if(iter%20000==0)
    {
      if(fabs(old_size-calc_en())<1e-16)
      {
        status =GSL_SUCCESS; 
      }
      else
      {
        old_size=calc_en();
      }
    }
  }
  while (status == GSL_CONTINUE && iter < n_gen);
  set_val_gsl(s->x);
  for(int i=0; i<2+n; i++)
  {
    cout<<"i="<<i<<" "<< gsl_vector_get(s->x,i)<<endl;
  }
  
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "GSL minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1gsl: solution :\n" << endl;
  cout<<"momnets"<<endl;
  for(int i=0; i<n; i++)
  {
    cout<< m[i].transpose()<<endl;
  }
    std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<endl;
    }
    cout<<"Bfield"<<endl;
     cout <<B.transpose()<<" "<<get_beta()<<" "<<acos(get_beta())*180/M_PI<<endl;
     cout<<"Angle"<<endl;
     cout <<theta<<" "<<phi<<" "<<pfi<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//     cout <<B.transpose()<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].dot(B)<<" ";
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//        cout << acos(m[i].dot(B))*180/M_PI<<" ";
//     }
//     cout <<endl;
//     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

void calc::minimize_energy_gsl_smart_kink(int kink1)
{
  cout<<"!!!!!!!!!minimize_energy_gsl"<<endl;
  kink=kink1;
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n);
  
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (n);
  /* Initialize method and iterate */
  minex_func.n =n;
  minex_func.f = fn_gsl_smart_kink;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, n);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    
  iter = 0;
  old_size=1.0e100;
  gsl_vector_set_all (x, 0.0);
  gsl_vector_set_all (ss, 0.01);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  do
  {
    iter++;
    // cout <<iter<<endl;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status)
    {
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
    if(iter%20000==0)
    {
      if(fabs(old_size-calc_en())<1e-16)
      {
        status =GSL_SUCCESS; 
      }
      else
      {
        old_size=calc_en();
      }
    }
  }
  while (status == GSL_CONTINUE && iter < n_gen);
  set_val_gsl_smart_kink(s->x);
  for(int i=0; i<n; i++)
  {
    cout<<"i="<<i<<" "<< gsl_vector_get(s->x,i)<<endl;
  }
  
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "GSL minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1gsl: solution :\n" << endl;
  cout<<"momnets"<<endl;
  for(int i=0; i<n; i++)
  {
    cout<< m[i].transpose()<<endl;
  }
    std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<endl;
    }
    cout<<"Bfield"<<endl;
     cout <<B.transpose()<<" "<<get_beta()<<" "<<acos(get_beta())*180/M_PI<<" "<<acos(get_beta1())*180/M_PI<<" "<<acos(get_beta2())*180/M_PI<<endl;
     vector <double> betas=get_betas();
    for(int i=0; i<betas.size();i++)
    {
      cout <<i<<" " <<betas[i]<<endl;
    }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//     cout <<B.transpose()<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].dot(B)<<" ";
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//        cout << acos(m[i].dot(B))*180/M_PI<<" ";
//     }
//     cout <<endl;
//     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

vector <double> calc::get_betas()
{
  Eigen::Vector3d dr;
  vector <double> betas; betas.resize(n-1);
  for(int i=0; i<n-1; i++)
  {
    dr=r[i+1]-r[i];
    if((dr.cross(B))[1]>0)
    {
      betas[i]=acos(B.dot(dr)/dr.norm())/M_PI*180;
    }
    else
    {
      betas[i]=-acos(B.dot(dr)/dr.norm())/M_PI*180;
    }
    
  }
  return betas; 
}


void calc::minimize_energy_gsl_cube( double q1)
{
  cout<<"!!!!!!!!!minimize_energy_gsl"<<endl;
  q_par=q1;
  /*
  if(n!=2)
  {
    cout <<"This function works only for 2 particles!!!!"<<endl;abort();
  }*/
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n-1);
  
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (n-1);
  /* Initialize method and iterate */
  minex_func.n =n-1;
  minex_func.f = fn_gsl_cube;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, n-1);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    
  iter = 0;
  old_size=1.0e100;
  gsl_vector_set_all (x, 0.0);
  gsl_vector_set_all (ss, 0.01);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  do
  {
    iter++;
    // cout <<iter<<endl;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status)
    {
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
    if(iter%20000==0)
    {
      if(fabs(old_size-calc_en())<1e-16)
      {
        status =GSL_SUCCESS; 
      }
      else
      {
        old_size=calc_en();
      }
    }
  }
  while (status == GSL_CONTINUE && iter < n_gen);
  set_val_gsl_cube(s->x);
  for(int i=0; i<n-1; i++)
  {
    cout<<"i="<<i<<" "<< gsl_vector_get(s->x,i)<<endl;
  }
  
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "GSL minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1gsl: solution :\n" << endl;
  cout<<"momnets"<<endl;
  for(int i=0; i<n; i++)
  {
    cout<< m[i].transpose()<<endl;
  }
    std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<endl;
    }
    cout<<"Bfield"<<endl;
     cout <<B.transpose()<<" "<<get_beta()<<" "<<acos(get_beta())*180/M_PI<<" "<<acos(get_beta1())*180/M_PI<<" "<<acos(get_beta2())*180/M_PI<<endl;

                  //     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//     cout <<B.transpose()<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].dot(B)<<" ";
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//        cout << acos(m[i].dot(B))*180/M_PI<<" ";
//     }
//     cout <<endl;
//     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

void calc::minimize_energy_gsl_superball( double q1)
{
  cout<<"!!!!!!!!!minimize_energy_gsl"<<endl;
  q_par=q1;
  /*
  if(n!=2)
  {
    cout <<"This function works only for 2 particles!!!!"<<endl;abort();
  }*/
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n-1);
  
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (n-1);
  /* Initialize method and iterate */
  minex_func.n =n-1;
  minex_func.f = fn_gsl_superball;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, n-1);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    
  iter = 0;
  old_size=1.0e100;
  gsl_vector_set_all (x, 0.0);
  gsl_vector_set_all (ss, 0.01);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  do
  {
    iter++;
    // cout <<iter<<endl;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status)
    {
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
    if(iter%20000==0)
    {
      if(fabs(old_size-calc_en())<1e-16)
      {
        status =GSL_SUCCESS; 
      }
      else
      {
        old_size=calc_en();
      }
    }
  }
  while (status == GSL_CONTINUE && iter < n_gen);
  set_val_gsl_superball(s->x);
  for(int i=0; i<n-1; i++)
  {
    cout<<"i="<<i<<" "<< gsl_vector_get(s->x,i)<<endl;
  }
  
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "GSL minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1gsl: solution :\n" << endl;
  cout<<"momnets"<<endl;
  for(int i=0; i<n; i++)
  {
    cout<< m[i].transpose()<<endl;
  }
    std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<endl;
    }
    cout<<"Bfield"<<endl;
     cout <<B.transpose()<<" "<<get_beta()<<" "<<acos(get_beta())*180/M_PI<<" "<<acos(get_beta1())*180/M_PI<<" "<<acos(get_beta2())*180/M_PI<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//     cout <<B.transpose()<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].dot(B)<<" ";
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//        cout << acos(m[i].dot(B))*180/M_PI<<" ";
//     }
//     cout <<endl;
//     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

void calc::minimize_energy_gsl_smart(bool first1)
{
  cout<<"!!!!!!!!!minimize_energy_gsl"<<endl;
  first=first1;
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (n);
  
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (n);
  /* Initialize method and iterate */
  minex_func.n =n;
  minex_func.f = fn_gsl_smart;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, n);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    
  iter = 0;
  old_size=1.0e100;
  gsl_vector_set_all (x, 0.0);
  gsl_vector_set_all (ss, 0.01);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  do
  {
    iter++;
    // cout <<iter<<endl;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status)
    {
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
    if(iter%20000==0)
    {
      if(fabs(old_size-calc_en())<1e-16)
      {
        status =GSL_SUCCESS; 
      }
      else
      {
        old_size=calc_en();
      }
    }
  }
  while (status == GSL_CONTINUE && iter < n_gen);
  set_val_gsl_smart(s->x);
  for(int i=0; i<n; i++)
  {
    cout<<"i="<<i<<" "<< gsl_vector_get(s->x,i)<<endl;
  }
  
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "GSL minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1gsl: solution :\n" << endl;
  cout<<"momnets"<<endl;
  for(int i=0; i<n; i++)
  {
    cout<< m[i].transpose()<<endl;
  }
    std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<endl;
    }
    cout<<"Bfield"<<endl;
     cout <<B.transpose()<<" "<<get_beta()<<" "<<acos(get_beta())*180/M_PI<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].transpose()<<endl;
//     }
//     cout <<endl;
//     cout <<B.transpose()<<endl;
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<m[i].dot(B)<<" ";
//     }
//     cout <<endl;
//     for(int i=0; i<n; i++)
//     {
//        cout << acos(m[i].dot(B))*180/M_PI<<" ";
//     }
//     cout <<endl;
//     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}



// void calc::minimize_energy_gsl(bool first1)
// {
//   first=first1;
//   int n_gen=1000000; double tol=1e-18;
//   cout << "start to initialize"<<endl;;
//   const gsl_multimin_fminimizer_type *T = 
//   gsl_multimin_fminimizer_nmsimplex2rand;
//   gsl_multimin_fminimizer *s = NULL;
//   gsl_vector *ss, *x;
//   gsl_multimin_function minex_func;
//   double tmp;
//   size_t iter = 0;
//   int status;
//   double size, old_size=1.0e100;
//   x = gsl_vector_alloc (2+n-1);
//   gsl_vector_set_all (x, 0.0);
//   for(int j=1; j<n; j++)
//   {
//       gsl_vector_set(x,(j+1),0);
//   }
//   
//   
//   /* Set initial step sizes to 1 */
//   ss = gsl_vector_alloc (2+n-1);
//   gsl_vector_set_all (ss, 0.01);
//   /* Initialize method and iterate */
//   minex_func.n =2+n-1;
//   minex_func.f = fn_gsl;
//   minex_func.params = this;
//   s = gsl_multimin_fminimizer_alloc (T, 2+n-1);
//   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
//   add_case=2;
//   
//   
//   
//   
//   
//    
//     std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
// 
//    
//   do
//   {
//     iter++;
//     // cout <<iter<<endl;
//     status = gsl_multimin_fminimizer_iterate(s);
//     
//     if (status)
//     {
//       break;
//     }
// 
//     size = gsl_multimin_fminimizer_size (s);
//     status = gsl_multimin_test_size (size, tol);
// 
//     if (status == GSL_SUCCESS)
//     {
//       printf ("converged to minimum at\n");
//     }
//     if(iter%20000==0)
//     {
//       if(fabs(old_size-calc_en())<1e-16)
//       {
//         status =GSL_SUCCESS; 
//         //cout<<old_size<<" "<<calc_en()<<endl;
//       }
//       else
//       {
//         //cout<<old_size<<" "<<calc_en()<<endl;
//         old_size=calc_en();
//       }
//     }
// //     if(iter%10==0)
// //     {
// //       cout <<iter<<" f() ="<<s->fval<<" "<<size<<endl;
// // //          printf ("%5d f() = %7.10f size = %.10f\n", 
// // //                   iter,
// // //                   s->fval, size);
// //     }
//   }
//   while (status == GSL_CONTINUE && iter < n_gen);
//    
//     cout<<"!!!!!!!!!!!!done"<<endl;
//     
//     //success = optim::de(x,fn,this);
//   cout << "\n0gsl: solution :\n" << endl;
//   for(int i=0; i<2+n-1; i++)
//   {
//     cout<<"i="<<i<<" "<< gsl_vector_get(s->x,i)<<endl;
//   }
//   set_val_gsl(s->x);
//   std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
//   std::chrono::duration<double> elapsed_seconds = end-start;
//   std::cout << "GSL minimization completed successfully.\n"
//                   << "elapsed time: " << elapsed_seconds.count() << "s\n";
//   cout << "\n1gsl: solution :\n" << endl;
//   cout<<"momnets"<<endl;
//   for(int i=0; i<n; i++)
//   {
//     cout<< m[i].transpose()<<endl;
//   }
//     std::cout<<"E="<<calc_en()<<" "<<endl;// << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
//     cout<<"rot_point="<<rot_point<<" add_case="<<add_case<<" cords:"<<endl;
//     for(int i=0; i<n; i++)
//     {
//         cout <<r[i].transpose()<<endl;
//     }
// //     cout <<endl;
// //     for(int i=0; i<n; i++)
// //     {
// //         cout <<m[i].transpose()<<endl;
// //     }
// //     cout <<endl;
// //     cout <<B.transpose()<<endl;
// //     cout <<endl;
// //     for(int i=0; i<n; i++)
// //     {
// //         cout <<m[i].dot(B)<<" ";
// //     }
// //     cout <<endl;
// //     for(int i=0; i<n; i++)
// //     {
// //        cout << acos(m[i].dot(B))*180/M_PI<<" ";
// //     }
// //     cout <<endl;
// //     cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
//   
//   gsl_vector_free(x);
//   gsl_vector_free(ss);
//   gsl_multimin_fminimizer_free (s);
// }



void calc::improve_energy_gsl()
{
  int n_gen=1000000; double tol=1e-18;
  cout << "start to initialize"<<endl;;
  const gsl_multimin_fminimizer_type *T = 
  gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  double tmp;
  size_t iter = 0;
  int status;
  double size, old_size=1.0e100;
  x = gsl_vector_alloc (3*n-1);
  gsl_vector_set_all (x, 0.0);
  for(int i=0; i<3*n-1; i++)
  {
    gsl_vector_set(x,i,vec_x[i]);
  }
  
  
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (3*n-1);
  gsl_vector_set_all (ss, 0.01);
  /* Initialize method and iterate */
  minex_func.n = 3*n-1;
  minex_func.f = fn_gsl;
  minex_func.params = this;
  s = gsl_multimin_fminimizer_alloc (T, 3*n-1);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  
  
  
  
  
   
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

   
  do
  {
    iter++;
    // cout <<iter<<endl;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status)
    {
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
    if(iter%20000==0)
    {
      if(fabs(old_size-calc_en())<1e-16)
      {
        status =GSL_SUCCESS; 
        //cout<<old_size<<" "<<calc_en()<<endl;
      }
      else
      {
        //cout<<old_size<<" "<<calc_en()<<endl;
        old_size=calc_en();
      }
    }
//     if(iter%10==0)
//     {
//       cout <<iter<<" f() ="<<s->fval<<" "<<size<<endl;
// //          printf ("%5d f() = %7.10f size = %.10f\n", 
// //                   iter,
// //                   s->fval, size);
//     }
  }
  while (status == GSL_CONTINUE && iter < n_gen);
   
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
  cout << "\n0gsl: solution :\n" << endl;
  for(int i=0; i<3*n-1; i++)
  {
    cout<< gsl_vector_get(s->x,i)<<endl;
  }
  set_val_gsl(s->x);
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "De minimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
  cout << "\n1de: solution :\n" << endl;
  for(int i=0; i<3*n-1; i++)
  {
    cout<< gsl_vector_get(s->x,i)<<endl;
  }
  cout << "\n1gsl: solution :\n" << endl;
  for(int i=0; i<n; i++)
  {
    cout<< r[i].transpose()<<endl;
  }
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

void calc::minimize_energy_shit3(int kink1)
{
  kink=kink1;
  cout << "start to initialize"<<endl;;
  arma::vec x = arma::ones(2*(n-1),1) - 1.0; // (2,2)
  //cout <<"A"<<endl;
  for(int i=1; i<n; i++)
  {
    x(2*i-2)=0.3892793079172506;
    x(2*i-1)=0.3892793079172506;
    m[i]=m[0];
  }
  B=m[0];
  //cout <<"A"<<endl;
  set_val_shift3(x);
  cout <<"Energy="<<calc_en()/*<<" "<<r[n-2].transpose()*/<<endl;

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    optim::algo_settings_t settings_1;
 
  settings_1.pso_center_particle = false;
//     settings_1.vals_bound = true;
//     arma::vec x_1 = arma::zeros(3*n,1);
//     arma::vec x_2 = arma::zeros(3*n,1);
//     for(int i=0; i<n-1; i++)
//     {
//       x_1(2*i)=-0.5*(n-1);
//       x_2(2*i)=0.5*(n-1);
//       x_1(2*i+1)=-0.5*(n-1);
//       x_2(2*i+1)=0.5*(n-1);
//     }
//     for(int i=2*n-2; i<3*n-1; i++)
//     {
//       x_1(i)=-M_PI;
//       x_2(i)=M_PI;
//     }
//      settings_1.lower_bounds = x_1;
//      settings_1.upper_bounds = x_2;
    //cout<<settings_1.pso_n_pop<<" "<<settings_1.pso_n_gen<<endl; abort();
    //cout<<settings_1.de_n_pop<<" "<<settings_1.de_n_gen<<" "<<settings_1.de_check_freq<<" "<<settings_1.de_mutation_method<<endl; abort();
    settings_1.de_n_pop=150;
    settings_1.de_n_gen=10000;
    //settings_1.de_check_freq=1000;
    settings_1.de_mutation_method=1;
    
//     settings_1.pso_n_pop = 4000;
//     settings_1.pso_n_gen = 4000;
    
    //bool success = optim::pso(x,fn,this,settings_1);
    bool success = optim::de(x,fn_shift3,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0;x(1)=0; x(2)=M_PI; x(3)=0.5*M_PI;x(4)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    //bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val_shift3(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
//     for(int i=2*n-2; i<3*n-1;i++)
//     {  
//       x(i)=fmod(x(i), 1);
//     }
    cout << "\n1de: solution :\n";//
    for(int i=0; i<x.n_elem; i++)
    { 
      cout<< x(i)<<endl;
    } cout<<endl;//
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy_de_smart()
{
  arma::vec x;
  cout << "start to initialize"<<endl;;
  if(abs(B_mag)<1e-6)
  {
    x = arma::ones(2*n-2,1) - 1.0; // (2,2)
  }
  else
  {
    x = arma::ones(2*n,1) - 1.0; // (2,2)
  }
  for(int j=1; j<n; j++)
  {
      x(j-1)=r[j][0];
  }
//   for(int i=n-1; i<2*n-2; i++)
//   {
//       x(i)=0;
//   }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    optim::algo_settings_t settings_1;
 
  settings_1.pso_center_particle = false;
//     settings_1.vals_bound = true;
//     arma::vec x_1 = arma::zeros(3*n,1);
//     arma::vec x_2 = arma::zeros(3*n,1);
//     for(int i=0; i<n-1; i++)
//     {
//       x_1(2*i)=-0.5*(n-1);
//       x_2(2*i)=0.5*(n-1);
//       x_1(2*i+1)=-0.5*(n-1);
//       x_2(2*i+1)=0.5*(n-1);
//     }
//     for(int i=2*n-2; i<3*n-1; i++)
//     {
//       x_1(i)=-M_PI;
//       x_2(i)=M_PI;
//     }
//      settings_1.lower_bounds = x_1;
//      settings_1.upper_bounds = x_2;
    //cout<<settings_1.pso_n_pop<<" "<<settings_1.pso_n_gen<<endl; abort();
    //cout<<settings_1.de_n_pop<<" "<<settings_1.de_n_gen<<" "<<settings_1.de_check_freq<<" "<<settings_1.de_mutation_method<<endl; abort();
    
    settings_1.de_n_pop=150;
    settings_1.de_n_gen=/*640000*/8000;
    //settings_1.de_check_freq=1000;
    settings_1.de_mutation_method=1;
    
    //bool success = optim::pso(x,fn_smart,this,settings_1);
    bool success = optim::de(x,fn_smart,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0;x(1)=0; x(2)=M_PI; x(3)=0.5*M_PI;x(4)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    //bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val_smart(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    if(abs(B_mag)<1e-6)
    {
      for(int i=n-1; i<2*n-2;i++)
      {  
        x(i)=fmod(x(i), 1);
      }
    }
    else
    {
      for(int i=n-1; i<2*n;i++)
      {  
        x(i)=fmod(x(i), 1);
      }
    }
    cout << "\n1de: solution :\n";
    for(int i=0; i<x.n_elem; i++)
    { 
      cout<< x(i)<<endl;
    } cout<<endl;
    
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<" "<<(m[0].dot(B))<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy_de_smart_smart( bool first)
{
  arma::vec x;
  
  cout << "start to initialize"<<endl;;
  int xsize;
  
  if (first)
  {
    if(abs(B_mag)<1e-6)
    {
      if(n%2==0)
      {
        xsize=n/2-1; // (2,2)
      }
      else
      {
        xsize=(n-1)/2;
      }
    }
    else
    {
      if(n%2==0)
      {
        xsize=n/2+1;
      }
      else
      {
        xsize=(n+3)/2;
      }
    }
    x = arma::ones(xsize,1) - 1.0;
    if(n%2==0)
    {
      for(int j=1; j<n/2; j++)
      {
          x(j-1)=r[j][0];
          r[n-j-1][0]=r[j][0];
      }
      r[n-1][0]=r[0][0];r[n-1][1]=r[0][1];
    }
    else
    {
      int n2=(n-1)/2;
      for(int j=1; j<=n2; j++)
      {
          x(j-1)=r[j][0];
      }
      for(int j=1; j<=n2; j++)
      {
        r[n2+j][0]=2*r[n2][0]-r[n2-j][0];
      }
    }
  }
  else
  {
    if(abs(B_mag)<1e-6)
    {
      xsize=1;
    }
    else
    {
      xsize=3;
    }
    x = arma::ones(xsize,1) - 1.0;
  }
  if (first)
  {
    for(int j=1; j<n; j++)
    {
      if(j%2==0)
      {
        m[j]=m[0];
      }
      else
      {
        m[j][0]=-m[0][0];m[j][1]=-m[0][1];m[j][2]=m[0][2];
      }
    }
  }
  else
  {
    for(int j=1; j<n; j++)
    {
      m[j]=m[0];
    }
  }
//   for(int i=n-1; i<2*n-2; i++)
//   {
//       x(i)=0;
//   }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    optim::algo_settings_t settings_1;
 
  settings_1.pso_center_particle = false;
    
    settings_1.de_n_pop=150;
    settings_1.de_n_gen=4000;
    //settings_1.de_check_freq=1000;
    settings_1.de_mutation_method=1;
    
    bool success;
    if(first)
    {
      success = optim::de(x,fn_smart_smart,this,settings_1);
    }
    else
    {
      success = optim::de(x,fn_smart_smart_second,this,settings_1);
    }
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0;x(1)=0; x(2)=M_PI; x(3)=0.5*M_PI;x(4)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    //bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    if(first)
    {
      set_val_smart_smart(x);
    }
    else
    {
      set_val_smart_smart_second(x);
    }
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    if(first)
    {
      for(int i=x.n_elem-2; i<x.n_elem;i++)
      {  
        x(i)=fmod(x(i), 1);
      }
    }
    cout << "\n1de: solution :\n";
    vec_smart_smart_x.resize(x.n_elem);
    for(int i=0; i<x.n_elem; i++)
    { 
      vec_smart_smart_x[i]=x(i);
      cout<< x(i)<<endl;
    } cout<<endl;
    
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<B.transpose()<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<" "<<(m[0].dot(B))<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy_pso()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
    for(int j=1; j<n; j++)
    {
      for(int i=0; i<2; i++)
      {
        x(2*j-2+i)=r[j][i];
      }
    }
    for(int i=2*n-2; i<3*n-1; i++)
    {
        x(i)=0;
    }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    optim::algo_settings_t settings_1;
 
  settings_1.pso_center_particle = false;
//     settings_1.vals_bound = true;
//     arma::vec x_1 = arma::zeros(3*n,1);
//     arma::vec x_2 = arma::zeros(3*n,1);
//     for(int i=0; i<n-1; i++)
//     {
//       x_1(2*i)=-0.5*(n-1);
//       x_2(2*i)=0.5*(n-1);
//       x_1(2*i+1)=-0.5*(n-1);
//       x_2(2*i+1)=0.5*(n-1);
//     }
//     for(int i=2*n-2; i<3*n-1; i++)
//     {
//       x_1(i)=-M_PI;
//       x_2(i)=M_PI;
//     }
//      settings_1.lower_bounds = x_1;
//      settings_1.upper_bounds = x_2;
    //cout<<settings_1.pso_n_pop<<" "<<settings_1.pso_n_gen<<endl; abort();
    //cout<<settings_1.de_n_pop<<" "<<settings_1.de_n_gen<<" "<<settings_1.de_check_freq<<" "<<settings_1.de_mutation_method<<endl; abort();
    
     settings_1.pso_n_pop = 4000;
     settings_1.pso_n_gen = 4000;
    
    bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0;x(1)=0; x(2)=M_PI; x(3)=0.5*M_PI;x(4)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    //bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    for(int i=2*n-2; i<3*n-1;i++)
    {  
      x(i)=fmod(x(i), 2*pi);
    }
    cout << "\n1de: solution :\n" << x << endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy_pso_smart()
{
  arma::vec x;
  cout << "start to initialize"<<endl;;
  if(abs(B_mag)<1e-6)
  {
    x = arma::ones(2*n-2,1) - 1.0; // (2,2)
  }
  else
  {
    x = arma::ones(2*n,1) - 1.0; // (2,2)
  }
  for(int j=1; j<n; j++)
  {
      x(j-1)=r[j][0];
  }
//   for(int i=n-1; i<2*n-2; i++)
//   {
//       x(i)=0;
//   }

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    optim::algo_settings_t settings_1;
 
  settings_1.pso_center_particle = false;
//     settings_1.vals_bound = true;
//     arma::vec x_1 = arma::zeros(3*n,1);
//     arma::vec x_2 = arma::zeros(3*n,1);
//     for(int i=0; i<n-1; i++)
//     {
//       x_1(2*i)=-0.5*(n-1);
//       x_2(2*i)=0.5*(n-1);
//       x_1(2*i+1)=-0.5*(n-1);
//       x_2(2*i+1)=0.5*(n-1);
//     }
//     for(int i=2*n-2; i<3*n-1; i++)
//     {
//       x_1(i)=-M_PI;
//       x_2(i)=M_PI;
//     }
//      settings_1.lower_bounds = x_1;
//      settings_1.upper_bounds = x_2;
    //cout<<settings_1.pso_n_pop<<" "<<settings_1.pso_n_gen<<endl; abort();
    //cout<<settings_1.de_n_pop<<" "<<settings_1.de_n_gen<<" "<<settings_1.de_check_freq<<" "<<settings_1.de_mutation_method<<endl; abort();
    
     settings_1.pso_n_pop = 4000;
     settings_1.pso_n_gen = 4000;
    
    bool success = optim::pso(x,fn_smart,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    //x(0)=0;x(1)=0; x(2)=M_PI; x(3)=0.5*M_PI;x(4)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    //bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val_smart(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    const double pi = arma::datum::pi;
    for(int i=n-1; i<2*n-2;i++)
    {  
      x(i)=fmod(x(i), 1);
    }
    cout << "\n1de: solution :\n" << x << endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy_first()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
 
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
   
    //bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    //x(0)=0.286;x(1)=0.286; x(2)=0; x(3)=0;x(4)=0;
    for(int i=0; i<2*n-2; i++)
    {
      x(i)=0;
    }
    for(int i=2*n-2; i<3*n-3;i++)
    {
      if(i%2==0)
      {
        x(i)=0.5;
      }
      else
      {
        x(i)=0;
      }
    }
    x(3*n-3)=0.25;x(3*n-2)=0.25;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cout << "1de: Ackley test completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cout << "1de: Ackley test completed unsuccessfully." << std::endl;
    }
    cout << "\n1de: solution :\n" << x << endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}

void calc::minimize_energy_second()
{
  cout << "start to initialize"<<endl;;
    arma::vec x = arma::ones(3*n-1,1) + 1.0; // (2,2)
 
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
   
    //bool success = optim::pso(x,fn,this,settings_1);
    //bool success = optim::de(x,fn,this,settings_1);
    x(0)=0.286;x(1)=0.286; //x(2)=0; x(3)=0;x(4)=0;
//     for(int i=0; i<2*n-2; i++)
//     {
//       x(i)=0;
//     }
    for(int i=1; i<n-1; i++)
    {
      x(2*i)=i*x(0);
      x(2*i+1)=i*x(0);
    }
    for(int i=2*n-2; i<3*n-3;i++)
    {
      x(i)=0;
    }
    x(3*n-3)=acos(m[0][0])/(2*M_PI);x(3*n-2)=acos(m[0][1]/sin(2*M_PI*x(3*n-3)))/(2*M_PI);
    //x(3*n-3)=0.5*M_PI;x(3*n-2)=0.5*M_PI;
    //x(0)=0.28689;x(1)=x(0); x(2)=0.0;x(3)=0;x(4)=0;
    bool success = optim::bfgs(x,fn_with_grad,this);
    cout<<"!!!!!!!!!!!!done"<<endl;
    
    //success = optim::de(x,fn,this);
    set_val(x);
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
    cout << "\n1de: solution :\n" << x << endl;
    std::cout <<calc_en()<<" " << acos(m[0].dot(m[1]))*180/M_PI<<" "<<(r10.dot(B)/r10.norm())<<" "<<acos(min((r10.dot(B)/r10.norm()),1.0))*180/M_PI<<endl;
    for(int i=0; i<n; i++)
    {
        cout <<r[i].transpose()<<" ";
    }
    cout <<endl;
    for(int i=0; i<n; i++)
    {
        cout <<m[i].transpose()<<" ";
    }
    cout <<endl;
    cout <<r10.transpose()<<" "<<B.transpose()<<" "<<endl;    
}


double calc::ackley_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
  calc* objfn_data = reinterpret_cast<calc*>(opt_data);
  cout <<"!!!!! a="<<objfn_data->xx<<endl;
    const double x = vals_inp(0);
    const double y = vals_inp(1);
    const double pi = arma::datum::pi;

    double obj_val = -20*std::exp( -0.2*std::sqrt(0.5*(x*x + y*y)) ) - std::exp( 0.5*(std::cos(2*pi*x) + std::cos(2*pi*y)) ) + 22.718282L;

    //

    return obj_val;
}

std::__cxx11::string calc::get_data()
{
  ostringstream os;
  os<<B_mag<<" "<<calc_en()<<" " << (m[0].dot(m[1]))<<" "<<(r10.dot(B))/(r10.norm())<<" ";
  os<<r[0](0)<<" "<<r[0](1)<<" "<<r[0](2)<<" ";
  os<<r[1](0)<<" "<<r[1](1)<<" "<<r[1](2)<<" ";
  os<<m[0](0)<<" "<<m[0](1)<<" "<<m[0](2)<<" ";
  os<<m[1](0)<<" "<<m[1](1)<<" "<<m[1](2)<<" ";
  os<<r10(0)<<" "<<r10(1)<<" "<<r10(2)<<" ";
  os<<B(0)<<" "<<B(1)<<" "<<B(2)/*<<" "*/;
  return os.str();
  
}

double calc::get_offset()
{
  return r[1][0]; 
}

double calc::get_maxoffset()
{
  if(n%2==0)
  {
    return r[n/2][0];
  }
  else
  {
    return r[(n-1)][0];
  }
}

double calc::get_angle()
{
  return m[0].dot(B); 
}
double calc::get_rangle()
{
  Eigen::Vector3d M=Eigen::Vector3d::Zero();
  for(int i=0;i<n; i++)
  {
    M+=m[i];
  }
  return M.dot(r10)/r10.norm()/M.norm(); 
}

double calc::get_Bangle()
{
  Eigen::Vector3d M=Eigen::Vector3d::Zero();
  for(int i=0;i<n; i++)
  {
    M+=m[i];
  }
  return M.dot(B)/M.norm(); 
}


double calc::orintation()
{
  return pow(r10.dot(m[0]),2)/pow(r10.norm(),5);
}

void calc::shift( double a)
{
  r[1][0]=a;
  r[1][1]=a;
  r10=r[n-1]-r[0];
}

void calc::setB(double B_mag1)
{
  B_mag=B_mag1;
}

void calc::fill_vec_x(const arma::vec& vals_inp)
{
  for(int i=0; i<3*n-1; i++)
  {
    vec_x[i]=vals_inp(i);
  }
}


double calc::get_beta()
{
  return B.dot(r10)/r10.norm(); 
}

double calc::get_beta1()
{
  Eigen::Vector3d dr=r[n-kink-1]-r[0];
  return B.dot(dr)/dr.norm(); 
}

double calc::get_beta2()
{
  if(kink==0)
  {
    return 1;
  }
  Eigen::Vector3d dr=r[n-1]-r[n-kink-1];
  return B.dot(dr)/dr.norm(); 
}
