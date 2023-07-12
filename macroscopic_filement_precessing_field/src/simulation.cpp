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

#include "simulation.h"


simulation::simulation(size_t Ns1, double Cm1, double H11, double om1, double torq1):Ns(Ns1), Np(Ns1+1),Nv(3*(Np)), alfa(Np*Np*Np), mu(Np) ,  Cm(Cm1),H1(H11), om(om1), torq(torq1), length(Ns), max_rel_error(0.0), npart(64), x(Nv)
{
  sys=new Integrator(Ns, Cm, H1, om );
  sys1=new Integrator1(Ns, Cm, H1, om, torq );
  x= Eigen::VectorXd::Zero(Nv);
  Omega3= Eigen::VectorXd::Zero(Np);
  init();
}
simulation::simulation(size_t Ns1, double Cm1, double H11, double om1, double torq1, Eigen::VectorXd x_):Ns(Ns1), Np(Ns1+1),Nv(3*(Np)), alfa(Np*Np*Np), mu(Np) ,  Cm(Cm1),H1(H11), om(om1), torq(torq1), length(Ns), max_rel_error(0.0), npart(64), x(x_)
{
  sys=new Integrator(Ns, Cm, H1, om );
  sys1=new Integrator1(Ns, Cm, H1, om, torq );
  Omega3= Eigen::VectorXd::Zero(Np);
  x_temp=x;
  renorm(x_temp,x);
}

simulation::~simulation()
{
  delete sys1;
  delete sys;
}

void simulation::init()
{
  cout<<"Start init"<<endl; 
  Eigen::VectorXd x_coord = Eigen::VectorXd::LinSpaced(Np, -0.005, 0.005);
  Eigen::VectorXd y_coord = Eigen::VectorXd::LinSpaced(Np, -0.5, 0.5);
  Eigen::VectorXd z_coord = Eigen::VectorXd::LinSpaced(Np, 0.0, 0.0);
  for(size_t i=0; i<Np; i++)
  {
    z_coord[i] = 0.01*sin(2*M_PI*y_coord[i]) + 0.001*pow(cosh(y_coord[i]),4) + 0.001*y_coord[i];
  }
//   for(size_t i=0; i<Np; i++)
//   {
//     x_coord[i]=cos(0.1*M_PI*y_coord[i]);
//     y_coord[i]=sin(0.1*M_PI*y_coord[i]);
//   }
  double h=1.0/Ns;
  double dd=1.0-(0.005*M_PI)*(0.005*M_PI);

  //x_coord[0]=0.0;y_coord[0]=0.0;z_coord[0]=0.0;
//   for(size_t i=1; i<Np; i++)
//   {
//     y_coord[i]=0.01*sin(M_PI*h*i);
//     x_coord[i] =x_coord[i-1]+sqrt(h*h-(y_coord[i]-y_coord[i-1])*(y_coord[i]-y_coord[i-1]));
//
//     z_coord[i] = 0.0;
//   }
  
  
  for(size_t i=0; i<Np; i++)
  {
    x[3*i]=x_coord[i];
    x[3*i+1]=y_coord[i];
    x[3*i+2]=z_coord[i];
  }
  //cout<<"Before"<<endl;
  //cout<<x_coord.transpose()<<endl<<endl;
  x_temp=x;
  //cout<<x.transpose()<<endl<<endl<<endl;
  renorm(x_temp,x);
  // cout<<x.transpose()<<endl;
  // abort();
  //cout<<"After"<<endl;
  //cout<<x.transpose()<<endl<<endl;
//    for(size_t i=0; i<Ns; i++)
//    {
//      cout <<i<<" "<<(x(3*(i+1))-x(3*i))*(x(3*(i+1))-x(3*i))+(x(3*(i+1)+1)-x(3*i+1))*(x(3*(i+1)+1)-x(3*i+1))+(x(3*(i+1)+2)-x(3*i+2))*(x(3*(i+1)+2)-x(3*i+2))<<endl;;
//    }
  cout<<"End init"<<endl; 
}

// void simulation::init()
// {
//   cout<<"Start init"<<endl; 
//   Eigen::VectorXd x_coord = Eigen::VectorXd::LinSpaced(Np, 0., 1.0);
//   Eigen::VectorXd y_coord = Eigen::VectorXd::LinSpaced(Np, 0., 0.0);
//   Eigen::VectorXd z_coord = Eigen::VectorXd::LinSpaced(Np, 0., 0.0);
//   for(size_t i=0; i<Np; i++)
//   {
//     x[3*i]=0.01-0.01*sin(M_PI*x_coord[i]);
//     x[3*i+1]=x_coord[i];
//     x[3*i+2]=z_coord[i];
//   }
//   x_temp=x;
//   renorm(x_temp,x);
// //    for(size_t i=0; i<Ns; i++)
// //    {
// //      cout <<i<<" "<<(x(3*(i+1))-x(3*i))*(x(3*(i+1))-x(3*i))+(x(3*(i+1)+1)-x(3*i+1))*(x(3*(i+1)+1)-x(3*i+1))+(x(3*(i+1)+2)-x(3*i+2))*(x(3*(i+1)+2)-x(3*i+2))<<endl;;
// //    }
//   cout<<"End init"<<endl; 
// }



void simulation::propagate(double dt1, double tmax1, double abs_tol,double rel_tol )
{
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_rk_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
//   ostringstream os1;
//   os1<<"data_Np_"<<Np<<"_gamma_"<<gamma<<"_lambda_"<<lambda;
//   folder=os1.str();
//   os1.str("");
//   os1<<"mkdir "<<folder;
//   system(os1.str().c_str());
    dt=dt1; t_max=tmax1;
//   typedef runge_kutta_fehlberg78< Eigen::VectorXd >  error_stepper_type;
  typedef runge_kutta_dopri5<  Eigen::VectorXd> error_stepper_type;
//   //typedef runge_kutta_cash_karp54<  Eigen::VectorXd> error_stepper_type;
  typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
  error_stepper_type rk;
 
  double ax = 1.0;
  double adxdt = 1.0;
  t=0;
  int i=1;
  size_t acpeted_steps=0;
  size_t ignored_steps=0;
  size_t limiter_steps=0;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
 
  controlled_stepper_type cs(
          default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
  controlled_step_result res;
  //int res;
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  while(t<t_max)
  {
    if (dt>1e-15)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      res= cs.try_step(*sys, x, t,dt);
      //res=0;rk.do_step(*sys, x, t,dt);t+=dt;
//        cout<<x<<endl;
//        abort();
      if(res==0)
      {
        acpeted_steps++;
        x_temp=x;
        renorm(x_temp,x);
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//         data[0].push_back(t);
//         for(int ii=0; ii<Np;ii++)
//         {
//           data[ii+1].push_back(x[ii]);
//           data[ii+Np+1].push_back(x[ii+Np]);
//           data[ii+2*Np+1].push_back(x[ii+2*Np]);
//         }
        i++;
      }
      else
      {
        //cout<< "reduced timest tp dt="<<dt<<endl;
        ignored_steps++;
      }

      if(i%100==0)
      {
        accept_step();
        cout << t << '\t' << dt <<" "<<t_max<<" "<<max_rel_error<< '\n';
      }
    }
    else
    {
      dt=2e-10;
      rk.do_step(*sys, x, t,dt);
      cout<<x<<endl;
      abort();
      t+=dt;
      limiter_steps++;
      i++;
      
      cout <<t<<" "<<dt<<endl;
    }
  }
  print_data(1e-4,t_max);
//   cout << t << '\t' << dt << '\n';
//   foo= std::chrono::system_clock::now();
//   diff1 = foo - now;
//   cout <<"Simulation too "<<diff1.count()<<" seconds"<<endl;
//   cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
//   print_pars();
// //   print_data(10.0,t_max);
//   print_data_sigle_file(0.1,t_max);
}



void simulation::propagate_CN(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_CN(x,x_temp,t, dt);
      renorm(x_temp,x);
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%100==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<< '\n';
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }

   print_data(1e-4,t_max);
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_CN_proper(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_CNp_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_CN_propper(x,x_temp,t, dt);
      renorm(x_temp,x);
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%100==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<< '\n';
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }

   print_data(1e-4,t_max);
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_Euler(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_Euler_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_Euler(x,x_temp,t, dt);
      renorm(x_temp,x);
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%100==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<< '\n';
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }

   print_data(1e-4,t_max);
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_Euler_impl_expl(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_EulerI_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_EulerI_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.); 
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_Euler_impl_expl(x,x_temp,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x=x_temp;
      renorm(x_temp,x);
       //calc_dist(x);
       //const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      // cout<<"!! " /*<<*min1*Ns<<" "<<*max1*Ns<<" "*/<<1-*min*Ns<< " "<<*max*Ns-1<<endl;
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%print==0)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_ARS_232(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_ARS232_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_ARS232_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.); 
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs1_IMEX_ARS_232(dt);
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_IMEX_ARS_232(x,x_temp,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      //x=x_temp;
      renorm(x_temp,x);
       //calc_dist(x);
       //const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      // cout<<"!! " /*<<*min1*Ns<<" "<<*max1*Ns<<" "*/<<1-*min*Ns<< " "<<*max*Ns-1<<endl;
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%print==0)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_IMEX_BDF2(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF2_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF2_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
    
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs2_IMEX_BDF2(dt);
  Eigen::VectorXd x_old=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF2(x, dt);
  int j=0;
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF2(x,x_old,x_temp,k1,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x_old=x;
      x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      
      if(t>j*1e-5)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
        j++;
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
   fclose(file);
//   print_data_sigle_file(0.1,t_max);
}


// void simulation::propagate_VSIMEX_BDF3(double dt1, double tmax1, double tol)
// {
//   
//   t=0; dt=dt1, t_max=tmax1;
//   cout <<"Startng propagation! x.size()="<<x.size()<<endl;
//   ostringstream os1;
//   os1<<"data_filament_VSIMEX_BDF3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol;
//   folder=os1.str();
//   os1.str("");
//   os1<<"mkdir "<<folder;
//   system(os1.str().c_str());
//   os1.str("");
//   os1<<"error_filament_VSIMEX_BDF3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol<<".txt";
//   FILE *file;
//   file=NULL;
//   file=fopen(os1.str().c_str(),"w");
//   int print=(int) max(( 1e-5/dt),1.);
//   
//   vec_t.push_back(t);
//   data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//   now= std::chrono::system_clock::now();
// //   data[0].push_back(t);
// //   for(int ii=0; ii<Np;ii++)
// //   {
// //     data[ii+1].push_back(x[ii]);
// //     data[ii+Np+1].push_back(x[ii+Np]);
// //     data[ii+2*Np+1].push_back(x[ii+2*Np]);
// //   }
//   int i=1;
//   sys1->fill_lhs2_VSIMEX_BDF3(dt,dt,dt);
//   Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x;
//   Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
//   Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
//   int j=0;
//   while(t<t_max)
//   {
//     if (dt>1e-14)
//     {
//       if(t+dt>t_max)
//       {
//         break;
//         dt=t_max-t;
//       }
//       sys1->do_VSIMEX_BDF3(x,x_old,x_old_old,x_temp,k1, k2,t, dt);
//       calc_dist(x_temp);
//       const auto [min, max] = std::minmax_element(length.begin(), length.end());
//       x_old_old=x_old;
//       x_old=x;
//       x=x_temp;
//       renorm(x_temp,x);
// //       calc_dist(x);
// //       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
//       t+=dt;
//       vec_t.push_back(t);
//       data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//       
//       if(t>j*1e-5)
//       {
//         double a=1-*min*Ns, b=*max*Ns-1;
//         cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
//         fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
//         j++;
//       }
//       i++;
//     }
//     else
//     {
//       cout <<"Doing something wrong!"<<endl;abort();
//     
//     }
//   }
//   dt=dt1;
// 
//    print_data(1e-6,t_max);
//    anal_data();
//    print_pars();
//    fclose(file);
//    cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
// //   print_data_sigle_file(0.1,t_max);
// }


void simulation::propagate_VSIMEX_BDF3(double dt1, double tmax1, double tol)
{
  double tol_max=tol/640;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_VSIMEX_BDF3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_VSIMEX_BDF3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  Eigen::VectorXd k1, k2;
  Eigen::VectorXd x_old, x_old_old;
  double dt_old, dt_old_old;
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  
  
  x_old_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  dt_old_old=dt;
  
  x_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  dt_old=dt;
  k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt_old);
  k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-dt_old-dt_old_old);
  
  int changed_timestep=0;
  int increase=1000;
  int j_store=0;
  sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
  int j=0;
  Eigen::MatrixXd aa;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_VSIMEX_BDF3(x,x_old,x_old_old,x_temp,k1, k2, t, dt);
      if (accept_step(tol))
      {
        if(t>j*1e-5/*true*/)
        {
          j++;
          cout << t << '\t' << dt <<" "<<dt_old<<" "<<dt_old_old<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
        }
        if(changed_timestep<=0)
        {
          dt_old_old=dt_old;
          dt_old=dt;
          sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
        }
        x_old_old=x_old;
        x_old=x;   
        renorm(x_temp,x);
        t+=dt;
        if(t>j_store*1e-6/*true*/){
          vec_t.push_back(t);
          data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          fill_0mega3();
          data1.push_back(Omega3);
          j_store++;
        }
        changed_timestep++;
        if(max_rel_error<tol_max&&changed_timestep>0)
        {
           if(increase>1000)
           {
            //increase=0;
            dt*=2;
            changed_timestep=-3;
            k1*=2; k2*=2;
            sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
            cout <<"Increasing time step dt="<<dt<<" "<<dt_old<<" "<<max_rel_error<<endl;
           }
           else
           {
             increase++;
           }

            //sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1,t, dt, dt_old);
        }
        else
        {
          //increase=0;
        }

      i++;
      }
      else
      {
        changed_timestep=-3;
//         aa=data[data.size()-2].transpose();
//         x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
//         aa=data[data.size()-3].transpose();
//         x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt/2, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt/2, t-dt-dt_old);
        dt/=2;
        cout<<" reducing timestep dt="<<dt<< " "<<max_rel_error<<endl;
        sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
        increase=0;
        //abort();
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-4,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_VSIMEX_BDF3_DD(double dt1, double tmax1, double tol)
{
  double tol_max=tol/640;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_VSIMEX_BDF3_DD_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_VSIMEX_BDF3_DD_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  Eigen::VectorXd k1, k2;
  Eigen::VectorXd x_old, x_old_old, x_old_old_old;
  double dt_old, dt_old_old;
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  
  x_old_old_old=x;
  x_old_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl_DD(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  //vec_t.push_back(t);
  //data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  dt_old_old=dt;
  
  x_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl_DD(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  //vec_t.push_back(t);
  //data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  dt_old=dt;
  k1=sys1->fill_K1_IMEX_BDF3_DD(x_old, dt, t-dt_old);
  k2=sys1->fill_K1_IMEX_BDF3_DD(x_old_old, dt, t-dt_old-dt_old_old);
  
  int changed_timestep=0;
  int increase=1000;
  sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
  int j=0;
  int j_store=0;
  Eigen::MatrixXd aa;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_VSIMEX_BDF3_DD(x,x_old,x_old_old,x_temp,k1, k2, t, dt);
      if (accept_step(tol))
      {
        if(t>j*1e-4/*true*/)
        {
          j++;
          cout << t << '\t' << dt <<" "<<dt_old<<" "<<dt_old_old<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
        }
        if(changed_timestep<=0)
        {
          dt_old_old=dt_old;
          dt_old=dt;
          sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
        }
        x_old_old_old=x_old_old;
        x_old_old=x_old;
        x_old=x;   
        renorm(x_temp,x);
        t+=dt;
        if(t>j_store*1e-6/*true*/){
          vec_t.push_back(t);
          data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          fill_0mega3();
          data1.push_back(Omega3);
          j_store++;
        }
        changed_timestep++;
        if(max_rel_error<tol_max&&changed_timestep>0)
        {
           if(increase>1000)
           {
            //increase=0;
            dt*=2;
            changed_timestep=-3;
            k1*=2; k2*=2;
            sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
            cout <<"Increasing time step dt="<<dt<<" "<<dt_old<<" "<<max_rel_error<<endl;
           }
           else
           {
             increase++;
           }
            
            //sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1,t, dt, dt_old);
        }
        else
        {
          //increase=0;
        }
      
      i++;
      }
      else
      {
        changed_timestep=-3;
        //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
        //cout<<data[data.size()-5]-aaa<<endl;abort();
        k1=sys1->fill_K1_IMEX_BDF3_DD(x_old, dt/2, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3_DD(x_old_old, dt/2, t-dt-dt_old);
        dt/=2;
        cout<<" reducing timestep dt="<<dt<< " "<<max_rel_error<<endl;
        sys1->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
        increase=0;
        //abort();
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-5,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_VSIMEX_BDF4(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_VSIMEX_BDF4_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_VSIMEX_BDF4_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  Eigen::VectorXd k1, k2, k3;
  Eigen::VectorXd x_old, x_old_old, x_old_old2;
  double dt_old, dt_old_old, dt_old_old2;
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  x_old_old2=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  dt_old_old2=dt;
  
  x_old_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  dt_old_old=dt;
  
  x_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  dt_old=dt;
  k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt_old);
  k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-dt_old-dt_old_old);
  k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-dt_old-dt_old_old-dt_old_old2);
  
  int changed_timestep=0;
  sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
  int j=0;
  Eigen::MatrixXd aa;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_VSIMEX_BDF4(x,x_old,x_old_old,x_old_old2,x_temp,k1, k2, k3, t, dt);
      if (accept_step(1e-10))
      {
        if(t>j*1e-5/*true*/)
        {
          j++;
          cout << t << '\t' << dt <<" "<<dt_old<<" "<<dt_old_old<<" "<<dt_old_old2<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
        }
        if(changed_timestep<=0)
        {
          dt_old_old2=dt_old_old;
          dt_old_old=dt_old;
          dt_old=dt;
          sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
        }
        x_old_old2=x_old_old;
        x_old_old=x_old;
        x_old=x;   
        renorm(x_temp,x);
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        fill_0mega3();
        data1.push_back(Omega3);
        changed_timestep++;
        if(max_rel_error<5e-12&&changed_timestep>0)
        {
          dt*=2;
          changed_timestep=-3;
          k1*=2; k2*=2; k3*=2;
          sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
          cout <<"Increasing time step dt="<<dt<<" "<<dt_old<<endl;
          
          //sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1,t, dt, dt_old);
        }
      
      i++;
      }
      else
      {
        changed_timestep=-3;
        aa=data[data.size()-2].transpose();
        x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
        aa=data[data.size()-3].transpose();
        x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
        aa=data[data.size()-4].transpose();
        x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt/2, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt/2, t-dt-dt_old);
        k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt/2,t-dt-dt_old-dt_old_old);
        dt/=2;
        cout<<" reducing timestep dt="<<dt<< " "<<max_rel_error<<endl;
        sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
        //abort();
      }
    }
//       calc_dist(x_temp);
//       const auto [min, max] = std::minmax_element(length.begin(), length.end());
//       x_old_old2=x_old_old;
//       x_old_old=x_old;
//       x_old=x;
//       x=x_temp;
//       renorm(x_temp,x);
// //       calc_dist(x);
// //       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
//       t+=dt;
//       vec_t.push_back(t);
//       data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//       
//       if(i%print==0)
//       {
//         double a=1-*min*Ns, b=*max*Ns-1;
//         cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
//         fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
//       }
//       i++;
//     }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-4,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_VSIMEX_BDF4_adaptive(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_VSIMEX_BDF4a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_VSIMEX_BDF4a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  Eigen::VectorXd x_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<1e-10)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  int changed_timestep=0;
  double dt_old=dt, dt_old_old=dt, dt_old_old2=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
  Eigen::VectorXd x_old_old=x_old; Eigen::VectorXd x_old_old2=x_old_old; 
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x_old, dt);
  Eigen::VectorXd k2=k1;
  Eigen::VectorXd k3=k2;
  int j=0;
  Eigen::MatrixXd aa;
  int rejected=0;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_VSIMEX_BDF4(x,x_old,x_old_old,x_old_old2,x_temp,k1, k2, k3, t, dt);
      if (accept_step(1e-10))
      {
        if(t>j*1e-5/*true*/)
        {
          j++;
          cout << t << '\t' << dt <<" "<<dt_old<<" "<<dt_old_old<<" "<<dt_old_old2<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
        }
        if(max_rel_error<5e-11)
        {
          x_old_old2=x_old_old;
          x_old_old=x_old;
          x_old=x;   
          renorm(x_temp,x);
          t+=dt;
          vec_t.push_back(t);
          data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          dt_old_old2=dt_old_old;
          dt_old_old=dt_old;
          dt_old=dt;
          dt*=min( 0.9*pow(max_rel_error/1e-10 , -0.25 ) , 5. );
          k1*=dt/dt_old;k2*=dt/dt_old;k3*=dt/dt_old;
          sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
          //cout <<"Increasing time step dt="<<dt<<" "<<dt_old<<endl;
        }
        else if(changed_timestep<=0)
        {
          dt_old_old2=dt_old_old;
          dt_old_old=dt_old;
          dt_old=dt;
          sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
          x_old_old2=x_old_old;
          x_old_old=x_old;
          x_old=x;   
          renorm(x_temp,x);
          t+=dt;
          vec_t.push_back(t);
          data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          changed_timestep++;
        }
        else
        {
          x_old_old2=x_old_old;
          x_old_old=x_old;
          x_old=x;   
          renorm(x_temp,x);
          t+=dt;
          vec_t.push_back(t);
          data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          changed_timestep++;
        }
      i++;
      }
      else
      {
        rejected++;
        changed_timestep=-3;
        aa=data[data.size()-2].transpose();
        x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
        aa=data[data.size()-3].transpose();
        x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
        aa=data[data.size()-4].transpose();
        x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
        double dt_new= dt*max( 0.9*pow(max_rel_error/1e-10 , -0.25 ) , 0.2 );
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt_new, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt_new, t-dt-dt_old);
        k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt_new,t-dt-dt_old-dt_old_old);
        //cout<<" reducing timestep dt="<<dt_new<< " "<<dt<< " "<<max_rel_error<<endl;
        dt=dt_new;
        sys1->fill_lhs1_VSIMEX_BDF4(dt, dt_old, dt_old_old, dt_old_old2);
        //abort();
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;
  cout<<" finished: rejected="<<rejected<<" accepted="<<i<<" ratio="<<(rejected/((double) i))<<endl;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}

// void simulation::propagate_VSIMEX_BDF2(double dt1, double tmax1)
// {
//   
//   t=0; dt=dt1, t_max=tmax1;
//   double dt_old=dt;
//   cout <<"Startng propagation! x.size()="<<x.size()<<endl;
//   ostringstream os1;
//   os1<<"data_filament_VSIMEX_BDF2_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
//   folder=os1.str();
//   os1.str("");
//   os1<<"mkdir "<<folder;
//   system(os1.str().c_str());
//   os1.str("");
//   os1<<"error_filament_VSIMEX_BDF2_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
//   FILE *file;
//   file=NULL;
//   file=fopen(os1.str().c_str(),"w");
//   int print=(int) max(( 1e-5/dt),1.);
//     
//   vec_t.push_back(t);
//   data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//   now= std::chrono::system_clock::now();
// //   data[0].push_back(t);
// //   for(int ii=0; ii<Np;ii++)
// //   {
// //     data[ii+1].push_back(x[ii]);
// //     data[ii+Np+1].push_back(x[ii+Np]);
// //     data[ii+2*Np+1].push_back(x[ii+2*Np]);
// //   }
//   int i=1;
//   sys1->fill_lhs_Euler_impl_expl(dt);
//   Eigen::VectorXd x_old=x;
//   while(true)
//   {
//     sys1->do_Euler_impl_expl(x,x_temp,t, dt);
//     accept_step();
//     cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
//     if(max_rel_error<1e-10)
//     {
//       renorm(x_temp,x);
//       break;
//     }
//     dt/=2;
//     sys1->fill_lhs_Euler_impl_expl(dt);
//   }
//   Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF2(x_old, dt);
//   t+=dt;
//   vec_t.push_back(t);
//   data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//   sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
//   //abort();
//   int j=0;
//   while(t<t_max)
//   {
//     if (dt>1e-14)
//     {
//       if(t+dt>t_max)
//       {
//         break;
//         dt=t_max-t;
//       }
//       sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1,t, dt, dt_old);
//       if (accept_step(1e-6))
//       {
//         if(fabs(1-dt_old/dt)>1e-10)
//         {
//           dt_old=dt;
//           sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
//         }
//         x_old=x;
//         renorm(x_temp,x);
//       
// //       calc_dist(x);
// //       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
//         t+=dt;
//         vec_t.push_back(t);
//         data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
//         if(t>j*1e-5)
//         {
//           j++;
//           cout << t << '\t' << dt <<" "<<dt_old<<" "<<max_rel_error<<endl;
//           fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
//         }
//         if(max_rel_error<1e-8)
//         {
//           dt*=2;
//           sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
//           cout <<"Increasing time step dt="<<dt<<" "<<dt_old<<endl;
//           k1*=2;
//           //sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1,t, dt, dt_old);
//         }
//       
//       i++;
//       }
//       else
//       {
//         dt/=2;
//         cout<<" reducing timestep dt="<<dt<< " "<<max_rel_error<<endl;
//         sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
//         abort();
//       }
//     }
//     else
//     {
//       cout <<"Doing something wrong!"<<endl;abort();
//     
//     }
//   }
//   dt=dt1;
// 
//    print_data(1e-6,t_max);
//    anal_data();
//    print_pars();
//    cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//    fclose(file);
// //   print_data_sigle_file(0.1,t_max);
// }



void simulation::propagate_VSIMEX_BDF2(double dt1, double tmax1, double tol)
{
  double max_tol=tol/6;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_VSIMEX_BDF2_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_VSIMEX_BDF2_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_tol_"<<tol<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  Eigen::VectorXd k1;
  Eigen::VectorXd x_old;
  double dt_old;
  int i=1;
  
  sys1->fill_lhs_Euler_impl_expl(dt);
  x_old=x;
  while(true)
  {
    sys1->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<tol)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys1->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  dt_old=dt;
  k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt_old);
  
  int changed_timestep=0;
  sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
  int j=0;
  Eigen::MatrixXd aa;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1, t, dt);
      if (accept_step(tol))
      {
        if(t>j*1e-5/*true*/)
        {
          j++;
          cout << t << '\t' << dt <<" "<<dt_old<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
        }
        if(changed_timestep<=0)
        {
          dt_old=dt;
          sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
        }
        x_old=x;   
        renorm(x_temp,x);
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        changed_timestep++;
        if(max_rel_error<max_tol&&changed_timestep>100)
        {
          dt*=2;
          changed_timestep=-1;
          k1*=2;
          sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
          cout <<"Increasing time step dt="<<dt<<" "<<dt_old<<endl;
          
          //sys1->do_VSIMEX_BDF2(x,x_old,x_temp,k1,t, dt, dt_old);
        }
      
      i++;
      }
      else
      {
        changed_timestep=-1;
        aa=data[data.size()-2].transpose();
        x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt/2, t-dt);
        dt/=2;
        cout<<" reducing timestep dt="<<dt<< " "<<max_rel_error<<endl;
        sys1->fill_lhs1_IMEX_BDF2(dt, dt_old);
        //abort();
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}



void simulation::propagate_IMEX_BDF3(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  fill_0mega3();
  data1.push_back(Omega3);
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs2_IMEX_BDF3(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  int j=0;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF3(x,x_old,x_old_old,x_temp,k1, k2,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x_old_old=x_old;
      x_old=x;
      x=x_temp;
//      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      fill_0mega3();
      data1.push_back(Omega3);
      if(t>j*1e-5)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
        j++;
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-5,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_IMEX_BDF5(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF5_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF5_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs2_IMEX_BDF5(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x; Eigen::VectorXd x_old_old2=x; Eigen::VectorXd x_old_old3=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k3=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k4=sys1->fill_K1_IMEX_BDF3(x, dt);
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF5(x,x_old,x_old_old,x_old_old2, x_old_old3,x_temp,k1, k2, k3, k4,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x_old_old3=x_old_old2;
      x_old_old2=x_old_old;
      x_old_old=x_old;
      x_old=x;
      x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%print==0)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}



void simulation::propagate_IMEX_TVB5(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_TVB5_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_TVB5_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=/*(int) max(( 1e-5/dt),1.)*/1;
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs2_IMEX_TVB5(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x; Eigen::VectorXd x_old_old2=x; Eigen::VectorXd x_old_old3=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k3=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k4=sys1->fill_K1_IMEX_BDF3(x, dt);
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_TVB5(x,x_old,x_old_old,x_old_old2, x_old_old3,x_temp,k1, k2, k3, k4,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x_old_old3=x_old_old2;
      x_old_old2=x_old_old;
      x_old_old=x_old;
      x_old=x;
      x=x_temp;
      if(t<1e-4)
      {
        renorm(x_temp,x);
      }
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%print==0)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_IMEX_BDF4(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF4_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF4_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs2_IMEX_BDF4(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x; Eigen::VectorXd x_old_old2=x; 
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k3=sys1->fill_K1_IMEX_BDF3(x, dt);
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF4(x,x_old,x_old_old,x_old_old2,x_temp,k1, k2, k3, t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x_old_old2=x_old_old;
      x_old_old=x_old;
      x_old=x;
      x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%print==0)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}




void simulation::propagate_IMEX_TVB3(double dt1, double tmax1)
{
  
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_TVB3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_TVB3_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs2_IMEX_TVB3(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  int j=0;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_TVB3(x,x_old,x_old_old,x_temp,k1, k2,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x_old_old=x_old;
      x_old=x;
      x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(t>j*1e-5)//if(i%print==0)
      {
        double a=1-*min*Ns, b=*max*Ns-1;
        cout << t << '\t' << dt <<" "<<t_max<<" "<<a<< " "<<b<<endl;
        fprintf(file,"%.7e %.14e %.14e\n",t, dt, (a<b)?b:a );
        j++;
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;
  fclose(file);

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_IMEX_TVB5_adaptive(double dt1, double tmax1)
{
  int j=0;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_TVB5a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_TVB5a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int failed_to_reduce_dt=0;
  int sucsessfull_in_row=0;
  int i=1;
  sys1->fill_lhs2_IMEX_TVB5(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x; Eigen::VectorXd x_old_old2=x; Eigen::VectorXd x_old_old3=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k3=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k4=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::MatrixXd aa;
  bool res=true;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_TVB5(x,x_old,x_old_old,x_old_old2, x_old_old3,x_temp,k1, k2, k3, k4,t, dt);
      if(i>0)
      {
        res=accept_step();        
      }
      else
      {
        res=accept_step_beginning();
      }
      if(res)
      {
        failed_to_reduce_dt=0;
        x_old_old3=x_old_old2;
        x_old_old2=x_old_old;
        x_old_old=x_old;
        x_old=x;
        renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
        if(t>j*1e-5)
        {
          cout << t << '\t' << dt <<" "<<t_max<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
          j++;
        }
        i++;
        if(max_rel_error<5.e-11)
        {
          sucsessfull_in_row++;
        }
        else
        {
          sucsessfull_in_row=0;
        }
        if(sucsessfull_in_row>100)
        {
        //  cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt)<<endl;abort();
          sucsessfull_in_row=0;
          aa=data[data.size()-3].transpose();
          x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
          aa=data[data.size()-5].transpose();
          x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
          aa=data[data.size()-7].transpose();
          x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          aa=data[data.size()-9].transpose();
          x_old_old3=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
         // cout<<"!!!!!!! Increasing dt"<<endl;
          //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
          dt*=2;
          sys1->fill_lhs2_IMEX_TVB5(dt);
          //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
          k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
          k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
          k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
          k4=sys1->fill_K1_IMEX_BDF3(x_old_old3, dt,t-4*dt);
        }
      }
      else
      {
        sucsessfull_in_row=0;
        //cout <<"i="<<i<<" failed_to_reduce_dt="<<failed_to_reduce_dt<<" max_rel_error="<<max_rel_error<<endl;
        x_old_old3=x;
        npart/=2;
        dt/=npart;
        while(true)
        {
          sys1->fill_lhs_Euler_impl_expl(dt);
          sys1->do_Euler_impl_expl(x, x_temp, t, dt);
          if(accept_step())
          {
            renorm(x_temp, x);
            break;
          }
          npart*=2;
          dt/=2;
        //  cout <<"Reducing dt="<<dt<<" "<<max_rel_error<<endl;
        }
        sys1->fill_lhs2_IMEX_TVB5(dt);
        t+=dt;
        x_old_old2=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        x_old_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        x_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
        k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
        k4=sys1->fill_K1_IMEX_BDF3(x_old_old3, dt,t-4*dt);
        while(true)
        {
          for (int i=0; i<8; i++)
          {
            sys1->do_IMEX_TVB5(x,x_old,x_old_old,x_old_old2, x_old_old3,x_temp,k1, k2, k3, k4,t, dt);
            x_old_old3=x_old_old2;
            x_old_old2=x_old_old;
            x_old_old=x_old;
            x_old=x;
            accept_step();
            renorm(x_temp,x);
            t+=dt;
            vec_t.push_back(t);
            data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          }
          if(max_rel_error<1e-11)
          {
            aa=data[data.size()-3].transpose();
            x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
            aa=data[data.size()-5].transpose();
            x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
            //cout<<data[data.size()-5]-aaa<<endl;abort();
            aa=data[data.size()-7].transpose();
            x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            aa=data[data.size()-9].transpose();
            x_old_old3=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
            dt*=2;
            // cout<<"!!!!!!! Increasing dt="<<dt<<endl;
            sys1->fill_lhs2_IMEX_TVB5(dt);
            //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
            k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
            k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
            k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
            k4=sys1->fill_K1_IMEX_BDF3(x_old_old3, dt,t-4*dt);
          }
          else
          {
            break;
          }
        }
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_IMEX_BDF5_adaptive(double dt1, double tmax1)
{
  int j=0;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF5a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF5a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int failed_to_reduce_dt=0;
  int sucsessfull_in_row=0;
  int i=1;
  sys1->fill_lhs2_IMEX_BDF5(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x; Eigen::VectorXd x_old_old2=x; Eigen::VectorXd x_old_old3=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k3=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k4=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::MatrixXd aa;
  bool res=true;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF5(x,x_old,x_old_old,x_old_old2, x_old_old3,x_temp,k1, k2, k3, k4,t, dt);
      if(i>0)
      {
        res=accept_step();        
      }
      else
      {
        res=accept_step_beginning();
      }
      if(res)
      {
        failed_to_reduce_dt=0;
        x_old_old3=x_old_old2;
        x_old_old2=x_old_old;
        x_old_old=x_old;
        x_old=x;
        renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
        if(t>j*1e-5)
        {
          cout << t << '\t' << dt <<" "<<t_max<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
          j++;
        }
        i++;
        if(max_rel_error<5.e-11)
        {
          sucsessfull_in_row++;
        }
        else
        {
          sucsessfull_in_row=0;
        }
        if(sucsessfull_in_row>100)
        {
        //  cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt)<<endl;abort();
          sucsessfull_in_row=0;
          aa=data[data.size()-3].transpose();
          x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
          aa=data[data.size()-5].transpose();
          x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
          aa=data[data.size()-7].transpose();
          x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          aa=data[data.size()-9].transpose();
          x_old_old3=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
         // cout<<"!!!!!!! Increasing dt"<<endl;
          //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
          dt*=2;
          sys1->fill_lhs2_IMEX_BDF5(dt);
          //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
          k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
          k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
          k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
          k4=sys1->fill_K1_IMEX_BDF3(x_old_old3, dt,t-4*dt);
        }
      }
      else
      {
        sucsessfull_in_row=0;
        //cout <<"i="<<i<<" failed_to_reduce_dt="<<failed_to_reduce_dt<<" max_rel_error="<<max_rel_error<<endl;
        x_old_old3=x;
        npart/=2;
        dt/=npart;
        while(true)
        {
          sys1->fill_lhs_Euler_impl_expl(dt);
          sys1->do_Euler_impl_expl(x, x_temp, t, dt);
          if(accept_step())
          {
            renorm(x_temp, x);
            break;
          }
          npart*=2;
          dt/=2;
        //  cout <<"Reducing dt="<<dt<<" "<<max_rel_error<<endl;
        }
        sys1->fill_lhs2_IMEX_BDF5(dt);
        t+=dt;
        x_old_old2=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        x_old_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        x_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
        k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
        k4=sys1->fill_K1_IMEX_BDF3(x_old_old3, dt,t-4*dt);
        while(true)
        {
          for (int i=0; i<8; i++)
          {
            sys1->do_IMEX_BDF5(x,x_old,x_old_old,x_old_old2, x_old_old3,x_temp,k1, k2, k3, k4,t, dt);
            x_old_old3=x_old_old2;
            x_old_old2=x_old_old;
            x_old_old=x_old;
            x_old=x;
            accept_step();
            renorm(x_temp,x);
            t+=dt;
            vec_t.push_back(t);
            data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          }
          if(max_rel_error<1e-11)
          {
            aa=data[data.size()-3].transpose();
            x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
            aa=data[data.size()-5].transpose();
            x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
            //cout<<data[data.size()-5]-aaa<<endl;abort();
            aa=data[data.size()-7].transpose();
            x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            aa=data[data.size()-9].transpose();
            x_old_old3=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
            dt*=2;
            // cout<<"!!!!!!! Increasing dt="<<dt<<endl;
            sys1->fill_lhs2_IMEX_BDF5(dt);
            //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
            k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
            k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
            k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
            k4=sys1->fill_K1_IMEX_BDF3(x_old_old3, dt,t-4*dt);
          }
          else
          {
            break;
          }
        }
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_IMEX_BDF4_adaptive(double dt1, double tmax1)
{
  int j=0;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF4a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF4a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int failed_to_reduce_dt=0;
  int sucsessfull_in_row=0;
  int i=1;
  sys1->fill_lhs2_IMEX_BDF4(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x; Eigen::VectorXd x_old_old2=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k3=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::MatrixXd aa;
  bool res=true;
  while(t<t_max)
  {
    if (dt>1e-14)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF4(x,x_old,x_old_old,x_old_old2,x_temp,k1, k2, k3,t, dt);
      if(i>0)
      {
        res=accept_step();        
      }
      else
      {
        res=accept_step_beginning();
      }
      if(res)
      {
        failed_to_reduce_dt=0;
        x_old_old2=x_old_old;
        x_old_old=x_old;
        x_old=x;
        renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
        if(t>j*1e-5)
        {
          cout << t << '\t' << dt <<" "<<t_max<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
          j++;
        }
        i++;
        if(max_rel_error<5.e-11)
        {
          sucsessfull_in_row++;
        }
        else
        {
          sucsessfull_in_row=0;
        }
        if(sucsessfull_in_row>100)
        {
        //  cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt)<<endl;abort();
          sucsessfull_in_row=0;
          aa=data[data.size()-3].transpose();
          x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
          aa=data[data.size()-5].transpose();
          x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
          aa=data[data.size()-7].transpose();
          x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
         // cout<<"!!!!!!! Increasing dt"<<endl;
          //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
          dt*=2;
          sys1->fill_lhs2_IMEX_BDF4(dt);
          //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
          k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
          k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
          k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
        }
      }
      else
      {
        sucsessfull_in_row=0;
        //cout <<"i="<<i<<" failed_to_reduce_dt="<<failed_to_reduce_dt<<" max_rel_error="<<max_rel_error<<endl;
        x_old_old2=x;
        npart/=2;
        dt/=npart;
        while(true)
        {
          sys1->fill_lhs_Euler_impl_expl(dt);
          sys1->do_Euler_impl_expl(x, x_temp, t, dt);
          if(accept_step())
          {
            renorm(x_temp, x);
            break;
          }
          npart*=2;
          dt/=2;
        //  cout <<"Reducing dt="<<dt<<" "<<max_rel_error<<endl;
        }
        sys1->fill_lhs2_IMEX_BDF4(dt);
        t+=dt;
        x_old_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        x_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
        k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
        while(true)
        {
          for (int i=0; i<6; i++)
          {
            sys1->do_IMEX_BDF4(x,x_old,x_old_old,x_old_old2,x_temp,k1, k2, k3,t, dt);
            x_old_old2=x_old_old;
            x_old_old=x_old;
            x_old=x;
            accept_step();
            renorm(x_temp,x);
            t+=dt;
            vec_t.push_back(t);
            data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          }
          if(max_rel_error<1e-11)
          {
            aa=data[data.size()-3].transpose();
            x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
            aa=data[data.size()-5].transpose();
            x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
            //cout<<data[data.size()-5]-aaa<<endl;abort();
            aa=data[data.size()-7].transpose();
            x_old_old2=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
            dt*=2;
            // cout<<"!!!!!!! Increasing dt="<<dt<<endl;
            sys1->fill_lhs2_IMEX_BDF4(dt);
            //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
            k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
            k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
            k3=sys1->fill_K1_IMEX_BDF3(x_old_old2, dt,t-3*dt);
          }
          else
          {
            break;
          }
        }
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_IMEX_BDF3_adaptive(double dt1, double tmax1)
{
  int j=0;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_IMEX_BDF3a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  os1<<"error_filament_IMEX_BDF3a_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt<<".txt";
  FILE *file;
  file=NULL;
  file=fopen(os1.str().c_str(),"w");
  int print=(int) max(( 1e-5/dt),1.);
  
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int failed_to_reduce_dt=0;
  int sucsessfull_in_row=0;
  int i=1;
  sys1->fill_lhs2_IMEX_BDF3(dt);
  Eigen::VectorXd x_old=x; Eigen::VectorXd x_old_old=x;
  Eigen::VectorXd k1=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::VectorXd k2=sys1->fill_K1_IMEX_BDF3(x, dt);
  Eigen::MatrixXd aa;
  bool res=true;
  while(t<t_max)
  {
    if (dt>1e-24)
    {
      if(t+dt>t_max)
      {
        break;
        dt=t_max-t;
      }
      sys1->do_IMEX_BDF3(x,x_old,x_old_old,x_temp,k1, k2,t, dt);
      if(i>0)
      {
        res=accept_step();        
      }
      else
      {
        res=accept_step_beginning();
      }
      if(res)
      {
        failed_to_reduce_dt=0;
        x_old_old=x_old;
        x_old=x;
        renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
        if(t>j*1e-5)
        {
          cout << t << '\t' << dt <<" "<<t_max<<" "<<max_rel_error<<endl;
          fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
          j++;
        }
        i++;
        if(max_rel_error<5.e-11)
        {
          sucsessfull_in_row++;
        }
        else
        {
          sucsessfull_in_row=0;
        }
        if(sucsessfull_in_row>100)
        {
        //  cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt)<<endl;abort();
          sucsessfull_in_row=0;
          aa=data[data.size()-3].transpose();
          x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
          aa=data[data.size()-5].transpose();
          x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
          //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
          //cout<<data[data.size()-5]-aaa<<endl;abort();
         // cout<<"!!!!!!! Increasing dt"<<endl;
          //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
          dt*=2;
          sys1->fill_lhs2_IMEX_BDF3(dt);
          //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
          k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
          k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
        }
      }
      else
      {
        sucsessfull_in_row=0;
        cout <<"i="<<i<<" failed_to_reduce_dt="<<failed_to_reduce_dt<<" max_rel_error="<<max_rel_error<<endl;
        x_old_old=x;
        npart/=2;
        dt/=npart;
        while(true)
        {
          sys1->fill_lhs_Euler_impl_expl(dt);
          sys1->do_Euler_impl_expl(x, x_temp, t, dt);
          if(accept_step())
          {
            renorm(x_temp, x);
            break;
          }
          npart*=2;
          dt/=2;
          cout <<"Reducing dt="<<dt<<" "<<max_rel_error<<endl;
        }
        sys1->fill_lhs2_IMEX_BDF4(dt);
        t+=dt;
        x_old=x;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        sys1->do_Euler_impl_expl(x, x_temp, t, dt);renorm(x_temp, x);
        t+=dt;
        vec_t.push_back(t);
        data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
        k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
        k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
        while(true)
        {
          for (int i=0; i<4; i++)
          {
            sys1->do_IMEX_BDF3(x,x_old,x_old_old,x_temp,k1, k2,t, dt);
            x_old_old=x_old;
            x_old=x;
            accept_step();
            renorm(x_temp,x);
            t+=dt;
            vec_t.push_back(t);
            data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
          }
          if(max_rel_error<1e-11)
          {
            aa=data[data.size()-3].transpose();
            x_old=(Eigen::Map<Eigen::VectorXd>((aa).data(), Nv));
            aa=data[data.size()-5].transpose();
            x_old_old=(Eigen::Map<Eigen::VectorXd>(aa.data(), Nv));
            //Eigen::MatrixXd aaa=((Eigen::Map<Eigen::MatrixXd> (x_old_old3.data(), 3,Np)).transpose());
            //cout<<data[data.size()-5]-aaa<<endl;abort();
            //cout<<k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-2*dt)<<endl;abort();
            dt*=2;
             cout<<"!!!!!!! Increasing dt="<<dt<<endl;
            sys1->fill_lhs2_IMEX_BDF3(dt);
            //cout<<2*k2-sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt)<<endl<<endl;;abort();
            k1=sys1->fill_K1_IMEX_BDF3(x_old, dt, t-dt);
            k2=sys1->fill_K1_IMEX_BDF3(x_old_old, dt, t-2*dt);
          }
          else
          {
            break;
          }
        }
      }
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}


void simulation::propagate_rk4_impl_expl(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_rk4I_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  while(t<t_max)
  {
    if (dt>1e-10)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_rk4_impl_expl(x,x_temp,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
      x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      cout <</**min1*Ns<<" "<<*max1*Ns<<" "<<*/*min*Ns<< " "<<*max*Ns<<endl;
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%100==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<< '\n';
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_Euler1_impl_expl(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_EulerI1_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  while(t<t_max)
  {
    if (dt>1e-15)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      if(i>1000)
      {
        sys1->do_Euler1_impl_expl(x,x_temp,t, dt);
      }
      else
      {
        sys1->do_Euler_impl_expl(x,x_temp,t, dt);
      }
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
//      x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%100==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<<" "<<1-*min*Ns<< " "<<*max*Ns-1<<endl;
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
//   print_data_sigle_file(0.1,t_max);
}

void simulation::propagate_midp_impl_expl(double dt1, double tmax1)
{
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_midp_Ns_"<<Ns<<"_omega_"<<om<<"_Cm_"<<Cm<<"_H1_"<<H1<<"_Torq_"<<torq<<"_dt_"<<dt;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  vec_t.push_back(t);
  data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  int i=1;
  sys1->fill_lhs_Euler_impl_expl(dt);
  sys1->fill_lhs1_Euler_impl_expl(dt);
  sys1->fill_lhs2_Euler_impl_expl(dt);
  while(t<t_max)
  {
    if (dt>1e-15)
    {
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      sys1->do_midpoint(x,x_temp,t, dt);
      calc_dist(x_temp);
      const auto [min, max] = std::minmax_element(length.begin(), length.end());
     // x=x_temp;
      renorm(x_temp,x);
//       calc_dist(x);
//       const auto [min1, max1] = std::minmax_element(length.begin(), length.end());
      t+=dt;
      vec_t.push_back(t);
      data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
      
      if(i%100==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<< '\n';
        cout <</**min1*Ns<<" "<<*max1*Ns<<" "<<*/1-*min*Ns<< " "<<*max*Ns-1<<endl;
      }
      i++;
    }
    else
    {
      cout <<"Doing something wrong!"<<endl;abort();
    
    }
  }
  dt=dt1;

   print_data(1e-6,t_max);
   anal_data();
   print_pars();
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
}

void simulation::anal_data()
{
  double vel=0.0, omega=0.0, R=0.0, N_rot=0.0;
  bool helix=true;
  ostringstream os1, os2;
  os1<<folder<<"/ranges.txt";
  Eigen::Vector3d maxVal, minVal, mean;
  dvec x_max, x_min,y_max, y_min,z_max, z_min, radius; 
  dvec helx, hely, helz,hel_angle,hel_angle1; helx.resize(Np);hely.resize(Np);helz.resize(Np);hel_angle.resize(Np);
  vector<Eigen::Vector3d> center, velocity;
  velocity.resize(data.size()-1);
  radius.resize(data.size());
  center.reserve(data.size());
  x_max.reserve(data.size()); x_min.reserve(data.size());
  y_max.reserve(data.size()); y_min.reserve(data.size());
  z_max.reserve(data.size()); z_min.reserve(data.size());
  for(size_t i=0; i<data.size(); i++)
  {
    maxVal = data[i].colwise().maxCoeff();
    minVal = data[i].colwise().minCoeff();
    center.push_back(data[i].colwise().mean());
    x_max.push_back(maxVal(0));x_min.push_back(minVal(0));
    y_max.push_back(maxVal(1));y_min.push_back(minVal(1));
    z_max.push_back(maxVal(2));z_min.push_back(minVal(2));
  }
  for(size_t i=0; i<data.size(); i++){
    for(size_t j=0; j<Np; j++){
      radius[i]+=sqrt((data[i](j,0)-center[i][0])*(data[i](j,0)-center[i][0])+(data[i](j,1)-center[i][1])*(data[i](j,1)-center[i][1])+(data[i](j,2)-center[i][2])*(data[i](j,2)-center[i][2]));
    }
    radius[i]/=Np;
  }
  ofstream file, file1;
  file.open(os1.str().c_str());
  std::vector<double>::iterator min, max;
  if(file.is_open())
  {
    min = std::min_element(x_min.begin(), x_min.end());
    file<<"xmin "<<*min<<endl;
    max = std::max_element(x_max.begin(), x_max.end());
    file<<"xmax "<<*max<<endl;
    min = std::min_element(y_min.begin(), y_min.end());
    file<<"ymin "<<*min<<endl;
    max = std::max_element(y_max.begin(), y_max.end());
    file<<"ymax "<<*max<<endl;
    min = std::min_element(z_min.begin(), z_min.end());
    file<<"zmin "<<*min<<endl;
    max = std::max_element(z_max.begin(), z_max.end());
    file<<"zmax "<<*max<<endl;
    file.close();
  }
  
  os1.clear();
  os2<<folder<<"/vel.txt";
  file1.open(os2.str().c_str());
  file1.precision(10);
  if(file1.is_open()){
    for(size_t i=0; i<data.size()-1; i++)
    {
      velocity[i]=(center[i+1]-center[i])/(vec_t[i+1]-vec_t[i]);
      file1<<vec_t[i+1]<<" "<<velocity[i][0]<<" "<<velocity[i][1]<<" "<<velocity[i][2]<<" "<<radius[i+1]<<endl;
    }
    file1.close();
  }
  else{
    cout<<"Can not create file"<<endl;abort();
  }
  for(size_t i=0;i<Np; i++){
    helx[i]=data[data.size()-1](i,0)-center[data.size()-1][0];
    hely[i]=data[data.size()-1](i,1)-center[data.size()-1][1];
    helz[i]=data[data.size()-1](i,2)-center[data.size()-1][2];
    hel_angle[i]=atan2( hely[i], helz[i]);
  }
  hel_angle1=hel_angle;
  double mult=0.0;
  double total2pi =0;
  double twopi=2*M_PI;
  for(size_t i=0;i<Np-1; i++){
    if(fabs(hel_angle[i+1]-hel_angle[i]+total2pi)>fabs(hel_angle[i+1]-hel_angle[i]+total2pi+twopi))
    {
      total2pi+=twopi;
    }
    else{
      if(fabs(hel_angle[i+1]-hel_angle[i]+total2pi)>fabs(hel_angle[i+1]-hel_angle[i]+total2pi-twopi)){
        total2pi-=twopi;
      }
    }
    hel_angle[i+1]+=total2pi;
        
  }
  os2=ostringstream();
  os2<<folder<<"/shape.txt";
  file1.open(os2.str().c_str());
  file1.precision(10);
  if(file1.is_open()){
    for(size_t i=0; i<Np; i++)
    {
      file1<<helx[i]<<" "<<hely[i]<<" "<<helz[i]<<" "<<hel_angle[i]<<" "<<hel_angle1[i]<<endl;
      x[3*i]=helx[i]; x[3*i+1]=hely[i]; x[3*i+2]=helz[i];
    }
    file1.close();
  }
  else{
    cout<<"Can not create file "<<os2.str()<<endl;abort();
  }
  vel=velocity[data.size()-2][0];
  if(fabs(vel)>1.0)
  {
    int j=ceil(0.9*data.size());
    if(fabs(velocity[j][0]-vel)>0.001*fabs(vel))
    {
      helix=false;
      cout<<"Here2!!"<<endl;
    }
    else{
      for(size_t k=j+1;k<data.size()-3; k++ ){
          if(fabs(velocity[k][0]-vel)>0.001*fabs(vel))
          {
            helix=false;
            cout<<"Here3!!"<<endl;
          }
      }
    }
  }
  else{
    helix=false;
    cout<<"Here4!! "<<vel <<endl;
  }
  os2=ostringstream();
  if(helix){
    R=radius[data.size()-1];
    N_rot=(hel_angle[Np-1]-hel_angle[0])/(2*M_PI);
    size_t j=data.size()-3;
    if(velocity[data.size()-2][0]>velocity[data.size()-1][0])
    {
      
      while((velocity[j][0]>velocity[data.size()-1][0])&&(velocity[j+1][0]<velocity[j][0]))
      {
        j--;
      }
    }
    else{
      while((velocity[j][0]>velocity[data.size()-1][0])&&(velocity[j+1][0]<velocity[j][0]))
      {
        j--;
      }
    }
    omega=2*(M_PI)/(vec_t[data.size()-1]-vec_t[j-1]);
    
  }
  os2<<helix<<" "<<vel<<" "<<omega<<" "<<R<<" "<<N_rot<<endl;
  cout<<"The helix par: "<<os2.str();
  
}


string simulation::anal_data1()
{
  cout<<"starting anal_data1()"<<endl;
  double vel=0.0, omega=0.0, R=0.0, N_rot=0.0;
  bool helix=true;
  ostringstream os1, os2;
  os1<<folder<<"/ranges.txt";
  Eigen::Vector3d maxVal, minVal, mean;
  dvec x_max, x_min,y_max, y_min,z_max, z_min, radius; 
  dvec helx, hely, helz,hel_angle,hel_angle1, omega3, curv; helx.resize(Np);hely.resize(Np); helz.resize(Np);hel_angle.resize(Np);omega3.resize(Np);curv.resize(Np);
  vector<Eigen::Vector3d> center, velocity;
  velocity.resize(data.size()-1);
  radius.resize(data.size());
  center.reserve(data.size());
  x_max.reserve(data.size()); x_min.reserve(data.size());
  y_max.reserve(data.size()); y_min.reserve(data.size());
  z_max.reserve(data.size()); z_min.reserve(data.size());
  for(size_t i=0; i<data.size(); i++)
  {
    maxVal = data[i].colwise().maxCoeff();
    minVal = data[i].colwise().minCoeff();
    center.push_back(data[i].colwise().mean());
    x_max.push_back(maxVal(0));x_min.push_back(minVal(0));
    y_max.push_back(maxVal(1));y_min.push_back(minVal(1));
    z_max.push_back(maxVal(2));z_min.push_back(minVal(2));
  }
  for(size_t i=0; i<data.size(); i++){
    for(size_t j=0; j<Np; j++){
      radius[i]+=sqrt((data[i](j,0)-center[i][0])*(data[i](j,0)-center[i][0])+(data[i](j,1)-center[i][1])*(data[i](j,1)-center[i][1])+(data[i](j,2)-center[i][2])*(data[i](j,2)-center[i][2]));
    }
    radius[i]/=Np;
  }
  ofstream file, file1;
  file.open(os1.str().c_str());
  std::vector<double>::iterator min, max;
  if(file.is_open())
  {
    min = std::min_element(x_min.begin(), x_min.end());
    file<<"xmin "<<*min<<endl;
    max = std::max_element(x_max.begin(), x_max.end());
    file<<"xmax "<<*max<<endl;
    min = std::min_element(y_min.begin(), y_min.end());
    file<<"ymin "<<*min<<endl;
    max = std::max_element(y_max.begin(), y_max.end());
    file<<"ymax "<<*max<<endl;
    min = std::min_element(z_min.begin(), z_min.end());
    file<<"zmin "<<*min<<endl;
    max = std::max_element(z_max.begin(), z_max.end());
    file<<"zmax "<<*max<<endl;
    file.close();
  }
  
  os1.clear();
  os2<<folder<<"/vel.txt";
  file1.open(os2.str().c_str());
  file1.precision(10);
  if(file1.is_open()){
    for(size_t i=0; i<data.size()-1; i++)
    {
      velocity[i]=(center[i+1]-center[i])/(vec_t[i+1]-vec_t[i]);
      file1<<vec_t[i+1]<<" "<<velocity[i][0]<<" "<<velocity[i][1]<<" "<<velocity[i][2]<<" "<<radius[i+1]<<endl;
    }
    file1.close();
  }
  else{
    cout<<"Can not create file"<<endl;abort();
  }
  for(size_t i=0;i<Np; i++){
    omega3[i]=data1[data.size()-1](i);
    helx[i]=data[data.size()-1](i,0)-center[data.size()-1][0];
    hely[i]=data[data.size()-1](i,1)-center[data.size()-1][1];
    helz[i]=data[data.size()-1](i,2)-center[data.size()-1][2];
    hel_angle[i]=atan2( hely[i], helz[i]);
  }
  hel_angle1=hel_angle;
  double mult=0.0;
  double total2pi =0;
  double twopi=2*M_PI;
  for(size_t i=0;i<Np-1; i++){
    if(fabs(hel_angle[i+1]-hel_angle[i]+total2pi)>fabs(hel_angle[i+1]-hel_angle[i]+total2pi+twopi))
    {
      total2pi+=twopi;
    }
    else{
      if(fabs(hel_angle[i+1]-hel_angle[i]+total2pi)>fabs(hel_angle[i+1]-hel_angle[i]+total2pi-twopi)){
        total2pi-=twopi;
      }
    }
    hel_angle[i+1]+=total2pi;
        
  }
  double ll=1.0/Ns;ll*=ll;
  for(size_t i=1;i<Np-1; i++){
    double c2=(helx[i+1]-helx[i-1])*(helx[i+1]-helx[i-1])+(hely[i+1]-hely[i-1])*(hely[i+1]-hely[i-1])+(helz[i+1]-helz[i-1])*(helz[i+1]-helz[i-1]);
   cout<<i<<" "<<c2<<endl;
    double c=sqrt(c2);
    //curv[i]=sqrt(4-(2-c2)*(2-c2))/c;
    curv[i]=sqrt(std::max(4*ll-c2,0.0))/ll;
  }
  const auto [mincurv, maxcurv] = std::minmax_element(curv.begin(), curv.end());


  os2=ostringstream();
  os2<<folder<<"/shape.txt";
  file1.open(os2.str().c_str());
  file1.precision(10);
  if(file1.is_open()){
    for(size_t i=0; i<Np; i++)
    {
      file1<<helx[i]<<" "<<hely[i]<<" "<<helz[i]<<" "<<omega3[i]<<" "<<hel_angle[i]<<" "<<hel_angle1[i]<<" "<<curv[i]<<endl;
    }
    file1.close();
  }
  else{
    cout<<"Can not create file "<<os2.str()<<endl;abort();
  }
  vel=velocity[data.size()-2][0];
  if(fabs(vel)>1.0)
  {
    int j=ceil(0.9*data.size());
    if(fabs(velocity[j][0]-vel)>0.001*fabs(vel))
    {
      helix=false;
      cout<<"Here2!!"<<endl;
    }
    else{
      for(size_t k=j+1;k<data.size()-3; k++ ){
          if(fabs(velocity[k][0]-vel)>0.001*fabs(vel))
          {
            helix=false;
            cout<<"Here3!!"<<endl;
          }
      }
    }
  }
  else{
    helix=false;
    cout<<"Here4!! "<<vel <<endl;
  }
  os2=ostringstream();
  size_t j=0;
  if(helix){
    R=radius[data.size()-1];
    N_rot=(hel_angle[Np-1]-hel_angle[0])/(2*M_PI);
    j=velocity.size()-3;
    if(velocity[velocity.size()-2][1]>velocity[velocity.size()-1][1]){      
      while(!((velocity[j][1]>velocity[velocity.size()-1][1]&&velocity[j+1][1]<velocity[velocity.size()-1][1]))&&j>1){
        j--;
      }
    }
    else{
      cout<<"LOW"<<endl;
      while(!((velocity[j][1]>velocity[velocity.size()-1][1]&&velocity[j+1][1]<velocity[velocity.size()-1][1]))&&j>1){
        j--;
      }
    }
    double ax=fabs(velocity[j+1][1]-velocity[velocity.size()-1][1])/fabs(velocity[j+1][1]-velocity[j][1]);
    double ay=fabs(velocity[j][1]-velocity[velocity.size()-1][1])/fabs(velocity[j+1][1]-velocity[j][1]);
    omega=2*(M_PI)/(vec_t[data.size()-1]-(ax*vec_t[j]+ay*vec_t[j+1]));
    
  }
  min = std::min_element(x_min.begin(), x_min.end());
  max = std::max_element(x_max.begin(), x_max.end());
  os2<<helix<<" "<<vel<<" "<<omega<<" "<<R<<" "<<N_rot<<" "<<*maxcurv<<" "<<data.size()-1-j<<" "<<j<<" "<<vec_t[j]<<" "<<velocity[velocity.size()-1][1]<<" "<<velocity[j][1]<<fabs(*min-*max) <<endl;
  cout<<"The helix par: "<<os2.str();
  cout<<*maxcurv<<endl;
  //abort();
  return os2.str();
  
}

double simulation::dist(const  Eigen::VectorXd &x , int i)
{
  return sqrt((x(3*(i+1))-x(3*i))*(x(3*(i+1))-x(3*i))+(x(3*(i+1)+1)-x(3*i+1))*(x(3*(i+1)+1)-x(3*i+1))+(x(3*(i+1)+2)-x(3*i+2))*(x(3*(i+1)+2)-x(3*i+2)));
}

void simulation::calc_dist(const  Eigen::VectorXd &x)
{
  for(size_t i=0; i<Ns; i++)
  {
    length[i]=dist(x,i);
  }
}


void simulation::renorm(const  Eigen::VectorXd &x1, Eigen::VectorXd &x2)
{
  double att, renorm, alpha1, alpha2;
  x2=x1;
//   for (size_t i=0; i<Ns; i++)
//   {
//     att=dist(x2,i);
//     renorm=1.-1.0/(Np)/att;
//     alpha1=(Np-(i+1))*renorm/Np;
//     alpha2=(i+1)*renorm/Np;
//     for(size_t j=0; j<=i; j++)
//     {
//       x2(3*j)+=alpha1*(x2(3*(i+1))-x2(3*i));
//       x2(3*j+1)+=alpha1*(x2(3*(i+1)+1)-x2(3*i+1));
//       x2(3*j+2)+=alpha1*(x2(3*(i+1)+2)-x2(3*i+2));
//     }
//     for(size_t j=i+1; j<Np; j++)
//     {
//       x2(3*j)-=alpha2*(x2(3*(i+1))-x2(3*i));
//       x2(3*j+1)-=alpha2*(x2(3*(i+1)+1)-x2(3*i+1));
//       x2(3*j+2)-=alpha2*(x2(3*(i+1)+2)-x2(3*i+2));
//     }
//   }
   for (size_t i=0; i<Ns; i++){
     att=dist(x1,i);
     renorm=1.-1.0/(Ns)/att;
     alpha1=(Np-(i+1))*renorm/Np;
     alpha2=(i+1)*renorm/Np;
     for(size_t j=0; j<=i; j++){
       x2(3*j)+=alpha1*(x1(3*(i+1))-x1(3*i));
       x2(3*j+1)+=alpha1*(x1(3*(i+1)+1)-x1(3*i+1));
       x2(3*j+2)+=alpha1*(x1(3*(i+1)+2)-x1(3*i+2));
     }
     for(size_t j=i+1; j<Np; j++){
       x2(3*j)-=alpha2*(x1(3*(i+1))-x1(3*i));
       x2(3*j+1)-=alpha2*(x1(3*(i+1)+1)-x1(3*i+1));
       x2(3*j+2)-=alpha2*(x1(3*(i+1)+2)-x1(3*i+2));
     }
   }
}

void simulation::print_data1( double step, double tmax1)
{

  FILE *file;int j=0;
  file=NULL;
  ostringstream os;
  os.str("");
  os<<folder<<"/d_"<<j<<".txt";
  file=fopen(os.str().c_str(),"w");
  if(file!=NULL)
  {
    for(int ii=0; ii<Np;ii++)
    {
      fprintf(file,"%.14e %.14e %.14e %.14e\n",data[j](ii,0), data[j](ii,1), data[j](ii,2), data1[j](ii));
    };
    fclose(file);
  }
  j++;
  double ax,ay;
  for(int k=1; k<vec_t.size()&&vec_t[k]<tmax1;k++)
  {
      if(vec_t[k]>step*j)
      {
        //cout <<data[0][k]<<" "<<step*j<<endl;
        ax=fabs(vec_t[k]-step*j)/fabs(vec_t[k]-vec_t[k-1]);
        ay=fabs(vec_t[k-1]-step*j)/fabs(vec_t[k]-vec_t[k-1]);
        os.str("");
        os<<folder<<"/d_"<<j<<".txt";
        file=fopen(os.str().c_str(),"w");
        if(file!=NULL)
        {
          for(int ii=0; ii<Np;ii++)
          {
            fprintf(file,"%.14e %.14e %.14e %.14e\n",ay*data[k](ii,0)+ax*data[k-1](ii,0), ay*data[k](ii,1)+ax*data[k-1](ii,1),ay*data[k](ii,2)+ax*data[k-1](ii,2),ay*data1[k](ii)+ax*data1[k-1](ii) );
          };
          fclose(file);
        }
        j++;
        k--;
      }
  }
}


void simulation::print_data( double step, double tmax1)
{
  
  FILE *file;int j=0;
  file=NULL;
  ostringstream os;
  os.str("");
  os<<folder<<"/d_"<<j<<".txt";
  file=fopen(os.str().c_str(),"w");
  if(file!=NULL)
  {
    for(int ii=0; ii<Np;ii++)
    {
      fprintf(file,"%.14e %.14e %.14e\n",data[j](ii,0), data[j](ii,1), data[j](ii,2));
    };
    fclose(file);
  }
  j++;
  double ax,ay;
  for(int k=1; k<vec_t.size()&&vec_t[k]<tmax1;k++)
  {
      if(vec_t[k]>step*j)
      {
        //cout <<data[0][k]<<" "<<step*j<<endl;
        ax=fabs(vec_t[k]-step*j)/fabs(vec_t[k]-vec_t[k-1]);
        ay=fabs(vec_t[k-1]-step*j)/fabs(vec_t[k]-vec_t[k-1]);
        os.str("");
        os<<folder<<"/d_"<<j<<".txt";
        file=fopen(os.str().c_str(),"w");
        if(file!=NULL)
        {
          for(int ii=0; ii<Np;ii++)
          {
            fprintf(file,"%.14e %.14e %.14e\n",ay*data[k](ii,0)+ax*data[k-1](ii,0), ay*data[k](ii,1)+ax*data[k-1](ii,1),ay*data[k](ii,2)+ax*data[k-1](ii,2) );
          };
          fclose(file);
        }
        j++;
        k--;
      }
  }
}

void simulation::print_pars()
{
  cout<<"Printing simulation parameters."<<endl;
    ostringstream os1;
  os1<<"Pars.txt";
  ofstream file;
  file.open(os1.str().c_str());
  if(file.is_open())
  {
    file <<"dt "<<dt<<endl;
    file <<"om "<<om<<endl;
    file <<"Cm "<<Cm<<endl;
    file <<"Ns "<<Ns<<endl;
    file <<"H1 "<<H1<<endl;
    file <<"Torq "<<torq<<endl;
    cout<<"Printing simulation parameters1."<<endl;
    file.close();
    cout<<"Printing simulation parameters1."<<endl;
  }
  else
  {
    cout <<"Serious error: Can not open file"<<endl;
  }
  cout<<"Parameters printed1"<<endl;
  os1.str("");
  os1<<folder<<"/Pars.txt";
  cout<<"Parameters printed2"<<endl;
  file.open(os1.str().c_str());
  if(file.is_open())
  {
    file <<"dt "<<dt<<endl;
    file <<"om "<<om<<endl;
    file <<"Cm "<<Cm<<endl;
    file <<"Ns "<<Ns<<endl;
    file <<"H1 "<<H1<<endl;
    file <<"Torq "<<torq<<endl;
    file.close();
  }
  cout<<"Parameters printed"<<endl;
}

bool simulation::accept_step(double err)
{  
  calc_dist(x_temp);
  const auto [min, max] = std::minmax_element(length.begin(), length.end());
  double a=1-*min*Ns, b=*max*Ns-1;
  max_rel_error=std::max(a,b);
  if(max_rel_error>err)
  {
    return false;
  }
  else
  {
    return true;
  }
}

// double simulation::max_error()
// {  
//   calc_dist(x_temp);
//   const auto [min, max] = std::minmax_element(length.begin(), length.end());
//   double a=1-*min*Ns, b=*max*Ns-1;
//   return std::max(a,b);
// }

bool simulation::accept_step_beginning()
{  
  calc_dist(x_temp);
  const auto [min, max] = std::minmax_element(length.begin(), length.end());
  double a=1-*min*Ns, b=*max*Ns-1;
  max_rel_error=std::max(a,b);
  if(max_rel_error>1)
  {
    cout<<"Programm crashed do to too big initial timestep!!!!!! Reduce it!!!!!"<<endl;
    abort();
    return false;
  }
  else
  {
    return true;
  }
}

void simulation::fill_0mega3()
{
  double c1=x[0]-x[3*Ns];
  //double c1=x1[0]-x1[3*Ns];
  double c2=-0.5*(x[0]+x[3*Ns]);
  //double c2=-0.5*(x1[0]+x1[3*Ns]);
  double h=1.0/Ns;
  for(size_t i=0; i<Np; i++){
    Omega3[i]=-torq*(x[3*i]+c1*(-0.5+h*i)+c2);
  }
}

