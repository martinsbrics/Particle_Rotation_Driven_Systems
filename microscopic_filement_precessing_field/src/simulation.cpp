// SPDX-FileCopyrightText: 2023 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#include "simulation.h"


simulation::simulation(int N_, double omega_B_, double beta_,  bool B_off, Eigen::VectorXd x_):Np(N_),Ns(N_-1), Nv(3*N_),
omega_B(omega_B_),beta(beta_), length(Ns), max_rel_error(0.0), x(x_)
{
  sys=new Integrator(Np, omega_B, beta,  B_off );
  x_temp=x;
  renorm(x_temp,x);
}

simulation::~simulation()
{
  delete sys;
}

void simulation::renorm(const  Eigen::VectorXd &x1, Eigen::VectorXd &x2)
{
    double att, renorm, alpha1, alpha2;
    x2=x1;
    for (size_t i=0; i<Ns; i++){
        att=dist(x1,i);
        renorm=1.-1.05002/att;
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

bool simulation::accept_step(double err)
{
  calc_dist(x_temp);
  const auto [min, max] = std::minmax_element(length.begin(), length.end());
  double a=1-*min/1.05002, b=*max/1.05002-1;
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



void simulation::propagate_VSIMEX_BDF3_DD(double dt1, double tmax1, double tol)
{
  double tol_max=tol/640;
  t=0; dt=dt1, t_max=tmax1;
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
  ostringstream os1;
  os1<<"data_filament_VSIMEX_BDF3_DD_Ns_"<<Ns<<"_omega_B_"<<omega_B<<"_beta_"<<beta<<"_tol_"<<tol;
  folder=os1.str();
  os1.str("");
  os1<<"mkdir "<<folder;
  system(os1.str().c_str());
  os1.str("");
  FILE *file;
  file=NULL;
  int print=(int) max(( 1e-5/dt),1.);

  //vec_t.push_back(t);
  //data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  now= std::chrono::system_clock::now();
//   data[0].push_back(t);
//   for(int ii=0; ii<Np;ii++)
//   {
//     data[ii+1].push_back(x[ii]);
//     data[ii+Np+1].push_back(x[ii+Np]);
//     data[ii+2*Np+1].push_back(x[ii+2*Np]);
//   }
  Eigen::VectorXd x_old, x_old_old;
  double dt_old, dt_old_old;
  int i=1;
  sys->fill_lhs_Euler_impl_expl(dt);

  x_old_old=x;
  while(true)
  {
    sys->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<0.1*tol)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys->fill_lhs_Euler_impl_expl(dt);
    if(dt<1e-13)
    {
      abort();
    }
  }
  t+=dt;
  //vec_t.push_back(t);
  //data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  dt_old_old=dt;

  x_old=x;
  while(true)
  {
    sys->do_Euler_impl_expl(x,x_temp,t, dt);
    accept_step();
    cout <<"dt="<<dt<<" error="<<max_rel_error<<endl;
    if(max_rel_error<0.1*tol)
    {
      renorm(x_temp,x);
      break;
    }
    dt/=2;
    sys->fill_lhs_Euler_impl_expl(dt);
  }
  t+=dt;
  //vec_t.push_back(t);
  //data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
  dt_old=dt;
  sys->fill_K1_IMEX_BDF3(x_old_old, dt, t-dt_old-dt_old_old);
  sys->K2=sys->K1;
  sys->fill_K1_IMEX_BDF3(x_old, dt, t-dt_old);
  cout<<"calculated inital steps!!"<<endl;
  int changed_timestep=0;
  int increase=1000;
  sys->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
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
      sys->do_VSIMEX_BDF3(x,x_old,x_old_old,x_temp, t, dt);
      if (accept_step(tol))
      {
        if(t>j*1e-4/*true*/)
        {
          j++;
          cout << t << '\t' << dt <<" "<<dt_old<<" "<<dt_old_old<<" "<<max_rel_error<<endl;
          //fprintf(file,"%.7e %.14e %.14e\n",t, dt, max_rel_error );
        }
        if(changed_timestep<=0)
        {
          dt_old_old=dt_old;
          dt_old=dt;
          sys->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
        }
        x_old_old=x_old;
        x_old=x;
        renorm(x_temp,x);
        t+=dt;
        if(t>j_store*1e-5/*true*/){
          vec_t.push_back(t);
          data.push_back((Eigen::Map<Eigen::MatrixXd> (x.data(), 3,Np)).transpose());
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
            sys->K1*=2; sys->K2*=2;
            sys->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
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
        sys->fill_K1_IMEX_BDF3(x_old_old, dt/2, t-dt-dt_old);
        sys->K2=sys->K1;
        sys->fill_K1_IMEX_BDF3(x_old, dt/2, t-dt);
        dt/=2;
        cout<<" reducing timestep dt="<<dt<< " "<<max_rel_error<<endl;
        sys->fill_lhs2_VSIMEX_BDF3(dt, dt_old, dt_old_old);
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
   //fclose(file);
   cout<<"!!!!!!!!!!!!!!!!!dt=" <<dt<<endl;
//   print_data_sigle_file(0.1,t_max);
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
    file <<"om "<<omega_B<<endl;
    file <<"al "<<beta<<endl;
    file <<"Ns "<<Ns<<endl;
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
    file <<"om "<<omega_B<<endl;
    file <<"al "<<beta<<endl;
    file <<"Ns "<<Ns<<endl;
    file.close();
  }
  cout<<"Parameters printed"<<endl;
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
            fprintf(file,"%.14e %.14e %.14e\n",ay*data[k](ii,0)+ax*data[k-1](ii,0), ay*data[k](ii,1)+ax*data[k-1](ii,1),ay*data[k](ii,2)+ax*data[k-1](ii,2));
          };
          fclose(file);
        }
        j++;
        k--;
      }
  }
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
