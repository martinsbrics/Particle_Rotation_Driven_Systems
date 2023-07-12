
using namespace std;
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream> 
#include <Eigen/Core>
#include "integrator.h"
#include "eigen_algebra.hpp"
using namespace boost::numeric::odeint;
namespace boost {
namespace numeric {
namespace odeint {

template<typename B,int S1,int S2,int O, int M1, int M2>
struct algebra_dispatcher< Eigen::Matrix<B,S1,S2,O,M1,M2> >
{
    typedef vector_space_algebra algebra_type;
};

}}}


bool propagate(Eigen::VectorXd &x, Eigen::VectorXd &dxdt, vector<dvec>& data,  Integrator *sys, double dt1, double tmax1, double abs_tol,double rel_tol )
{
  cout <<"Startng propagation! x.size()="<<x.size()<<endl;
//   ostringstream os1;
//   os1<<"1data_Np_"<<Np<<"_gamma_"<<gamma<<"_lambda_"<<lambda;
//   folder=os1.str();
//   os1.str("");
//   os1<<"mkdir "<<folder;
//   system(os1.str().c_str());
  double dt=dt1; double t_max=tmax1;
  //typedef runge_kutta_fehlberg78<Eigen::VectorXd> error_stepper_type;
  typedef runge_kutta_dopri5< Eigen::VectorXd> error_stepper_type;
  //typedef runge_kutta_cash_karp54< Eigen::VectorXd> error_stepper_type;
  typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
  error_stepper_type rk;
  std::chrono::time_point<std::chrono::system_clock> foo;
  std::chrono::duration<double> diff1;
  std::chrono::time_point<std::chrono::system_clock> now;
  double xx=0.0;
  double xy=0.0;
  double ax = 1.0;
  double adxdt = 1.0;
  double t=0;
  size_t i=1;
  size_t acpeted_steps=0;
  size_t ignored_steps=0;
  size_t limiter_steps=0;
  size_t N=x.size()/6;
  dvec a,b;a.resize(N-1);b.resize(N-1);

  controlled_stepper_type cs(
          default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
  controlled_step_result res;
  now= std::chrono::system_clock::now();
  data[0].push_back(t);
  for(int ii=0; ii<x.size();ii++)
  {
    data[ii+1].push_back(x[ii]);
  }
  data[1+x.size()].push_back(sys->get_Bx(t));
  data[2+x.size()].push_back(sys->get_By(t));
  data[3+x.size()].push_back(sys->get_phi(t)-x[2]-sys->get_mom_angle());
  if(N>1){
    xx=(x[6]-x[0])*cos(x[6+2])+(x[6+1]-x[1])*sin(x[6+2]);
    xy=-(x[6]-x[0])*sin(x[6+2])+(x[6+1]-x[1])*cos(x[6+2]);
  }
  data[4+x.size()].push_back(xx);
  //data[5+x.size()].push_back(std::sqrt(xx*xx+xy*xy));
  data[5+x.size()].push_back(xy);
  for(int ii=0; ii<2;ii++)
//   {
//     data[ii*4+1+x.size()].push_back(sys->get_HD_force(ii)[0]);
//     data[ii*4+2+x.size()].push_back(sys->get_HD_force(ii)[1]);
//     data[ii*4+3+x.size()].push_back(sys->get_steric_force(ii)[0]);
//     data[ii*4+4+x.size()].push_back(sys->get_steric_force(ii)[1]);
//   }
  while(t<t_max)
  {
    if (dt>1e-10)
    {
//       if(dt>dt1)
//       {
//         dt=dt1;
//       }
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      res= cs.try_step(*sys, x, t,dt);
      if(res==0)
      {
        acpeted_steps++;
        data[0].push_back(t);
        for(int ii=0; ii<x.size();ii++)
        {
          data[ii+1].push_back(x[ii]);
        }
        data[1+x.size()].push_back(sys->get_Bx(t));
        data[2+x.size()].push_back(sys->get_By(t));
        data[3+x.size()].push_back(sys->get_phi(t)-x[5]-sys->get_mom_angle());
        if(N>1){
          xx=(x[12]-x[0])*cos(x[12+2])+(x[12+1]-x[1])*sin(x[12+2]);
          xy=-(x[12]-x[0])*sin(x[12+2])+(x[12+1]-x[1])*cos(x[12+2]);
        }
        data[4+x.size()].push_back(xx);
        //data[5+x.size()].push_back(std::sqrt(xx*xx+xy*xy));
        data[5+x.size()].push_back(xy);
//         for(int ii=0; ii<2;ii++)
//         {
//           data[ii*4+1+x.size()].push_back(sys->get_HD_force(ii)[0]);
//           data[ii*4+2+x.size()].push_back(sys->get_HD_force(ii)[1]);
//           data[ii*4+3+x.size()].push_back(sys->get_steric_force(ii)[0]);
//           data[ii*4+4+x.size()].push_back(sys->get_steric_force(ii)[1]);
//         }
        i++;
      }
      else
      {
        ignored_steps++;
      }

      if(i%100000==0)
      //if(i%1==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<<" "<< i/5000000.0<<"\% of max\n";
        if(N>1){
          for(size_t ix=1; ix<N; ix++){
          a[ix-1]=(x[6*ix]-x[6*(ix-1)])*cos(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*sin(x[6*ix+2]);
          b[ix-1]=-(x[6*ix]-x[6*(ix-1)])*sin(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*cos(x[6*ix+2]);
          }
          cout<<*max_element(a.begin(), a.end())<<" "<<*min_element(a.begin(), a.end())<<" "<<*max_element(b.begin(), b.end())<<" "<<*min_element(b.begin(), b.end())<<endl;
        }
        //cout<<sys->get_Mag_torque(x)<<" "<<sys->get_Steric_torque()<<" "<<sys->get_steric_force(0).transpose() <<endl;
        if(i>500000000){
          return false;
        }
      }
    }
    else
    {
      dt=2e-10;
      rk.do_step(*sys, x, t,dt);
      limiter_steps++;
      i++;
      if(limiter_steps>1000)
      {
        return false;
      }
      //abort();
      //cout <<t<<endl;
    }
    
  }
  cout << t << '\t' << dt << '\n';
  foo= std::chrono::system_clock::now();
  diff1 = foo - now;
  cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
  cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}


bool propagate_breakcheck(Eigen::VectorXd &x, Eigen::VectorXd &dxdt, vector<dvec>& data,  Integrator *sys, double dt1, double tmax1, double abs_tol,double rel_tol )
{
  cout <<"Startng propagate_breakchec! omega="<<sys->get_omega()<<endl;
//   ostringstream os1;
//   os1<<"1data_Np_"<<Np<<"_gamma_"<<gamma<<"_lambda_"<<lambda;
//   folder=os1.str();
//   os1.str("");
//   os1<<"mkdir "<<folder;
//   system(os1.str().c_str());
  double dt=dt1; double t_max=tmax1;
  typedef runge_kutta_fehlberg78<Eigen::VectorXd> error_stepper_type;
  //typedef runge_kutta_dopri5< Eigen::VectorXd> error_stepper_type;
  //typedef runge_kutta_cash_karp54< Eigen::VectorXd> error_stepper_type;
  typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
  error_stepper_type rk;
  std::chrono::time_point<std::chrono::system_clock> foo;
  std::chrono::duration<double> diff1;
  std::chrono::time_point<std::chrono::system_clock> now;
  double xx=0.0;
  double xy=0.0;
  double ax = 1.0;
  double adxdt = 1.0;
  double t=0;
  size_t i=1;
  size_t acpeted_steps=0;
  size_t ignored_steps=0;
  size_t limiter_steps=0;
  size_t N=x.size()/6;
  dvec a,b;a.resize(N-1);b.resize(N-1);

  controlled_stepper_type cs(
          default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
  controlled_step_result res;
  now= std::chrono::system_clock::now();
 
  while(t<t_max)
  {
    if (dt>1e-10)
    {
//       if(dt>dt1)
//       {
//         dt=dt1;
//       }
      if(t+dt>t_max)
      {
        dt=t_max-t;
      }
      res= cs.try_step(*sys, x, t,dt);
      if(res==0)
      {
        acpeted_steps++;
        if(N>1){
          for(size_t ix=1; ix<N; ix++){
          a[ix-1]=(x[6*ix]-x[6*(ix-1)])*cos(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*sin(x[6*ix+2]);
          b[ix-1]=-(x[6*ix]-x[6*(ix-1)])*sin(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*cos(x[6*ix+2]);
          }
          //cout<<*max_element(a.begin(), a.end())<<" "<<*min_element(a.begin(), a.end())<<" "<<*max_element(b.begin(), b.end())<<" "<<*min_element(b.begin(), b.end())<<endl;
        }
        if(max(*max_element(a.begin(), a.end()),fabs(*min_element(a.begin(), a.end())))>M_SQRT2||max(*max_element(b.begin(), b.end()),fabs(*min_element(b.begin(), b.end())))>M_SQRT2){
          cout<<"break!!!!!!!!"<<endl;
          //abort();
          return false;
        }
        i++;
      }
      else
      {
        ignored_steps++;
      }

      if(i%100000==0)
      //if(i%1==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<<" "<< i/5000000.0<<"\% of max\n";
        if(N>1){
          for(size_t ix=1; ix<N; ix++){
            a[ix-1]=(x[6*ix]-x[6*(ix-1)])*cos(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*sin(x[6*ix+2]);
            b[ix-1]=-(x[6*ix]-x[6*(ix-1)])*sin(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*cos(x[6*ix+2]);
          }
          cout<<*max_element(a.begin(), a.end())<<" "<<*min_element(a.begin(), a.end())<<" "<<*max_element(b.begin(), b.end())<<" "<<*min_element(b.begin(), b.end())<<endl;
        }
        //cout<<sys->get_Mag_torque(x)<<" "<<sys->get_Steric_torque()<<" "<<sys->get_steric_force(0).transpose() <<endl;
        if(i>500000000){
          return false;
        }
      }
    }
    else
    {
      dt=2e-10;
      rk.do_step(*sys, x, t,dt);
      limiter_steps++;
      i++;
      if(limiter_steps>1000)
      {
        return false;
      }
      //abort();
      //cout <<t<<endl;
    }
    
  }
  cout << t << '\t' << dt << '\n';
  foo= std::chrono::system_clock::now();
  diff1 = foo - now;
  cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
  cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}

void print_data_sigle_file( double step, double tmax1, vector<dvec>& data, string folder )
{
  cout << "Start print_data single file"<<endl;
  FILE *file;int j=0;
  file=NULL;
  ostringstream os;
  os.str("");
  os<<folder<<".txt";
  file=fopen(os.str().c_str(),"w");
  if(file!=NULL)
  {
    fprintf(file,"%.4f",data[0][j]-data[0][0]);
    for(int ii=0; ii<data.size()-1;ii++)
    {
        fprintf(file," %.14e",data[ii+1][j]);
    }
    fprintf(file,"\n");
    double ax,ay;
    j++;
    for(int k=1; k<data[0].size()&&data[0][k]<tmax1;k++)
    {
      if(data[0][k]>data[0][0]+step*j)
      {
        //cout <<data[0][k]<<" "<<step*j<<endl;
        ax=fabs(data[0][k]-data[0][0]-step*j)/fabs(data[0][k]-data[0][k-1]);
        ay=fabs(data[0][k-1]-data[0][0]-step*j)/fabs(data[0][k]-data[0][k-1]);
        
        fprintf(file,"%.4f",step*j);
        for(int ii=0; ii<data.size()-1;ii++)
        {
          fprintf(file," %.14e",ay*data[ii+1][k]+ax*data[ii+1][k-1]);
        }
          fprintf(file,"\n");
        
      j++;
      k--;
      }
    }
    fclose(file);
  }
}

dvec analyze_data( vector<dvec>& data )
{
  dvec res(6);
  size_t size=data[0].size()-1;
  size_t Ndata=data.size()-1;
  cout<<"start analyzing"<<" "<<size<<endl;
  double t_end=data[0][size];
  double tmid=t_end*0.5;
  size_t i, j, mid; i=0; j=size; mid=0;
  double period=-1.0;
  double Wn_num=-1.0;
  while (i < j-1) {
    mid = (i + j) / 2;
    cout <<i<<" "<<j<<" "<<mid<<" "<<data[0][mid]<<endl;
    if(tmid < data[0][mid]) { 
      j = mid; 
    }     
    else {
      i = mid;  
    } 
  }
  double t_int=t_end-data[0][mid];
  cout<<tmid<<" "<<data[0][i]<<" "<<i<<" " <<size<<endl;
  
  auto mnmx = minmax_element(data[Ndata-1].begin()+mid, data[Ndata-1].end());
  double min_a=*mnmx.first;
  double max_a=*mnmx.second;
  bool period_from_a=true;
  double aver_a=0;
  double mean_a=0.5*(min_a+max_a);
  if(max_a-min_a<1e-5)
  {
    aver_a=0.5*(min_a+max_a);
    period_from_a=false;
  }
  else
  {
    aver_a+=data[Ndata-1][mid]*(data[0][mid+1]-data[0][mid]);
    aver_a+=data[Ndata-1][size]*(data[0][size]-data[0][size-1]);
    for(size_t i=mid+1; i<size; i++){
      aver_a+=data[Ndata-1][i]*(data[0][i+1]-data[0][i-1]);
    }
    aver_a*=0.5/t_int;
  }
  cout<<" min a="<<min_a<<" max a="<<max_a<<" aver_a="<<aver_a<<" mean_a="<<mean_a<<endl;
  i=mid;
  dvec per, Wn;
  //cout<<"here"<<endl;
  double nper=fabs((data[Ndata-2][size]- data[Ndata-2][mid])/(2*M_PI));
  double approx_per;
  if(nper<1.0){
    approx_per=0.0;
  }
  else{
    approx_per=t_int/nper;
  }
  //cout<<"here"<<endl;
  if(period_from_a){
    while(i<size){
      if(data[Ndata-1][i]>mean_a&&data[Ndata-1][i-1]<mean_a){
        break;
      }
      else{
        i++;
      }
    }
    double ax=fabs(mean_a-data[Ndata-1][i])/fabs(data[Ndata-1][i]-data[Ndata-1][i-1]);
    double ay=fabs(mean_a-data[Ndata-1][i-1])/fabs(data[Ndata-1][i]-data[Ndata-1][i-1]);
    double Period_start=ay*data[0][i]+ax*data[0][i-1];
    double Wn_start=ay*data[3][i]+ax*data[3][i-1];
    for(size_t j=0; i<size-2; j++){
      i+=2;
      while(i<size){
        if(data[Ndata-1][i]>mean_a&&data[Ndata-1][i-1]<mean_a){
          break;
        }
        else{
          i++;
        }
      }
      size_t i_end=i;
      if(i<size){
        ax=fabs(mean_a-data[Ndata-1][i])/fabs(data[Ndata-1][i]-data[Ndata-1][i-1]);
        ay=fabs(mean_a-data[Ndata-1][i-1])/fabs(data[Ndata-1][i]-data[Ndata-1][i-1]);
        period=ay*data[0][i]+ax*data[0][i-1];
        Wn_num=ay*data[3][i]+ax*data[3][i-1];
        per.push_back(period-Period_start);
        Wn.push_back((Wn_num-Wn_start)/(2*M_PI));
        Period_start=period;
        Wn_start=Wn_num;
      }
    }
    i=0;
    period=0;
    Wn_num=0;
    if(max(fabs(min_a),fabs(max_a))<M_SQRT2){
      if(nper>1.0){
        while(i<per.size())
        {
          period+=per[i];
          Wn_num+=Wn[i];
          i++;
        }
        period/=per.size();
        Wn_num/=per.size();
      }
    }
  }
  else{
    if(fabs(data[Ndata-2][mid])>2*M_PI)
    {
      double Period_start=data[0][mid];
      double lag_start=fabs(data[Ndata-2][mid]);
      double lag_end=lag_start+(2*M_PI);
      double ax=0.0;
      double ay=0.0;
      double Wn_start=data[3][mid];
      for(size_t j=0; i<size-2; j++){
        i+=2;
        while(i<size){
          if(fabs(data[Ndata-2][i])>lag_end&&fabs(data[Ndata-2][i-1])<lag_end){
            break;
          }
          else{
            i++;
          }
        }
        size_t i_end=i;
        if(i<size){
          ax=fabs(lag_end-fabs(data[Ndata-2][i]))/fabs(data[Ndata-2][i]-data[Ndata-2][i-1]);
          ay=fabs(lag_end-fabs(data[Ndata-2][i-1]))/fabs(data[Ndata-2][i]-data[Ndata-2][i-1]);
          period=ay*data[0][i]+ax*data[0][i-1];
          Wn_num=ay*data[3][i]+ax*data[3][i-1];
          per.push_back(period-Period_start);
          Wn.push_back((Wn_num-Wn_start)/(2*M_PI));
          Period_start=period;
          Wn_start=Wn_num;
          lag_end+=2*M_PI;
        }
      }
      i=0;
      period=0;
      Wn_num=0;
      if(max(fabs(min_a),fabs(max_a))<M_SQRT2)
      {
//         cout<<"!!!!!!!!!!!!! "<<max(fabs(min_a),fabs(max_a))<<endl;
        if(nper>1.0){
          while(i<per.size())
          {
            period+=per[i];
            Wn_num+=Wn[i];
            i++;
          }
          period/=per.size();
          Wn_num/=per.size();
        }
        else{
         period=-1.0; 
         Wn_num=-1.0;
        }
      }
    }
    else{
      period=-1.0;
    }
  }
  cout.precision(10);
  cout<<"period from a is "<<period<<" "<<approx_per<<" "<<nper<<" "<<2*M_PI/period<<" wn="<<Wn_num<<" per.size()="<<per.size()<<endl;
  for(size_t i=0; i<per.size(); i++){
    cout<<per[i]<<" "<<Wn[i]<<endl;
  }
  res[0]=min_a;
  res[1]=max_a;
  res[2]=aver_a;
  res[3]=mean_a;
  res[4]=period;
  res[5]=Wn_num;
  return(res);
}


// void wafunction_complex_1d::print_cacl_fftw_force(string prefix, dvec &force, double tmaxs, double time_interval, int option)
// {
//    fftw_complex *out, *in;
//   fftw_plan p;
//   int size=force.size();
//   
//   cvec local_atocor,en_stectr;
//   if (size>tmaxs/time_interval)
//   {
//     size=tmaxs/time_interval;
//   }
//   if(size%2==1)
//   {
//     size-=1;
//   }
//   cout <<"size="<<size<< " "<<tmaxs/dt<<endl;
//   local_atocor.resize(size); en_stectr.reserve(size);
//   for(int i=0; i<size; i++)
//   {
//     local_atocor[i]=force[i];      
//   }
//     int sz2=size/2;
//   if( option==0)
//   {
//     Hann(local_atocor, tmaxs);
//   }
//   if( option==2)
//   {
//     sin8_window(local_atocor, tmaxs);
//   }
// 
//   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (size));
//   p = fftw_plan_dft_1d(size,  reinterpret_cast<fftw_complex*>(&local_atocor[0]), out,FFTW_FORWARD, FFTW_ESTIMATE);
//   
//   fftw_execute(p); /* repeat as needed */
//   for(int i=0; i<sz2; i++)
//   {
//     en_stectr.push_back(sqrt(0.5/M_PI)*time_interval*dcmplx(out[i+sz2][0], out[i+sz2][1]));
//   }
//   for(int i=0; i<sz2; i++)
//   {
//     en_stectr.push_back(sqrt(0.5/M_PI)*time_interval*dcmplx(out[i][0], out[i][1]));
//   }
//   fftw_destroy_plan(p);
//   fftw_free(out);
//   
//   FILE *file;
//   string fname=prefix;
//   if( option==0)
//   {
//     fname+="_Hann";
//   }
//   if( option==0)
//   {
//     fname+="_sin8";
//   }
//   fname+="_tmax_";
//   fname+=int2string(tmaxs);
//   fname+="_";
//   fname+=ham_opt.Mylaser->get_laser_par();
//   fname+="_";
//   fname+=grid_info();
//   fname+=".txt";
//   char *filename = str2char(fname);
//   file=fopen(filename,"w");
//   double offset=en_stectr.size()/2.0/tmaxs*(2*M_PI);
//   cout <<"offset="<<offset<< " dw="<<1./tmaxs*(2*M_PI)<<endl;
//   double frek;
//   for(int i=0; i<en_stectr.size(); i++)
//   {
//     frek=i*1./tmaxs*(2*M_PI)-offset;
// //     if(fabs(frek)<fmax)
// //     {
//       fprintf(file,"%.14e %.14e %.14e %.14e\n",frek, en_stectr[i].real(),en_stectr[i].imag(),force[i]);
// //    }
//   }
//   fclose(file);
//   delete [] filename;
//   filename=NULL;
// }


// void print_data_w_phi6( double step, double tmax1)
// {
//   dvec xlocal(2*Np);
//   cout << "Start print_data"<<endl;
//   FILE *file, *file1;int j=0;
//   file=NULL;file1=NULL;
//   ostringstream os;
//   os.str("");
//   os<<folder<<"/d_"<<j<<".txt";
//   file=fopen(os.str().c_str(),"w");
//   os.str("");
//   os<<"Phi6_"<<folder<<".txt";
//   file1=fopen(os.str().c_str(),"w");
//   if(file1!=NULL)
//   {
//     if(file!=NULL)
//     {
//       for(int ii=0; ii<Np;ii++)
//       {
//         xlocal[2*ii]=data[ii+1][j];  xlocal[2*ii+1]=data[ii+Np+1][j];
//       };
//       fprintf(file1,"%.6f %.14e %.14e\n",step*j, phi6(xlocal), phi6_local(xlocal));
//       for(int ii=0; ii<Np;ii++)
//       {
//         fprintf(file,"%.14e %.14e %.14e\n",data[ii+1][j], data[ii+Np+1][j], data[ii+2*Np+1][j]);
//       };
//       fclose(file);
//     /*  for(int i=0; i<Np; i++)
//     {
//       cout <<data[i+1][j]<<" "<<data[i+Np+1][j]<<" "<<data[i+2*Np+1][j]<<endl;
//     }
//     abort();*/
//     }
//     j++;
//     double ax,ay;
//     cout<<data[0].size()<<" "<<data[0][0]<<" tmax="<< tmax1<<endl;
//     for(int k=1; k<data[0].size()&&data[0][k]<tmax1;k++)
//     {
//       if(data[0][k]>data[0][0]+step*j)
//       {
//         cout <<j<< " "<<k <<" "<<step*j<<endl;
//         ax=fabs(data[0][k]-data[0][0]-step*j)/fabs(data[0][k]-data[0][k-1]);
//         ay=fabs(data[0][k-1]-data[0][0]-step*j)/fabs(data[0][k]-data[0][k-1]);
//         os.str("");
//         os<<folder<<"/d_"<<j<<".txt";
//         file=fopen(os.str().c_str(),"w");
//         for(int ii=0; ii<Np;ii++)
//         {
//           xlocal[2*ii]=ay*data[ii+1][k]+ax*data[ii+1][k-1];  
//           xlocal[2*ii+1]=ay*data[Np+ii+1][k]+ax*data[Np+ii+1][k-1];
//         };
//         fprintf(file1,"%.6f %.14e %.14e\n",step*j, phi6(xlocal), phi6_local(xlocal));
//         if(file!=NULL)
//         {
//           for(int ii=0; ii<Np;ii++)
//           {
//             fprintf(file,"%.14e %.14e %.14e\n",ay*data[ii+1][k]+ax*data[ii+1][k-1], ay*data[ii+Np+1][k]+ax*data[ii+Np+1][k-1],ay*data[ii+2*Np+1][k]+ax*data[ii+2*Np+1][k-1]);
//             //cout <<ay*data[ii+1][k]+ax*data[ii+1][k-1]<<" "<<ay*data[ii+Np+1][k]+ax*data[ii+Np+1][k-1]<<" "<<ay*data[ii+2*Np+1][k]+ax*data[ii+2*Np+1][k-1]<<endl;
//           };
//           fclose(file);
//         }
//         j++;
//         k--;
//       }
//     }
//     fclose(file1);
//     file1=NULL;
//   }
// }


int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    Eigen::VectorXd x; Eigen::VectorXd dxdt;
    
    
// // // //single write 
//    vector <double> angle;
//     vector <Eigen::Vector2d> pos;
//    vector <dvec> data;
//     int N=10;//number of particles;
//     double alpha=0.94;
//     double komega_kv=1.84;
//      double shiftx=0.35579;
//      double shifty=0.984809;
// //    double shifty=-0.35579;
// //    double shiftx=-0.984809;
//     double omega=-0.05804/*1.15470053838*/;
//     x=Eigen::VectorXd::Zero(6*N);dxdt=Eigen::VectorXd::Zero(6*N);data.resize(6*N+6); ;
//     if(N>1){
//       dvec shift(2*(N-1));
//       for(size_t i=1; i<N; i++){
//         shift[2*(i-1)]=shiftx*i;
//         shift[2*(i-1)+1]=shifty*i;
//       }
//   //     //kink
//   //     shift[0]=shiftx;
//   //     shift[1]=shifty;
//   //     for(size_t i=2; i<N/2+1; i++){
//   //       shift[2*(i-1)]=shift[2*(i-2)]+shifty;
//   //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shiftx;
//   //     }
//   //     for(size_t i=N/2+1; i<N; i++){
//   //       shift[2*(i-1)]=shift[2*(i-2)]+shiftx;
//   //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shifty;
//   //     }
//       shiftx=shift[2*(N-2)]/2;
//       shifty=shift[2*(N-2)+1]/2;
//       cout<<shiftx<<" "<<shifty<<" "<<endl;
//       pos.push_back(Eigen::Vector2d(-shiftx,-shifty));
//       angle.push_back(0.0);
//       for(size_t i=1; i<N; i++){
//         pos.push_back(Eigen::Vector2d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1]));
//         angle.push_back(0.0);
//       }
//     }
//     else{
//       pos.push_back(Eigen::Vector2d(0.0,0.0));
//       angle.push_back(0.0);
//     }
//     double phi0 = atan(M_SQRT1_2) * 180 / M_PI;
//     Eigen::Vector2d mom (1., tan((12+phi0)*M_PI/180.));
//     mom.normalize();
//     cout<<mom<<endl;
//     //abort();
//     Integrator *sys;
//     double T=2*M_PI/fabs(omega);
//     double t_end=100*T/*+asin(fabs(omega))/(2*M_PI)*T*//*50.0*/;
//     sys=new Integrator(N,omega,mom,alpha,komega_kv,0.0, false);
//     //sys=new Integrator();
//     for(size_t i=0; i<N; i++){
//       x[6*i]=pos[i][0];
//       x[6*i+1]=pos[i][1];
//       x[6*i+2]=angle[i];
//     }
//     cout<<x<<endl;
//     cout<<"Start ptopagation"<<endl;
//     propagate(x, dxdt, data, sys, 0.01,t_end,1e-14, 1e-14);
//     dvec res=analyze_data( data);
//     //abort();
// //    sys->fill_steric_foreces_and_torques();
//     cout <<"finished propagation x[7]="<<x[7]<<endl;
//     cout<< x<<endl;
//     for(size_t i=1; i<N; i++){
//       double a=(x[6*i]-x[6*(i-1)])*cos(x[6*i+2])+(x[6*i+1]-x[6*(i-1)+1])*sin(x[6*i+2]);
//       double b=-(x[6*i]-x[6*(i-1)])*sin(x[6*i+2])+(x[6*i+1]-x[6*(i-1)+1])*cos(x[6*i+2]);
//       cout<<a<<" "<<b<<" "<<0.5*sqrt(a*a+b*b)<<endl;
//     }
//     double lag=sys->get_phi(t_end)-x[2]-sys->get_mom_angle();
//     cout<<" "<<lag<<endl;
//     sys->print_forces(x, dxdt,t_end);
//     delete sys;
//     //string fname="x";
//     stringstream fname;
//     fname<<"x_"<<N<<"shift_alpha_"<<alpha<<"_komega_kv_"<<komega_kv<<"omega_"<<omega/*<<".txt"*/;
// //     file.open(fname.str().c_str());
//     print_data_sigle_file(0.01*T,100000*T, data, fname.str());
//     cout<<"omega="<<omega<<endl;
//     for(int i=0; i<res.size(); i++){
//     cout<<" "<<res[i];
//     }
//     cout<<endl;
//     
  
//Wn_num scan    
    std::cout<<" here"<<std::endl;
    int N=1;//number of particles;
    ofstream file;
    dvec res;
    double lag=0.0;
    double shiftx=0.35;
    double shifty=0.999;
    double alpha=0.94;
    double komega_kv=1.84;
    double phi0 = atan(M_SQRT1_2) * 180 / M_PI;
    Eigen::Vector3d mom (1., tan((12+phi0)*M_PI/180.),0.0);
    mom.normalize();
    std::cout<<" here"<<std::endl;
    stringstream fname;
    fname<<N<<"Wn_alpha_"<<alpha<<"_komega_kv_"<<komega_kv<<".txt";
    file.open(fname.str().c_str());
    double omega=1.05;
    bool stop=true;
    x=Eigen::VectorXd::Zero(12*N);dxdt=Eigen::VectorXd::Zero(12*N);
    
//     for(size_t i=0; i<N; i++){
//       x[6*i]=pos[i][0];
//       x[6*i+1]=pos[i][1];
//       x[6*i+2]=angle[i];
//     }
    if(file.is_open()){
     // for(size_t i=0; omega>-1.9999&&stop; i++)
      {
        vector <Eigen::Vector3d> pos, angle;
       // omega-=0.1;
        if(N>1){
          dvec shift(2*(N-1));
          for(size_t i=1; i<N; i++){
            shift[2*(i-1)]=shiftx*i;
            shift[2*(i-1)+1]=shifty*i;
          }
          shiftx=shift[2*(N-2)]/2;
          shifty=shift[2*(N-2)+1]/2;
          cout<<shiftx<<" "<<shifty<<" "<<endl;
          pos.push_back(Eigen::Vector3d(-shiftx,-shifty,0.0));
          angle.push_back(Eigen::Vector3d(0.0,0.0,0.0));
          for(size_t i=1; i<N; i++){
            pos.push_back(Eigen::Vector3d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1],0.0));
            angle.push_back(Eigen::Vector3d(0.0,0.0,0.0));
          }
        }
        else{
          pos.push_back(Eigen::Vector3d(0.0,0.0,0.0));
          angle.push_back(Eigen::Vector3d(0.2,0.1,0.3));
        }
        for(size_t i=0; i<N; i++){
          x[12*i]=pos[i][0];
          x[12*i+1]=pos[i][1];
          x[12*i+2]=pos[i][2];
          x[12*i+3]=angle[i][0];
          x[12*i+4]=angle[i][1];
          x[12*i+5]=angle[i][2];
        }
        
    //     for(size_t i=0; i<N; i++){
    //       pos.push_back(Eigen::Vector2d(0.28*i,0.999*i));
    //       angle.push_back(0.0);
    //     }
        //    Eigen::Vector2d mom (0.,1.);
        Integrator *sys;
        double T=2*M_PI/fabs(omega);
          std::cout<<" here"<<std::endl;
        sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
          std::cout<<"initalizeing here1"<<std::endl;
        //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
        //sys=new Integrator();
        
        cout<<x<<endl;
        cout<<"Start ptopagation"<<endl;
        vector <dvec> data;data.resize(12*N+6);
        stop=propagate(x, dxdt, data, sys, 0.01,1000*T/*max(100*T,270.)*/,1e-13, 1e-13);
        res=analyze_data( data);
        print_data_sigle_file(0.01,1000*T, data,"test");
    //    sys->fill_steric_foreces_and_torques();
//        cout <<"finished propagation x[7]="<<x[7]<<endl;
        //cout<< x<<endl;
//         for(size_t i=1; i<N; i++){
//           cout<<(x[6*i]-x[0])*cos(x[6*i+2])+(x[6*i+1]-x[1])*sin(x[6*i+2])<<" "<<-(x[6*i]-x[0])*sin(x[6*i+2])+(x[6*i+1]-x[1])*cos(x[6*i+2])<<endl;
//         }
//         if(shift[0]>1||shift[0]<-1||shift[1]>1||shift[1]<0)
//         {
//           stop=false;
//         }
//        stringstream fname;
//        fname<<"x"<<i;
//         if(i>40)
//         {
//           print_data_sigle_file(0.01*T,100000*T, data, fname.str());
//         }
        if(stop){
          file<<omega;
          for(int i=0; i<res.size(); i++)
          {
            file<<" "<<res[i];
          }
          file<<endl;
        }
        delete sys;
      }
      
      file.close();
    }
    
        
// //omega crit solid;
//   int N=3;//number of particles;
//   double komega_kv=1.84;
//   stringstream fname;
//   bool posangle=false;
//   if(posangle){
//     fname<<N<<"omega_crit_solid_komega_kv_"<<komega_kv<<"_pos1.txt";
//   }
//   else{
//     fname<<N<<"omega_crit_solid_komega_kv_"<<komega_kv<<"_neg3.txt";
//   }
//   ofstream file;
//   file.open(fname.str().c_str());
//   //for(double alpha=0.09; alpha<5; alpha+=0.05){
//   for(double alpha=0.94; alpha<0.94001; alpha+=0.01){
//     
//     dvec res;
//     double lag=0.0;
//     double shiftx=0.35;
//     double shifty=0.999;
//     //double alpha=0.94;
//     double phi0 = atan(M_SQRT1_2) * 180 / M_PI;
//     Eigen::Vector2d mom (1., tan((12+phi0)*M_PI/180.));
//     mom.normalize();
//     double omega;
//     double omega_min;
//     if(posangle){
//       omega_min=0.6;
//     }
//     else{
//       omega_min=-0.6;
//     }
//     omega=omega_min;
//     bool stop=true;
//     x=Eigen::VectorXd::Zero(6*N);dxdt=Eigen::VectorXd::Zero(6*N);
//     dvec shift(2*(N-1));
//     for(size_t i=1; i<N; i++){
//       shift[2*(i-1)]=shiftx*i;
//       shift[2*(i-1)+1]=shifty*i;
//     }
//     vector <Eigen::Vector2d> pos;
//     vector <double> angle;
//     pos.push_back(Eigen::Vector2d(-shiftx,-shifty));
//     angle.push_back(0.0);
//     for(size_t i=1; i<N; i++){
//       pos.push_back(Eigen::Vector2d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1]));
//       angle.push_back(0.0);
//     }
//     double T;
//     Integrator *sys;
//     vector <dvec> data;
//     while(true)
//     {
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear(); data.resize(6*N+6);
//       stop=propagate(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       res=analyze_data( data);
//       if(fabs(res[4]+1)<1e-13){
//         omega_min=omega;
//         break;
//       }
//       else{
//         omega/=2;
//       }
//     }
//     double omega_max;
//     if(posangle){
//       omega_max=2;
//     }
//     else{
//       omega_max=-2;
//     }
//     omega=omega_max;
//     while(true)
//     {
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear(); data.resize(6*N+6);
//       stop=propagate(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       res=analyze_data( data);
//       if(fabs(res[4]+1)<1e-13){
//         omega*=2;
//       }
//       else{
//         omega_max=omega;
//         break;
//       }
//     }
// //     //kink
// //     shift[0]=shiftx;
// //     shift[1]=shifty;
// //     for(size_t i=2; i<N/2+1; i++){
// //       shift[2*(i-1)]=shift[2*(i-2)]+shifty;
// //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shiftx;
// //     }
// //     for(size_t i=N/2+1; i<N; i++){
// //       shift[2*(i-1)]=shift[2*(i-2)]+shiftx;
// //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shifty;
// //     }
// //     shiftx=shift[2*(N-2)]/2;
// //     shifty=shift[2*(N-2)+1]/2;
//     cout<<"omega_min="<<omega_min<<" omega_max="<<omega_max<<endl;
//     for(size_t i=0; fabs(omega_min-omega_max)>1e-5&&stop; i++){
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       omega=(omega_max+omega_min)*0.5;
//       cout<<"!!!!!!!!!!!omega="<<omega<<endl;
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear();data.resize(6*N+6);
//       stop=propagate(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       res=analyze_data( data);
//       if (fabs(res[4]+1)<1e-13){
//         omega_min=omega;        
//       }
//       else{
//         omega_max=omega;
//       }
//       
//   //    sys->fill_steric_foreces_and_torques();
//       delete sys;
//     }
//     cout.precision(14);
//     cout <<"critical omega="<<(omega_max+omega_min)*0.5<<endl;;
//     file<<alpha<<" "<<omega<<endl;
// }
// file.close();    
//     
    
// //omega crit solid nuber of particles;
//   //int N=2;//number of particles;
//   double komega_kv=1.84;
//   double alpha=0.94;
//   stringstream fname;
//   bool posangle=false;
//   if(posangle){
//     fname<<"var_N_omega_crit_solid_alpha_"<<alpha<<"_komega_kv_"<<komega_kv<<"_pos1.txt";
//   }
//   else{
//     fname<<"var_N_omega_crit_solid_alpha_"<<alpha<<"_komega_kv_"<<komega_kv<<"_neg1.txt";
//   }
//   ofstream file;
//   file.open(fname.str().c_str());
//   //for(double alpha=0.09; alpha<5; alpha+=0.05){
//   for(int N=4; N<20; N++){
//     dvec res;
//     double lag=0.0;
//     double shiftx=0.35;
//     double shifty=0.999;
//     //double alpha=0.94;
//     double phi0 = atan(M_SQRT1_2) * 180 / M_PI;
//     Eigen::Vector2d mom (1., tan((12+phi0)*M_PI/180.));
//     mom.normalize();
//     double omega;
//     double omega_min;
//     if(posangle){
//       omega_min=0.2;
//     }
//     else{
//       omega_min=-0.2;
//     }
//     omega=omega_min;
//     bool stop=true;
//     x=Eigen::VectorXd::Zero(6*N);dxdt=Eigen::VectorXd::Zero(6*N);
//     dvec shift(2*(N-1));
//     for(size_t i=1; i<N; i++){
//       shift[2*(i-1)]=shiftx*i;
//       shift[2*(i-1)+1]=shifty*i;
//     }
//     vector <Eigen::Vector2d> pos;
//     vector <double> angle;
//     pos.push_back(Eigen::Vector2d(-shiftx,-shifty));
//     angle.push_back(0.0);
//     for(size_t i=1; i<N; i++){
//       pos.push_back(Eigen::Vector2d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1]));
//       angle.push_back(0.0);
//     }
//     double T;
//     Integrator *sys;
//     vector <dvec> data;
//     while(true)
//     {
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear(); data.resize(6*N+6);
//       stop=propagate_breakcheck(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       //cout<<stop<<" "<<false<<endl;
//       //abort();
//       if(stop){
//         omega_min=omega;
//         break;
//       }
//       else{
//         omega/=2;
//       }
//     }
//     cout<<"!!!!!!!!!!!!Found-lower bound"<<endl;
//     double omega_max;
//     if(posangle){
//       omega_max=0.3;
//     }
//     else{
//       if(N==4){
//         omega_max=-0.333;
//       }
//       else
//       {
//         omega_max=-0.28;
//       }
//     }
//     omega=omega_max;
//     while(true)
//     {
//       cout<<"omega="<<omega<<endl;
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear(); data.resize(6*N+6);
//       stop=propagate_breakcheck(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       cout<<stop<<" "<<false<<endl;
//       //abort();
//       if(stop){
//         omega*=1.05;
//       }
//       else{
//         omega_max=omega;
//         break;
//       }
//     }
// //     //kink
// //     shift[0]=shiftx;
// //     shift[1]=shifty;
// //     for(size_t i=2; i<N/2+1; i++){
// //       shift[2*(i-1)]=shift[2*(i-2)]+shifty;
// //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shiftx;
// //     }
// //     for(size_t i=N/2+1; i<N; i++){
// //       shift[2*(i-1)]=shift[2*(i-2)]+shiftx;
// //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shifty;
// //     }
// //     shiftx=shift[2*(N-2)]/2;
// //     shifty=shift[2*(N-2)+1]/2;
//     cout<<"!!!!!!!!!!!!Found-upper bound"<<endl;
//     cout<<"omega_min="<<omega_min<<" omega_max="<<omega_max<<endl;
//     for(size_t i=0; fabs(omega_min-omega_max)>1e-5; i++){
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       omega=(omega_max+omega_min)*0.5;
//       cout<<"!!!!!!!!!!!omega="<<omega<<endl;
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear();data.resize(6*N+6);
//       stop=propagate_breakcheck(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-13, 1e-13);
//       if(stop){
//         omega_min=omega;        
//       }
//       else{
//         omega_max=omega;
//       }
//       
//   //    sys->fill_steric_foreces_and_torques();
//       delete sys;
//     }
//     cout.precision(14);
//     cout <<"critical omega="<<(omega_max+omega_min)*0.5<<endl;;
//     file<<N<<" "<<omega<<endl;
// }
// file.close();
    
// //omega crit back-forth;
//   int N=2;//number of particles;
//   double komega_kv=1.84;
//   stringstream fname;
//   bool posangle=false;
//   if(posangle){
//     fname<<N<<"omega_crit_back-forth_komega_kv_"<<komega_kv<<"pos1.txt";
//   }
//   else{
//     fname<<N<<"omega_crit_back-forth_komega_kv_"<<komega_kv<<"neg1.txt";
//   }
//   ofstream file;
//   file.open(fname.str().c_str());
//   for(double alpha=0.5; alpha<0.65; alpha+=0.01){
//     
//     dvec res;
//     double lag=0.0;
//     double shiftx=0.35;
//     double shifty=0.999;
//     //double alpha=0.94;
//     double phi0 = atan(M_SQRT1_2) * 180 / M_PI;
//     Eigen::Vector2d mom (1., tan((12+phi0)*M_PI/180.));
//     mom.normalize();
//     double omega;
//     double omega_min;
//       if(posangle){
//         omega_min=0.6;
//       }
//       else{
//         omega_min=-0.6;
//       }
//     omega=omega_min;
//     bool stop=true;
//     x=Eigen::VectorXd::Zero(6*N);dxdt=Eigen::VectorXd::Zero(6*N);
//     dvec shift(2*(N-1));
//     for(size_t i=1; i<N; i++){
//       shift[2*(i-1)]=shiftx*i;
//       shift[2*(i-1)+1]=shifty*i;
//     }
//     vector <Eigen::Vector2d> pos;
//     vector <double> angle;
//     pos.push_back(Eigen::Vector2d(-shiftx,-shifty));
//     angle.push_back(0.0);
//     for(size_t i=1; i<N; i++){
//       pos.push_back(Eigen::Vector2d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1]));
//       angle.push_back(0.0);
//     }
//     double T;
//     Integrator *sys;
//     vector <dvec> data;
//     while(true)
//     {
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear(); data.resize(6*N+6);
//       stop=propagate(x, dxdt, data, sys, 0.01,1000*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       res=analyze_data( data);
//       if(res[4]<1e-13){
//         omega_min=omega;
//         break;
//       }
//       else{
//         omega/=2;
//       }
//     }
//     double omega_max;
//     if(posangle){
//       omega_max=2;
//     }
//     else{
//       omega_max=-2;
//     }
//     omega=omega_max;
//     while(true)
//     {
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear(); data.resize(6*N+6);
//       stop=propagate(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       res=analyze_data( data);
//       if(res[4]<1e-13){
//         omega*=2;
//       }
//       else{
//         omega_max=omega;
//         break;
//       }
//     }
// //     //kink
// //     shift[0]=shiftx;
// //     shift[1]=shifty;
// //     for(size_t i=2; i<N/2+1; i++){
// //       shift[2*(i-1)]=shift[2*(i-2)]+shifty;
// //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shiftx;
// //     }
// //     for(size_t i=N/2+1; i<N; i++){
// //       shift[2*(i-1)]=shift[2*(i-2)]+shiftx;
// //       shift[2*(i-1)+1]=shift[2*(i-2)+1]+shifty;
// //     }
// //     shiftx=shift[2*(N-2)]/2;
// //     shifty=shift[2*(N-2)+1]/2;
//     cout<<"omega_min="<<omega_min<<" omega_max="<<omega_max<<endl;
//     for(size_t i=0; fabs(omega_min-omega_max)>1e-5&&stop; i++){
//       for(size_t i=0; i<N; i++){
//         x[6*i]=pos[i][0];
//         x[6*i+1]=pos[i][1];
//         x[6*i+2]=angle[i];
//       }
//       omega=(omega_max+omega_min)*0.5;
//       cout<<"!!!!!!!!!!!omega="<<omega<<endl;
//       T=2*M_PI/fabs(omega);
//       sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//       //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//       //sys=new Integrator();
//       //cout<<"Start ptopagation"<<endl;
//       data.clear();data.resize(6*N+6);
//       stop=propagate(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-10, 1e-10);
//       res=analyze_data( data);
//       if (res[4]<1e-13){
//         omega_min=omega;        
//       }
//       else{
//         omega_max=omega;
//       }
//       
//   //    sys->fill_steric_foreces_and_torques();
//       delete sys;
//     }
//     cout.precision(14);
//     cout <<"critical omega="<<(omega_max+omega_min)*0.5<<endl;;
//     file<<alpha<<" "<<omega<<endl;
// }
// file.close();   
    
// ////scan omega
// //  for(double alpha=0.94; alpha>0.05; alpha-=0.05)
// //  {
//     int N=6;//number of particles;
//     ofstream file;
//     dvec res;
//     double lag=0.0;
//     double shiftx=0.35579;
//     double shifty=0.984809;
//     double alpha=.94;
//     double komega_kv=1.84;
//     double phi0 = atan(M_SQRT1_2) * 180 / M_PI;
//     Eigen::Vector2d mom (1., tan((12+phi0)*M_PI/180.));
//     mom.normalize();
//     dvec shift(2*(N-1));
//     for(size_t i=1; i<N; i++){
//       shift[2*(i-1)]=shiftx*i;
//       shift[2*(i-1)+1]=shifty*i;
//     }
//     stringstream fname;
//     fname<<N<<"shift_alpha_"<<alpha<<"_komega_kv_"<<komega_kv<<".txt";
//     file.open(fname.str().c_str());
//     double omega=-0.34;
//     bool stop=true;
//     x=Eigen::VectorXd::Zero(6*N);dxdt=Eigen::VectorXd::Zero(6*N);
//     
// //     for(size_t i=0; i<N; i++){
// //       x[6*i]=pos[i][0];
// //       x[6*i+1]=pos[i][1];
// //       x[6*i+2]=angle[i];
// //     }
// //    omega=-0.71/*+0.0002*/;
//     if(file.is_open()){
//       for(size_t i=0; omega>-1.0&&stop; i++){
//         vector <Eigen::Vector2d> pos;
//         vector <double> angle;
//         omega-=0.01;
//         shiftx=shift[2*(N-2)]/2;
//         shifty=shift[2*(N-2)+1]/2;
//         cout<<shiftx<<" "<<shifty<<" "<<endl;
//         pos.push_back(Eigen::Vector2d(-shiftx,-shifty));
//         angle.push_back(0.0);
//         for(size_t i=1; i<N; i++){
//           pos.push_back(Eigen::Vector2d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1]));
//           angle.push_back(0.0);
//         }
//         for(size_t i=0; i<N; i++){
//           x[6*i]=pos[i][0];
//           x[6*i+1]=pos[i][1];
//           x[6*i+2]=angle[i];
//         }
//         
//     //     for(size_t i=0; i<N; i++){
//     //       pos.push_back(Eigen::Vector2d(0.28*i,0.999*i));
//     //       angle.push_back(0.0);
//     //     }
//         //    Eigen::Vector2d mom (0.,1.);
//         Integrator *sys;
//         double T=2*M_PI/fabs(omega);
//         sys=new Integrator(N, omega,mom, alpha, komega_kv, lag, true);
//         //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//         //sys=new Integrator();
//         
//         cout<<x<<endl;
//         cout<<"Start ptopagation"<<endl;
//         vector <dvec> data;data.resize(6*N+6);
//         stop=propagate(x, dxdt, data, sys, 0.01,/*100*T*/min(100*T,1e7),1e-14, 1e-14);
//         res=analyze_data( data);
//     //    sys->fill_steric_foreces_and_torques();
//         cout <<"finished propagation x[7]="<<x[7]<<endl;
//         cout<< x<<endl;
//         dvec shift(2*(N-1));
//         for(size_t i=1; i<N; i++){
//           cout<<(x[6*i]-x[0])*cos(x[6*i+2])+(x[6*i+1]-x[1])*sin(x[6*i+2])<<" "<<-(x[6*i]-x[0])*sin(x[6*i+2])+(x[6*i+1]-x[1])*cos(x[6*i+2])<<endl;
//           //shift[2*(i-1)]=(x[6*i]-x[0])*cos(x[6*i+2])+(x[6*i+1]-x[1])*sin(x[6*i+2]);
//          //shift[2*(i-1)+1]=-(x[6*i]-x[0])*sin(x[6*i+2])+(x[6*i+1]-x[1])*cos(x[6*i+2]);
//         }
//         for(int i=0; i<res.size(); i++){
//         cout<<" "<<res[i];
//         }
//         cout<<endl;
// //         if(i>1)
// //         {
// //           abort();
// //         }
// //         if(shift[0]>1||shift[0]<-1||shift[1]>1||shift[1]<0)
// //         {
// //           for(size_t i=1; i<N; i++){
// //           shift[2*(i-1)]=shiftx*i;
// //           shift[2*(i-1)+1]=shifty*i;
// //           }
// //         }
// //        stringstream fname;
// //        fname<<"x"<<i;
// //         if(i>40)
// //         {
// //           print_data_sigle_file(0.01*T,100000*T, data, fname.str());
// //         }
//         if(stop){
//           file<<omega;
//           file<<" "<<shift[0];file<<" "<<shift[1];
//           for(int i=2; i<shift.size(); i++)
//           {
//             file<<" "<<shift[i]-shift[i-2];
//           }
//           lag=sys->get_phi(100*T)-x[2]-std::atan2(mom[1],mom[0]);
//           file<<" "<<lag;
//           for(int i=0; i<res.size(); i++)
//           {
//             file<<" "<<res[i];
//           }
//           file<<endl;
//         }
// //         else{
// //           cout <<"Why!!!!!!"<<endl;
// //           abort();
// //         }
//         delete sys;
//       }
//       lag=0.0;
//       shiftx=0.35;
//       shifty=0.999;
//       for(size_t i=1; i<N; i++){
//         shift[2*(i-1)]=shiftx*i;
//         shift[2*(i-1)+1]=shifty*i;
//       }
//       omega=0.00;
//       stop=true;
//       for(size_t i=0; omega<0.9999&&stop; i++){
//         vector <Eigen::Vector2d> pos;
//         vector <double> angle;
//         omega+=0.01;
//         shiftx=0.0;
//         shifty=0.0;
//         for(size_t i=1; i<N; i++){
//           shiftx+=shift[2*(i-1)];
//           shifty+=shift[2*(i-1)+1];
//         }
//         shiftx/=2;
//         shifty/=2;
//         pos.push_back(Eigen::Vector2d(-shiftx,-shifty));
//         angle.push_back(0.0);
//         for(size_t i=1; i<N; i++){
//           pos.push_back(Eigen::Vector2d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1]));
//           angle.push_back(0.0);
//         }
//         for(size_t i=0; i<N; i++){
//           x[6*i]=pos[i][0];
//           x[6*i+1]=pos[i][1];
//           x[6*i+2]=angle[i];
//         }
//         
//     //     for(size_t i=0; i<N; i++){
//     //       pos.push_back(Eigen::Vector2d(0.28*i,0.999*i));
//     //       angle.push_back(0.0);
//     //     }
//         //    Eigen::Vector2d mom (0.,1.);
//         Integrator *sys;
//         double T=2*M_PI/fabs(omega);
//         sys=new Integrator(N, omega,mom, alpha, komega_kv,lag, true);
//         //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//         //sys=new Integrator();
//         
//         cout<<x<<endl;
//         cout<<"Start ptopagation"<<endl;
//         vector <dvec> data;data.resize(6*N+6);
//         stop=propagate(x, dxdt, data, sys, 0.01,100*T/*max(100*T,270.)*/,1e-14, 1e-14);
//         res=analyze_data( data);
//     //    sys->fill_steric_foreces_and_torques();
//         cout <<"finished propagation x[7]="<<x[7]<<endl;
//         cout<< x<<endl;
//         dvec shift(2*(N-1));
//         for(size_t i=1; i<N; i++){
//           cout<<(x[6*i]-x[0])*cos(x[6*i+2])+(x[6*i+1]-x[1])*sin(x[6*i+2])<<" "<<-(x[6*i]-x[0])*sin(x[6*i+2])+(x[6*i+1]-x[1])*cos(x[6*i+2])<<endl;
//           //shift[2*(i-1)]=(x[6*i]-x[0])*cos(x[6*i+2])+(x[6*i+1]-x[1])*sin(x[6*i+2]);
//           //shift[2*(i-1)+1]=-(x[6*i]-x[0])*sin(x[6*i+2])+(x[6*i+1]-x[1])*cos(x[6*i+2]);
//         }
//         if(shift[0]>1||shift[0]<-1||shift[1]>1||shift[1]<0)
//         {
//           for(size_t i=1; i<N; i++){
//           shift[2*(i-1)]=shiftx*i;
//           shift[2*(i-1)+1]=shifty*i;
//           }
//         }
// //        stringstream fname;
// //        fname<<"x"<<i;
// //         if(i>40)
// //         {
// //           print_data_sigle_file(0.01*T,100000*T, data, fname.str());
// //         }
//         if(stop){
//           file<<omega;
//           file<<" "<<shift[0];file<<" "<<shift[1];
//           for(int i=2; i<shift.size(); i++){
//             file<<" "<<shift[i]-shift[i-2];
//           }
//           lag=sys->get_phi(100*T)-x[2]-std::atan2(mom[1],mom[0]);
//           file<<" "<<lag;
//           for(int i=0; i<res.size(); i++){
//             file<<" "<<res[i];
//           }
//           file<<endl;
//         }
//         delete sys;
//       }
//       
//       file.close();
//     }
// //  }
    
}






