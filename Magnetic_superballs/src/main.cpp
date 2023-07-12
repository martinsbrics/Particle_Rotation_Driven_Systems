
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


bool propagate(Eigen::VectorXd &x, Eigen::VectorXd &dxdt, vector<dvec>& data,  Integrator *sys, double dt1,double t0, double tmax1, double abs_tol,double rel_tol )
{
  cout <<"Startng propagation! x.size()="<<x.size()<<" "<<tmax1<<" "<<dt1<<" " <<endl;
  Eigen::Vector3d mom;
  Eigen::Quaterniond quat;
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
  double t=t0; double appended_t=t;
  size_t i=1;
  size_t acpeted_steps=0;
  size_t ignored_steps=0;
  size_t limiter_steps=0;
  size_t N=x.size()/7;
  dvec a,b;a.resize(N-1);b.resize(N-1);

  controlled_stepper_type cs(
          default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
  controlled_step_result res;
  now= std::chrono::system_clock::now();
  sys->append_data_quat(x, data, t);
  //sys->debug(x, t);
  //abort();
  while(t<t_max)
  {
    //cout<<"t="<<t<<endl;
    if(dt>0.01)
    {
      dt=0.01;
    }
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
      if(res==0&&(t-appended_t>1.0e-3))
      {
        acpeted_steps++;
        sys->append_data_quat(x, data, t);
        //sys->debug(x, t);
        //abort();
        appended_t=t;
        i++;
      }
      else
      {
        ignored_steps++;
      }

      //if(i%100000==0)
      if(i%100000==0)
      //if(i%1==0)
      {
        cout << t << '\t' << dt <<" "<<t_max<<" "<< i/5000000.0<<"\% of max\n";
        if(N>1){
          for(size_t ix=1; ix<N; ix++){
          a[ix-1]=(x[6*ix]-x[6*(ix-1)])*cos(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*sin(x[6*ix+2]);
          b[ix-1]=-(x[6*ix]-x[6*(ix-1)])*sin(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*cos(x[6*ix+2]);
          }
          Eigen::Vector3d avec={x[7]-x[0],x[8]-x[1],x[9]-x[2]};
          //cout<<*max_element(a.begin(), a.end())<<" "<<*min_element(a.begin(), a.end())<<" "<<*max_element(b.begin(), b.end())<<" "<<*min_element(b.begin(), b.end())<<endl;
          cout<<avec.transpose()<<" "<<avec.norm()<<endl;
        }
        //cout<<sys->get_Mag_torque(x)<<" "<<sys->get_Steric_torque()<<" "<<sys->get_steric_force(0).transpose() <<endl;
        if(i>500000000){
          return false;
        }
      }
    }
    else
    {
      cout<<"failing at t="<<t<<" "<<dt<<endl;
      cout <<"!!!!!!!!!!!"<<x.transpose()<<endl;
      //cs.try_step(*sys, x, t,dt);
      sys->debug(x, t);
      abort();
      dt=2e-8;
      rk.do_step(*sys, x, t,dt);
      limiter_steps++;
      i++;
      if(limiter_steps>1000)
      {
        cout<<"!!!!!!!!!!!More than 1000 limiter steps at t="<<t<<endl;
        return false;
      }
      //cout <<t<<endl;
    }
    
  }
  cout << t << '\t' << dt << '\n';
  foo= std::chrono::system_clock::now();
  diff1 = foo - now;
  cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
  cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
  sys->test_do_something();
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}

bool propagate_virtual(Eigen::VectorXd &x, Eigen::VectorXd &dxdt, vector<dvec>& data,  Integrator *sys, double dt1,double t0, double tmax1, double abs_tol,double rel_tol )
{
  cout <<"Startng propagation! x.size()="<<x.size()<<" "<<tmax1<<" "<<dt1<<" " <<endl;
  Eigen::Vector3d mom;
  Eigen::Quaterniond quat;
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
  double t=t0; double appended_t=t;
  size_t i=1;
  size_t acpeted_steps=0;
  size_t ignored_steps=0;
  size_t limiter_steps=0;
  size_t N=x.size()/7;
  dvec a,b;a.resize(N-1);b.resize(N-1);

  controlled_stepper_type cs(
          default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
  controlled_step_result res;
  now= std::chrono::system_clock::now();
  sys->append_data_virtual(x, data, t);
  while(t<t_max)
  {
    if(dt>0.01)
    {
      dt=0.01;
    }
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
      if(res==0&&(t-appended_t>1.0e-3))
      {
        acpeted_steps++;
        sys->append_data_virtual(x, data, t);
        //sys->debug(x, t);
        //abort();
        appended_t=t;
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
        cout<<"!!!!!!!!!!!More than 1000 limiter steps"<<endl;
        return false;
      }
      //cout <<t<<endl;
    }
    
  }
  cout << t << '\t' << dt << '\n';
  foo= std::chrono::system_clock::now();
  diff1 = foo - now;
  cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
  cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
  sys->test_do_something();
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}


bool propagate_silent(Eigen::VectorXd &x, Eigen::VectorXd &dxdt,  Integrator *sys, double dt1,double t0, double tmax1, double abs_tol,double rel_tol )
{
  cout <<"Startng propagation! x.size()="<<x.size()<<" "<<tmax1<<" "<<dt1<<" " <<endl;
  Eigen::Vector3d mom;
  Eigen::Quaterniond quat;
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
  double t=t0; double appended_t=t;
  size_t i=1;
  size_t acpeted_steps=0;
  size_t ignored_steps=0;
  size_t limiter_steps=0;
  size_t N=x.size()/7;
  dvec a,b;a.resize(N-1);b.resize(N-1);

  controlled_stepper_type cs(
          default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
  controlled_step_result res;
  now= std::chrono::system_clock::now();
  //sys->append_data(x, data, t);
  while(t<t_max)
  {
    if(dt>0.01)
    {
      dt=0.01;
    }
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
      if(res==0&&(t-appended_t>1.0e-3))
      {
        acpeted_steps++;
        //sys->append_data(x, data, t);
        appended_t=t;
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
        cout<<"!!!!!!!!!!!More than 1000 limiter steps"<<endl;
        return false;
      }
      //cout <<t<<endl;
    }
    
  }
  cout << t << '\t' << dt << '\n';
  foo= std::chrono::system_clock::now();
  diff1 = foo - now;
  cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
  cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
  sys->test_do_something();
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}

void print_data_sigle_file( double step, double tmax1, vector<dvec>& data, string folder )
{
  cout << "Start print_data single file "<<step<<" "<<tmax1<<" "<<data[0].size()<<endl;
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
        fprintf(file," %.4e",data[ii+1][j]);
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
          fprintf(file," %.4e",ay*data[ii+1][k]+ax*data[ii+1][k-1]);
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
  dvec res(10);
  size_t size=data[0].size()-1;
  size_t Ndata=data.size()-3;
  size_t posphi=data.size()-1;
  double total2pi =0;
  double twopi=2*M_PI;
  for(size_t i=0;i<size; i++){
    //cout <<i<<endl;
    if(fabs(data[posphi][i+1]-data[posphi][i]+total2pi)>fabs(data[posphi][i+1]-data[posphi][i]+total2pi+twopi))
    {
      total2pi+=twopi;
    }
    else{
      if(fabs(data[posphi][i+1]-data[posphi][i]+total2pi)>fabs(data[posphi][i+1]-data[posphi][i]+total2pi-twopi)){
        total2pi+=twopi;
      }
    }
    data[posphi][i+1]+=total2pi;
  }

  auto mnmx = minmax_element(data[Ndata].begin(), data[Ndata].end());
  double min_a=*mnmx.first;
  double max_a=*mnmx.second;
  double aver_a=0.5*(min_a+max_a);
  Ndata=data.size()-2;
  mnmx = minmax_element(data[Ndata].begin(), data[Ndata].end());
  double min_alpha=*mnmx.first;
  double max_alpha=*mnmx.second;
  double aver_alpha=0.5*(min_alpha+max_alpha);
  Ndata=data.size()-1;
  mnmx = minmax_element(data[Ndata].begin(), data[Ndata].end());
  double min_beta=*mnmx.first;
  double max_beta=*mnmx.second;
  double aver_beta=0.5*(min_beta+max_beta);
  double not_conv;
  if(max_a-min_a>1e-5){
    not_conv=-1.0;
  }
  else{
    not_conv=1.0;
  }
  res[0]=not_conv;
  res[1]=aver_a;
  res[2]=min_a;
  res[3]=max_a;
  res[4]=aver_alpha;
  res[5]=min_alpha;
  res[6]=max_alpha;
  res[7]=aver_beta;
  res[8]=min_beta;
  res[9]=max_beta;
  
  cout<<"!!!!!!!"<<endl;
  cout<<res[0]<<" "<<res[1]<<" "<<res[2]<<" "<<res[3]<< endl;
  
  return(res);
}

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    Eigen::VectorXd x; Eigen::VectorXd dxdt;
    
    
    std::cout<<" here"<<std::endl;
    int N=2;//number of particles;
    ofstream file;
    dvec res;
    double lag=0.0;
    //R22 = (0.964357, -0.350553, 0.0);
    double shiftx=0.964357;
    double shifty=-0.350553;
    double alpha=0.94;
    double komega_kv=1.84;
    Eigen::Vector3d mom (2.,0.0,0.0);
    mom.normalize();
    std::cout<<" here"<<std::endl;
    stringstream fname;
    
    double omega=-0.95;
    lag=0.0;
    double angle1=0.01;
    //double angle1=-(acos(1/omega))/(2. * M_PI);
    int Nvirtual=92;// 8 20 44 68 92
    double force_prefactor= 0.001;
    bool stop=true;
    x=Eigen::VectorXd::Zero(7*N);dxdt=Eigen::VectorXd::Zero(7*N);
    fname<<"DGravity_"<<N<<"omega_"<<omega<<".txt";
    file.open(fname.str().c_str());
    
//     for(size_t i=0; i<N; i++){
//       x[6*i]=pos[i][0];
//       x[6*i+1]=pos[i][1];
//       x[6*i+2]=angle[i];
//     }
    if(file.is_open()){
      //for(size_t i=0; force_prefactor<2.5&&stop; i++)
      {
        vector <Eigen::Vector3d> pos;
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
          //angle.push_back(Eigen::Vector3d(0.0,0.0,0.0));
          for(size_t i=1; i<N; i++){
            pos.push_back(Eigen::Vector3d(-shiftx+shift[2*(i-1)],-shifty+shift[2*(i-1)+1],0.0));
            //angle.push_back(Eigen::Vector3d(0.0,0.0,0.0));
          }
        }
        else{
          pos.push_back(Eigen::Vector3d(0.0,0.0,0.0));
          //angle.push_back(Eigen::Vector3d(0.2,0.1,0.3));
        }
        Eigen::Vector4d quat={0.,0.,0.,1.};
        for(size_t i=0; i<N; i++){
          x.segment<3>(7*i)=pos[i];
          x.segment<4>(7*i+3)=quat;
        }
        
    //     for(size_t i=0; i<N; i++){
    //       pos.push_back(Eigen::Vector2d(0.28*i,0.999*i));
    //       angle.push_back(0.0);
    //     }
        //    Eigen::Vector2d mom (0.,1.);
        Integrator *sys;
        double T;
        if(fabs(omega)>0.001){
        T=2*M_PI/fabs(omega);
        }
        else{
          T=1.;
        }
          std::cout<<" here"<<std::endl;
        sys=new Integrator(N, omega,mom, alpha, komega_kv, lag,force_prefactor, angle1, Nvirtual, true);
          std::cout<<"initalizeing here1"<<std::endl;
          //sys->print_mathematica(x,"mathematica_pos");
          //abort();
        //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
        //sys=new Integrator();
        
        cout<<x<<endl;
        cout<<"Start ptopagation"<<endl;
        //vector <dvec> data;data.resize(38*N+1);
        //vector <dvec> data;data.resize((6+4*44)*N+1);
        vector <dvec> data;data.resize(13*N+3);
        stop=propagate_silent(x, dxdt, sys, 0.01,0.0,100*T/*max(100*T,270.)*/,1e-12, 1e-12);
//        for(size_t i=0; i<data.size(); i++){
//        data[i].clear();
//        }
//        cout<<data.size()<<" "<<data[0].size()<<" !!!!!!!!!!!!!!"<<endl;;
        stop=propagate(x, dxdt, data, sys, 0.01,0,50*T/*max(100*T,270.)*/,1e-12, 1e-12);
        //stop=propagate_virtual(x, dxdt, data, sys, 0.01,0.0,10*T/*max(100*T,270.)*/,1e-12, 1e-12);
//         cout<<data[0].size()<<" !!!!!!!!!!!!!!"<<endl;;
        res=analyze_data( data);
        ostringstream os1;
        os1<<"axtest_N_"<<N<<"_Nv_"<<Nvirtual<<"_omega_"<<omega<<"_F_"<<force_prefactor<<"_alpha_"<<alpha;
        print_data_sigle_file(0.1,50*T, data,os1.str());
        //sys->print_mathematica(x,"mathematica_pos");
    //    sys->fill_steric_foreces_and_torques();
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
          file<<force_prefactor;
          for(int i=0; i<res.size(); i++)
          {
            file<<" "<<res[i];
          }
          file<<endl;
          force_prefactor+=0.01;
        }
        delete sys;
      }
      
      file.close();
    }
    
    
}






