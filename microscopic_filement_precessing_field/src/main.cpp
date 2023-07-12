
#include "integrator.h"
#include "simulation.h"




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
  //typedef runge_kutta_dopri5< Eigen::VectorXd> error_stepper_type;
  typedef runge_kutta_cash_karp54< Eigen::VectorXd> error_stepper_type;
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
  sys->append_data(x, data, t);
  //sys->debug(x, t);
  //abort();
  while(t<t_max)
  {
    
    //cout<<"t="<<t<<endl;
    if(dt>0.01)
    {
      dt=0.01;
    }
    if (dt>1e-12)
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
        sys->append_data(x, data, t);
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
      if((i+ ignored_steps) %10000==0)
      {
        cout<<i+ ignored_steps<<" " << t << '\t' << dt <<" "<<t_max<<" "<< i/5000000.0<<"\% of max\n";
        //cout<<sys->get_Mag_torque(x)<<" "<<sys->get_Steric_torque()<<" "<<sys->get_steric_force(0).transpose() <<endl;
        if(i>500000000){
          return false;
        }
      }
    }
    else
    {
      //cout<<"failing at t="<<t<<" "<<dt<<endl;
      //cout <<"!!!!!!!!!!!"<<x.transpose()<<endl;
      //cs.try_step(*sys, x, t,dt);
      //sys->debug(x, t);
      //abort();
      dt=1e-10;
      rk.do_step(*sys, x, t,dt);
      t+=dt;
      limiter_steps++;
      i++;
      if(limiter_steps>100000)
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
//  sys->test_do_something();
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}


bool propagate_rk4(Eigen::VectorXd &x, Eigen::VectorXd &dxdt, vector<dvec>& data,  Integrator *sys, double dt1,double t0, double tmax1)
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
  //typedef runge_kutta_dopri5< Eigen::VectorXd> error_stepper_type;
  typedef runge_kutta_cash_karp54< Eigen::VectorXd> error_stepper_type;
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

  now= std::chrono::system_clock::now();
  sys->append_data(x, data, t);
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
      rk.do_step(*sys, x, t,dt);
      t+=dt;
      acpeted_steps++;
      sys->append_data(x, data, t);
      //sys->debug(x, t);
      //abort();
      appended_t=t;
      i++;
      //if(i%100000==0)
      if((i+ ignored_steps) %10000==0)
      {
        cout<<i+ ignored_steps<<" " << t << '\t' << dt <<" "<<t_max<<" "<< i/5000000.0<<"\% of max\n";
        //cout<<sys->get_Mag_torque(x)<<" "<<sys->get_Steric_torque()<<" "<<sys->get_steric_force(0).transpose() <<endl;
        if(i>500000000){
          return false;
        }
      }
    }
  }
  cout << t << '\t' << dt << '\n';
  foo= std::chrono::system_clock::now();
  diff1 = foo - now;
  cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
  cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
//  sys->test_do_something();
//  print_pars();
//   print_data(10.0,t_max);
//  print_data_sigle_file(0.1,t_max);
  return true;
}

// bool propagate_silent(Eigen::VectorXd &x, Eigen::VectorXd &dxdt,  Integrator *sys, double dt1,double t0, double tmax1, double abs_tol,double rel_tol )
// {
//   cout <<"Startng propagation! x.size()="<<x.size()<<" "<<tmax1<<" "<<dt1<<" " <<endl;
//   Eigen::Vector3d mom;
//   Eigen::Quaterniond quat;
// //   ostringstream os1;
// //   os1<<"1data_Np_"<<Np<<"_gamma_"<<gamma<<"_lambda_"<<lambda;
// //   folder=os1.str();
// //   os1.str("");
// //   os1<<"mkdir "<<folder;
// //   system(os1.str().c_str());
//   double dt=dt1; double t_max=tmax1;
//   //typedef runge_kutta_fehlberg78<Eigen::VectorXd> error_stepper_type;
//   typedef runge_kutta_dopri5< Eigen::VectorXd> error_stepper_type;
//   //typedef runge_kutta_cash_karp54< Eigen::VectorXd> error_stepper_type;
//   typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
//   error_stepper_type rk;
//   std::chrono::time_point<std::chrono::system_clock> foo;
//   std::chrono::duration<double> diff1;
//   std::chrono::time_point<std::chrono::system_clock> now;
//   double xx=0.0;
//   double xy=0.0;
//   double ax = 1.0;
//   double adxdt = 1.0;
//   double t=t0; double appended_t=t;
//   size_t i=1;
//   size_t acpeted_steps=0;
//   size_t ignored_steps=0;
//   size_t limiter_steps=0;
//   size_t N=x.size()/7;
//   dvec a,b;a.resize(N-1);b.resize(N-1);
// 
//   controlled_stepper_type cs(
//           default_error_checker< double , vector_space_algebra , default_operations >( abs_tol, rel_tol, ax, adxdt) );;
//   controlled_step_result res;
//   now= std::chrono::system_clock::now();
//   //sys->append_data(x, data, t);
//   while(t<t_max)
//   {
//     if(dt>0.01)
//     {
//       dt=0.01;
//     }
//     if (dt>1e-10)
//     {
// //       if(dt>dt1)
// //       {
// //         dt=dt1;
// //       }
//       if(t+dt>t_max)
//       {
//         dt=t_max-t;
//       }
//       res= cs.try_step(*sys, x, t,dt);
//       if(res==0&&(t-appended_t>1.0e-3))
//       {
//         acpeted_steps++;
//         //sys->append_data(x, data, t);
//         appended_t=t;
//         i++;
//       }
//       else
//       {
//         ignored_steps++;
//       }
// 
//       if(i%100000==0)
//       //if(i%1==0)
//       {
//         cout << t << '\t' << dt <<" "<<t_max<<" "<< i/5000000.0<<"\% of max\n";
//         if(N>1){
//           for(size_t ix=1; ix<N; ix++){
//           a[ix-1]=(x[6*ix]-x[6*(ix-1)])*cos(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*sin(x[6*ix+2]);
//           b[ix-1]=-(x[6*ix]-x[6*(ix-1)])*sin(x[6*ix+2])+(x[6*ix+1]-x[6*(ix-1)+1])*cos(x[6*ix+2]);
//           }
//           cout<<*max_element(a.begin(), a.end())<<" "<<*min_element(a.begin(), a.end())<<" "<<*max_element(b.begin(), b.end())<<" "<<*min_element(b.begin(), b.end())<<endl;
//         }
//         //cout<<sys->get_Mag_torque(x)<<" "<<sys->get_Steric_torque()<<" "<<sys->get_steric_force(0).transpose() <<endl;
//         if(i>500000000){
//           return false;
//         }
//       }
//     }
//     else
//     {
//       dt=2e-10;
//       rk.do_step(*sys, x, t,dt);
//       limiter_steps++;
//       i++;
//       if(limiter_steps>1000)
//       {
//         cout<<"!!!!!!!!!!!More than 1000 limiter steps"<<endl;
//         return false;
//       }
//       //cout <<t<<endl;
//     }
//     
//   }
//   cout << t << '\t' << dt << '\n';
//   foo= std::chrono::system_clock::now();
//   diff1 = foo - now;
//   cout <<"Simulation took "<<diff1.count()<<" seconds"<<endl;
//   cout <<" acpeted_steps="<<acpeted_steps<<" ignored_steps="<<ignored_steps<<" limiter_steps="<<limiter_steps<<endl;
//   sys->test_do_something();
// //  print_pars();
// //   print_data(10.0,t_max);
// //  print_data_sigle_file(0.1,t_max);
//   return true;
// }
// 
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
    int N=60;//number of particles;
    ofstream file;
    dvec res;
    stringstream fname;
    
    double omega=100;
    bool stop=true;
    double beta=M_PI*62./180;
    x=Eigen::VectorXd::Zero(3*N);dxdt=Eigen::VectorXd::Zero(3*N);
    fname<<"DGravity_"<<N<<"omega_"<<omega<<".txt";
    //file.open(fname.str().c_str());
    
//     for(size_t i=0; i<N; i++){
//       x[6*i]=pos[i][0];
//       x[6*i+1]=pos[i][1];
//       x[6*i+2]=angle[i];
//     }
    //if(file.is_open())
    //{
      //for(size_t i=0; force_prefactor<2.5&&stop; i++)
      {
        vector <Eigen::Vector3d> pos;pos.resize(N);
        for(size_t i=0; i<N; i++){
          pos[i]={i*1.04958, 0.255-0.03044*i, 0.0};
        }
        for(size_t i=0; i<N-1; i++){
           cout<<"i="<<i<<" lengh="<<(pos[i+1]-pos[i]).norm()<<endl;
        }
        for(size_t i=0; i<N; i++){
          x.segment<3>(3*i)=pos[i];
        }
        
    //     for(size_t i=0; i<N; i++){
    //       pos.push_back(Eigen::Vector2d(0.28*i,0.999*i));
    //       angle.push_back(0.0);
    //     }
        //    Eigen::Vector2d mom (0.,1.);
        double T;
        double t_max=5.;
        if(fabs(omega)>0.001){
        T=2*M_PI/fabs(omega);
        }
        else{
          T=1.;
        }
        std::cout<<" here"<<std::endl;
        simulation *sim;
        sim=new simulation(N, omega,beta, true,x);
        double dt=1e-6;
        sim->propagate_VSIMEX_BDF3_DD(dt/10*2.5, 5000*1e-3, 1e-4);
//         Integrator *sys;
//         sys=new Integrator(N, omega,beta, true);
//           std::cout<<"initalizeing here1"<<std::endl;
//           //sys->print_mathematica(x,"mathematica_pos");
//           //abort();
//         //sys=new Integrator(N, 1,pos,angle, 1,omega,mom);
//         //sys=new Integrator();
//
//         cout<<x<<endl;
//         cout<<"Start ptopagation"<<endl;
//         //vector <dvec> data;data.resize(38*N+1);
//         //vector <dvec> data;data.resize((6+4*44)*N+1);
//         vector <dvec> data;data.resize(3*N+1);
//         //stop=propagate_silent(x, dxdt, sys, 0.01,0.0,t_max/*max(100*T,270.)*/,1e-12, 1e-12);
// //        for(size_t i=0; i<data.size(); i++){
// //        data[i].clear();
// //        }
// //        cout<<data.size()<<" "<<data[0].size()<<" !!!!!!!!!!!!!!"<<endl;;
//         stop=propagate(x, dxdt, data, sys, 0.01,0.0,t_max/*max(100*T,270.)*/,1e-5, 1e-5);
//         //stop=propagate_rk4(x, dxdt, data, sys, 1.0e-5,0.0,t_max/*max(100*T,270.)*/);
//         //stop=propagate_virtual(x, dxdt, data, sys, 0.01,0.0,10*T/*max(100*T,270.)*/,1e-12, 1e-12);
// //         cout<<data[0].size()<<" !!!!!!!!!!!!!!"<<endl;;
//         res=analyze_data( data);
//         ostringstream os1;
//         os1<<"filament_N_"<<N<<"_beta_"<<beta<<"_omega_"<<omega;
//         print_data_sigle_file(0.1,t_max, data,os1.str());
//         //sys->print_mathematica(x,"mathematica_pos");
//     //    sys->fill_steric_foreces_and_torques();
//         //cout<< x<<endl;
// //         for(size_t i=1; i<N; i++){
// //           cout<<(x[6*i]-x[0])*cos(x[6*i+2])+(x[6*i+1]-x[1])*sin(x[6*i+2])<<" "<<-(x[6*i]-x[0])*sin(x[6*i+2])+(x[6*i+1]-x[1])*cos(x[6*i+2])<<endl;
// //         }
// //         if(shift[0]>1||shift[0]<-1||shift[1]>1||shift[1]<0)
// //         {
// //           stop=false;
// //         }
// //        stringstream fname;
// //        fname<<"x"<<i;
// //         if(i>40)
// //         {
// //           print_data_sigle_file(0.01*T,100000*T, data, fname.str());
// //         }
//         if(stop){
// //           file<<force_prefactor;
// //           for(int i=0; i<res.size(); i++)
// //           {
// //             file<<" "<<res[i];
// //           }
// //           file<<endl;
// //           force_prefactor+=0.01;
//         }
//         delete sys;
//       }
      
      //file.close();
    }
    
    
}






