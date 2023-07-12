#include "simulation.h"


void read_file(string fname, int N, Eigen::VectorXd &x)
{
  cout<<"Start read file:"<<fname<<endl;
  x.resize(3*N);
  string line;
  ifstream myfile (fname.c_str());
  double junk, junk1;
  double x1,y1,z1;
  if (myfile.is_open()){
    for(size_t i=0; i<N; i++){
      //cout<<i<<endl;
      myfile>>x1>>y1>>z1>>junk>>junk1;  
      x[3*i]=x1;x[3*i+1]=y1;x[3*i+2]=z1; 
      //cout <<x1<<" "<<y1<<" "<<z1<<endl;
    }
    myfile.close();
  }
  else{
    cout<<"can not open file:"<<fname<<endl;
    abort();
  }
  cout<<"End read file:"<<fname<<endl;
}

int main()
{
  time_t start,end;
  std::time (&start);
  double diff;
  int i;
  double t_max=100000;
  double extra_t=1000;
  int iter_per_step=10;
  simulation *sim;
  double dt=1e-6;
  double t=0;
  double y=0.; // degree
  double x=0.;
//   double Cm=-100;
  double Cm=20;
  double H1=0.0;
  double om=200*16;
  double torq=150;
  Eigen::VectorXd x_old;
  cout.precision(15);
  //vector<int> sizes={11,12,12,13,14,15,16,17,18,19,20,30,40,50,60,70, 80,90,100,150,200,300,400,500,600,700,800,900,1000};
  //vector<int> sizes={25};
   dvec torqs{100};
//    for(int i=800; i<801; i++){torqs.push_back(i);
//    }
  //dvec torqs={150};
  //dvec torqs={200,190,180,170,160};
  //vector<int> sizes={50};
  //double H1=2.5;
  //vector<double> cm={-10, -20, -30,-40,-50,-60,-70,-80,-90,-110,-120,-130,-140,-150,-160,-170,-180,-190,-200  }; 
  //vector<double> torqx={0.9, 0.95, 1, 1.05, 1.25, 1.5, 1.75, 2., 3. };
  vector<double> cm={-50}; 
  //vector<double> torqx={1.5};
  vector<double> torqx1={10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,130,140,150,160,170,180,190,200,210,220,230, 240, 250, 260};
  //vector<double> torqx1={70,60,50,40,30};
  
  
  vector<std::pair<double, double> > sims;
  for(int i=0; i<cm.size(); i++){
    for(int j=0; j<torqx1.size(); j++){
      //sims.push_back(std::make_pair(cm[i],fabs(cm[i])*torqx[j]));
      sims.push_back(std::make_pair(cm[i],torqx1[j]));
    }
  }
   
  string data;
  ofstream file;
  file.open("helix_scan_cm_-50_N_100.dat");
  if(file.is_open()){
    //for(size_t i=0; i<sizes.size(); i++){
    for(size_t i=0; i< sims.size(); i++){  
      //int size=sizes[i];
      int size=100;
      torq=sims[i].second;
      Cm=sims[i].first;
      if(i>1000){
        sim=new simulation(size,Cm, H1, om, torq,x_old);
      }
      else{
        //read_file("../tesrt/data_filament_VSIMEX_BDF3_DD_Ns_300_omega_3200_Cm_-100_H1_1.8_Torq_200_tol_1e-08/shape.txt",size+1, x_old);
        //sim=new simulation(size,Cm, H1, om, torq,x_old);
        sim=new simulation(size,Cm, H1, om, torq);
      }
  //sim->propagate(1e-12, 0.01, 1e-8,1e-8);
  //sim->propagate_CN(dt, 0.01);
  //sim->propagate_CN_proper(dt, 0.01);
  //sim->propagate_Euler(dt/10, 0.01);
  //sim->propagate_midp_impl_expl(1e-9, 1e-3);
  //sim->propagate_Euler1_impl_expl(1e-7, 1e-3);
  //sim->propagate_Euler_impl_expl(1e-7, 1e-3);
  //sim->propagate_Euler_impl_expl(dt/**5*//10*2.5, 5e-3);
  //sim->propagate_VSIMEX_BDF2(dt, 5e-3, 1e-10);
 // sim->propagate_IMEX_BDF2(dt/10*2.5, 5e-3);
  //sim->propagate_IMEX_BDF3(dt/10*2.5, 5e-3);
      sim->propagate_VSIMEX_BDF3_DD(dt/10*2.5, 1000*1e-3, 1e-7);
      //sim->propagate_VSIMEX_BDF3(dt/10*2.5, 100*1e-3, 1e-8);
      //sim->propagate_IMEX_BDF3(1e-6, 10e-3);
      //sim->propagate_VSIMEX_BDF4(dt/10*2.5, 100*1e-3);
      data=sim->anal_data1();
      x_old=sim->get_x();
      //file<<sizes[i]<<" "<<data<<std::flush;;
      file<<Cm<<" "<<torq<<" "<<data<<std::flush;
  //sim->propagate_ARS_232(dt/500, 5e-3);
  //sim->propagate_IMEX_TVB3(dt/10*2, 5e-3);
//was  //sim->propagate_IMEX_BDF4(1e-9, 5e-3);
  //sim->propagate_VSIMEX_BDF4(dt/17, 5e-5);
  //sim->propagate_VSIMEX_BDF4_adaptive(dt/20, 1e-3);
  //sim->propagate_IMEX_TVB5(dt/17, 5e-3);
  //sim->propagate_IMEX_BDF5(dt/70, 5e-3);
  //sim->propagate_IMEX_TVB5_adaptive(dt, 1e-3);
  //sim->propagate_IMEX_BDF5_adaptive(dt, 1e-3);
  //sim->propagate_IMEX_BDF4_adaptive(dt, 1e-3);
  //sim->propagate_IMEX_BDF3_adaptive(dt, 1e-3);
  //sim->propagate_IMEX_TVB3(dt, 1e-3);
   //sim->propagate_rk4_impl_expl(dt, 1e-3);
      delete sim;
    }
    file.close();
  }
  time (&end);
  diff = difftime (end,start);  
  cout<<"finished calculation in " <<diff<<" seconds."<<endl;
    cout << "Hello world!" << endl;
    return 0;
}
