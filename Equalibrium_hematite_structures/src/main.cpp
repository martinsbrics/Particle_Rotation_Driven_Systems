

#include "calculate.h"
#include "calc.h"

//
// Ackley function

double ackley_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    const double x = vals_inp(0);
    const double y = vals_inp(1);
    const double pi = arma::datum::pi;

    double obj_val = -20*std::exp( -0.2*std::sqrt(0.5*(x*x + y*y)) ) - std::exp( 0.5*(std::cos(2*pi*x) + std::cos(2*pi*y)) ) + 22.718282L;

    //

    return obj_val;
}

int main()
{
    // initial values:
//     arma::vec x = arma::ones(2,1) + 1.0; // (2,2)
// 
//     //
// 
//     std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
// 
//     bool success = optim::de(x,ackley_fn,nullptr);
// 
//     std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end-start;
// 
//     if (success) {
//         std::cout << "de: Ackley test completed successfully.\n"
//                   << "elapsed time: " << elapsed_seconds.count() << "s\n";
//     } else {
//         std::cout << "de: Ackley test completed unsuccessfully." << std::endl;
//     }
// 
//     arma::cout << "\nde: solution to Ackley test:\n" << x << arma::endl;
//     

// //TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
//   
//    calc *b; 
// //   ofstream file;
// //   file.open("a_scan.txt");
//   //for(double b_crit=2; b_crit<2.0000000000001; b_crit=b_crit+0.001)
//   {
//     double b_crit=1e-5;
//     double en1=0, en2=0;
//     double angle=-35/*12*/;
//     b=new calc (2, b_crit ,angle);
//   //b->minimize_energy_pso();
//   //b->minimize_energy();
//   b->minimize_energy_gsl(false);
//   b->minimize_energy_gsl_cube(1.5);
//     //b->minimize_energy_gsl_smart(false);
//     //b->minimize_energy_gsl_smart_kink(1);
//    // b->minimize_energy_gsl_superball(1.5);
//     
//     //b->minimize_energy_gsl_new();
//      //b->minimize_energy_gsl_new_kink(2);
// //     en1=b->calc_en();
// //     b->minimize_energy_gsl_new_kink(3);
// //     en2=b->calc_en();
//     //b->minimize_energy_gsl_new_kink_smart();
//     //b->minimize_energy_smart_smart_gsl(false);
//   //b->improve_energy_second_cout();
//  // b->improve_energy_first();
//   //b->minimize_energy_first();
// //   
// //     if(en2<en1)
// //     {
// //       cout<<"For angle="<<angle<<"and Bcrit="<<b_crit<<" kinks are favorable. En diff="<< en1-en2<<" "<<-(en1-en2)/en2<<endl;
// //     }
// //     else
// //     {
// //       cout<<"For angle="<<angle<<"and Bcrit="<<b_crit<<" straight chains are favorable. En diff="<< en1-en2<<" "<<(en1-en2)/en2<<endl;
// //     }
// //     
// //   file<<b_crit<<" "<<b->r[1][0]<<" "<<b->r[1][1]<<" "<<b->calc_en()<<endl;
// //   delete b;
//    }
// //   file.close();
  
  
//   calc *b; 
//   ofstream file;
//   b=new calc (3, 1.7 ,12);
//   //b->minimize_energy_pso();
//   b->minimize_energy_de_smart_smart(false);
//   //b->minimize_energy_gsl();
//  //
//   //b->improve_energy_first();
//   //b->minimize_energy_first();
//   cout <<"Energy="<<b->calc_en()<<endl;
//   //b->minimize_energy_de_smart_smart();
//   b->minimize_energy_smart_smart_gsl(false);
//   //b->improve_energy_smart_smart_gsl();
//   //swap(b->r[1][0],b->r[1][1]);swap(b->r[2][0],b->r[2][1]);
//   
//   cout <<"Energy="<<b->calc_en()<<endl;
//   delete b;  
  
// //Angle scan 3 particles
//   calc *b;
//   double ener1,ener2,ener3,ener4,a1,a2,b1,b2, Bcrit, r1, r2,r3,r4, ori, b_m_angle, r_m_angle, b_m_angle2, r_m_angle2;
//   ofstream file;file.open("Angle_scan_3particles_2d.txt"); file.precision(14);
//   for (int angle=0; angle<145; angle++)
//   {
//     b=new calc (3, 4 ,angle);
//     b->minimize_energy_smart_smart_gsl(false);;
//     r_m_angle=b->get_rangle();
//     r1=b->get_offset();
//     file<<angle<<" "<<r_m_angle<<" "<<acos(r_m_angle)<<" "<<r1 <<endl;
//   }
//   file.close();  
// 

  
/*//offset scan   
  calc *b;
  double ener1,ener2,ener3,ener4,a1,a2,b1,b2, Bcrit, r1, r2,r3,r4, ori, b_m_angle, r_m_angle, b_m_angle2, r_m_angle2;
  ofstream file;file.open("offset.txt"); file.precision(14);
  for (int np=2; np<5; np++)
  {
    b=new calc (np, 4 ,12);
    b->minimize_energy_smart_smart_gsl(false);;
    r_m_angle=b->get_rangle();
    r1=b->get_offset();
    file<<np<<" "<<r1<<" "<<r_m_angle<<" "<<acos(r_m_angle)*180/M_PI <<endl;
    delete b;
  }
  file.close();  */  
  
//     //calculate a(0.01);
//       double en1, en2, en3;
//       calc b(7,1.0,-0.5);
//       //b.minimize_energy_de_smart();
// //      b.minimize_energy_de_smart_smart(false);
//       b.minimize_energy_shit3(0);
//       en2=b.calc_en();
//       b.minimize_energy_shit3(1);
//       en1=b.calc_en();
//       b.minimize_energy_shit3(3);
//       en3=b.calc_en();
//       //b.minimize_energy_de_smart();
//        cout <<en2<<" "<<en1<<" "<<en3<<endl;
//        cout<< en1-en2<<" "<<en3-en1<<endl;
//       //b.minimize_energy_de_smart_smart(false);
//       //b.minimize_energy_pso_smart();
//       //b.minimize_energy_pso();
//       //b.minimize_energy_first();
//       //b.minimize_energy_second();
//      // b.improve_energy_second();
//       //  b.improve_energy_first();
// //        b.setB(0.1);;
// //        b.improve_energy_first();`` 
//       //b.improve_energy_second();
// //       en1=b.calc_en();
// //       b.minimize_energy_second();
// //       en2=b.calc_en();
// //       cout<< en1-en2<<endl;
//       
       
       
       
       
       
       
       
       
       
/*      b.minimize_energy_second();
      b.shift(0.5);
      ofstream file;file.open("Dipole_Check_second_12_Deg_0.9_shift.txt");file.precision(16);
      cout <<"Energy="<<b.calc_en()<<endl;
      for(int i=1; i<40; i++)
      {
        //cout <<"Energy"<<i<<"="<<b.calc_en_ncube(i)<<endl;
         file<<i<<" "<<b.calc_en_ncube(i)<<endl;
      }
      file.close()*/;
//     calc *b;
//     double Bmag=0.3;
//     ofstream file;file.open("B_scan_2_cubes_0_deg.txt");
//     for(int i=0; i<300; i++)
//     {
//       b=new calc (2, Bmag,0);
//       file<<b->get_data()<<endl;
//       delete b;
//       Bmag+=0.001;
//     }
//     file.close();

// //offset scan no field
//   calc *c;
//   ofstream file;file.open("Cube_possition_no_field.txt"); file.precision(14);
//   for (int cubes=2; cubes<20; cubes+=2)
//   {
//     c=new calc (cubes, 0.0 ,12);
//     c->minimize_energy_de_smart_smart();
//     c->improve_energy_first();
//     file<<cubes;
//     for(int i=0;i<cubes; i++)
//     {
//       file <<" "<< c->r[i][0];
//     }
//     file<<endl;
//     delete c;
//   }
//   file.close();      
      
// //offset scan no field odd
//   calc *c;
//   ofstream file;file.open("Cube_possition_no_field_odd.txt"); file.precision(14);
//   for (int cubes=3; cubes<40; cubes+=2)
//   {
//     c=new calc (cubes, 0.0 ,12);
//     c->minimize_energy_de_smart_smart();
//     c->improve_energy_first();
//     file<<cubes<<" "<<c->get_maxoffset()/*c->get_angle()*/;
//     for(int i=0;i<cubes; i++)
//     {
//       file <<" "<< c->r[i][0];
//     }
//     file<<endl;
//     delete c;
//   }
//   file.close();      
            


      
// //critical filed scan  particles
//   calc *b;
//   double ener1,ener2,ener3,ener4,a1,a2,b1,b2, Bcrit, r1, r2,r3,r4, ori, b_m_angle, r_m_angle, b_m_angle2, r_m_angle2;
//   ofstream file;file.open("Critical_field_particles.txt"); file.precision(14);
//   for (int np=3; np<20; np++)
//   {
//     b=new calc (np, 0.2 ,12);
//     b->minimize_energy_de_smart_smart();;
//     b->improve_energy_first();
//     ener1=b->calc_en();
//     r1=b->get_offset();
//     b_m_angle=b->get_angle();
//     r_m_angle=b->get_rangle();
//     delete b;
//     b=new calc (np, 0.3 ,12);
//     b->minimize_energy_de_smart_smart();;
//     b->improve_energy_first();
//     ener2=b->calc_en();
//     r2=b->get_offset();
//     b_m_angle2=b->get_angle();
//     r_m_angle2=b->get_rangle();
//     delete b;
//     b=new calc (np, 10 ,12);
//     b->minimize_energy_de_smart_smart(false);;
//     b->improve_energy_second();
//     ener3=b->calc_en();
//     r3=b->get_offset();
//     delete b;
//     b=new calc (np, 11 ,12);
//     b->minimize_energy_de_smart_smart(false);;
//     b->improve_energy_second();
//     //b->minimize_energy();
//     ener4=b->calc_en();
//     r4=b->get_offset();
//     ori=b->orintation();
//     delete b;
//     b1=10*(ener1-ener2);
//     a1=ener1+0.2*b1;
//     b2=(ener3-ener4);
//     a2=ener3+10*b2;
//     Bcrit=(a1-a2)/(b1-b2);
//     file<<np<<" "<<Bcrit<<" "<<r4<<" "<<b_m_angle<<" "<<r_m_angle<<" "<<b_m_angle2<<" "<<r_m_angle2<<" "<<a1<<" "<<b1<<" "<<a2<<" "<<b2<<" "<<r1<<" "<<r2<<" "<<r3<<" "<<a1-(ener2+0.3*b1)<<" "<<a2-(ener4+11*b2)<<" "<<ori<<" "<<ener1<<" "<<ener2<<" "<<ener3<<" "<<ener4 <<endl;
//   }
//   file.close();      
  

//Angle scan
  calc *b;
  int npart=2;
  double ener1,ener2,ener3,ener4,a1,a2,b1,b2, Bcrit, r1, r2,r3,r4, ori, beta, beta1, beta1_sb, beta2_sb, beta3_sb, beta4_sb, beta5_sb, beta6_sb;
  ofstream file;file.open("Angle_scan_2part_2d.txt"); file.precision(14);
  for (double angle=-35; angle<55; angle+=0.1)
  {
    b=new calc (npart, 0.2 ,angle);
    b->minimize_energy_smart_smart_gsl(true);
    ener1=b->calc_en();
    r1=b->get_offset();
     beta1=acos(b->get_beta())*180/M_PI;
    delete b;
    b=new calc (npart, 0.3 ,angle);
    b->minimize_energy_smart_smart_gsl(true);
    ener2=b->calc_en();
    r2=b->get_offset();
    delete b;
    b=new calc (npart, 11 ,angle);
    b->minimize_energy_smart_smart_gsl(false);
    ener3=b->calc_en();
    r3=b->get_offset();
    delete b;
    b=new calc (npart, 12 ,angle);
    b->minimize_energy_smart_smart_gsl(false);
    ener4=b->calc_en();
    r4=b->get_offset();
    ori=b->orintation();
    b->minimize_energy_gsl_cube(2.0);
    beta=acos(b->get_beta())*180/M_PI;
    b->minimize_energy_gsl_superball(2.0);
    beta1_sb=acos(b->get_beta())*180/M_PI;
    b->minimize_energy_gsl_superball(1.9);
    beta2_sb=acos(b->get_beta())*180/M_PI;
    b->minimize_energy_gsl_superball(1.8);
    beta3_sb=acos(b->get_beta())*180/M_PI;
    b->minimize_energy_gsl_superball(1.7);
    beta4_sb=acos(b->get_beta())*180/M_PI;
    b->minimize_energy_gsl_superball(1.6);
    beta5_sb=acos(b->get_beta())*180/M_PI;
    b->minimize_energy_gsl_superball(1.5);
    beta6_sb=acos(b->get_beta())*180/M_PI;
    //cout<<beta<<endl;abort();
    delete b;
    b1=10*(ener1-ener2);
    a1=ener1+0.2*b1;
    b2=(ener3-ener4);
    a2=ener3+11*b2;
    Bcrit=(a1-a2)/(b1-b2);
    file<<angle<<" "<<Bcrit<<" "<<r4<<" "<<beta<<" "<<beta1_sb<<" "<<beta2_sb<<" "<<beta3_sb<<" "<<beta4_sb<<" "<<beta5_sb<<" "<<beta6_sb<<" "<<beta1<<" "<<a1<<" "<<b1<<" "<<a2<<" "<<b2<<" "<<r1<<" "<<r2<<" "<<r3<<" "<<a1-(ener2+0.3*b1)<<" "<<a2-(ener4+2*b2)<<" "<<ori<<" "<<ener1<<" "<<ener2<<" "<<ener3<<" "<<ener4 <<endl;
  }
  file.close();
//   
      
      
// //Angle scan
//   calc *b;
//   double ener1,ener2,ener3,ener4,a1,a2,b1,b2, Bcrit, r1, r2,r3,r4, ori, beta, beta1;
//   ofstream file;file.open("Angle_scan_12deg_2d.txt"); file.precision(14);
//   for (int np=2; np<20; np++)
//   {
//     b=new calc (np, 0.2 ,12);
//     b->minimize_energy_smart_smart_gsl();
//     ener1=b->calc_en();
//     r1=b->get_offset();
//     beta1=acos(b->get_beta())*180/M_PI;
//     delete b;
//     b=new calc (np, 0.3 ,12);
//     b->minimize_energy_smart_smart_gsl();
//     ener2=b->calc_en();
//     r2=b->get_offset();
//     delete b;
//     b=new calc (np, 4 ,12);
//     b->minimize_energy_smart_smart_gsl(false);
//     ener3=b->calc_en();
//     r3=b->get_offset();
//     delete b;
//     b=new calc (np, 6 ,12);
//     b->minimize_energy_smart_smart_gsl(false);
//     ener4=b->calc_en();
//     r4=b->get_offset();
//     ori=b->orintation();
//     beta=acos(b->get_beta())*180/M_PI;
//     delete b;
//     b1=10*(ener1-ener2);
//     a1=ener1+0.2*b1;
//     b2=0.5*(ener3-ener4);
//     a2=ener3+4*b2;
//     Bcrit=(a1-a2)/(b1-b2);
//     file<<np<<" "<<Bcrit<<" "<<r4<<" "<<beta<<" "<<beta1<<" "<<a1<<" "<<b1<<" "<<a2<<" "<<b2<<" "<<r1<<" "<<r2<<" "<<r3<<" "<<a1-(ener2+0.3*b1)<<" "<<a2-(ener4+2*b2)<<" "<<ori<<" "<<ener1<<" "<<ener2<<" "<<ener3<<" "<<ener4 <<endl;
//   }
//   file.close();     



  
    return 0;
}
