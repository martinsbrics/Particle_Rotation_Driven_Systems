This is modified espresso version which allows to perform simulations with paramagnetc particles and for paramagnetic and ferromagnetic particles in precessing magnetic field.

compare to ordinary espresso user has to specity 3 aditional variables. PARAMAGNETIC_DIPOLES=false for ferromagnetic and treu for paramagnetic particles.  EXTERNAL_MAGNETIC_ROTATING_FIELD=false for static field, otherwise precessing field is used with magnitude Bmax, frequency omega, phase and pitch. The last 5 parameters of External_magnetic_Field_par  are for paramagnetic particles: magnetic_field_prefactor-> \mu_0/(4\pi), magnetization_prefactor-> dipole magnitude= magnetization_prefactor* tatal magnetic field; finite_time_relaxation, if <0, then it is assumed that finite time ralaxation is ised; finite_relaxation_time parameter ->ratio of new to the old magnetic moment dirrection 

system.integrator.External_magnetic_Field_par=(Bmax,omega, phase, pitch, magnetic_field_prefactor, magnetization_prefactor, finite_time_relaxation, finite_relaxation_time parameter ,uned_variable)
system.integrator.PARAMAGNETIC_DIPOLES=True
system.integrator.EXTERNAL_MAGNETIC_ROTATING_FIELD=True

All newly added function start with prefix martins_ and modified section of code start "//MB start"  and end with //MB end

the two scripts for simulating paramagenet and feromagnetic two particles with cubic shape are provided.
