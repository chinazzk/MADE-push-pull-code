# Slim Fast Input Script
# MADE Site


#set k      Fbm

file mkdir "madesite_slim_Fbm"
cd         "madesite_slim_Fbm"

#run id and other info
set parameter_file     "slim.inject.par"
set run_name           "Facies MADE Site Injection Portion"
set log_file           "slog.inject.txt"


#grid and domain info
set nx                  128
set ny                  128
set nz                  640
set dx                  0.390625
set dy                  0.390625
set dz                  0.016796875

#input files for head, perm, and porosity
set v_type              calc_trans
set kx_file             "../../../parflow_Fbm/MADE.out.perm_x.pfb"
set ky_file             "../../../parflow_Fbm/MADE.out.perm_y.pfb"
set kz_file             "../../../parflow_Fbm/MADE.out.perm_z.pfb"

set press              1 

set head_list_file     "../facies_inject_press2.txt"
set time_file          "../inject_time.txt"

set saturated           "no"
set vga_file            "../../../../vga_homog.pfb"
set vgn_file            "../../../../vgn_homog.pfb"
set sres_file           "../../../../sres_homog.pfb"

set phi_type             "constant"
set phi                  0.104

#particle number and numerics
set npmax               10000000
set give_up             100000
set epsilon             1e-10

#number of constituents
set num_constituents    1

#radioactive decay
set half_life(1)        0.0

#Sorption
set sorption_type(1)    "constant"
set sorption_value(1)   1.0

#Attachement
set attach_type(1)      "constant"
set attachment_value(1) 0.0
set detach_type(1)      "constant"
set detachment_value(1) 0.0

#Concentrations
set write_concentrations    "no"
set concentration_filename  "conc.Fbm"
set concentration_header(1) "Tracer"

set write_out_particles      "no"
set particle_out_file        "parttrace.txt"
set write_out_end_particles  "yes"
set particle_out_end_file    "endInjectionparticles.txt"
set write_out_moments        "no"
set moment_out_file          "moments.txt"

set temp_averaging       "no"

#mode of operation
set simulation_mode      "forward"
set part_split           "no"
set min_conc             1e-20
set fast_kin             "no"

#Dispersion and diffusion
set alpha_l              0.289
set alpha_t              0.071082
set diffusivity          3.33e-6

#Timing Info
set num_increments       31
set time_increment       1
set vel_nskip            0

#Wells
set number_wells             0
set write_well_breakthrough  "yes"
set well_overwrite           "yes"
set well_out_file(1)         "well.injection.txt"

set well_breakthrough_dt 1.0
set number_well_steps    31

set well_x_location(1)    25.1953125
set well_y_location(1)    25.1953125
set well_screen_top(1)    8.5
set well_screen_bottom(1) 2.4

set well_pumping_rate(1)  0.00

#Plane Information
set number_planes        0
set plane_out_file(1)   "noplane.txt"

#Well IC
set slim_initial_condition_type(1) "well"
set well_ic_x(1)             25.1953125
set well_ic_y(1)             25.1953125
set well_ic_bot(1)           2.4
set well_ic_top(1)           8.5
set well_ic_pump_rate(1)     0.00
set well_ic_num_part_per(1)  3
set well_ic_mass_per_part(1) 0.03072857
set well_ic_rad(1)           0.2
# 44649 total particles
# Total injected mass = 1372 g
# 41 particle injection time steps (4.1 hours) hardwired in slim.exe
#tfADE
set immobile_pause "no"

#Slim Output
source $env(SLIM_DIR)/tcl/slim.run.tcl

cd ".."