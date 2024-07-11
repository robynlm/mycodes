#!/bin/bash

#scp config_files/sciama.ini munozr@login1.sciama.icg.port.ac.uk:code/Cactus/simfactory/mdb/machines/sciama.ini
#scp config_files/sciama.cfg munozr@login1.sciama.icg.port.ac.uk:code/Cactus/simfactory/mdb/optionlists/sciama.cfg
#scp config_files/sciama.run munozr@login1.sciama.icg.port.ac.uk:code/Cactus/simfactory/mdb/runscripts/sciama.run
#scp config_files/sciama.sub munozr@login1.sciama.icg.port.ac.uk:code/Cactus/simfactory/mdb/submitscripts/sciama.sub

#scp -r code/orig/Cactus/arrangements/ExternalLibraries/ munozr@login1.sciama.icg.port.ac.uk:code/Cactus/arrangements/

#scp -r code/Cactus/arrangements/CTThorns/CT_Dust/ munozr@login1.sciama.icg.port.ac.uk:code/Cactus/arrangements/CTThorns/
#scp -r OtherThorns/ctthorns/CT_Dust munozr@login1.sciama.icg.port.ac.uk:code/Cactus/arrangements/CTThorns/
#scp OtherThorns/ctthorns/par/eds.par munozr@login1.sciama.icg.port.ac.uk:code/Cactus/par/

#scp -r workingThorns/ICPertFLRW/ munozr@login1.sciama.icg.port.ac.uk:.

#scp workingPar/pflrw.par munozr@login1.sciama.icg.port.ac.uk:code/Cactus/par/
#scp workingPar/eds.par munozr@login1.sciama.icg.port.ac.uk:code/Cactus/par/

#scp ../OtherThorns/dust.par munozr@login1.sciama.icg.port.ac.uk:code/2019_10/Cactus/par/

#scp ../code/Cactus/par/pflrw_rhoEq.par munozr@login1.sciama.icg.port.ac.uk:code/2019_10/Cactus/par/pflrw.par

#scp -r /home/robynm/ET/pflrwcodes/workingThorns/ctthorn_m_modified/CT_Dust munozr@login1.sciama.icg.port.ac.uk:.

#mkdir /home/robynm/simulations/MPI/testmpi2
#scp munozr@login1.sciama.icg.port.ac.uk:/mnt/lustre/munozr/runs/testmpi2/output-0000/pflrw/all_iterations/testmpi2_it_000000.hdf5 /home/robynm/simulations/MPI/testmpi2/
#mkdir /home/robynm/simulations/MPI/testmpi4

#mkdir /home/robynm/simulations/pflrw_mpi3
#scp munozr@login1.sciama.icg.port.ac.uk:code/pflrwcodes/data_analysis_codes/exe/runSlice.sh data_analysis_codes/exe/
 
#SIM="pflrw_d3e2_L1821_t1_N64_LCDM"
#scp -r munozr@login1.sciama.icg.port.ac.uk:simulations/$SIM/Time_dt.csv /home/robynm/simulations/$SIM/

#mkdir /home/robynm/simulations/pflrw_d5e2_L50_z500_N16_LCDM
#mkdir /home/robynm/simulations/pflrw_d5e2_L50_z500_N16_LCDM/output-0000
#mkdir /home/robynm/simulations/pflrw_d5e2_L50_z500_N16_LCDM/output-0000/pflrw_L50
#mkdir /home/robynm/simulations/pflrw_d5e2_L50_z500_N16_LCDM/output-0000/pflrw_L50

#for FILE in admbase-curv.file_0.h5 admbase-curv.file_1.h5 admbase-curv.file_2.h5 admbase-curv.file_3.h5 admbase-curv.file_4.h5
#do
#  scp munozr@login2.sciama.icg.port.ac.uk:/mnt/lustre/munozr/runs/pflrw_d5e2_L50_z500_N16_LCDM/output-0000/pflrw_L50/$FILE ~/simulations/pflrw_d5e2_L50_z500_N16_LCDM/output-0000/pflrw_L50/$FILE
#done

#scp munozr@login2.sciama.icg.port.ac.uk:/mnt/lustre/munozr/runs/pflrw_d3e2_L1821_t1_N32_LCDM/output-0000/pflrw_L1821/ct_dust-ct_rho.average.asc ~/simulations/pflrw_d3e2_L1821_t1_N32_LCDM/output-0000/pflrw_L1821/ct_dust-ct_rho.average.asc

#$(seq -f "%06g" 0 5000 48700)

#scp munozr@login2.sciama.icg.port.ac.uk:~/mycode/Notebooks/Dust_Spin/Sim_with_GRHydro.png 

#scp munozr@login2.sciama.icg.port.ac.uk:~/code/pflrwcodes/ebweyl/Plot_EBWeyl_BianchiII_slicing_error.ipynb ~/MyCodes/ebweyl/Plot_EBWeyl_BianchiII_slicing_error.ipynb
#$(seq -f "%06g" 20100 100 24300)

#for IT in $(seq -f "%06g" 0 100 12200)
#do
#   scp munozr@login2.sciama.icg.port.ac.uk:/users/munozr/simulations/pflrw_d3e2_L1821_t1_N32_LCDM/invar_diag_FD6/invar_diag_$IT.hdf5 ~/simulations/pflrw_d3e2_L1821_t1_N32_LCDM/output-0000/pflrw_L1821/invar_diag_FD6/invar_diag_$IT.hdf5 
#done

#for IT in $(seq -f "%06g" 0 200 24300)
#do
#   scp munozr@login2.sciama.icg.port.ac.uk:/users/munozr/simulations/pflrw_d3e2_L1821_t1_N64_LCDM/invar_diag_FD6/invar_diag_$IT.hdf5 ~/simulations/pflrw_d3e2_L1821_t1_N64_LCDM/output-0000/pflrw_L1821/invar_diag_FD6/invar_diag_$IT.hdf5 
#done

#for IT in $(seq -f "%06g" 0 400 49000)
#do
#   scp munozr@login2.sciama.icg.port.ac.uk:/users/munozr/simulations/pflrw_d3e2_L1821_t1_N128_LCDM/invar_diag_FD6/invar_diag_$IT.hdf5 ~/simulations/pflrw_d3e2_L1821_t1_N128_LCDM/output-0000/pflrw_L1821/invar_diag_FD6/invar_diag_$IT.hdf5 
#done

#scp 


#for SIM in pflrw_d3e2_L1821_t1_N128_LCDM pflrw_d3e2_L1821_t1_N64_LCDM pflrw_d3e2_L1821_t1_N32_LCDM
#do
   # Uni laptop
   #HOMEPATH=/home/robynm/simulations
   # My mac
   #HOMEPATH=/Users/robynmunoz/simulations

   #mkdir $HOMEPATH/$SIM
   #mkdir $HOMEPATH/$SIM/output-0000
   #mkdir $HOMEPATH/$SIM/output-0000/$SIM
   #mkdir $HOMEPATH/$SIM/output-0000/$SIM/all_iterations
   #mkdir $HOMEPATH/$SIM/output-0000/$SIM/DataSlice
   
   #SCIAMAPATH=munozr@login2.sciama.icg.port.ac.uk:/mnt/lustre2/ET_sims/$SIM/output-0000
   #SCIAMAPATH=munozr@login2.sciama.icg.port.ac.uk:/users/munozr/simulations
   #HOMEPATH=$HOMEPATH/$SIM/output-0000
   
   #scp $SCIAMAPATH/$SIM.par $HOMEPATH/$SIM.par
   #scp $SCIAMAPATH/$SIM.out $HOMEPATH/$SIM.out
   #scp $SCIAMAPATH/$SIM.err $HOMEPATH/$SIM.err
   #scp $SCIAMAPATH/$SIM/Time_*.csv $HOMEPATH/$SIM/Time_*.csv
   #scp $SCIAMAPATH/$SIM/h5_data.csv $HOMEPATH/$SIM/h5_data.csv
   #scp $SCIAMAPATH/$SIM/constraints.csv $HOMEPATH/$SIM/constraints.csv
   #scp $SCIAMAPATH/$SIM/mass_radius_evo.csv $HOMEPATH/$SIM/mass_radius_evo.csv

   #for NODE in $(seq -f "%01g" 6 1 8)
   #do
   #   scp $SCIAMAPATH/$SIM/ct_dust-ct_rho.file_$NODE.h5 $HOMEPATH/$SIM/ct_dust-ct_rho.file_$NODE.h5
   #done

   #scp $SCIAMAPATH/$SIM/all_iterations/pflrw_d3e2_L1821_t1_N128_LCDM_it_000000.hdf5 $HOMEPATH/$SIM/all_iterations/pflrw_d3e2_L1821_t1_N128_LCDM_it_000000.hdf5
#done

for SIM in pflrw_d3e2_L1821_t1_N64_LCDM_withlapse_1plog pflrw_d3e2_L1821_t1_N32_LCDM_withlapse_1plog  pflrw_d3e2_L1821_t1_N16_LCDM_withlapse_1plog pflrw_d3e2_L1821_t1_N64_LCDM_withlapse_1plog_pfix pflrw_d3e2_L1821_t1_N32_LCDM_withlapse_1plog_pfix  pflrw_d3e2_L1821_t1_N16_LCDM_withlapse_1plog_pfix
do
   #mkdir ~/simulations/$SIM
   #mkdir ~/simulations/$SIM/output-0000
   #mkdir ~/simulations/$SIM/output-0000/$SIM
   #mkdir ~/simulations/$SIM/output-0000/$SIM/videos

   scp munozr@login2.sciama.icg.port.ac.uk:/mnt/lustre2/ET_sims/$SIM/output-0000/$SIM.par ~/simulations/$SIM/output-0000/$SIM.par
   #scp munozr@login2.sciama.icg.port.ac.uk:/mnt/lustre2/ET_sims/$SIM/output-0000/$SIM/h5_data.csv ~/simulations/$SIM/output-0000/$SIM/h5_data.csv
   #scp munozr@login2.sciama.icg.port.ac.uk:/mnt/lustre2/ET_sims/$SIM/output-0000/$SIM/constraints.csv ~/simulations/$SIM/output-0000/$SIM/constraints.csv
done





