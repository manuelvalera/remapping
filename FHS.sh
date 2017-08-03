#!/bin/tcsh
set sourceMods = "/home/$USER/remapping/SourceMods2/"
set caze       = "fhs-spinup"
setenv src "physgrid-new" 

/home/mvalera/src/$src/cime/scripts/create_newcase --case /scratch/cluster/mvalera/$caze --compset FHS94 --res ne30_ne30_mg16 --walltime 32:00:00 --q verylong --pecount 192 --compiler nag --run-unsupported
#/home/mvalera/src/physgrid/cime/scripts/create_newcase --case /scratch/cluster/mvalera/$caze --compset FHS94 --res ne30_ne30 --walltime 32:00:00 --q overnight --pecount 384 --compiler nag --run-unsupported

cd /scratch/cluster/mvalera/$caze

./xmlchange NTHRDS=1
./xmlchange STOP_OPTION=ndays,STOP_N=450
./xmlchange DOUT_S=FALSE
./xmlchange CAM_CONFIG_OPTS="-phys held_suarez -nadv_tt 5 -analytic_ic"
#./xmlchange DEBUG=TRUE
./case.setup

if ( -d $sourceMods) then
echo "Copying source modes from $sourceMods to case directory /scratch/cluster/$USER/$caze/SourceMods/src.cam/"
foreach file ($sourceMods/*.F90)
  ln -s $file /scratch/cluster/$USER/$caze/SourceMods/src.cam/
  echo "Creating symbolic link $file -> /scratch/cluster/$USER/$caze/SourceMods/src.cam/"
end
foreach file ($sourceMods/*.h)
  ln -s $file /scratch/cluster/$USER/$caze/SourceMods/src.cam/
  echo "Creating symbolic link $file -> /scratch/cluster/$USER/$caze/SourceMods/src.cam/"
end

echo "se_statefreq       = 144                                                    ">> user_nl_cam
echo "adding energy diagnostics to fincl"
echo "nhtfrq= 0,-6,0,0    ">> user_nl_cam
echo "fincl2 = 'U:I','V:I','T:I'" >> user_nl_cam
echo "fincl3 = 'SE_dED','KE_dED', ">> user_nl_cam
echo "           'SE_dAF','KE_dAF', ">> user_nl_cam
echo "           'SE_dBD','KE_dBD', ">> user_nl_cam
echo "           'SE_dAD','KE_dAD', ">> user_nl_cam
echo "           'SE_dAR','KE_dAR', ">> user_nl_cam
echo "           'SE_dBF','KE_dBF', ">> user_nl_cam
echo "           'SE_dBH','KE_dBH', ">> user_nl_cam
echo "           'SE_dCH','KE_dCH', ">> user_nl_cam
echo "           'SE_dAH','KE_dAH'  ">> user_nl_cam
echo "fincl4 = 'TT_UN','TT_LW','TT_MD','TT_HI','TTRMD'" >>  user_nl_cam
echo "analytic_ic_type='held_suarez_1994'">> user_nl_cam
#echo "ncdata = '/scratch/cluster/mvalera/FHS94-1200d/run/FHS94.cam.i.0002-01-01-00000.nc' " 
#start at year 2
echo "ndens              = 2,2,1,2                         ">> user_nl_cam
echo "interpolate_output = .true.,.true.,.false.,.true." >> user_nl_cam
#./case.build --clean-all
#./case.build --clean
./case.build
./case.submit
