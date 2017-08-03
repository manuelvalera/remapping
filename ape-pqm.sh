#!/bin/tcsh
set sourceMods = "/home/$USER/remapping/SourceMods2/"
set caze       = "ape-spinup"
setenv src "physgrid-new" 

/home/mvalera/src/$src/cime/scripts/create_newcase --case /scratch/cluster/mvalera/$caze --compset QPC4  --res ne30_ne30_mg16 --walltime 32:00:00 --q overnight --pecount 384 --compiler nag --run-unsupported

cd /scratch/cluster/mvalera/$caze

./xmlchange NTHRDS=1
./xmlchange STOP_OPTION=ndays,STOP_N=200
#./xmlchange CAM_CONFIG_OPTS="-nadv_tt 5 -ocn aquaplanet -aquaplanet"
./xmlchange --append CAM_CONFIG_OPTS="-nadv_tt 5"
./xmlchange DOUT_S=FALSE
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

echo "adding energy diagnostics to fincl"
echo "fincl3 = 'WV_pBF','WL_pBF','WI_pBF','SE_pBF','KE_pBF', ">> user_nl_cam
echo "           'WV_pBP','WL_pBP','WI_pBP','SE_pBP','KE_pBP', ">> user_nl_cam
echo "           'WV_pAP','WL_pAP','WI_pAP','SE_pAP','KE_pAP', ">> user_nl_cam
echo "           'WV_pAM','WL_pAM','WI_pAM','SE_pAM','KE_pAM', ">> user_nl_cam
echo "           'WV_dED','WL_dED','WI_dED','SE_dED','KE_dED', ">> user_nl_cam
echo "           'WV_dAF','WL_dAF','WI_dAF','SE_dAF','KE_dAF', ">> user_nl_cam
echo "           'WV_dBD','WL_dBD','WI_dBD','SE_dBD','KE_dBD', ">> user_nl_cam
echo "           'WV_dAD','WL_dAD','WI_dAD','SE_dAD','KE_dAD', ">> user_nl_cam
echo "           'WV_dAR','WL_dAR','WI_dAR','SE_dAR','KE_dAR', ">> user_nl_cam
echo "           'WV_dBF','WL_dBF','WI_dBF','SE_dBF','KE_dBF', ">> user_nl_cam
echo "           'WV_dBH','WL_dBH','WI_dBH','SE_dBH','KE_dBH', ">> user_nl_cam
echo "           'WV_dCH','WL_dCH','WI_dCH','SE_dCH','KE_dCH', ">> user_nl_cam
echo "           'WV_dAH','WL_dAH','WI_dAH','SE_dAH','KE_dAH'  ">> user_nl_cam

echo "fincl4 = 'TT_UN','TT_LW','TT_MD','TT_HI','TTRMD'" >>  user_nl_cam
echo "ndens              = 2,2,1,2                         ">> user_nl_cam
echo "interpolate_output = .true.,.true.,.false.,.true." >> user_nl_cam
#./case.build --clean-all
#./case.build --clean
./case.build
./case.submit
