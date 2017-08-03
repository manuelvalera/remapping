#!/bin/tcsh

if ( "$#argv" != 2) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 is run case"
  echo "  -arg 2 is history file number (e.g., h0)"
  exit
endif
set n = 1
set case = "$argv[$n]"
set n = 2
set hn = "$argv[$n]"
if (`hostname` == "hobart.cgd.ucar.edu") then
  set data_dir = "/scratch/cluster/$USER/"
else
  set data_dir = "/glade/scratch/$USER/"
endif
set files = `ls $data_dir/$case-1200d/run/months78/$case.cam.$hn.*`
#set files = `ls $data_dir/$case-0ld/run/$case.cam.$hn.*`
#set files = `ls $data_dir/$case/run/months78/$case.cam.$hn.*`
echo $files

ncrcat -O -v T $files {$case}.total.$hn.nc
ncwa -O -v T -a lev {$case}.total.$hn.nc {$case}.var.$hn.nc
ncbo -O -v T {$case}.total.$hn.nc {$case}.var.$hn.nc {$case}.var.$hn.nc
ncra -O -y avgsqr {$case}.var.$hn.nc {$case}.var.$hn.nc

unset case
unset hn
unset files
unset data_dir

#It is possible to use a combination of these operations to compute the variance and standard deviation of a field stored in a single file or across multiple files. The procedure to compute the temporal standard deviation of the surface pressure at all points in a single file in.nc involves three steps.

#ncwa -O -v prs_sfc -a time in.nc out.nc
#ncbo -O -v prs_sfc in.nc out.nc out.nc 
#ncra -O -y rmssdn out.nc out.nc

#First construct the temporal mean of prs_sfc in the file out.nc. Next overwrite out.nc with the anomaly (deviation from the mean). Finally overwrite out.nc with the root-mean-square of itself. Note the use of ‘-y rmssdn’ (rather than ‘-y rms’) in the final step. This ensures the standard deviation is correctly normalized by one fewer than the number of time samples. The procedure to compute the variance is identical except for the use of ‘-y avgsqr’ instead of ‘-y rmssdn’ in the final step.

#ncwa -O -v T -a lev FHS-1200-phi.cam.h1.0001-01-31-00000.nc out.nc
#ncbo -O -v T FHS-1200-phi.cam.h1.0001-01-31-00000.nc out.nc out.nc
#ncra -O -y avgsqr out.nc out.nc
