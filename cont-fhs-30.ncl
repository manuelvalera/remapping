; Script to make the final plots for the SiParCS internship, Summer 2017,
; Manuel Valera. 
; Based on script examples from the ncl manual & resources.
;-----------------------------------------------------------------
;  NCL User Guide Example:   NUG_contour_dashpattern.ncl
;
;  KMF 30.10.14
;-----------------------------------------------------------------
; These load commands are not required in NCL versions 6.2.0 and later.
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;---- read the data and define variable reference var
  diri  = "/scratch/cluster/mvalera/fhs-15-pcm/run/"
  ;  fili  = "fhs-30-avg-pcm.h3.nc"
  fili  = "fhs-15-pcm.ave.h3.nc" 
  if (fileexists(diri+fili)) then
	  file1 = addfile(diri+fili, "r")
	  else
	  print("")
	  print("You don't have the necessary data for this script. You can download it from:â€‹")
	  print("")
      	  print("http://www.ncl.ucar.edu/Document/Manuals/NCL_User_Guide/Data/"+fili)
          print("")
     	  print("or use the wget command:")
	  print("")
          print("wget http://www.ncl.ucar.edu/Document/Manuals/NCL_User_Guide/Data/"+fili)
          print("")
          exit
  end if

  var   = file1->TT_HI(0,:,:,30)

  ;---- define the workstation (plot output type and name)
  wks = gsn_open_wks("png","fhs-30-pcm")

  ;---- set resources
  res                       =  True
  res@gsnMaximize           =  True
  res@cnMonoLineDashPattern =  True
;  res@cnLineThicknessF      =  3
  res@trYReverse	    = True
  res@tiMainString          = "PCM - High Tracer day 30 - Held Suarez"   ;-- title string
  res@tiYAxisString       = "Pressure (hPa)"
  res@tiXAxisString       = "latitude"
  res@tiMainFontHeightF     =  0.02
;  res@cnLevelSpacingF = 4.0                      ; contour interval
  res@cnFillOn        = True                    ; turn on color
  res@cnFillPalette   = "BlAqGrYeOrReVi200"     ; set color map
  res@cnLineLabelsOn = False   
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels    = (/0.,0.0001,0.0005,0.001,0.002,0.003,0.004,0.005,0.01,0.02,0.03,0.04,0.05,0.07,0.09,0.1,0.2,0.3,0.4,0.5,1. /)   ; set levels


;e res@cnLevelSelectionMode = "ManualLevels" 
;  res@cnLevelSpacingF = 0.05 
;  res@cnMinLevelValF = 0 
;  res@cnMaxLevelValF = 0.9 

  ;---- draw the contour map
 ;  plot = gsn_csm_contour_map(wks, var, res) 
  plot = gsn_csm_hov(wks, var, res)

;---- read the data and define variable reference var
    diri  = "/scratch/cluster/mvalera/fhs-15-ppm-mom/run/"
    fili  = "fhs-15-ppm-mom.ave.h3.nc"

    if (fileexists(diri+fili)) then
       file1 = addfile(diri+fili, "r")
    else
       print("")
       print("You don't have the necessary data for this script. You can download it from:â€‹")
       print("")
       print("http://www.ncl.ucar.edu/Document/Manuals/NCL_User_Guide/Data/"+fili)
       print("")
       print("or use the wget command:")
       print("")
       print("wget http://www.ncl.ucar.edu/Document/Manuals/NCL_User_Guide/Data/"+fili)
       print("")
       exit
     end if

     var   = file1->TT_HI(0,:,:,30)

     ;---- define the workstation (plot output type and name)
     wks = gsn_open_wks("png","fhs-30-ppm")

     ;---- set resources
     res                       =  True
     res@gsnMaximize           =  True
     res@cnMonoLineDashPattern =  True
    ;  res@cnLineThicknessF      =  3
     res@trYReverse            = True
     res@tiMainString          = "PPM - High Tracer day 30 - Held Suarez"   ;-- title string
     res@tiYAxisString       = "Pressure (hPa)"
     res@tiXAxisString       = "latitude"
     res@tiMainFontHeightF     =  0.02
    ;  res@cnLevelSpacingF = 4.0                      ; contour interval
     res@cnFillOn        = True                    ; turn on color
     res@cnFillPalette   = "BlAqGrYeOrReVi200"     ; set color map
     res@cnLineLabelsOn = False
     res@cnLevelSelectionMode = "ExplicitLevels"
     ;     res@cnLevels    = (/0.,0.0005,0.005,0.01,0.02,0.03,0.04,0.05,0.07,0.09,0.1,0.2,0.3,0.4,0.5,1. /)   ; set levels

      plot = gsn_csm_hov(wks, var, res)


; 0------- figure 3 - pqm 30:      
      diri  = "/scratch/cluster/mvalera/fhs-15-pqm-ih4ih3/run/"
      fili  = "fhs-15-pqm-ih4ih3.ave.h3.nc"

      if (fileexists(diri+fili)) then
         file1 = addfile(diri+fili, "r")
      else
	print(" ERROR - PQM 30 FILE NOT FOUND")
      end if
      var   = file1->TT_HI(0,:,:,30)
      wks = gsn_open_wks("png","fhs-30-pqm")

      res@tiMainString          = "PQM - High Tracer day 30 - Held Suarez"   ;-- title string
      plot = gsn_csm_hov(wks, var, res)

; 0------- figure 3 - pqm 30:      
   diri  = "/scratch/cluster/mvalera/fhs-15-ppm-mom/run/"
   fili  = "avediff.nc"

   if (fileexists(diri+fili)) then
      file1 = addfile(diri+fili, "r")
					       else
      print(" ERROR - PQM-PPM-DIFF 30 FILE NOT FOUND")
   end if

   var   = file1->TT_HI(0,:,:,30)
   wks = gsn_open_wks("png","diff-fhs-30-pqm-ppm")

   res@tiMainString          = "PQM vs PPM Diffs - High Tracer day 30 - Held Suarez"   ;-- title string
   plot = gsn_csm_hov(wks, var, res)
   
  end

