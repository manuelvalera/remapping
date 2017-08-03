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
  diri  = "/scratch/cluster/mvalera/ape-pcm/run/"
;  fili  = "ape-pcm.ave.15.h3.nc" 
  fili = "ape-pcm.cam.h3.0001-01-01-00000.nc"
  nt = 15
  if (fileexists(diri+fili)) then
	  file1 = addfile(diri+fili, "r")
	  else
	  print("")
	  print("You don't have the necessary data for this script. MISSING:"+fili)
	  print("")
          exit
  end if


  tmp  = file1->TT_HI
  ;  tmp  = file1->TT_MD
  ;  tmp  = file1->TT_LW 
  tmp_zonal = dim_avg_Wrap( tmp(time|:, lev|:, lat|:, lon|:) )
  var  = tmp_zonal(nt,:,:)


  ;---- define the workstation (plot output type and name)
  wks = gsn_open_wks("png","ape-15-pcm")

  ;---- set resources
  res                       =  True
  res@gsnMaximize           =  True
  res@cnMonoLineDashPattern =  True
;  res@cnLineThicknessF      =  3
  res@trYReverse	    = True
  res@tiMainString          = "PCM - High Tracer day 15 - Aquaplanet"   ;-- title string
  res@tiYAxisString       = "Pressure (hPa)"
  res@tiXAxisString       = "latitude"
  res@tiMainFontHeightF     =  0.02
;  res@cnLevelSpacingF = 4.0                      ; contour interval
  res@cnFillOn        = True                    ; turn on color
  res@cnFillPalette   = "BlAqGrYeOrReVi200"     ; set color map
  res@cnLineLabelsOn = False   
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels    = (/0.,0.0001,0.0005,0.001,0.002,0.003,0.004,0.005,0.01,0.02,0.03,0.04,0.05,0.07,0.09,0.1,0.2,0.3,0.4,0.5,1. /)   ; set levels

  ;---- draw the contour map
 ;  plot = gsn_csm_contour_map(wks, var, res) 
  plot = gsn_csm_hov(wks, var, res)

;---- read the data and define variable reference var
    diri  = "/scratch/cluster/mvalera/ape-ppm-mom/run/"
   ; fili  = "ape-ppm-mom.ave.15.h3.nc"
    fili = "ape-ppm-mom.cam.h3.0001-01-01-00000.nc"
    if (fileexists(diri+fili)) then
       file1 = addfile(diri+fili, "r")
    else
       print("")
       print("You don't have the necessary data for this script. MISSING: "+fili)
       print("")
       exit
     end if

  tmp  = file1->TT_HI
  ;  tmp  = file1->TT_MD
  ;  tmp  = file1->TT_LW 
  tmp_zonal = dim_avg_Wrap( tmp(time|:, lev|:, lat|:, lon|:) )
  var  = tmp_zonal(nt,:,:)


     ;---- define the workstation (plot output type and name)
     wks = gsn_open_wks("png","ape-15-ppm")

     ;---- set resources
     res                       =  True
     res@gsnMaximize           =  True
     res@cnMonoLineDashPattern =  True
    ;  res@cnLineThicknessF      =  3
     res@trYReverse            = True
     res@tiMainString          = "PPM - High Tracer day 15 - Aquaplanet"   ;-- title string
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
      diri  = "/scratch/cluster/mvalera/ape-pqm-ih4ih3/run/"
  ;    fili  = "ape-pqm-ih4ih3.ave.15.h3.nc"
  fili = "ape-pqm-ih4ih3.cam.h3.0001-01-01-00000.nc"

      if (fileexists(diri+fili)) then
         file1 = addfile(diri+fili, "r")
      else
	print(" ERROR - PQM 15 FILE NOT FOUND")
	exit
      end if
  
    tmp  = file1->TT_HI
    ;  tmp  = file1->TT_MD
    ;  tmp  = file1->TT_LW 
    tmp_zonal = dim_avg_Wrap( tmp(time|:, lev|:, lat|:, lon|:) )
      var  = tmp_zonal(nt,:,:)


      wks = gsn_open_wks("png","ape-15-pqm")

      res@tiMainString          = "PQM - High Tracer day 15 - Aquaplanet"   ;-- title string
      plot = gsn_csm_hov(wks, var, res)
  
  end

