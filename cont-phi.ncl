;
;
;
begin
;---- read the data and define variable reference var
  diri  = "/scratch/cluster/mvalera/FHS-1200-phi-ovn619/run/"
    ;  fili  = "fhs-30-avg-pcm.h3.nc"
    fili  = "FHS-1200-phi-ovn619.cam.h1.0004-04-01-00000.nc"
      if (fileexists(diri+fili)) then
              file1 = addfile(diri+fili, "r")
      else
              print("")
              print("You don't have the necessary data for this script. MISSING:"+fili)
              print("")
              exit
      end if

      var   = file1->T(0,1,:,:)
      wks = gsn_open_wks("png","fhs-phi-3D")

       ;---- set resources
      res                       =  True
      res@gsnMaximize           =  True
      res@cnMonoLineDashPattern =  True
    ;  res@cnLineThicknessF      =  3
      res@trYReverse            = True
      res@tiMainString          = "Geopotential Artifacts, Temperature field - Held Suarez"   ;-- title string
 ;     res@tiYAxisString       = "latitude"
 ;     res@tiXAxisString       = "longitude"
      res@tiMainFontHeightF     =  0.02
	      ;  res@cnLevelSpacingF = 4.0                      ; contour interval
      res@cnFillOn        = True                    ; turn on color
      res@cnFillPalette   = "BlAqGrYeOrReVi200"     ; set color map
      res@cnLineLabelsOn = False


       plot = gsn_csm_contour_map(wks, var, res)


      diri  = "/scratch/cluster/mvalera/FHS-spinup/run/"
      ;  fili  = "fhs-30-avg-pcm.h3.nc"
      fili  = "FHS-spinup.cam.h1.0001-06-30-00000.nc"
      if (fileexists(diri+fili)) then
	            file1 = addfile(diri+fili, "r")
      else
	      print("")
              print("You don't have the necessary data for this script. MISSING:"+fili)
	      print("")
	      exit
      end if
      var   = file1->T(0,1,:,:)
      wks = gsn_open_wks("png","fhs-nophi-3D")
     ;---- set resources
     res                       =  True
     res@gsnMaximize           =  True
     res@cnMonoLineDashPattern =  True
    ;  res@cnLineThicknessF      =  3
     res@trYReverse            = True
     res@tiMainString          = "Original Temperature Field - Held Suarez"   ;-- title string
;     res@tiYAxisString       = "latitude"
;     res@tiXAxisString       = "longitude"
     res@tiMainFontHeightF     =  0.02
    ;  res@cnLevelSpacingF = 4.0                      ; contour interval
     res@cnFillOn        = True                    ; turn on color
     res@cnFillPalette   = "BlAqGrYeOrReVi200"     ; set color map
     res@cnLineLabelsOn = False
     plot = gsn_csm_contour_map(wks, var, res)

     diri  = "/scratch/cluster/mvalera/FHS-1200-phi-ovn618/run/months78/"
     ;  fili  = "fhs-30-avg-pcm.h3.nc"
;     fili  = "FHS-1200-phi-ovn618.var.h1.nc"
     fili  = "diff.nc"

     if (fileexists(diri+fili)) then
                     file1 = addfile(diri+fili, "r")
     else
	             print("")
                     print("You don't have the necessary data for this script. MISSING:"+fili)
                     print("")
                     exit
     end if

     var   = file1->T(0,0,:,:)
     wks = gsn_open_wks("png","fhs-var-3D")

     res@tiMainString          = "Variance difference in Temperature field  - Held Suarez"   ;-- title 
     plot = gsn_csm_contour_map(wks, var, res)


  end
