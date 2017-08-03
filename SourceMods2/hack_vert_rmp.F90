module hack_vert_rmp
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: nlev, np, qsize, nc
  use element_mod,            only: element_t
  use hybvcoord_mod,          only: hvcoord_t
  implicit none

  logical, parameter, public :: lhack_vert_rmp=.false.
  !
  public :: get_levels
  public :: overwrite_state
  public :: write_data
  public :: write_data_TE,write_data_TE_eul,TE_probe
  public :: diagnostic,diagnostic_eul
  public :: trapezoid_integration

contains

  subroutine get_levels(option,pint1,pint2,t0)
    integer,intent(in) :: option
    real(KIND=r8),intent(inout),dimension(31) :: pint1,pint2   
    real(KIND=r8),dimension(31) :: p11,p21
    real(KIND=r8),dimension(47) :: p12,p22
    real(KIND=r8),dimension(61) :: p13,p23
    real(KIND=r8),dimension(84) :: p14,p24
    real(KIND=r8),dimension(67) :: p15,p25     

    real(KIND=r8),intent(inout),dimension(30) :: t0
    real(KIND=r8),dimension(30) :: t1
    real(KIND=r8),dimension(46) :: t2
    real(KIND=r8),dimension(60) :: t3
    real(KIND=r8),dimension(83) :: t4
    real(KIND=r8),dimension(66) :: t5

    !allocate(p11(nlev+1))    !TODO: Find a way to assign this data dinamycally
    !allocate(p21(nlev+1))
    !allocate(t1(nlev))

    if(option==1) then !CAM5 levels
           ! coordinate 1
      data p11 /2.2552, 5.0317, 10.158, 18.555, 30.669, 45.867, 63.323, &
                  80.701, 94.941, 111.69, 131.40, 154.59, 181.86, 213.95, &
                  251.70, 296.12, 348.37, 409.84, 482.15, 567.22, 652.33, &
                  730.45, 796.36, 845.35, 873.72, 900.32, 924.96, 947.43, &
                  967.54, 985.11, 1000./

        ! coordinate 2
      data p21 /2.2552, 4.0317, 12.158, 16.555, 28.669, 41.867, 60.323, &
                  84.701, 99.941, 118.69, 128.40, 150.59, 175.86, 223.95, &
                  261.70, 286.12, 338.37, 399.84, 472.15, 547.22, 632.33, &
                  750.45, 786.36, 825.35, 883.72, 908.32, 920.96, 940.43, &
                  960.54, 988.11, 1000./

        ! temperature profile under coordinate 1
      data t1    /2.625231323242188e+02, 2.422447662353516e+02,2.311667327880859e+02, &
                  2.167325744628906e+02, 2.105011749267578e+02,1.955296630859375e+02, &
                  1.822942810058594e+02, 1.818112792968750e+02,1.885818328857422e+02, &
                  1.960807952880859e+02, 2.042872161865234e+02,2.128411560058594e+02, &
                  2.217830505371094e+02, 2.308533477783203e+02,2.400057525634766e+02, &
                  2.491661071777344e+02, 2.572124023437500e+02,2.653801269531250e+02, &
                  2.727942810058594e+02, 2.780478820800781e+02,2.837588806152344e+02, &
                  2.879168090820312e+02, 2.906744995117188e+02,2.918422851562500e+02, &
                  2.924194641113281e+02, 2.934298095703125e+02,2.948598327636719e+02, &
                  2.965856018066406e+02, 2.982542114257812e+02,2.996733398437500e+02/

      pint1 = p11
      pint2 = p21
      t0    = t1

    else if(option==2)then !CAML46 levels
   
       ! coordinate 1
      data p12 /0.27964, 0.48322, 0.75879, 1.0827, 1.4069, 1.8189, 2.3398,2.9951, 3.8147, 4.8345, &
       6.0964, 7.6494, 9.5501, 11.864, 14.666, 18.038, 22.076, 26.883, 32.574,39.273,&
       47.115, 56.241, 66.8,   80.701, 94.941, 111.69, 131.4,  154.59, 181.86,213.95,&
        251.7, 296.12, 348.37, 409.84, 482.15, 567.22, 652.33, 730.45, 796.36,845.35,&
       873.72, 900.32, 924.96, 947.43, 967.54, 985.11, 1000./

        ! coordinate 2
      data p22 /0.27964, 0.58322, 0.85879, 1.1827, 1.6069, 1.9989, 2.4398,2.8951, 3.9147, 4.7345, &
       5.8964, 6.6494, 8.5501, 10.864, 15.666, 17.038, 20.076, 25.883, 30.574,38.273,&
       45.115, 53.241, 62.8,   81.701, 91.941, 107.69, 125.4,  148.59, 189.86,219.95,&
        245.7, 286.12, 338.37, 419.84, 472.15, 547.22, 632.33, 720.45, 786.36,855.35,&
       883.72, 905.32, 934.96, 940.43, 955.54, 975.11, 1000./

        ! temperature profile under coordinate 1
      data t2    /250.76, 251.96, 250.88, 248.57, 245.59, 242.38, 238.73,235.58, 236.07, 235.54, &
       233.47, 227.96, 223.38, 220.49, 217.79, 215.81, 215.28, 214.45, 212.4,210.07,&
       209.18, 209.66, 211.19, 214.42, 216.26, 217.84, 217.55, 215.8,  214.63,215.32,&
       219.42, 226.11, 233.85, 242.21, 251.98, 259.42, 263.56, 267.21, 270.26,272.46,&
       273.8,  275.34, 277.03, 278.75, 280.34, 281.73/

    else if(option==3)then !CAML60 levels
        ! coordinate 1
      data p13 /1.8341, 3.1692, 5.2058, 8.1295, 12.069, 17.033, 22.854,29.151, 35.35, 40.751,&
        44.66, 47.82, 51.203,  54.826, 58.705, 62.859, 67.306, 72.069, 77.168,82.628,&
       88.474, 94.733, 101.44, 108.61, 116.3,  124.53, 133.34, 142.77, 152.87,163.69,&
       175.27, 187.67, 200.95, 215.17, 230.39, 246.69, 264.15, 282.83, 302.85,324.27,&
       347.22, 371.78, 398.09, 426.25, 456.41, 488.7,  523.28, 560.31, 599.95,642.4,&
       687.85, 736.52, 788.63, 844.42, 872.98, 899.76, 924.56, 947.18, 967.42,985.11,&
         1000./

        ! coordinate 2
      data p23 /1.8341, 2.3692, 4.2058, 6.1295, 10.069, 14.033, 20.854,26.151, 33.35, 42.751,&
        45.66, 48.82, 50.203,  55.826, 59.705, 60.859, 65.306, 70.069, 79.168,85.628,&
       90.474, 97.733, 100.44, 110.61, 119.3,  127.53, 130.34, 140.77, 150.87,157.69,&
       170.27, 182.67, 195.95, 210.17, 222.39, 236.69, 252.15, 272.83, 309.85,314.27,&
       337.22, 361.78, 388.09, 416.25, 446.41, 490.7,  513.28, 565.31, 589.95,632.4,&
       677.85, 716.52, 768.63, 824.42, 862.98, 889.76, 914.56, 957.18, 969.42,989.11,&
         1000./

        ! temperature profile under coordinate 1
      data t3    /247.63, 240.61, 232.26, 228.98, 226.05, 222.72, 219.34,216.86, 213.74, 211.99,&
       209.79, 206.56, 202.72, 199.46, 196.96, 195.13, 193.93, 192.75,190.77,188.72,&
       187.61, 187.4,  188.26, 190.15, 192.6,  195.32, 198.26, 201.48,204.82,208.34,&
       211.97, 215.68, 219.47, 223.33, 227.26, 231.24, 235.25, 239.32,243.25,247.12,&
       250.77, 254.44, 257.9,  261.18, 264.23, 267.45, 270.88, 273.5,276.12,279.66,&
       282.85, 285.46, 287.96, 290.32, 291.77, 293.37, 294.97, 296.66,298.18,299.44/

    else if(option==4)then !CAML83 levels
        ! coordinate 1
      data p14 /0.24578, 0.42471, 0.6669, 0.95253, 1.2783, 1.6859, 2.1836,2.7773, 3.4695, 4.2589,&
       5.1417, 6.1115, 7.1607, 8.282, 9.4684, 10.715, 12.018, 13.376, 14.791,16.264,&
       17.799, 19.402, 21.08,  22.84, 24.689, 26.639, 28.696, 30.873, 33.179,35.626,&
       38.226, 40.992, 43.936, 47.072, 50.416, 53.983, 57.789, 61.853, 66.192,70.81,&
       75.727, 80.985, 86.608, 92.622, 99.054, 105.93, 113.29, 121.15, 129.57,138.56,&
       148.19, 158.48, 169.48, 181.25, 193.83, 207.29, 221.69, 237.08, 253.54,271.15,&
       289.98, 310.11, 331.65, 354.68, 379.3, 405.64, 433.81,  463.93, 496.15,530.6,&
       567.44, 606.85, 648.98, 694.05, 742.24, 796.39, 845.37, 873.73, 900.34,924.97,&
       947.44, 967.54, 985.11, 1000./

        ! coordinate 2
      data p24 /0.24578, 0.32471, 0.7669, 0.85253, 1.0783, 1.3859, 2.0836,2.9773, 3.1695, 4.0589,&
       5.5417, 6.4115, 7.0607, 8.482, 9.0684, 11.715, 12.518, 13.876, 15.791,16.964,&
       18.799, 20.402, 21.58,  23.84, 24.989, 27.639, 29.696, 31.873, 33.979,36.626,&
       39.226, 40.592, 44.936, 46.072, 51.416, 54.983, 56.789, 63.853, 67.192,72.81,&
       77.727, 82.985, 85.608, 90.622, 97.054, 103.93, 116.29, 124.15, 127.57,140.56,&
       151.19, 159.48, 166.48, 185.25, 197.83, 212.29, 226.69, 239.08, 250.54,273.15,&
       285.98, 302.11, 321.65, 344.68, 369.3, 395.64, 423.81,  453.93, 486.15,520.6,&
       577.44, 600.85, 640.98, 699.05, 732.24, 786.39, 825.37, 863.73, 910.34,934.97,&
       957.44, 969.54, 989.11, 1000./

        ! temperature profile under coordinate 1
      data t4    /251.35, 261.83, 262.18, 261.98, 261.1, 256.67, 250.63, 244.67,239.84, 235.77, &
       231.95, 229.25, 226.85, 224.52, 222.28, 220.41, 218.93, 217.39, 215.78,214.66,&
       214.22, 213.8,  213.8,  213.79, 213.66, 213.43, 213.2,  213.3,  213.4,213.6,&
       213.95, 214.33, 214.08, 213.58, 212.98, 211.91, 210.76, 210.09, 209.8,209.53,&
       209.82, 210.14, 210.05, 209.71, 209.41, 209.34, 209.28, 209.46, 209.76,210.06,&
       210.27, 210.49, 211.17, 212.01, 213.03, 214.34, 215.75, 218.16, 220.78,223.96,&
       227.63, 231.43, 234.94, 238.69, 242.51, 246.55, 250.38, 253.78, 257.41,260.83,&
       264.49, 267.25, 269.49, 271.98, 274.85, 275.8,  275.65, 276.82, 278.41,279.88,&
       281.22, 282.29, 282.97/


    else if(option==5)then ! WACMM levels

        ! coordinate 1
      data p15 /4.5005e-06, 7.4201e-06, 1.2234e-05, 2.017e-05, 3.3255e-05,5.4827e-05, 9.0398e-05, 0.00014904, 0.00024572, &
                  0.00040512, &
                  0.00066794, 0.0011013, 0.0018157, 0.0029935, 0.004963,0.0081507, 0.013477, 0.022319, 0.036796, 0.060665, &
                  0.099157, 0.15739, 0.23885, 0.3452, 0.47514, 0.63181, 0.82916,1.0827, 1.4069, 1.8188, &
                  2.3398, 2.995, 3.8147, 4.8344, 6.0964, 7.6494, 9.5501, 11.864,14.665, 18.038, &
                  22.076, 26.883, 32.573, 39.273, 47.114, 56.24, 66.8, 78.949,92.366, 108.66, &
                  127.84, 150.39, 176.93, 208.15, 244.88, 288.09, 338.92,398.72, 469.07, 551.84, &
                  649.21, 744.38, 831.02, 903.3, 956.0, 985.11, 1000.0/

        ! coordinate 2
      data p25 /4.5005e-06, 8.0e-06, 1.5234e-05, 2.417e-05, 4.3255e-05,6.4827e-05, 8.0398e-05, 0.00018904, 0.00028572, &
                  0.00042512, &
                  0.00086794, 0.0013013, 0.0017157, 0.0028935, 0.004663,0.0076507, 0.012477, 0.021319, 0.033796, 0.058665, &
                  0.095157, 0.13739, 0.21885, 0.3052, 0.44514, 0.60181, 0.80916,1.1827, 1.6069, 1.9188, &
                  2.1398, 2.895, 3.6147, 4.6344, 6.1964, 7.2494, 9.1501, 10.864,13.665, 17.038, &
                  23.076, 27.883, 31.573, 38.273, 48.114, 55.24, 67.8, 74.949,90.366, 105.66, &
                  123.84, 152.39, 172.93, 205.15, 248.88, 290.09, 335.92,393.72, 460.07, 555.84, &
                  660.21, 749.38, 838.02, 910.3, 950.0, 980.11, 1000.0/

        ! temperature profile under coordinate 1
      data t5    /675.28, 601.1, 512.79, 420.03, 338.77, 287.33, 241.95, 217.02,202.61, 198.48, &
       196.08, 196.68, 198.07, 198.92,       199.71,       202.58,       205.17,213.47,       218.88,       229.18,&
       238.36, 250.23, 258.78, 264.16,       265.99,       266.31,       265.14,263.07,       262.69,       258.98,&
       254.01, 250.6,  246.45, 239.42,       234.52,       232.41,       231.07,228.65,       226.42,       223.94,&
       220.65, 218.62, 215.07, 209.61,        204.8,       198.95,       188.04,188.99,       192.39,       195.71,&
       201.52, 209.14, 218.24, 227.86,       237.54,       246.87,       255.34,263.13,       270.15,        275.3,&
       282.88, 286.95, 291.98, 295.81,       298.64,       300.25/


    endif

    !deallocate(p11)
    !deallocate(p21)
    !deallocate(t1)


  end subroutine get_levels
   



  subroutine overwrite_state(elem)
    type (element_t), intent(inout):: elem
  end subroutine overwrite_state

  subroutine write_data(pint1,t0,v,filenum)
    use spmd_utils,       only: masterproc
    integer           :: unitn,k
    integer,intent(in):: filenum
    !character(len=256):: filename = 'test.dat' 
    character(len=256):: filename  

    real(KIND=r8), dimension(31),intent(in) :: pint1
    real(KIND=r8), dimension(31) :: p_k
    real(KIND=r8), dimension(30),intent(in) :: t0,v

    if (masterproc) then

      p_k = 0.0d0 
      unitn = 8
   
      !writing 2-column file for temp state 
      write (filename, '("p_k_T_", I0.3, ".dat")' )  filenum    
      open(unitn, file=trim(filename), status='replace' )
      do k=1,size(pint1)-1
        p_k(k) = (1/2.)*(pint1(k+1)+pint1(k))
        write(unitn,*) p_k(k),t0(k)
      end do
      close(unitn)
    
      p_k = 0.0d0
      
      !writing 2-column file for vel state 
      write (filename, '("p_k_u_", I0.3, ".dat")' )  filenum   
      open(unitn, file=trim(filename), status='replace' )
      do k=1,size(pint1)-1
        p_k(k) = (1/2.)*(pint1(k+1)+pint1(k))
        write(unitn,*) p_k(k),v(k)
      end do
      close(unitn)

    end if
  end subroutine write_data

  subroutine write_data_TE(pint1,ttmp,filenum)
    use spmd_utils,       only: masterproc
    integer           :: unitn,k
    integer,intent(in):: filenum
    character(len=256):: filename

    real(KIND=r8), dimension(31),intent(in) :: pint1
    real(KIND=r8), dimension(31) :: p_k
    real(KIND=r8), dimension(30),intent(in) :: ttmp

    if (masterproc) then

      p_k = 0.0d0
      unitn = 8

      !writing 2-column file for temp state 
      write (filename, '("p_k_TE_", I0.3, ".dat")' )  filenum
      open(unitn, file=trim(filename), status='replace' )
      do k=1,size(pint1)-1
             p_k(k) = (1/2.)*(pint1(k+1)+pint1(k))
             write(unitn,*) p_k(k),ttmp(k)
      end do
      close(unitn)

    end if
  end subroutine write_data_TE

  subroutine write_data_TE_eul(pint1,ttmp,filenum)
    use spmd_utils,       only: masterproc
    integer           :: unitn,k
    integer,intent(in):: filenum
    character(len=256):: filename

    real(KIND=r8), dimension(31),intent(in) :: pint1
    real(KIND=r8), dimension(31) :: p_k
    real(KIND=r8), dimension(30),intent(in) :: ttmp

    if (masterproc) then

      p_k = 0.0d0
      unitn = 8

      !writing 2-column file for temp state 
      write (filename, '("p_k_TE_eul", I0.3, ".dat")' )  filenum
      open(unitn, file=trim(filename), status='replace' )
      do k=1,size(pint1)-1
             p_k(k) = (1/2.)*(pint1(k+1)+pint1(k))
             write(unitn,*) p_k(k),ttmp(k),sum(ttmp(:))
      end do
      close(unitn)

    end if
  end subroutine write_data_TE_eul


  subroutine diagnostic_eul(dpu,dpt,dpk,dpphi,filenum)
    use spmd_utils,       only: masterproc
    integer           :: unitn,k,sz,dpk_k
    integer,intent(in):: filenum
    character(len=256):: filename
    real(KIND=r8), dimension(np,np,nlev,1),intent(in)::dpu,dpk,dpt
    real(KIND=r8), intent(in)     :: dpphi
    real(KIND=r8)                 :: dpphi_0
    real(KIND=r8)      :: delta_phi_0,m_0,phi_max,phi_min,l_inf,l_2,I_n,I_d,end_val
    real(KIND=r8),dimension(nlev)  :: dpu_0,dpt_0,dpk_0
    real(KIND=r8),  parameter :: PI=4*atan(1.0)

    if (masterproc) then

      end_val = 2*PI
      unitn = 8

      !Treating state 0:
      if (filenum == 1)then

        !dpu_phi_max:
        delta_phi_0 = MAXVAL(dpu(1,1,:,1)) - MINVAL(dpu(1,1,:,1))
        m_0 =  MAXVAL(dpu(1,1,:,1))
        write (filename, '("dpu_phi_max_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpu_phi_min:
        m_0 =  MINVAL(dpu(1,1,:,1))
        write (filename, '("dpu_phi_min_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpu_l_inf, dpu_l_2:
        write (filename, '("dpu_l_inf_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpu(1,1,:,1)
        close(unitn)

        write (filename, '("dpu_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpu_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpu(1,1,:,1))
        close(unitn)


      !dpt
        !dpt_phi_max:
        delta_phi_0 = MAXVAL(dpt(1,1,:,1)) - MINVAL(dpt(1,1,:,1))
        m_0 =  MAXVAL(dpt(1,1,:,1))
        write (filename, '("dpt_phi_max_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpt_phi_min:
        m_0 =  MINVAL(dpt(1,1,:,1))
        write (filename, '("dpt_phi_min_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpt_l_inf:
        write (filename, '("dpt_l_inf_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpt(1,1,:,1)
        close(unitn)


        write (filename, '("dpt_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpt_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpt(1,1,:,1))
       close(unitn)
      !dpk
        !dpk_phi_max:
        delta_phi_0 = MAXVAL(dpk(1,1,:,1)) - MINVAL(dpk(1,1,:,1))
        m_0 =  MAXVAL(dpk(1,1,:,1))
        write (filename, '("dpk_phi_max_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpk_phi_min:
        m_0 =  MINVAL(dpk(1,1,:,1))
        write (filename, '("dpk_phi_min_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpk_l_inf:
        write (filename, '("dpk_l_inf_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpk(1,1,:,1)
        close(unitn)


        write (filename, '("dpk_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpk_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpk(1,1,:,1))
        close(unitn)

      !dpk + dpt + dpphi

        write (filename, '("gamma_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpt(1,1,:,1)) + SUM(dpk(1,1,:,1)) + dpphi,SUM(dpt(1,1,:,1)),SUM(dpk(1,1,:,1)),dpphi
        close(unitn)


      !state /= 0
      else
      !dpu 
        !dpu_phi_max:
        filename = "dpu_phi_max_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        if (delta_phi_0/=0) then
                phi_max = (MAXVAL(dpu(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_max = 1.0d0
        endif

        !dpu_phi_min:
        filename = "dpu_phi_min_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        if (delta_phi_0/=0) then
                phi_min = (MINVAL(dpu(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_min = 1.0d0
        end if

        !dpu_l_inf:
        filename = "dpu_l_inf_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpu_0                         !read state 0    
        close(unitn)

         if (sum(dpu_0(:))/=0) then
                l_inf = MAXVAL(ABS(dpu(1,1,:,1) - dpu_0)) / MAXVAL(ABS(dpu_0))
        else
                l_inf = 1.0d0
        end if

        !dpu_l_2:

        call trapezoid_integration((dpu(1,1,:,1)-dpu_0)**2,end_val,I_n)
        call trapezoid_integration(dpu_0**2,end_val,I_d)

        I_n = I_n/(4*PI)
        I_d = I_d/(4*PI)

        if(I_d/=0)then
                l_2 = SQRT( I_n / I_d )
        else
                l_2 = 1.0d0
        end if

        write (filename, '("dpu_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) phi_max,phi_min,l_inf,l_2
        close(unitn)

        write (filename, '("dpu_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpu(1,1,:,1))
        close(unitn)


      !dpt
        !dpt_phi_max:
        filename = "dpt_phi_max_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        phi_max = (MAXVAL(dpt(1,1,:,1)) - m_0) / delta_phi_0

        !dpt_phi_min:
        filename = "dpt_phi_min_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        phi_min = (MINVAL(dpt(1,1,:,1)) - m_0) / delta_phi_0

        !dpt_l_inf:
        filename = "dpt_l_inf_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpt_0                         !read state 0    
        close(unitn)

        l_inf = MAXVAL(ABS(dpt(1,1,:,1) - dpt_0)) / MAXVAL(ABS(dpt_0))

        !dpt_l_2:

        call trapezoid_integration((dpt(1,1,:,1)-dpt_0)**2,end_val,I_n)
        call trapezoid_integration(dpt_0**2,end_val,I_d)

        I_n = I_n/(4*PI)
        I_d = I_d/(4*PI)

        l_2 = SQRT( I_n / I_d )

        write (filename, '("dpt_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) phi_max,phi_min,l_inf,l_2
        close(unitn)

        write (filename, '("dpt_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpt(1,1,:,1))
        close(unitn)


       !dpk
        !dpk_phi_max:
        filename = "dpk_phi_max_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        if(delta_phi_0/=0)then
                phi_max = (MAXVAL(dpk(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_max = 1.0d0
        end if

        !dpk_phi_min:
        filename = "dpk_phi_min_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        if (delta_phi_0/=0) then
                phi_min = (MINVAL(dpk(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_min = 1.0d0
        end if

        !dpk_l_inf:
        filename = "dpk_l_inf_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpk_0                         !read state 0    
        close(unitn)

        if (sum(dpk_0(:))/=0) then
                l_inf = MAXVAL(ABS(dpk(1,1,:,1) - dpk_0))  / MAXVAL(ABS(dpk_0))
        else
                l_inf = 1.0d0
        end if


        !dpk_l_2:

        call trapezoid_integration((dpk(1,1,:,1)-dpk_0)**2,end_val,I_n)
        call trapezoid_integration(dpk_0**2,end_val,I_d)

        I_n = I_n/(4*PI)
        I_d = I_d/(4*PI)

        if(I_d/=0)then
                l_2 = SQRT( I_n / I_d )
        else
                l_2 = 1.0d0
        end if

        write (filename, '("dpk_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) phi_max,phi_min,l_inf,l_2
        close(unitn)

        write (filename, '("dpk_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpk(1,1,:,1))
        close(unitn)
     

      !dpt + dpk + dpphi
        write (filename, '("gamma_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpt(1,1,:,1)) + SUM(dpk(1,1,:,1)) + dpphi,SUM(dpt(1,1,:,1)),SUM(dpk(1,1,:,1)),dpphi
        close(unitn)



      end if

      !dpphi - gotta do it separately...
      if(filenum==1)then

        !dpphi_l_inf, dpphi_l_2:
        write (filename, '("dpphi_l_inf_0_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpphi
        close(unitn)

        write (filename, '("dpphi_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpphi_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpphi
        close(unitn)

      else if(filenum>1)then
        !dpphi_l_inf:
        filename = "dpphi_l_inf_0_eul.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpphi_0                         !read state 0    
        close(unitn)

        if(dpphi_0/=0)then
                l_inf = ABS(dpphi - dpphi_0) / ABS(dpphi_0)
        else
                l_inf = 0.0d0
        end if
        !dpphi_l_2:

        I_n = ((dpphi-dpphi_0)**2)/(2*PI)
        I_d = (dpphi_0**2)/(2*PI)

        if(I_d/=0)then
                l_2 = SQRT( I_n / I_d )
        else
                l_2 = 0.0d0
        end if

        write (filename, '("dpphi_errors_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) 0,0,l_inf,l_2
        close(unitn)


        write (filename, '("dpphi_sum_eul.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) dpphi
        close(unitn)


       end if

    end if

  end subroutine diagnostic_eul

  subroutine TE_probe(cpt,k,phi,dp,pint,system,filenum)
  use spmd_utils,       only: masterproc
  real(KIND=r8), dimension(nlev),intent(in) :: cpt,k,dp
  real(KIND=r8), dimension(nlev+1),intent(in) :: pint,phi
  real(KIND=r8) :: te
  integer :: i,unitn
  integer, intent(in) :: filenum, system ! 1 = lagrangian, 2 = eulerian
  character(len=256):: filename

  if (masterproc) then

      unitn = 8
      te = 0.0_r8

      !Treating state 0:
      if (filenum == 1)then

        if(system == 1)then
               write (filename, '("n_gamma_lag.dat")' ) 
        else
               write (filename, '("n_gamma_eul.dat")' ) 
        end if

        do i=1,nlev
               te = te + dp(i)*(cpt(i)+k(i)+( pint(i+1)*phi(i+1) - pint(i)*phi(i) )/dp(i) )
        end do

        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) te
        close(unitn)

      else
        if(system == 1)then
               write (filename, '("n_gamma_lag.dat")' )
        else
               write (filename, '("n_gamma_eul.dat")' ) 
        end if

        do i=1,nlev
               te = te + dp(i)*(cpt(i)+k(i)+( pint(i+1)*phi(i+1) - pint(i)*phi(i) )/dp(i) )
        end do

        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) te
        close(unitn)

      end if

  end if

  end subroutine TE_probe

  subroutine diagnostic(dpu,dpt,dpk,dpphi,filenum)
    use spmd_utils,       only: masterproc
    integer           :: unitn,k,sz,dpk_k
    integer,intent(in):: filenum
    character(len=256):: filename
    real(KIND=r8), dimension(np,np,nlev,1),intent(in)::dpu,dpk,dpt
    real(KIND=r8), intent(in)     :: dpphi
    real(KIND=r8)                 :: dpphi_0
    real(KIND=r8)      :: delta_phi_0,m_0,phi_max,phi_min,l_inf,l_2,I_n,I_d,end_val
    real(KIND=r8),dimension(nlev)  :: dpu_0,dpt_0,dpk_0
    real(KIND=r8),  parameter :: PI=4*atan(1.0)

    if (masterproc) then

      end_val = 2*PI
      unitn = 8

      !Treating state 0:
      if (filenum == 0)then
        
        !dpu_phi_max:
        delta_phi_0 = MAXVAL(dpu(1,1,:,1)) - MINVAL(dpu(1,1,:,1))
        m_0 =  MAXVAL(dpu(1,1,:,1))
        write (filename, '("dpu_phi_max_0.dat")' ) 
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpu_phi_min:
        m_0 =  MINVAL(dpu(1,1,:,1))
        write (filename, '("dpu_phi_min_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpu_l_inf, dpu_l_2:
        write (filename, '("dpu_l_inf_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpu(1,1,:,1)
        close(unitn)

        write (filename, '("dpu_errors.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpu_sum.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpu(1,1,:,1))
        close(unitn)


      !dpt
        !dpt_phi_max:
        delta_phi_0 = MAXVAL(dpt(1,1,:,1)) - MINVAL(dpt(1,1,:,1))
        m_0 =  MAXVAL(dpt(1,1,:,1))
        write (filename, '("dpt_phi_max_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpt_phi_min:
        m_0 =  MINVAL(dpt(1,1,:,1))
        write (filename, '("dpt_phi_min_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpt_l_inf:
        write (filename, '("dpt_l_inf_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpt(1,1,:,1)
        close(unitn)


        write (filename, '("dpt_errors.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpt_sum.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpt(1,1,:,1))
        close(unitn)
      !dpk
        !dpk_phi_max:
        delta_phi_0 = MAXVAL(dpk(1,1,:,1)) - MINVAL(dpk(1,1,:,1))
        m_0 =  MAXVAL(dpk(1,1,:,1))
        write (filename, '("dpk_phi_max_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpk_phi_min:
        m_0 =  MINVAL(dpk(1,1,:,1))
        write (filename, '("dpk_phi_min_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) delta_phi_0, m_0
        close(unitn)


        !dpk_l_inf:
        write (filename, '("dpk_l_inf_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpk(1,1,:,1)
        close(unitn)


        write (filename, '("dpk_errors.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpk_sum.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpk(1,1,:,1))
        close(unitn)

        write (filename, '("gamma.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) SUM(dpt(1,1,:,1)) + SUM(dpk(1,1,:,1)) + dpphi
        close(unitn)

      !state /= 0
      else
      !dpu 
        !dpu_phi_max:
        filename = "dpu_phi_max_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        if (delta_phi_0/=0) then
                phi_max = (MAXVAL(dpu(1,1,:,1)) - m_0) / delta_phi_0  
        else    
                phi_max = 1.0d0
        endif

        !dpu_phi_min:
        filename = "dpu_phi_min_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)
       
        if (delta_phi_0/=0) then
                phi_min = (MINVAL(dpu(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_min = 1.0d0
        end if

        !dpu_l_inf:
        filename = "dpu_l_inf_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpu_0                         !read state 0    
        close(unitn)

        if (sum(dpu_0(:))/=0) then
                l_inf = MAXVAL(ABS(dpu(1,1,:,1) - dpu_0)) / MAXVAL(ABS(dpu_0))
        else 
                l_inf = 1.0d0  
        end if      

        !dpu_l_2:

        call trapezoid_integration((dpu(1,1,:,1)-dpu_0)**2,end_val,I_n) 
        call trapezoid_integration(dpu_0**2,end_val,I_d) 

        I_n = I_n/(4*PI)
        I_d = I_d/(4*PI)

        if(I_d/=0)then
                l_2 = SQRT( I_n / I_d )
        else
                l_2 = 1.0d0
        end if

        write (filename, '("dpu_errors.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) phi_max,phi_min,l_inf,l_2
        close(unitn)

        write (filename, '("dpu_sum.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpu(1,1,:,1))
        close(unitn)


      !dpt
        !dpt_phi_max:
        filename = "dpt_phi_max_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        phi_max = (MAXVAL(dpt(1,1,:,1)) - m_0) / delta_phi_0

        !dpt_phi_min:
        filename = "dpt_phi_min_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        phi_min = (MINVAL(dpt(1,1,:,1)) - m_0) / delta_phi_0

        !dpt_l_inf:
        filename = "dpt_l_inf_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpt_0                         !read state 0    
        close(unitn)

        l_inf = MAXVAL(ABS(dpt(1,1,:,1) - dpt_0)) / MAXVAL(ABS(dpt_0))

        !dpt_l_2:

        call trapezoid_integration((dpt(1,1,:,1)-dpt_0)**2,end_val,I_n)
        call trapezoid_integration(dpt_0**2,end_val,I_d)

        I_n = I_n/(4*PI)
        I_d = I_d/(4*PI)

        l_2 = SQRT( I_n / I_d )

        write (filename, '("dpt_errors.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) phi_max,phi_min,l_inf,l_2
        close(unitn)

        write (filename, '("dpt_sum.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpt(1,1,:,1))
        close(unitn)


       !dpk
        !dpk_phi_max:
        filename = "dpk_phi_max_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)

        if(delta_phi_0/=0)then      
                phi_max = (MAXVAL(dpk(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_max = 1.0d0
        end if
    
        !dpk_phi_min:
        filename = "dpk_phi_min_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) delta_phi_0, m_0              !read state 0
        close(unitn)  
   
        if (delta_phi_0/=0) then 
                phi_min = (MINVAL(dpk(1,1,:,1)) - m_0) / delta_phi_0
        else
                phi_min = 1.0d0
        end if

        !dpk_l_inf:
        filename = "dpk_l_inf_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpk_0                         !read state 0    
        close(unitn)

        if (sum(dpk_0(:))/=0) then      
                l_inf = MAXVAL(ABS(dpk(1,1,:,1) - dpk_0))  / MAXVAL(ABS(dpk_0))
        else
                l_inf = 1.0d0
        end if
        

        !dpk_l_2:

        call trapezoid_integration((dpk(1,1,:,1)-dpk_0)**2,end_val,I_n)
        call trapezoid_integration(dpk_0**2,end_val,I_d)

        I_n = I_n/(4*PI)
        I_d = I_d/(4*PI)

        if(I_d/=0)then
                l_2 = SQRT( I_n / I_d )
        else
                l_2 = 1.0d0
        end if
      
        write (filename, '("dpk_errors.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) phi_max,phi_min,l_inf,l_2
        close(unitn)

        write (filename, '("dpk_sum.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpk(1,1,:,1))
        close(unitn)

        write (filename, '("gamma.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) SUM(dpt(1,1,:,1)) + SUM(dpk(1,1,:,1)) + dpphi
        close(unitn)

      
      end if

      !dpphi - gotta do it separately...
      if(filenum==1)then

        !dpphi_l_inf, dpphi_l_2:
        write (filename, '("dpphi_l_inf_0.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpphi
        close(unitn)

        write (filename, '("dpphi_errors.dat")' )
        open(unitn, file=trim(filename), status='replace')
        close(unitn)

        write (filename, '("dpphi_sum.dat")' )
        open(unitn, file=trim(filename), status='replace')
        write(unitn,*) dpphi
        close(unitn)

      else if(filenum>1)then
        !dpphi_l_inf:
        filename = "dpphi_l_inf_0.dat"
        open(unitn, file=trim(filename),status='old')
        read(unitn,*) dpphi_0                         !read state 0    
        close(unitn)

        if(dpphi_0/=0)then
                l_inf = ABS(dpphi - dpphi_0) / ABS(dpphi_0)
        else
                l_inf = 0.0d0
        end if
        !dpphi_l_2:

        I_n = ((dpphi-dpphi_0)**2)/(2*PI)
        I_d = (dpphi_0**2)/(2*PI)

        if(I_d/=0)then
                l_2 = SQRT( I_n / I_d )
        else
                l_2 = 0.0d0
        end if

        write (filename, '("dpphi_errors.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) 0,0,l_inf,l_2
        close(unitn)

        write (filename, '("dpphi_sum.dat")' )
        open(unitn, file=trim(filename), status='old',position='append')
        write(unitn,*) dpphi
        close(unitn)

       
       end if      
      
    end if
  end subroutine diagnostic  

  subroutine trapezoid_integration(array,end_val,integral)
  !Routine for numerical integration of an array using trapezoid rule,
  !adapted from http://f90in15minutes.wikidot.com/numerical-integration
      implicit none
     
      real(KIND=r8),intent(in),dimension(nlev) ::  array
      real(KIND=r8),intent(in) ::  end_val
      real(KIND=r8),intent(out) ::  integral
      real(KIND=r8) :: u,h
      integer :: i

      integral = 0.0

      do i=1,size(array)        
         if ((i.eq.1).or.(i.eq.size(array))) then
            integral = integral+array(i)
         else
            integral = integral+(2.0*array(i))
         end if
      end do

      h=end_val/size(array)
      integral = (h/2.0)*integral

  end subroutine trapezoid_integration

  subroutine diagnostic_old(dpu,dpt,dpk,filenum)
    use spmd_utils,       only: masterproc
    integer           :: unitn,k,sz,dpk_k
    integer,intent(in):: filenum
    character(len=256):: filename
    real(KIND=r8), dimension(np,np,nlev,1),intent(in)::dpu,dpk,dpt
    !real(KIND=r8), dimension(np,np,nlev,1),intent(in)::dpu,dpk
    
    
    if (masterproc) then

      unitn = 8
      sz = size(dpu(1,1,:,1))
      !writing 2-column file for temp state 
      write (filename, '("dpu_", I0.3, ".dat")' )  filenum
      open(unitn, file=trim(filename), status='replace' )
      do k=1,sz
        write(unitn,*) dpu(1,1,k,1)
      end do
      close(unitn)

      sz = size(dpt(1,1,:,1))
      write (filename, '("dpt_", I0.3, ".dat")' )  filenum
      open(unitn, file=trim(filename), status='replace' )
      do k=1,sz
        write(unitn,*) dpt(1,1,k,1)
      end do
      close(unitn)

      sz = size(dpk(1,1,:,1))
      write (filename, '("dpk_", I0.3, ".dat")' )  filenum
      open(unitn, file=trim(filename), status='replace' )
      do k=1,sz
        dpk_k = (1/2.)*dpk(1,1,k,1)
        write(unitn,*) dpk_k
      end do
      close(unitn)


    end if
  end subroutine diagnostic_old

  subroutine remap_TE_forward(ttmp,dp,dp_s,remap_te,filtered,ppm,pqm,phi_inc)
    real (kind=r8), dimension(np,np,nlev),intent(in)   :: dp,dp_s
!    real (kind=r8), dimension(np,np,nlev),intent(in)   :: dp_inv,dp_dry,dp_star_dry,dp_s_inv
!    real (kind=r8), dimension(np,np,nlev)  :: q_test,q_test_s
    real (kind=r8), dimension(np,np,nlev,2),intent(inout) :: ttmp
!    real (kind=r8), dimension(np,np,nlev+1),intent(in) :: phi
!    real(KIND=r8), dimension(nlev+1) :: pint1,pint2,dx1
!    real(KIND=r8), dimension(nlev) :: t0,lnp1,lnp2
    logical,intent(in)          :: remap_te,filtered,phi_inc,ppm,pqm
    real(KIND=r8) :: r_universal
    real(KIND=r8) :: dpphi,kappa
    real(KIND=r8),dimension(np,np,nlev) :: num_tmp,den_tmp
!    type(remapping_CS)                  :: CS
    real, dimension(nlev) :: dp_mom,ttmp_mom,dp_s_mom,ttmp_s_mom,qdp_tmp,qdp_s_tmp


!        ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + &
!                        elem(ie)%state%v(:,:,2,:,np1)**2)/2._r8 + &
!                        elem(ie)%state%t(:,:,:,np1)*cpair !Energy minus phi term
!        if(phi_inc)then
!                        do k=nlev,1,-1
!                           phi(:,:,k) = phi(:,:,k+1)+&
!                                        r_universal*elem(ie)%state%t(:,:,k,np1)*&
!                                        (log(pint1(k+1))-log(pint1(k)))
                           !if(k==nlev)then
                           !     dpphi = pint1(k)*phi(1,1,k) !saving
                           !     p*phi(surface)
                           !end if 
!                           ttmp(:,:,k,1)= ttmp(:,:,k,1) + (pint1(k+1)*phi(:,:,k+1) - &
!                                        pint1(k)*phi(:,:,k))/dp_star_moist(:,:,k) !Energy  
!                        end do

                     !   dpphi = pint1(1)*phi(1,1,1) !- dpphi !(k==1==top)
                     !   p*phi(top) - p*phi(surface)

                     !   if(masterproc)then
                     !     if(i==1)then
                     !           write (filename, '("phi_TE_lag.dat")' )  
                     !           open(unitn, file=trim(filename),
                     !           status='replace' )
                     !           write(unitn,*)
                     !           dpphi,SUM(ttmp(1,1,:,1)*dp_star_moist(1,1,:)),cpair
                     !           close(unitn)
                     !     else 
                     !           write (filename, '("phi_TE_lag.dat")' )  
                     !           open(unitn,file=trim(filename),status='old',position='append'
                     !           )
                     !           write(unitn,*)
                     !           dpphi,SUM(ttmp(1,1,:,1)*dp_star_moist(1,1,:) )
                     !           close(unitn)
                     !     end if
                     !   endif

!                end if

!                elem(ie)%state%Qdp(:,:,:,6,np1_qdp) = q_test(:,:,:)*dp_star_dry(:,:,:)

!                if(i==1)then
!                        call write_data_TE(pint1,ttmp(1,1,:,1),i)
!                end if

!                ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star_moist !E*dp_star_moist 

!                if(filtered)then
!                       call remap1(ttmp,np,1,1,1,dp_star_moist,dp_moist)
!                       !E_rmp*dp_moist
!                elseif(ppm)then
!                       call remap_Q_ppm(ttmp,np,1,1,1,dp_star_moist,dp_moist,1)
!                       call remap_Q_ppm(elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp),np,1,qsize,qsize,dp_star_dry,dp_dry,2)
!                elseif(pqm)then

!                       dp_mom = dp_moist(1,1,:)
!                       dp_s_mom = dp_star_moist(1,1,:)
!                       ttmp_mom = ttmp(1,1,:,1)*dp_inv(1,1,:)
!                       ttmp_s_mom = ttmp(1,1,:,1)*dp_s_inv(1,1,:)

!                       call remapping_core_h(CS,nlev,dp_s_mom,ttmp_s_mom,nlev,dp_mom,ttmp_mom)

!                       elem(ie)%state%t(1,1,:,np1) = ttmp_mom

!                else
                        !call remap1_nofilter(ttmp,np,2,dp_star_moist,dp_moist)   
!                end if

!                if(.not.pqm)then
!                        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)*dp_inv  !E_rmp (as t)
!                        q_test_s(:,:,:) = elem(ie)%state%Qdp(:,:,:,6,np1_qdp)/dp_dry(:,:,:)
!                end if
!                if(masterproc)then
!                          if(i==1)then
!                                write (filename, '("TE_eul.dat")' )
!                                open(unitn, file=trim(filename),status='replace')
!                                write(unitn,*) SUM(elem(ie)%state%t(1,1,:,np1)*dp_moist(1,1,:)),cpair
!                                close(unitn)
!                          else
!                                write (filename, '("TE_eul.dat")' )
!                                open(unitn,file=trim(filename),status='old',position='append')
!                                write(unitn,*) SUM(elem(ie)%state%t(1,1,:,np1)*dp_moist(1,1,:) )
!                                close(unitn)
!                          end if
!                        endif

  end subroutine remap_TE_forward


end module hack_vert_rmp
