!----OpenCFD-EC 2d--------------------------------------------------------------
!  Copyright by Li Xinliang, LHD, Institute of Mechanics, lixl@imech.ac.cn                 
!  A multi-block Finite Volume Method Navier-Stokes solver
!  Ref. J. Blazek's Book: Computational Fluid Dynamics: principle and application 
!-

!------------------------------------------------------------------------------------------
! 主程序 主程序 主程序
!-----------------------------------------------------------------------------------------
    program main
    use Global_Var
    use Com_ctrl
    implicit none
    integer  nMesh
    integer :: nm,mblock
    call DATE_AND_TIME(VALUES=initialtime)
    call MPI_init(ierr)
    call MPI_comm_rank(mpi_comm_world,myid,ierr)
    call MPI_comm_size(mpi_comm_world,nprocs,ierr)
    if (myid == 0)then
        write(*,*)" ========================================= "
        write(*,*)" ===   The initial calculation time:   === "
        write(*,"(1X,A5,I4,A2,5(I2,A2),I3,A4,A5)")" === ", &
                    & initialtime(1),"年",initialtime(2),"月",initialtime(3),"日", &
                    & initialtime(5),"时",initialtime(6),"分",initialtime(7),"秒", &
                    & initialtime(8),"毫秒"," === "
        write(*,*)" ========================================= "
    endif
    call read_parameter                     ! 读取流动参数及控制信息
    call Init                               ! 初始化，创建数据结构; 读入几何及物理信息

    if(Iflag_init .eq. 0) then
    
        call init_flow_zero                      ! 从零流场算起 

    else
       call init_flow_read
    endif

    if(myid==0)then 
        print*, " Initialize OK ......"
        print*, " Start ......"
    end if

    do while(Mesh(1)%Kstep .lt. 100000)     
        call   Hybrid_Time_advance_LU_SGS(1)
        if(mod(Mesh(1)%Kstep, Kstep_show).eq.0)  call DETFP_res
        
        if(myid==0)then
            if(mod(Mesh(1)%Kstep, Kstep_show).eq.0)    call output_Res_GKUA(1)         ! 打印残差 (最密网格)
            if(mod(Mesh(1)%Kstep,Kstep_save) .eq. 0)   call output (1)           ! 输出流场 (最密网格)  
        end if
        if(mod(Mesh(1)%Kstep+500, 1000) .eq. 0)  call backup_Macro_write(1)
        if(mod(Mesh(1)%Kstep, 1000) .eq. 0)  call backup_Macro_write(2)
   
    enddo
    


    call MPI_GROUP_FREE(new_group,IERR)
    call MPI_COMM_FREE(new_comm,IERR)
    call MPI_FINALIZE(IERR)
    stop
      
    end program 
    
!----------------------------------------------------------------------------------------------



!  -----------------------读取流动参数及控制变量----------------------
    subroutine read_parameter
    use Global_var
    use Com_ctrl
    implicit none
    integer:: i,k,kbody_tmp(6)
    character*80 fname
    
    fname="input/"//"Quadrature_Para.dat"
    open(99,FILE=trim(fname),STATUS='OLD')
    rewind(99)
    read(99,*)
    read(99,*)
    read(99,*)Quadrature_Mode
    if(Quadrature_Mode==GKUA_GH)then
        do i=1,3
            read(99,*)
        end do
        read(99,*) k
        if(k==1)then
            NVI=7;NVJ=7
        elseif(k==2)then
            NVI=8;NVJ=8
        elseif(k==3)then
            NVI=16;NVJ=16
        end if
        NVIt=2*NVI
        NVJt=2*NVJ
    elseif(Quadrature_Mode==GKUA_GL)then        
        do i=1,6
            read(99,*)
        end do
        
        read(99,*) DVX,NVI,VXdown,VXup
        read(99,*)
        read(99,*) DVY,NVJ,VYdown,VYup
        read(99,*)
        read(99,*) KGLN 
        NVIt=KGLN*NVI
        NVJt=KGLN*NVJ
    else if(Quadrature_Mode==GKUA_NC)then       
        do i=1,13
            read(99,*)
        end do
        read(99,*) DVX,NVI,VXdown,VXup
        read(99,*)
        read(99,*) DVY,NVJ,VYdown,VYup
        NVIt=NVI
        NVJt=NVJ
    end if
    close(99)
        
    open(99,file="input/"//"control.in")
    read(99,*)
    read(99,*)
    read(99,*) Mole_Mass, gamma, KHS
    read(99,*)
    read(99,*)  ALPgas, OMGgas, vis_ref ,  TREF
    read(99,*)
    read(99,*)  Sgas , vis_ref,TREF
    read(99,*)
    read(99,*)  DSLJ, TSLJ, vis_ref,TREF
    read(99,*)
    read(99,*) Ma, ND_inf, T_inf, Twall, P_outlet, AoA, Lref
   	read(99,*)
    read(99,*)
	read(99,*)If_viscous,Time_Method,Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction
    read(99,*)
    read(99,*)
    read(99,*) Iflag_init,Iflag_local_dt,dt_global,CFL,dtmax,dtmin,t_end,Kstep_save,Kstep_show
    close(99)
    NVAR=4  !!!不考虑湍流，求解变量4
!-----------------------------------------------------------------------------
    continue
    end subroutine read_parameter
    


