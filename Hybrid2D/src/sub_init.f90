! Copyright By Li Xinliang, Code by Li Xinliang
! 初始化：包括创建数据结构及赋初值
! 对于多重网格，根据上级网格的信息，创建各级网格
!---------------------------------------------------------------------------------
!------------------------------------------------------------------------------     
    subroutine init
    use Global_var
    use Com_ctrl
    implicit none
    real(sp) ,external::Omiga22
 !--------------------------------------------------------------------
    allocate( Mesh(1) )             ! 主数据结构： “网格” （其成员是“网格块”）
    eps=1.E-6
    x1s=1.-OMGgas
!---- PR is prlant number:Pr=Mu.Cp/k,which is equal to 2/3 for monatomic gas.
    PR=gamma/( 0.25*(9.*gamma-5.) )
    RCON=Boltz/(Matom*Mole_Mass)
    if(KHS==VHS_Mode)then
        vis_inf=vis_ref*(T_inf/TREF)**OMGgas
    else if(KHS==Sutherland_Mode)then
        vis_inf=vis_ref*(T_inf/TREF)**1.5*(Sgas+TREF)/(Sgas+T_inf)
    else if (KHS==LJ_Mode)then
        vis_inf=5.0*sqrt(PI*Boltz*Matom*Mole_Mass*T_inf)/(16.0*PI*DSLJ*DSLJ*Omiga22(T_inf/TSLJ))
    else
        if(myid==0)then
            write(*,*)"Viscosity Model error"
            call MPI_FINALIZE(IERR)
            stop
        end if
    end if 
 
    Re=ND_inf*Matom*Mole_Mass*Ma*sqrt(gamma*RCON*T_inf)*Lref/vis_inf
    AoA=AoA*PI/180.d0
    Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma) 
    Cp=Cv*gamma 
     
    if(myid==0)then
        write(*,*)"Re number of this case is ",Re,RCON
    end if
      
    if(KHS==VHS_Mode)then
        p_inf=ND_inf*Matom*Mole_Mass*T_inf*RCON
        MFP_inf=(4.0*ALPgas*(5.-2.*OMGgas)*(7.-2.*OMGgas))/(5.0*(ALPgas+1)*(ALPgas+2))*&
            &(vis_ref/(2.0*sqrt(Pi)*p_inf))*sqrt(2.0*RCON*T_inf)
    else if(KHS==Sutherland_Mode)then
        p_inf=ND_inf*Matom*Mole_Mass*T_inf*RCON
    else if (KHS==LJ_Mode)then
        vis_inf=5.0*sqrt(PI*Boltz*Matom*Mole_Mass*T_inf)/(16.0*PI*DSLJ*DSLJ*Omiga22(T_inf/TSLJ))
        p_inf=ND_inf*Matom*Mole_Mass*T_inf*RCON
        Dgas_inf=1.016*DSLJ*DSLJ*Omiga22(T_inf/TSLJ)
        MFP_inf=1.0/(sqrt(2.0)*Pi*Dgas_inf*ND_inf)
    else
        if(myid==0)then
            write(*,*)"Viscosity Model error"
            call MPI_FINALIZE(IERR)
            stop
        end if
    end if  
    
    Kn_inf=MFP_inf/Lref
    Rdof=(5.-3.*gamma)/(gamma-1.)
    if (myid == 0)then
        write(*,*)"KHS,Ma,Re_inf,Kn_inf,Rcon",KHS,Ma,Re,Kn_inf,Rcon
    endif
    cfc=Re/(Ma*sqrt(2.0*gamma))
    cpr=(1.-PR)/(5./2.)
!--------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------

    call Creat_Mesh            ! 创建最细的网格 (从网格文件Mesh2d.dat)
    call read_bcin              ! 读网格连接信息 (bc2d.in)
    call Update_coordinate_buffer_onemesh(1)           ! 利用连接信息，给出虚网格的坐标
    call set_control_para                   ! 设定各重网格上的控制信息（数值方法、通量技术、湍流模型、时间推进方式）
    call Vel_Phase_det
    call Alloc_GKUA_FRD

 

    end  subroutine init
    
!---------------------------------------------------------




!--------------------------------------------------------------------------------------
!   创建数据结构： 最细网格 （储存几何量及守恒变量）
    subroutine Creat_Mesh
    use Global_var
    use Com_ctrl
    implicit none
    integer,allocatable,dimension(:):: NI,NJ
    integer:: NB,m,nx,ny,i,j,k
    Type (Block_TYPE),pointer:: B
    real(sp) :: dx,dy,ztemp
    maxnode = 0
    if(myid==0)then
        print*, "-------------------------------------"
        print*, "read generic.dat"
    end if
    
    open(99,file="input/"//"generic.dat")
    read(99,*) NB   ! Block number
    Mesh(1)%Num_Block=NB
    Mesh(1)%Num_Cell=0
!    时间步、时间
    Mesh(1)%Kstep=0
    Mesh(1)%tt=0.d0
    allocate(Mesh(1)%Block(NB))
    allocate( NI(NB),NJ(NB) )
    do m=1,NB
        read(99,*) NI(m), NJ(m)
        maxnode = max(NI(m), NJ(m),maxnode)
        B => Mesh(1)%Block(m)
        B%Block_no=m
        B%nx=NI(m); B%ny=NJ(m) 
        nx=B%nx ; ny= B%ny
        Mesh(1)%Num_Cell=Mesh(1)%Num_Cell+(nx-1)*(ny-1)

!   申请内存  几何量: (x,y) 节点坐标； (x1,y1)网格中心坐标; s0 控制体体积； 
  
        allocate(B%vol(nx,ny),B%si(nx,ny),B%sj(nx,ny),B%ni1(nx,ny),B%ni2(nx,ny),B%nj1(nx,ny),B%nj2(nx,ny)) 
        allocate(B%x(1-LAP:nx+LAP,1-LAP:ny+LAP), B%y(1-LAP:nx+LAP,1-LAP:ny+LAP))   
        allocate(B%x1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%y1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
!   
        allocate(B%U(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))    !守恒变量 (4个流体变量)
        allocate(B%U1(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))    !守恒变量 (4个流体变量)
         allocate(B%U2(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))    !守恒变量 (4个流体变量)
         allocate(B%U3(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))    !守恒变量 (4个流体变量)
        allocate(B%U0(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Un(NVAR,0:nx,0:ny))             ! 上一时间步的值   
        allocate(B%Res(NVAR,1:nx-1,1:ny-1))        ! 残差
        allocate(B%Res_i(NVAR,1:nx-1,1:ny-1))        ! 残差
        allocate(B%dt(1:nx-1,1:ny-1))              ! 时间步长    
        allocate(B%Lci(nx,ny),B%Lcj(nx,ny),B%Lvi(nx,ny),B%Lvj(nx,ny))        ! 谱半径    
        allocate(B%dU(NVAR,nx,ny))                 ! =U(n+1)-U(n), LU-SGS中使用    
        
        allocate(B%TaoNSi_in(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TaoNSj_in(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Qi_in(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%Qj_in(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        
        allocate(B%TaoNSi_out(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TaoNSj_out(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Qi_out(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%Qj_out(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        
        allocate(B%TaoNS_out(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Q_out(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        
        allocate(B%TaoNSi_G(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TaoNSj_G(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Qi_G(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%Qj_G(2,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Qx(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%Qy(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Taoxx(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%Taoxy(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),&
            &B%Taoyy(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
        allocate(B%Amu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), B%Amu_t(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))




!   初始化	   
        B%x1(:,:)=0.d0; B%y1(:,:)=0.d0; B%vol(:,:)=0.d0   
         B%U0(1,:,:)=1.d0; B%U0(2,:,:)=0.d0; B%U0(3,:,:)=0.d0; B%U0(4,:,:)=1.d0
        B%U(1,:,:)=1.d0; B%U(2,:,:)=0.d0; B%U(3,:,:)=0.d0; B%U(4,:,:)=1.d0
        B%U1(1,:,:)=1.d0; B%U1(2,:,:)=0.d0; B%U1(3,:,:)=0.d0; B%U1(4,:,:)=1.d0
        B%U2(1,:,:)=1.d0; B%U2(2,:,:)=0.d0; B%U2(3,:,:)=0.d0; B%U2(4,:,:)=1.d0
        
        B%Un(1,:,:)=1.d0; B%Un(2,:,:)=0.d0; B%Un(3,:,:)=0.d0; B%Un(4,:,:)=0.d0
        B%Lci(:,:)=0.d0; B%Lcj(:,:)=0.d0; B%Lvi(:,:)=0.d0; B%Lvj(:,:)=0.d0
        B%dU(:,:,:)=0.d0
        B%Amu(:,:)=0.d0;  B%Amu_t(:,:)=0.d0
        ! B%DeltU=0.d0
        B%Res=0.d0
        B%dt=0.d0
        B%TaoNSi_in=0.d0;B%TaoNSj_in=0.d0
        B%TaoNSi_out=0.d0;B%TaoNSj_out=0.d0
        B%TaoNSi_G=0.d0;B%TaoNSj_G=0.d0
        B%Qi_in=0.d0;B%Qj_in=0.d0
        B%Qi_out=0.d0;B%Qj_out=0.d0
        B%Qi_G=0.d0;B%Qj_G=0.d0
        B%QX=0.0;B%QY=0.0;B%Taoxx=0.0;B%Taoxy=0.0;B%Taoyy=0.0   
        read(99,*) ((B%x(i,j),i=1,nx),j=1,ny),  ((B%y(i,j),i=1,nx),j=1,ny),((ztemp,i=1,nx),j=1,ny)
        B%x=B%x/Lref
        B%y=B%y/Lref
!   计算控制体的体积 
!   控制体中心坐标需要等获得Ghost Cell 区信息后计算       
        do j=1,ny-1
        do i=1,nx-1
            B%vol(i,j) =abs((B%x(i,j)-B%x(i+1,j+1))*(B%y(i+1,j)-B%y(i,j+1)) -  &
            (B%x(i+1,j)-B%x(i,j+1))*(B%y(i,j)-B%y(i+1,j+1)) )*0.5d0           ! 控制体体积（面积）
        enddo
        enddo

!     几何量 （边长，法方向）      
        do j=1,ny-1
        do i=1,nx
            dx=B%x(i,j+1)-B%x(i,j)
            dy=B%y(i,j+1)-B%y(i,j)
            B%si(i,j)=sqrt(dx*dx+dy*dy)
            B%ni1(i,j)=dy/B%si(i,j); B%ni2(i,j)=-dx/B%si(i,j)   ! normal vector at (i,j) or (I-1/2,J) 
        enddo
        enddo
        

        do j=1,ny
        do i=1,nx-1
            dx=B%x(i+1,j)-B%x(i,j)
            dy=B%y(i+1,j)-B%y(i,j)
            B%sj(i,j)=sqrt(dx*dx+dy*dy)     ! length 
            B%nj1(i,j)=-dy/B%sj(i,j); B%nj2(i,j)=dx/B%sj(i,j)     ! normal vector at i, j+1/2
        enddo
        enddo
   
    enddo
    nmt = maxnode+nad
    close(99) 
    deallocate(NI,NJ)

    if(myid==0)then
        print*, "read generic.dat OK"
    end if
    
    end subroutine Creat_Mesh


!----Mesh control message (bc2d.in)------------------------------------------
    subroutine read_bcin 
    use Global_Var
    use Com_ctrl
    implicit none
    integer:: NB,m,ksub,nx,ny,temp,ni_temp,nj_temp
    integer :: ibegin,iend,jbegin,jend,mode_check_ns,mode_check_es
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc
    if(myid==0)then
        print*, "read bc2d.in ......"
    end if
    mode_check_ns=0
    mode_check_es=0
    open(88,file="input/"//"generic.inp")
    read(88,*)
    read(88,*) NB
    
    if(NB .ne. Mesh(1)%Num_Block) then
        if(myid==0)then
            print*, "Error!  Block number in bc2d.in is not equal to that in Mesh2d.dat !"
        end if
        stop
    endif
    khalf=0;kall=1
    do m=1,NB
        B => Mesh(1)%Block(m)
        nx=B%nx ; ny= B%ny
        read(88,*) ni_temp,nj_temp
        if((ni_temp/=nx).or.(nj_temp/=ny))then
            if(myid==0)then
                write(*,*)"Error! demensions of blocks is not equal to that in mesh2d.dat"
            end if
            stop
        end if
        read(88,*) B%solver
        if(B%solver==NS_solver)then
            mode_check_ns=1
        else if (B%solver==GKUA)then
            mode_check_es=1
        endif
        read(88,*) B%subface   
        allocate (B%bc_msg(B%subface))
        do ksub=1, B%subface
            Bc => B%bc_msg(ksub)
            Bc%f_no=ksub
            read(88,*)   Bc%ist, Bc%iend, Bc%jst, Bc%jend,  Bc%neighb
            if(Bc%ist==Bc%iend.and.Bc%ist==1)then
                Bc%face=1
            else if(Bc%ist==Bc%iend.and.Bc%ist==nx)then
                Bc%face=3
            else if(Bc%jst==Bc%jend.and.Bc%jst==1)then
                Bc%face=2
            else if(Bc%jst==Bc%jend.and.Bc%jst==ny)then
                Bc%face=4
            end if
            Bc%neighb=INT(-1*Bc%neighb)
            if(Bc%neighb==BC_Symmetry_or_slidewall)then
                khalf=1;kall=0
            end if   
            if(Bc%neighb==1)then
                read(88,*) Bc%ist_neighb, Bc%iend_neighb, Bc%jst_neighb, Bc%jend_neighb,Bc%neighb
            end if
            ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
            if(( Bc%neighb== BC_Wall).or.(Bc%neighb== BC_Farfield))then
                temp=Bc%iend-Bc%ist+Bc%jend-Bc%jst
                allocate(BC%F_Wall(2,2))
                if(Bc%face .eq. 1 ) then 
                    allocate(BC%Swall( 19,1,jbegin:jend))
                    allocate(BC%Swall0(19,1,jbegin:jend))
                    allocate(BC%Qnwall(2,1,jbegin:jend))
                    allocate(BC%Qnwall0(2,1,jbegin:jend))
                else if(Bc%face .eq. 2 ) then 
                    allocate(BC%Swall( 19,ibegin:iend,1))
                    allocate(BC%Swall0(19,ibegin:iend,1))
                    allocate(BC%Qnwall(2,ibegin:iend,1))
                  allocate(BC%Qnwall0(2,ibegin:iend,1))
             else if(Bc%face .eq. 3 ) then 
                  allocate(BC%Swall( 19,1,jbegin:jend))
                  allocate(BC%Swall0(19,1,jbegin:jend))
                  allocate(BC%Qnwall(2,1,jbegin:jend))
                  allocate(BC%Qnwall0(2,1,jbegin:jend))
             else if(Bc%face .eq. 4 ) then 
                  allocate(BC%Swall( 19,ibegin:iend,1))
                  allocate(BC%Swall0(19,ibegin:iend,1))
                  allocate(BC%Qnwall(2,ibegin:iend,1))
                  allocate(BC%Qnwall0(2,ibegin:iend,1))
             end if

      
         BC%F_Wall=0.0;BC%Swall=0.0;BC%Swall0=0.0;BC%Qnwall=0.0;BC%Qnwall0=0.0
        
      end if
      
     enddo
     
  
     
    enddo
    close(88)
     
    call Get_neighb_info
    if(myid==0)then
        print*, "read bc2d.in OK"
    end if
    end  subroutine read_bcin

    subroutine  Alloc_GKUA_FRD
    use Global_Var
    use Com_ctrl
    implicit none
    integer:: NB,m,ksub,nx,ny,temp,ni_temp,nj_temp
    integer :: ibegin,iend,jbegin,jend
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc
    
    do m=1,Mesh(1)%Num_Block
        B => Mesh(1)%Block(m)
            
        nx=B%nx ; ny= B%ny
        if(B%solver==GKUA)then  !@CARDC
              nvt = max(nvit,nvjt)
              allocate(B%FRD(NRD,-3:nx+3,-3:ny+3,I1MAX,J1MAX)) 
               allocate(B%FRD0(NRD,-3:nx+3,-3:ny+3,I1MAX,J1MAX)) 
              allocate(B%FRD00(NRD,-3:nx+3,-3:ny+3,I1MAX,J1MAX))
              allocate(B%FRD000(NRD,-3:nx+3,-3:ny+3,I1MAX,J1MAX))
              allocate(B%ATSVX(-3:nx+3,-3:ny+3),B%ATSVY(-3:nx+3,-3:ny+3)) 
              allocate(B%DUFRD(NRD,-3:nx+3,-3:ny+3),B%RHSFRD(NRD,-3:nx+3,-3:ny+3)) 

              allocate(B%S_GKUA(Nsum,-3:nx+3,-3:ny+3),B%BNIRHS(6,-3:nx+3,-3:ny+3),B%S00_GKUA(Nsum,-3:nx+3,-3:ny+3))
               allocate(B%S1_GKUA(Nsum,-3:nx+3,-3:ny+3))
               allocate(B%HOT_GKUA(6,-3:nx+3,-3:ny+3),B%HOQ_GKUA(4,-3:nx+3,-3:ny+3)) 
              allocate(B%S0_GKUA(Nsum,-3:nx+3,-3:ny+3),B%CF00(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%CF(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%CF0(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
              allocate(B%vi(NVIt),B%vj(NVJt)) 
              allocate(B%wi(NVIt),B%wj(NVJt))
              allocate(B%FBAVX(nvt),B%FBAVY(nvt))
     
     

              allocate(B%QXES(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%QYES(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%TaoXXES(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TaoXYES(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%TAOYYES(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TAOZZES(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
              allocate(B%TaoXXES0(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TaoXYES0(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%TAOYYES0(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TAOZZES0(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%TLaXX(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TLaXY(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%TLaYY(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%TLaZZ(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
              allocate(B%GKUA_P(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%GKUA_RON(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%GKUA_T(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%GKUA_Tpt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%GKUA_U(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%GKUA_V(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%GKUA_Ppt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%GKUA_RONpt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              allocate(B%GKUA_Upt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%GKUA_Vpt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
              B%FRD=0.0;B%FRD00=0.0;B%FRD000=0.0;B%CF=0.0;B%vi=0.0;B%vj=0.0;B%wi=0.0;B%wj=0.0
              B%ATSVX=0.0;B%ATSVY=0.0;B%DUFRD=0.0;B%RHSFRD=0.0
              B%S_GKUA=0.0;B%S0_GKUA=0.0;B%S1_GKUA=0.0
              B%TaoXXES=0.0;B%TaoYYES=0.0;B%TaoZZES=0.0;B%TaoXYES=0.0
             
              B%S_GKUA=0.0;B%BNIRHS=0.0
              B%FBAVX=0.0;B%FBAVY=0.0
      
     end if
    end do
    

    end subroutine Alloc_GKUA_FRD
    

! 用来流初始化； 
    subroutine init_flow_zero
    use Global_var
    use Com_ctrl
    implicit none
    real(sp)  ,external::Omiga22
    real(sp) :: d0,u0,v0,p0,T0,tmp,miu_tmp,T_tmp,T_tmp2
    integer:: i,j,m,step,nMesh
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
	
!-------------------------------------------------------------------
    d0=1.d0; u0=1.d0*cos(AoA); v0=1.d0*sin(AoA) ; p0=1.d0/(gamma*Ma*Ma)
    Twall=Twall/T_inf

    T00=1.0
    R00=d0
    u00=sqrt(gamma/2.)*Ma*u0
    v00=sqrt(gamma/2.)*Ma*v0
      
    MP=>Mesh(1)   
    Init_mass=0.0
    do m=1,MP%Num_Block
        B => MP%Block(m)    
        do j=-1,B%ny+1
        do i=-1,B%nx+1
          
            B%U(1,i,j)=d0
            B%U(2,i,j)=d0*u0
            B%U(3,i,j)=d0*v0
            B%U(4,i,j)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0)
        enddo
        enddo
        if(B%solver==GKUA)then
            do j=-1,B%ny+1
            do i=-1,B%nx+1
                B%S_GKUA(1,i,j)=d0 !n 
                B%S_GKUA(2,i,j)=u0*sqrt(gamma/2.)*Ma !u
                B%S_GKUA(3,i,j)=v0*sqrt(gamma/2.)*Ma  !v
                B%S_GKUA(4,i,j)=T00  !Ttra
                if(KHS==VHS_Mode)then
                    if(x1s .le. 1.e-8) then
                        B%CF(i,j)=cfc*B%S_GKUA(1,i,j)
                    else
                        B%CF(i,j)=cfc*B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j) **x1s
                    endif
                else if(KHS==Sutherland_Mode)then 
                    Miu_tmp=(vis_ref/vis_inf)*(B%S_GKUA(4,i,j)*T_inf/Tref)**(3.0/2.0)*(Sgas+Tref)/(Sgas+B%S_GKUA(4,i,j)*T_inf)
                    B%CF(i,j)=cfc*B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j)/Miu_tmp
                else if (KHS==LJ_Mode)then
                    tmp=Omiga22(B%S_GKUA(4,i,j)*T_inf/TSLJ)/Omiga22(T_inf/TSLJ)
                    Miu_tmp=sqrt(B%S_GKUA(4,i,j))/tmp
                    B%CF(i,j)=cfc*B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j)/Miu_tmp
                else 
                    if(myid==0)then
                        write(*,*)"Viscosity Model error"
                        stop
                    end if
                end if
            enddo
            enddo 
            B%CF0= B%CF
            do j=1,B%ny-1
            do i=1,B%nx-1
                Init_mass=Init_mass+B%S_GKUA(1,i,j)*B%vol(i,j)
            end do
            end do
        end if
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call GKUA_sudufenbu
    call GET_SYM_NUM
    call Boundary_condition_onemesh(1)                   ! 边界条件 （设定Ghost Cell的值）     
    call update_buffer_onemesh(1)                        ! 同步各块的交界区
    end subroutine init_flow_zero 
    
 


    
    subroutine Vel_Phase_det
    use Global_var
    use Com_ctrl
    implicit none
    Type (Block_TYPE),pointer:: B
    Type (BC_MSG_TYPE),pointer:: Bc
    integer:: nMesh,mBlock,ksub
    integer :: inprocs,PES,MM,J,I,totblk
    integer :: ix0,iy0,iv,jv,iv1,jv1,memsize,id,id_x
    XDIM=0;YDIM=0
    totblk=Mesh(1)%Num_Block
    if(khalf.EQ.1)then
        if(mod(nvit*nvjt,nprocs).NE.0)then
            if(myid==0)write(*,*)"khalf.EQ.1, why 'mod(nvit*nvjt,nprocs).NE.0'?"
            stop
        endif
        YDIM=NINT(sqrt(nprocs*1.))-1
        do inprocs=1,nprocs-YDIM  !while(mod(nprocs,XDIM).EQ.0)
            YDIM=YDIM+1
            XDIM=nprocs/YDIM
            PES = XDIM*YDIM
           !write(*,*) 0,XDIM,YDIM
            if(mod(nprocs,YDIM).EQ.0 .and. &
           &  mod(nprocs,XDIM).EQ.0 .and. &
           &  mod(nvit,XDIM).EQ.0 .and. &
           &  mod(nvjt,YDIM).EQ.0 .and. &
           &      nprocs .EQ. PES)then
                 exit 
           else
                cycle 
           endif
        enddo
        if( mod(nvit,XDIM).NE.0 .or. mod(nvjt,YDIM).NE.0 )then
            XDIM=NINT(sqrt(nprocs*1.))-1
            do inprocs=1,nprocs-XDIM  
                XDIM=XDIM+1
                YDIM=nprocs/XDIM
                PES = XDIM*YDIM
          
                if(mod(nprocs,YDIM).EQ.0 .and. &
                    &  mod(nprocs,XDIM).EQ.0 .and. &
                    &  mod(nvit,XDIM).EQ.0 .and. &
                    &  mod(nvjt,YDIM).EQ.0 .and. &
                    &      nprocs .EQ. PES)then
                    exit 
                    else
                        cycle 
                    endif
            enddo
        endif
    else
        YDIM=NINT(sqrt(nprocs*1.))-1
        do inprocs=1,nprocs-YDIM  
            if(YDIM.GT.0 .and. mod(nprocs,YDIM).EQ.0)exit
            YDIM=YDIM+1
        enddo
        XDIM=nprocs/YDIM
    endif
    if( khalf.EQ.1 )then
        if((nprocs.NE.1) .and. (mod(nvit,XDIM).NE.0 .or. mod(nvjt,YDIM).NE.0 ))then
            if (myid == 0)then
                write(*,*)"nvit,XDIM,nvjt,YDIM",nvit,XDIM,nvjt,YDIM
                write(*,*)"Some bdtype is symmetry, please make sure mod(nvit,XDIM)=0 and mod(nvjt,YDIM)=0!!! "
            endif
            stop
        endif
    endif
    PES = XDIM*YDIM
    if(nprocs .ne. PES)then
        if (myid == 0) then
            print *,'ERROR NUMBER OF PROCESS'
        endif
        stop
    endif
    I1MAX = (nvit+(XDIM-1)) / XDIM
    J1MAX = (nvjt+(YDIM-1)) / YDIM
    nvt = max(nvit,nvjt)
    myid_x=mod(myid,XDIM)
    myid_y=myid/XDIM
    GKUA_IST=1
    MM=MOD(NVIt,XDIM)
    IF (myid_x.LT.MM) THEN
        GKUA_IEND=NVIt/XDIM+1
        GKUA_IPT=myid_x*GKUA_IEND
    ELSE
        GKUA_IEND=NVIt/XDIM
        GKUA_IPT=myid_x*GKUA_IEND+MM
    ENDIF
    GKUA_JST=1
    MM=MOD(NVJt,YDIM)
    IF (myid_y.LT.MM) THEN
        GKUA_JEND=NVJt/YDIM+1
        GKUA_JPT=myid_y*GKUA_JEND
    else
        GKUA_JEND=NVJt/YDIM
        GKUA_JPT=myid_y*GKUA_JEND+MM
    endif
    !!!寻找关于X轴对称的进程（对YDIM能否被2整除无要求）
    do j=1,YDIM
    do i=(j-1)*XDIM,j*XDIM-1
        if(myid==i)next= i+(YDIM-2*j+1)*XDIM !(YDIM-j)*XDIM+(YDIM-2*j+1)*XDIM
    enddo
    enddo
    nextwall=nprocs-1-myid
    !当mod(YDIM,2)==1时，myid在[CPUst,CPUed]之间时为本进程对称
    CPUst=(YDIM-1)/2*XDIM
    CPUed=(YDIM+1)/2*XDIM-1
    ntyptmp=myid*totblk*GKUA_IEND*GKUA_JEND
    allocate(rank1(0:XDIM-1,YDIM),rank2(YDIM))
 !****for 'subroution calmonitor' to print ***********
    ix0=nvit/2
    iy0=nvjt/2
    prt_id = -1
    DO jv=GKUA_JST,GKUA_JEND
        jv1=jv+GKUA_JPT
    DO iv=GKUA_IST,GKUA_IEND
        iv1=iv+GKUA_IPT
        if ((iv1.eq.ix0) .and. (jv1.eq.iy0)) then
            prt_id=myid
            prt_x=iv
            prt_y=jv
            print *,"iv,jv,prt_id=",iv,jv,prt_id
            goto 1000
        endif
    ENDDO
    ENDDO
1000     call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    memsize = 0
    do mBlock=1,Mesh(1)%Num_Block
        B => Mesh(1)%Block(mBlock)
        memsize = (B%nx+4)*(B%ny+4)*2*6*NRD/1024*I1MAX*J1MAX/1024 + memsize
    enddo 
    if (myid == 0) then
        print *,''
        print *,'FRD allocate memory:', memsize,' (MB) in every cpu.'
        print *,'myid=',myid,'nprocs=',nprocs
        print *,''
    endif
    if (myid.eq.0) then
        do i=0,XDIM-1
            j=1
            do id=0,nprocs-1
                id_x=mod(id,XDIM)
                if (id_x .eq. i) then
                    rank1(i,j)=id
                    j=j+1
                endif
            enddo
        enddo
    endif
        
    call MPI_bcast(rank1,XDIM*YDIM,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_COMM_GROUP(MPI_COMM_WORLD,orig_group,ierr)
    do i=0,XDIM-1
        do j=1,YDIM
            rank2(j)=rank1(i,j)
        enddo
        if (myid_x .eq. i) then
            call MPI_GROUP_INCL(orig_group,YDIM,rank2,new_group,ierr)
        endif
    enddo
    call MPI_COMM_CREATE(MPI_COMM_WORLD,new_group,new_comm,IERR)
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
     
      
    end subroutine Vel_Phase_det
    
    subroutine Get_neighb_info
    use Global_Var
    use Com_ctrl
    implicit none
    integer:: NB,m,ksub,nx,ny,temp,ni_temp,nj_temp,ksub_2,flag
    integer :: ibegin,iend,jbegin,jend,ibegin_neighb,iend_neighb,jbegin_neighb,jend_neighb
    Type (Block_TYPE),pointer:: B,B_neighb
    TYPE (BC_MSG_TYPE),pointer:: Bc,Bc_neighb
    
    
    do m=1,Mesh(1)%Num_Block
        B => Mesh(1)%Block(m)
        do ksub=1, B%subface
            Bc => B%bc_msg(ksub)
            if(Bc%neighb.gt.0)then
                temp=Bc%neighb
                ibegin_neighb=Bc%ist_neighb; iend_neighb=Bc%iend_neighb
                jbegin_neighb=Bc%jst_neighb; jend_neighb=Bc%jend_neighb  
                    B_neighb=>Mesh(1)%Block(temp)
                    flag=0
                    do ksub_2=1,B_neighb%subface
                        Bc_neighb => B_neighb%bc_msg(ksub_2)
                        ibegin=Bc_neighb%ist; iend=Bc_neighb%iend; jbegin=Bc_neighb%jst; jend=Bc_neighb%jend
                        if((ibegin==ibegin_neighb).and.(iend==iend_neighb).and.&
                            &(jbegin==jbegin_neighb).and.(jend==jend_neighb))then
                            Bc%subface=ksub_2
                            Bc%orient=2
                            flag=1
                            else if((ibegin==iend_neighb).and.(iend==ibegin_neighb).and.&
                            &(jbegin==jbegin_neighb).and.(jend==jend_neighb))then
                            Bc%subface=ksub_2
                            Bc%orient=0
                            flag=1
                            else if((ibegin==ibegin_neighb).and.(iend==iend_neighb).and.&
                            &(jbegin==jend_neighb).and.(jend==jbegin_neighb))then
                            Bc%subface=ksub_2
                            Bc%orient=0  
                            flag=1
                       end if
                       
                        
                    end do
                    if(flag==0)then
                        if(myid==0)then
                            write(*,*)"can not find  the correct boundary for inner bound. condition"
                            stop
                        end if
                    end if
                    
                    
            end if
            
        end do
        
    end do
    
    
    end subroutine Get_neighb_info
    


        !------------------------------------------------------------------    
!  设定各重网格上的控制信息
    subroutine set_control_para
    use Global_var
    implicit none
    integer nMesh
    TYPE (Mesh_TYPE),pointer:: MP
    MP=>Mesh(1)            
    MP%Iflag_turbulence_model=Iflag_turbulence_model
    MP%Iflag_Scheme=Iflag_Scheme
    MP%IFlag_flux=IFlag_flux
    MP%IFlag_Reconstruction=IFlag_Reconstruction
    MP%Nvar=Nvar   ! 变量（方程）数目（最密网格使用湍流模型，数目与Nvar相同）

    end subroutine set_control_para

!---------Continue boundary (inner boundary) -------------------------------
!    根据连接关系，将对应点的坐标写入缓冲区，并计算中心点坐标 (处理1套网格) 
!    两块网格交界区布置 LAP层虚网格
!    交换虚网格的几何信息
!    对于物理边界，采用外插的方法构造虚网格的几何信息 （利用高阶外插）
     
    subroutine Update_coordinate_buffer_onemesh(nMesh)
    use Global_Var
    use Com_ctrl
    implicit none
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B,B1
    Type (BC_MSG_TYPE),pointer:: Bc,Bc1
    real(sp) ,allocatable:: Ux(:,:,:)
    integer:: nMesh,i,j,k,k1,m,mBlock,ksub,n1,m_neighbour,msub,orient,Kflag_initial
    integer:: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2
    integer:: nx,ny,i1,i2,j1,j2
    character(len=30)::filename
 ! ------------------------------------------------------ 
    MP=>Mesh(nMesh)   ! 网格 （nMesh=3,2,1代表 粗、中、密网格） 

    do mBlock=1,MP%Num_Block
        B => MP%Block(mBlock)
        do  ksub=1,B%subface
            Bc=> B%bc_msg(ksub)
            ibegin2=Bc%ist; iend2=Bc%iend; jbegin2=Bc%jst; jend2=Bc%jend     ! write
            if(Bc%neighb .gt. 0 ) then   ! inner boundary
                n1=Bc%iend-Bc%ist+Bc%jend-Bc%jst+1  
                allocate(Ux(2,n1,LAP))
                m_neighbour=Bc%neighb;  msub= Bc%subface   
                B1 => Mesh(nMesh)%Block(m_neighbour)
                Bc1=> B1%bc_msg(msub)
                ibegin1=Bc1%ist; iend1=Bc1%iend; jbegin1=Bc1%jst; jend1=Bc1%jend   ! read
                do i=1,LAP
                    if(Bc1%face .eq. 1 ) then                     !  boundary  i-
                        Ux(1,:,i)=B1%x(ibegin1+i,jbegin1:jend1) ;  Ux(2,:,i)=B1%y(ibegin1+i,jbegin1:jend1); 
                    else if ( Bc1%face .eq. 2) then               !  boundary  j-
                        Ux(1,:,i)=B1%x(ibegin1:iend1,jbegin1+i);   Ux(2,:,i)=B1%y(ibegin1:iend1,jbegin1+i)
                    else if( Bc1%face .eq. 3) then                !  boundary  i+
                        Ux(1,:,i)=B1%x(iend1-i,jbegin1:jend1) ;    Ux(2,:,i)=B1%y(iend1-i,jbegin1:jend1);
                    else                                                   !  boundary  j+
                        Ux(1,:,i)=B1%x(ibegin1:iend1,jend1-i);     Ux(2,:,i)=B1%y(ibegin1:iend1,jend1-i)
                    endif
                enddo
                orient=Bc%orient
                do k=1,n1
                    if(orient .eq. 2) then
                        k1=k
                    else
                        k1=n1+1-k
                    endif
                    do i=1,LAP
                        if(Bc%face .eq. 1 ) then             !  boundary  i-
                            B%x(ibegin2-i,jbegin2+k-1)=Ux(1,k1,i) ; B%y(ibegin2-i,jbegin2+k-1)=Ux(2,k1,i)
                        else if(Bc%face .eq. 2 ) then        !  boundary  j-
                            B%x(ibegin2+k-1,jbegin2-i)=Ux(1,k1,i) ; B%y(ibegin2+k-1,jbegin2-i)=Ux(2,k1,i)
                        else if (Bc%face .eq. 3 ) then       !  boundary  i+
                            B%x(iend2+i,jbegin2+k-1)=Ux(1,k1,i) ; B%y(iend2+i,jbegin2+k-1)=Ux(2,k1,i)
                         else                                          !  boundary  j+
                            B%x(ibegin2+k-1,jend2+i)=Ux(1,k1,i) ; B%y(ibegin2+k-1,jend2+i)=Ux(2,k1,i)
                         endif
                    enddo
                enddo
                deallocate(Ux)
            else  ! not inner boundary
 ! 非内边界点，Ghost Cell的坐标采用外推方法获得 （有待改进）
                do i=1,LAP
                    if(Bc%face .eq. 1 ) then             !  boundary  i-
                        B%x(ibegin2-i,jbegin2:jend2)=2.d0*B%x(ibegin2,jbegin2:jend2)-B%x(ibegin2+i,jbegin2:jend2) 
                        B%y(ibegin2-i,jbegin2:jend2)=2.d0*B%y(ibegin2,jbegin2:jend2)-B%y(ibegin2+i,jbegin2:jend2)   
                    else if(Bc%face .eq. 2 ) then    !  boundary  j-
                        B%x(ibegin2:iend2,jbegin2-i)=2.d0*B%x(ibegin2:iend2,jbegin2)-B%x(ibegin2:iend2,jbegin2+i) 
                        B%y(ibegin2:iend2,jbegin2-i)=2.d0*B%y(ibegin2:iend2,jbegin2)-B%y(ibegin2:iend2,jbegin2+i)  
                    else if (Bc%face .eq. 3 ) then       !  boundary  i+
                        
                        B%x(iend2+i,jbegin2:jend2)=2.d0*B%x(iend2,jbegin2:jend2)-B%x(iend2-i,jbegin2:jend2) 
                        B%y(iend2+i,jbegin2:jend2)=2.d0*B%y(iend2,jbegin2:jend2)-B%y(iend2-i,jbegin2:jend2) 
                    else                                          !  boundary  j+
                        
                        B%x(ibegin2:iend2,jend2+i)=2.d0*B%x(ibegin2:iend2,jend2)-B%x(ibegin2:iend2,jend2-i)
                        B%y(ibegin2:iend2,jend2+i)=2.d0*B%y(ibegin2:iend2,jend2)-B%y(ibegin2:iend2,jend2-i)
                    endif
                    
                enddo  
            endif       
        enddo
!  Conner point     
        nx=B%nx; ny=B%ny
        do j=1,LAP
        do i=1,LAP
            i1=1-i; j1=1-j
            i2=nx+i;j2=ny+j 
            B%x(i1,j1)=B%x(i1,1)+B%x(1,j1)-B%x(1,1); B%y(i1,j1)=B%y(i1,1)+B%y(1,j1)-B%y(1,1) 
            B%x(i2,j2)=B%x(i2,ny)+B%x(nx,j2)-B%x(nx,ny); B%y(i2,j2)=B%y(i2,ny)+B%y(nx,j2)-B%y(nx,ny) 
            B%x(i1,j2)=B%x(i1,ny)+B%x(1,j2)-B%x(1,ny); B%y(i1,j2)=B%y(i1,ny)+B%y(1,j2)-B%y(1,ny)
            B%x(i2,j1)=B%x(i2,1)+B%x(nx,j1)-B%x(nx,1); B%y(i2,j1)=B%y(i2,1)+B%y(nx,j1)-B%y(nx,1)
        enddo
        enddo
	    

  ! Coordinates of the Cell center 
        do j=1-LAP,ny+LAP-1
        do i=1-LAP,nx+LAP-1
            B%x1(i,j)=(B%x(i,j)+B%x(i+1,j)+B%x(i,j+1)+B%x(i+1,j+1))*0.25
            B%y1(i,j)=(B%y(i,j)+B%y(i+1,j)+B%y(i,j+1)+B%y(i+1,j+1))*0.25
        enddo
        enddo

    enddo
 
    end  subroutine Update_coordinate_buffer_onemesh
  
    
    
