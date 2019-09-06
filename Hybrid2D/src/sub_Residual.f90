 !-----The core subroutines: Comput Residual: includes the inviscous and viscous flux ---------
 !     Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS. lixl@imech.ac.cn
 !     Ver 1.1
 !--------------------------------------------------------------------------------------------
 ! Since V0.96, Original value is used in Reconstruction
 ! Since V0.971, User can choose using original variables or conservative variables
 ! Since V0.99, User can choose using characteristic variables
 ! Ver 1.1 , 2011-3-28: multi-grid is supported, Code by Li Xinliang
 ! Ver 1.13, 2011-5-22: OMUSCL2 is adding
 !--------------------------------------------------------------------------------------------
   subroutine Residual(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   integer:: nMesh,mBlock
 
     MP=> Mesh(nMesh)
     B => MP%Block(mBlock)
  
   if(If_viscous .eq. 1) then   
    call get_viscous(nMesh,mBlock)         ! 计算层流粘性系数
    B%Amu_t(:,:)=0.d0                       ! 湍流粘性系数清零
!   湍流模型，计算湍流粘性系数Amu_t

   
   endif
    call Residual_FVM(nMesh,mBlock)
    end subroutine Residual
    

!------------------------------------------------------------------------------
! 计算残差 (残差=右端项=净流量); 是求解N-S方程的核心模块; 
! 湍流模型也包含在这个模块里
! Flux (= inviscous flux + viscous flux )

   subroutine Residual_FVM(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   
   implicit none
   real(sp) ,dimension(4):: UL,UR,UL1,UR1,UL2,UR2,UL3,UR3,QL,QR,Flux0
   real(sp) :: U0(4,4)
   integer:: nMesh,mBlock,Scheme,IFlux,Reconstruction
   integer:: nx,ny,i,j,m
   real(sp) :: minmod,dx,dy,si,ni1,ni2,sj,nj1,nj2
   real(sp) :: Diu,Div,Dju,Djv,DiT,DjT,ux,vx,Tx,uy,vy,Ty,t11,t12,t22,t33,E1,E2,qx,qy
   real(sp) :: Dix,Diy,Djx,Djy,Ds,Amu1,u1,v1,Amk1,tmp
   real(sp) :: dl,uul,vl,ppl,dr,uur,vr,ppr,pr1
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx ; ny=B%ny

   Scheme=MP%Iflag_Scheme         ! 数值格式 （不同网格上采用不同格式）
   IFlux=MP%Iflag_Flux            ! 通量技术
   Reconstruction=MP%IFlag_Reconstruction  ! 重构方式
  
   Fluxi=0.d0  ! 清零
   Fluxj=0.d0



!----- i- direction ----------------------------------------------------------------------------------

   do j=1,ny-1
   do i=1,nx
	 si=B%si(i,j)   ! 控制体边界长度
	 ni1=B%ni1(i,j); ni2=B%ni2(i,j)   ! 法方向
!-----无粘项 -----------Inviscous term-------------------------------------------------------------

!-----选择重构（插值）界面值的 方式  （使用原始变量、守恒变量或特征变量重构）

     if(Reconstruction .eq. Reconst_Original) then
       U0(:,1)=d(i-2:i+1,j) ; U0(:,2)=uu(i-2:i+1,j) ; U0(:,3)=v(i-2:i+1,j); U0(:,4)=p(i-2:i+1,j)
       call Reconstuction_original(U0,UL,UR,gamma,Scheme)
     else if (Reconstruction .eq. Reconst_Conservative) then
      do m=1,4
       U0(:,m)=B%U(m,i-2:i+1,j)
      enddo
      call Reconstuction_conservative(U0,UL,UR,gamma,Scheme)
     else
      do m=1,4
       U0(:,m)=B%U(m,i-2:i+1,j)
      enddo
      call Reconstuction_Characteristic(U0,UL,UR,gamma,Scheme)
     endif
   

!-------坐标旋转到垂直于界面 （法向-切向 坐标系）
      QL(1)=UL(1); QL(2)=UL(2)*ni1+UL(3)*ni2 ; QL(3)= -UL(2)*ni2+UL(3)*ni1 ; QL(4)=UL(4) ! 密度、压力、法向速度、切向速度 （左值）
      QR(1)=UR(1); QR(2)=UR(2)*ni1+UR(3)*ni2 ; QR(3)= -UR(2)*ni2+UR(3)*ni1 ; QR(4)=UR(4) ! 密度、压力、法向速度、切向速度 （右值）

!-------选择通量分裂方法----(Steger_Warming,Van Leer, Roe, AUSM, HLL, HLLC----
      if(IFlux .eq. Flux_Steger_Warming ) then
         call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)   
      else if (IFlux .eq. Flux_HLL ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,0)
      else if (IFlux .eq. Flux_HLLC ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,1)
      else if  (IFlux .eq. Flux_VanLeer ) then
         call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_AUSM ) then
         call Flux_Ausm_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_Roe ) then
         call Flux_Roe_1D(QL,QR,Flux0,gamma)   
      endif


!------------------------------------------------	  
!   算出通量 （变换回x-y坐标系）  
        Fluxi(1,i,j)=-Flux0(1)*si                            ! 质量通量
        Fluxi(2,i,j)=-(Flux0(2)*ni1-Flux0(3)*ni2)*si         ! x-方向动量通量
        Fluxi(3,i,j)=-(Flux0(2)*ni2+Flux0(3)*ni1)*si         ! y-方向动量通量
        Fluxi(4,i,j)=-Flux0(4)*si                            ! 能量通量
      
	  if (If_viscous .eq. 0) cycle   ! 如求解Euler方程（无粘），则无需粘性项计算

!--i-方向无粘通量计算结束，计算i-方向粘性通量
!--------------------------------------------------------------------------------------------------------- 
!----------- Viscous term---------------------------------------------------------------------------------  
! 扩散系数（粘性系数、热传导系数）= 界面两侧值的平均 （边界处采用单侧值)    
    if( i .eq. 1) then 
     Amu1=B%Amu(i,j) + B%Amu_t(i,j)                  ! 粘性系数 (层流+湍流), 界面上的值=两侧值的平均
     Amk1=Cp*(B%Amu(i,j)/Pr + B%Amu_t(i,j)/PrT )     ! 热传导系数
    else if (i .eq. nx ) then
      Amu1=B%Amu(i-1,j) + B%Amu_t(i-1,j)                  ! 粘性系数 (层流+湍流), 界面上的值=两侧值的平均
      Amk1=Cp*(B%Amu(i-1,j)/Pr + B%Amu_t(i-1,j)/PrT )     ! 热传导系数
    else
      Amu1=(B%Amu(i-1,j)+B%Amu(i,j) + B%Amu_t(i-1,j)+B%Amu_t(i,j) )*0.5d0                 ! 粘性系数 (层流+湍流), 界面上的值=两侧值的平均
      Amk1=Cp*((B%Amu(i-1,j)+B%Amu(i,j))/Pr + (B%Amu_t(i-1,j)+B%Amu_t(i,j))/PrT )*0.5d0   ! 热传导系数
    endif

!----Jocabian系数 （物理坐标对计算坐标的导数, 用于计算物理量的导数）
    Dix=B%x1(i,j)-B%x1(i-1,j)
    Diy=B%y1(i,j)-B%y1(i-1,j)
    Djx=(B%x1(i-1,j+1)+B%x1(i,j+1)-B%x1(i-1,j-1)-B%x1(i,j-1))*0.25d0
    Djy=(B%y1(i-1,j+1)+B%y1(i,j+1)-B%y1(i-1,j-1)-B%y1(i,j-1))*0.25d0
    Ds=1.d0/(Dix*Djy-Djx*Diy)
! 物理量对计算坐标的导数    
    Diu=uu(i,j)-uu(i-1,j)
    Div=v(i,j)-v(i-1,j)
    DiT=T(i,j)-T(i-1,j)
    Dju=(uu(i-1,j+1)+uu(i,j+1)-uu(i-1,j-1)-uu(i,j-1))*0.25d0
    Djv=(v(i-1,j+1)+v(i,j+1)-v(i-1,j-1)-v(i,j-1))*0.25d0
    DjT=(T(i-1,j+1)+T(i,j+1)-T(i-1,j-1)-T(i,j-1))*0.25d0
! 物理量对x,y坐标的导数
    ux=(Diu*Djy-Dju*Diy)*Ds
    vx=(Div*Djy-Djv*Diy)*Ds
    Tx=(DiT*Djy-DjT*Diy)*Ds
    uy=(-Diu*Djx+Dju*Dix)*Ds
    vy=(-Div*Djx+Djv*Dix)*Ds
    Ty=(-DiT*Djx+DjT*Dix)*Ds
!  粘性应力及能量通量
    u1=(uu(i,j)+uu(i-1,j))*0.5d0
    v1=(v(i,j)+v(i-1,j))*0.5d0
    
    t11=((4.d0/3.d0)*ux-(2.d0/3.d0)*vy)*Amu1
    t22=((4.d0/3.d0)*vy-(2.d0/3.d0)*ux)*Amu1
    t12=(uy+vx)*Amu1
    t33=Amu1*(-(2.d0/3.d0)*ux-(2.d0/3.d0)*vy)
    qx=Amk1*Tx
    qy=Amk1*Ty
    
    
    B%TaoNSi_out(1,i,j)=t11
    B%TaoNSi_out(2,i,j)=t22
    B%TaoNSi_out(3,i,j)=t12
    B%TaoNSi_out(4,i,j)=t33
    B%Qi_out(1,i,j)=qx
    B%Qi_out(2,i,j)=qy

     
    if((B%solver==GKUA))then
        t11=t11+B%TaoNSi_G(1,i,j)
        t22=t22+B%TaoNSi_G(2,i,j)
        t12=t12+B%TaoNSi_G(3,i,j)
        t33=t33+B%TaoNSi_G(4,i,j)
         
        qx=qx+B%Qi_G(1,i,j)
        qy=qy+B%Qi_G(2,i,j)
        !t11=B%TaoNSi_G(1,i,j)
        !t22=B%TaoNSi_G(2,i,j)
        !t12=B%TaoNSi_G(3,i,j)
        !t33=B%TaoNSi_G(4,i,j)
        ! 
        !qx=B%Qi_G(1,i,j)
        !qy=B%Qi_G(2,i,j)
        !

    end if 

      
    
    E1=u1*t11+v1*t12+qx
    E2=u1*t12+v1*t22+qy
! 添加粘性通量
    Fluxi(2,i,j)=Fluxi(2,i,j)   +(t11*ni1+t12*ni2)*si   
    Fluxi(3,i,j)=Fluxi(3,i,j)   +(t12*ni1+t22*ni2)*si
    Fluxi(4,i,j)=Fluxi(4,i,j)   +(E1*ni1+ E2*ni2)*si
  
   enddo
   enddo
   
!==================================================================================================================
!                                         j-方向的无粘及粘性通量    
!-----------------------------------------j- direction -------------------------------------------------------------	 
   do j=1,ny
   do i=1,nx-1
  
   ! 边长，法方向   
	sj=B%sj(i,j)
	nj1=B%nj1(i,j); nj2=B%nj2(i,j)

!-----Inviscous term --------------------------------------------
 
    if(Reconstruction .eq. Reconst_Original) then
       U0(:,1)=d(i,j-2:j+1) ; U0(:,2)=uu(i,j-2:j+1) ; U0(:,3)=v(i,j-2:j+1); U0(:,4)=p(i,j-2:j+1)
       call Reconstuction_original(U0,UL,UR,gamma,Iflag_Scheme)
    else if (Reconstruction .eq. Reconst_Conservative) then
      do m=1,4
       U0(:,m)=B%U(m,i,j-2:j+1)
      enddo
       call Reconstuction_conservative(U0,UL,UR,gamma,Iflag_Scheme)
    else
      do m=1,4
       U0(:,m)=B%U(m,i,j-2:j+1)
      enddo
      call Reconstuction_characteristic(U0,UL,UR,gamma,Iflag_Scheme)
    endif	   
    

       QL(1)=UL(1); QL(2)=UL(2)*nj1+UL(3)*nj2 ; QL(3)= -UL(2)*nj2+UL(3)*nj1 ; QL(4)=UL(4) !density, normal/tangitial velocity and pressure
       QR(1)=UR(1); QR(2)=UR(2)*nj1+UR(3)*nj2 ; QR(3)= -UR(2)*nj2+UR(3)*nj1 ; QR(4)=UR(4)
 
!-------选择通量分裂方法----(Steger_Warming,Van Leer, Roe, AUSM, HLL, HLLC)----
      if(IFlux .eq. Flux_Steger_Warming ) then
         call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)   
      else if (IFlux .eq. Flux_HLL ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,0)
      else if (IFlux .eq. Flux_HLLC ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,1)
      else if  (IFlux .eq. Flux_VanLeer ) then
         call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_AUSM ) then
         call Flux_Ausm_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_Roe ) then
         call Flux_Roe_1D(QL,QR,Flux0,gamma)   
      endif
     
!   无粘通量   	   
    Fluxj(1,i,j)=-Flux0(1)*sj 
    Fluxj(2,i,j)=-(Flux0(2)*nj1-Flux0(3)*nj2)*sj
    Fluxj(3,i,j)=-(Flux0(2)*nj2+Flux0(3)*nj1)*sj
    Fluxj(4,i,j)=-Flux0(4)*sj

    if (If_viscous .eq. 0) cycle  ! 如求解Euler方程（无粘），则无需粘性项计算
 !---------Viscous term -----------------------------------------------------------------------------
 ! 扩散系数（粘性系数、热传导系数）= 界面两侧值的平均 （边界处采用单侧值)    
    if( j .eq. 1) then
     Amu1=B%Amu(i,j)+B%Amu_t(i,j)
     Amk1=Cp*(B%Amu(i,j)/Pr +B%Amu_t(i,j)/PrT )   ! 热传导系数
    else if (j .eq. ny) then
      Amu1=B%Amu(i,j-1)+B%Amu_t(i,j-1)
      Amk1=Cp*(B%Amu(i,j-1)/Pr +B%Amu_t(i,j-1)/PrT )   ! 热传导系数
    else
     Amu1=(B%Amu(i,j)+B%Amu(i,j-1) +B%Amu_t(i,j)+B%Amu_t(i,j-1))*0.5d0
     Amk1=Cp*((B%Amu(i,j-1)+B%Amu(i,j))/Pr + (B%Amu_t(i,j-1)+B%Amu_t(i,j))/PrT )*0.5d0   ! 热传导系数
    endif
    
    Dix=(B%x1(i+1,j-1)+B%x1(i+1,j)-B%x1(i-1,j-1)-B%x1(i-1,j))*0.25d0
    Diy=(B%y1(i+1,j-1)+B%y1(i+1,j)-B%y1(i-1,j-1)-B%y1(i-1,j))*0.25d0
    Djx=B%x1(i,j)-B%x1(i,j-1)
    Djy=B%y1(i,j)-B%y1(i,j-1)
   
    Ds=1.d0/(Dix*Djy-Djx*Diy)

    Diu=(uu(i+1,j-1)+uu(i+1,j)-uu(i-1,j-1)-uu(i-1,j))*0.25d0
    Div=(v(i+1,j-1)+v(i+1,j)-v(i-1,j-1)-v(i-1,j))*0.25d0
    DiT=(T(i+1,j-1)+T(i+1,j)-T(i-1,j-1)-T(i-1,j))*0.25d0
    Dju=uu(i,j)-uu(i,j-1)
    Djv=v(i,j)-v(i,j-1)
    DjT=T(i,j)-T(i,j-1)
!     
    ux=(Diu*Djy-Dju*Diy)*Ds
    vx=(Div*Djy-Djv*Diy)*Ds
    Tx=(DiT*Djy-DjT*Diy)*Ds
    uy=(-Diu*Djx+Dju*Dix)*Ds
    vy=(-Div*Djx+Djv*Dix)*Ds
    Ty=(-DiT*Djx+DjT*Dix)*Ds
    t11=((4.d0/3.d0)*ux-(2.d0/3.d0)*vy)*Amu1
    t22=((4.d0/3.d0)*vy-(2.d0/3.d0)*ux)*Amu1
    t12=(uy+vx)*Amu1
    t33=Amu1*(-(2.d0/3.d0)*ux-(2.d0/3.d0)*vy)
    qx=Amk1*Tx
    qy=Amk1*Ty
    
    B%TaoNSj_out(1,i,j)=t11
    B%TaoNSj_out(2,i,j)=t22
    B%TaoNSj_out(3,i,j)=t12
    B%TaoNSj_out(4,i,j)=t33
    
    B%Qj_out(1,i,j)=qx  
    B%Qj_out(2,i,j)=qy
    
  
    if((B%solver==GKUA))then
         t11=t11+B%TaoNSj_G(1,i,j)
         t22=t22+B%TaoNSj_G(2,i,j)
         t12=t12+B%TaoNSj_G(3,i,j)
         t33=t33+B%TaoNSj_G(4,i,j)
         
         qx=qx+B%Qj_G(1,i,j)
         qy=qy+B%Qj_G(2,i,j)
   
         ! t11=B%TaoNSj_G(1,i,j)
         !t22=B%TaoNSj_G(2,i,j)
         !t12=B%TaoNSj_G(3,i,j)
         !t33=B%TaoNSj_G(4,i,j)
         !
         !qx=B%Qj_G(1,i,j)
         !qy=B%Qj_G(2,i,j)
         
    end if

 
    
    u1=(uu(i,j)+uu(i,j-1))*0.5d0
    v1=(v(i,j)+v(i,j-1))*0.5d0

    E1=u1*t11+v1*t12+qx
    E2=u1*t12+v1*t22+qy

    Fluxj(2,i,j)=Fluxj(2,i,j) +(t11*nj1+t12*nj2)*sj
    Fluxj(3,i,j)=Fluxj(3,i,j) +(t12*nj1+t22*nj2)*sj
    Fluxj(4,i,j)=Fluxj(4,i,j) +(E1*nj1+ E2*nj2)*sj

 
   enddo
   enddo
  
 !-------计算残差 （净流量）----------------------------------------------------  
 !                                    不含强迫函数 （多重网格时，需添加强迫函数）
    do j=1,ny-1
    do i=1,nx-1
    do m=1,4
      B%Res(m,i,j)= Fluxi(m,i+1,j)-Fluxi(m,i,j)+Fluxj(m,i,j+1)-Fluxj(m,i,j)
	enddo
	enddo
	enddo

  end  


      function minmod(a,b)
       use const_var ,only :sp,dop
      implicit none
      real(sp) :: a,b,minmod
      if(a*b .le. 0.d0) then
       minmod=0.d0
      else if (abs(a) .le. abs(b)) then
       minmod=a
      else
       minmod=b
      endif
      end
!--------------------------------------------------------------------------------


!-------------------------------------------------------------------
! 利用原始变量进行平均 （用于计算角点的值）
  subroutine U_average(U1,U2,U0,gamma)
  use const_var ,only :sp,dop
  implicit none
  real(sp) ,dimension(4):: U1,U2,U0
  real(sp) :: gamma,d1,uu1,v1,p1,d2,uu2,v2,p2,d0,uu0,v0,p0
    d1=U1(1); uu1=U1(2)/d1; v1=U1(3)/d1; p1=(U1(4)-(uu1*U1(2)+v1*U1(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure 
    d2=U2(1); uu2=U2(2)/d2; v2=U2(3)/d2; p2=(U2(4)-(uu2*U2(2)+v2*U2(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure 
    d0=(d1+d2)*0.5d0; uu0=(uu1+uu2)*0.5d0; v0=(v1+v2)*0.5d0; p0=(p1+p2)*0.5d0
    U0(1)=d0; U0(2)=d0*uu0; U0(3)=d0*v0; U0(4)=p0/(gamma-1.d0)+(d0*uu0*uu0+d0*v0*v0)*0.5d0
  end

!--------------------------------------------------------------------------
! 利用原始变量重构 
   subroutine Reconstuction_original(U0,UL,UR,gamma,Iflag_Scheme)
    use const_var ,only :sp,dop
     implicit none
     real(sp) :: U0(4,4),UL(4),UR(4),gamma
     integer:: Iflag_Scheme,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1    ;   m=1,4 for d,u,v,p
     do m=1,4
       call scheme_fP(UL(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
       call scheme_fm(UR(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
     enddo
   end

! 利用守恒变量重构
   subroutine Reconstuction_conservative(U0,UL,UR,gamma,Iflag_Scheme)
     use const_var ,only :sp,dop 
     implicit none
     real(sp) :: U0(4,4),UL(4),UR(4),QL(4),QR(4),gamma
     integer:: Iflag_Scheme,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1   ; m for the conservative variables U0(1,m)=d, U0(2,m)=d*u, ....
     do m=1,4
        call scheme_fP(QL(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
        call scheme_fm(QR(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
      enddo
       UL(1)=QL(1); UL(2)=QL(2)/UL(1); UL(3)=QL(3)/UL(1); UL(4)=(QL(4)-(UL(2)*QL(2)+UL(3)*QL(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
       UR(1)=QR(1); UR(2)=QR(2)/UR(1); UR(3)=QR(3)/UR(1); UR(4)=(QR(4)-(UR(2)*QR(2)+UR(3)*QR(3))*0.5d0)*(gamma-1.d0)  ! find a bug, removed
   end

!  利用特征变量重构
   subroutine Reconstuction_Characteristic(U0,UL,UR,gamma,Iflag_Scheme)
     use const_var ,only :sp,dop
     implicit none
     real(sp) :: U0(4,4),UL(4),UR(4),gamma
     real(sp) :: Uh(4),S(4,4),S1(4,4),V0(4,4),VL(4),VR(4),QL(4),QR(4)
     real(sp) :: V2,d1,u1,v1,p1,c1,tmp0,tmp1,tmp3,tmp5
     integer:: Iflag_Scheme,i,j,k,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1   ; m for the conservative variables U0(1,m)=d, U0(2,m)=d*u, ....
     Uh(:)=0.5d0*(U0(2,:)+U0(3,:))           ! conservative variables in the point I-1/2  (or i)
     d1=Uh(1); u1=Uh(2)/d1; v1=Uh(3)/d1; p1=(Uh(4)-(Uh(2)*u1+Uh(3)*v1)*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
     c1=sqrt(gamma*p1/d1)

	  V2=(u1*u1+v1*v1)*0.5d0
	  tmp1=(gamma-1.d0)/c1
      tmp3=(gamma-1.d0)/(c1*c1)
	  tmp5=1.d0/(2.d0*c1)
	  tmp0=1.d0/tmp3

!  A=S(-1)*LAMDA*S    see 《计算空气动力学》 158-159页   (with alfa=1, beta=0)

          S(1,1)=V2-tmp0;       S(1,2)=-u1 ;         S(1,3)=-v1 ;      S(1,4)=1.d0
          S(2,1)=-v1 ;          S(2,2)=0.d0 ;        S(2,3)=1.d0 ;     S(2,4)=0.d0 
          S(3,1)=-u1-V2*tmp1;   S(3,2)=1.d0+tmp1*u1; S(3,3)=tmp1*v1;   S(3,4)=-tmp1
          S(4,1)=-u1+V2*tmp1;   S(4,2)=1.d0-tmp1*u1; S(4,3)=-tmp1*v1;  S(4,4)=tmp1 
 
          S1(1,1)=-tmp3;    S1(1,2)=0.d0;   S1(1,3)=-tmp5 ;         S1(1,4)=tmp5
          S1(2,1)=-tmp3*u1; S1(2,2)=0.d0;   S1(2,3)=0.5d0-u1*tmp5 ; S1(2,4)=0.5d0+u1*tmp5
          S1(3,1)=-tmp3*v1; S1(3,2)=1.d0;   S1(3,3)=-v1*tmp5;       S1(3,4)=v1*tmp5
          S1(4,1)=-tmp3*V2; S1(4,2)=v1;     S1(4,3)=(c1*u1-V2-tmp0)*tmp5; S1(4,4)=(c1*u1+V2+tmp0)*tmp5
     
! V=SU      V(k)=S*U(k)
    do k=1,4
     
     do m=1,4
     V0(k,m)=0.d0
      do j=1,4
      V0(k,m)=V0(k,m)+S(m,j)*U0(k,j)
      enddo
     enddo
    enddo
     
     do m=1,4
       call scheme_fP(VL(m),V0(1,m),V0(2,m),V0(3,m),V0(4,m),Iflag_Scheme) 
       call scheme_fm(VR(m),V0(1,m),V0(2,m),V0(3,m),V0(4,m),Iflag_Scheme) 
     enddo

    do m=1,4
    QL(m)=0.d0; QR(m)=0.d0
    do j=1,4
    QL(m)=QL(m)+S1(m,j)*VL(j)
    QR(m)=QR(m)+S1(m,j)*VR(j)
    enddo
    enddo
     UL(1)=QL(1); UL(2)=QL(2)/UL(1); UL(3)=QL(3)/UL(1); UL(4)=(QL(4)-(UL(2)*QL(2)+UL(3)*QL(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
     UR(1)=QR(1); UR(2)=QR(2)/UR(1); UR(3)=QR(3)/UR(1); UR(4)=(QR(4)-(UR(2)*QR(2)+UR(3)*QR(3))*0.5d0)*(gamma-1.d0)  ! find a bug, removed

   end


!-----------------------------------------------------------
! 数值格式，构造UL=U(j+1/2,L); u1=u(j-1),u2=u(j),u3=u(j+1),u4=u(j+2)
    subroutine scheme_fP(uL,u1,u2,u3,u4,Iflag_Scheme)
     use   const_var
     
     implicit none
     real(sp) :: uL,u1,u2,u3,u4,minmod
     real(sp) :: IS1,IS2,a1,a2,w1,w2,up,um,s
     real(sp) :: k=1.d0,b=2.d0
	 real(sp) :: r1,r2,f1,f
     real(sp) ,parameter::k3=1.d0/3.d0,ep=1.d-6
     integer:: Iflag_Scheme
     if(Iflag_Scheme .eq. Scheme_UD1 ) then
       UL=u2                            ! 1阶迎风格式
	 else if(Iflag_Scheme .eq. Scheme_UD3 ) then
       UL=(-u1+5.d0*u2+2.d0*u3  )/6.d0   !UD 3nd order  3阶迎风格式
     else if(Iflag_Scheme .eq. Scheme_NND2 ) then
       UL=u2+0.5d0*minmod(u2-u1,u3-u2)   ! NND 2nd order  2阶NND
     else if(Iflag_Scheme .eq. Scheme_WENO3) then  ! WENO3 scheme 
       IS1=(u2-u1)**2 ;    IS2=(u3-u2)**2
       a1=1.d0/(3.d0*(1.d-6+IS1)**2) ;   a2=2.d0/(3.d0*(1.d-6+IS2)**2)
       w1=a1/(a1+a2) ;     w2=a2/(a1+a2) 
       UL=w1*(-u1/2.d0+3.d0*u2/2.d0)+w2*(u2/2.d0+u3/2.d0)   !3阶WENO
     else if(Iflag_Scheme .eq. Scheme_MUSCL2 ) then         !2阶MUSCL(Minmod限制器）
       UL=u2+0.25d0*((1.d0-k)*minmod(u2-u1,b*(u3-u2))+(1.d0+k)*minmod(u3-u2,b*(u2-u1))) ! MUSCL
     else if(Iflag_Scheme .eq. Scheme_MUSCL3 ) then          ! 3阶MUSCL (Van Albada限制器)
        up=u3-u2; um=u2-u1                                   ! 1阶 前差、后差
        s=(2.d0*up*um+ep)/(up*up+um*um+ep)                   ! Van Albada限制器 （光滑区，前差与后差接近，该值接近1）
        UL=u2+0.25d0*s*((1.d0-k3*s)*um+(1.d0+k3*s)*up)       ! 3阶MUSCL (光滑区逼近3阶迎风)
     else if(Iflag_Scheme .eq. Scheme_OMUSCL2 ) then         ! 2阶优化的MUSCL 格式 (Developed by Leng Yan)
        r1=(u2-u1+ep)/(u3-u2+ep); r2=(u3-u2+ep)/(u4-u3+ep)
		f1=0.8d0-0.175d0/r2+0.375d0*r1
		f=max(0.d0,min(2.d0,f1,2.d0*r1))
        UL=u2+0.5d0*f*(u3-u2)                               ! 2阶优化的MUSCL: OMUSCL2
	 endif
     end

! 数值格式，构造UR=U(j+1/2,R) ; u1=u(j-1),u2=u(j),u3=u(j+1),u4=u(j+2)
     subroutine scheme_fm(uR,u1,u2,u3,u4,Iflag_Scheme)
     use   const_var
     
     implicit none
     real(sp) :: uR,u1,u2,u3,u4,minmod,up,um,s
     real(sp) :: IS1,IS2,a1,a2,w1,w2
  	 real(sp) :: r1,r2,f1,f
     real(sp) :: k=1.d0/3.0,b=1.d0
     real(sp) ,parameter::k3=1.d0/3.d0,ep=1.d-6
     integer:: Iflag_Scheme
     if(Iflag_Scheme .eq. Scheme_UD1) then
	   UR=u3                              ! 1阶迎风
	 else if(Iflag_Scheme .eq. Scheme_UD3 ) then
       UR=(2.d0*u2+5.d0*u3-u4  )/6.d0     !UD 3nd order
     else if(Iflag_Scheme .eq.Scheme_NND2 ) then
        UR=u3-0.5d0*minmod(u3-u2,u4-u3)    ! NND 2nd order
     else if(Iflag_Scheme .eq. Scheme_WENO3 ) then   ! WENO3   
        IS1=(u3-u4)**2;   IS2=(u2-u3)**2
        a1=1.d0/(3.d0*(1.d-6+IS1)**2) ;   a2=2.d0/(3.d0*(1.d-6+IS2)**2)
        w1=a1/(a1+a2) ;     w2=a2/(a1+a2) 
        UR=w1*(-u4/2.d0+3.d0*u3/2.d0)+w2*(u3/2.d0+u2/2.d0)       !WENO 3nd order
     else if(Iflag_Scheme .eq. Scheme_MUSCL2 ) then
        uR=u3-0.25d0*((1.d0-k)*minmod(u4-u3,b*(u3-u2))+(1.d0+k)*minmod(u3-u2,b*(u4-u3)))   ! MUSCL
     else if(Iflag_Scheme .eq. Scheme_MUSCL3 ) then      ! 3阶MUSCL (Van Albada限制器)
         up=u4-u3 ; um=u3-u2                             !前差、后差
         s=(2.d0*up*um+ep)/(up*up+um*um+ep)
         UR=u3-0.25d0*s*((1.d0-k3*s)*up+(1.d0+k3*s)*um)
     else if(Iflag_Scheme .eq. Scheme_OMUSCL2 ) then         ! 2阶优化的MUSCL 格式 (Developed by Leng Yan)
        r1=(u4-u3+ep)/(u3-u2+ep); r2=(u3-u2+ep)/(u2-u1+ep)
		f1=0.8d0-0.175d0/r2+0.375d0*r1
		f=max(0.d0,min(2.d0,f1,2.d0*r1))
        UR=u3-0.5d0*f*(u3-u2)                               ! 2阶优化的MUSCL: OMUSCL2
    
     endif
     end
!------------------------------------------------------------------------------
 


!----------------------------------------------------------------
! 分子粘性系数的计算
   subroutine get_viscous(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   
   implicit none
   real(sp) :: Tsb,OMGgas
   integer:: nMesh,mBlock,i,j,nx,ny
   Type (Block_TYPE),pointer:: B

   B => Mesh(nMesh)%Block(mBlock)
   OMGgas=0.81
!  Ref_Amu_T0=288.15d0
!  T_inf （来流）参考温度
 ! Tsb=110.4d0/T_inf  !空气
  Tsb=148.0/T_inf  !Ar
     nx=B%nx; ny=B%ny
     do j=0,ny
     do i=0,nx
     !  B%Amu(i,j)=1.d0/Re*(1.d0+Tsb)*sqrt(T(i,j)**3)/(Tsb+T(i,j))   !sutherland equation!7-9
       B%Amu(i,j)=1.d0/Re*T(i,j)**OMGgas
 !     Amu(i,j)=1.d0/Re
     enddo
     enddo
  end

