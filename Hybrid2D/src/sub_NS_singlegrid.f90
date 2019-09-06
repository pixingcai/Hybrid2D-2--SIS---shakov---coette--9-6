
!---------------------------------------------------------------------------------------------
! 采用LU_SGS法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
 subroutine NS_Time_advance_LU_SGS(nMesh)
   use Global_var
   use const_var ,only :sp,dop
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(sp) :: du

    call Comput_Residual_one_mesh(nMesh)     ! 单重网格上计算残差
!    if(nMesh .ne. 1) call Add_force_function(nMesh)                !  添加强迫函数（多重网格的粗网格使用,目前LU-SGS不支持）

    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock)                          ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   时间推进 
      do j=3,ny-3
      do i=3,nx-3
      do m=1,MP%NVAR
        du=B%du(m,i,j)  
        B%U(m,i,j)=B%U(m,i,j)+du                      ! U(n+1)=U(n)+dU
      enddo
      enddo
      enddo
    enddo    

!!---------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------  
        !call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
!!      call update_buffer_onemesh(nMesh)              ! 同步各块的交界区



 end subroutine NS_Time_advance_LU_SGS


!-------------------------------------------------------------------------------------------------------   
! 计算残差（网格的全部块）
   Subroutine Comput_Residual_one_mesh(nMesh)
   use Global_Var
   use Flow_Var 
   use const_var ,only :sp,dop
   implicit none
   real(sp) :: Res
   integer:: nMesh,mBlock,nx,ny,i,j,m,Nvar1
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
!---------------------------------------------  
       MP=>Mesh(1)               ! 本层网格
        

!----------------------------------------------
   do mBlock=1,MP%Num_Block
      B => MP%Block(mBlock)                 !第nMesh 重网格的第mBlock块
      B%Res_max(:) =0.d0           ! 最大残差
      B%Res_rms(:) =0.d0           ! 均方根残差    
     nx=B%nx; ny=B%ny
!    分配临时变量
!     allocate(d(0:nx,0:ny),uu(0:nx,0:ny),v(0:nx,0:ny),T(0:nx,0:ny),cc(0:nx,0:ny),p(0:nx,0:ny))   ! Bug found, 2012-5-1
     allocate(d(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),uu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), &
              v(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),T(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),  &
              p(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), cc(nx-1,ny-1))
     allocate(Fluxi(4,nx,ny),Fluxj(4,nx,ny))


!------------------------------------------------------------------------------------     
	 call comput_duvtpckw(nMesh,mBlock)    ! 计算基本量 d,u,v,T,p,cc,Kt,Wt
!---------------------------------------------------------------------------------------
!  求解N-S方程的最核心模块  
!  计算一个网格块的残差（右端项） ; 第nMesh 重网格的第mBlock块  （湍流模型也在该块中计算）  

   call Residual (nMesh,mBlock)       
  
     
!--------------------------------------------------------------------------------------
!   统计最大残差和均方根残差 
    do j=3,ny-3
    do i=3,nx-3
! -------------------------------------------------------------------------------------------
!    时间推进
       do m=1,MP%Nvar
        Res=B%Res(m,i,j)

!--------------------------------------------------------------------------------------------------
        if(abs(Res) .gt. B%Res_max(m))then
            B%Res_max(m)=abs(Res)          ! 最大残差
            B%Resmax_location_i=i
            B%Resmax_location_j=j
        end if
        
        B%Res_rms(m)=B%Res_rms(m)+Res*Res                           ! 均方根残差      
!--------------------------------------------------------------------------------------------------       
	   enddo
    enddo
    enddo
            

!---------------------------------------------------------------------------------------- 
     call comput_Lvc(nMesh,mBlock)
	 call comput_dt(nMesh,mBlock)        ! 计算(当地) 时间步长
 
 !   删除掉临时变量    
	 deallocate(d,uu,v,T,cc,p,Fluxi,Fluxj)
 
   enddo    
   
   B%Res_rms(:)=sqrt(B%Res_rms(:)/(1.d0*Mesh(nMesh)%Num_Cell))   !    全部网格点的总均方根残差
  end Subroutine Comput_Residual_one_mesh



!---------------------------------------------------------
!  利用守恒变量，计算基本量 (d,u,v,T,p,c) 
!----------------------------------------------------------
    subroutine comput_duvtpckw(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
 
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B
	 integer nMesh,mBlock,nx,ny,i,j
     real(sp) :: ptmp
	 ptmp=1.d0/(gamma*Ma*Ma)
     MP=> mesh(nMesh)
     B => MP%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny
 
	  do j=1-LAP,ny+LAP-1
      do i=1-LAP,nx+LAP-1
        d(i,j)=B%U(1,i,j)
        uu(i,j)=B%U(2,i,j)/d(i,j)
        v(i,j)=B%U(3,i,j)/d(i,j)
        T(i,j)=(B%U(4,i,j)-0.5d0*d(i,j)*(uu(i,j)*uu(i,j)+v(i,j)*v(i,j)))/(Cv*d(i,j))
        p(i,j)=ptmp*d(i,j)*T(i,j)
      enddo
      enddo

      do j=1,ny-1
      do i=1,nx-1
        if(T(i,j) .lt. 0.1) then
          print*, "T <0.1 ! ", i,j,mBlock, T(i,j)
          stop
        endif
        cc(i,j)=sqrt(T(i,j))/Ma
      enddo
      enddo
      
    
    end subroutine comput_duvtpckw
!---------------------------------------------------------------------------------
! 计算（当地）时间步长
   subroutine comput_dt(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real(sp) ,parameter:: C0=1
     Type (Block_TYPE),pointer:: B
    
     B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny

     do j=1,ny-1
     do i=1,nx-1
      if(Iflag_local_dt .eq. 0) then   ! 全局时间步长
        B%dt(i,j)=dt_global   
       else                             ! 当地时间步长
         B%dt(i,j)=CFL*B%vol(i,j)/(B%Lci(i,j)+B%Lcj(i,j) +C0*(B%Lvi(i,j)+B%Lvj(i,j)) )     ! 局部时间步长
         if(B%dt(i,j) .gt. dtmax) B%dt(i,j)=dtmax
         if(B%dt(i,j) .lt. dtmin) B%dt(i,j)=dtmin
      endif
	 enddo
	 enddo
  end subroutine comput_dt

!----------------------------------------------------------
! 计算无粘性及粘性项的谱半径，在加速收敛技术（局部时间步长，隐格式，残差光顺等）中使用
   subroutine comput_Lvc(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real(sp) :: C0,si,sj,s0,ni1,ni2,nj1,nj2,uni,unj,tmp1
     Type (Block_TYPE),pointer:: B
    
	 C0=1.d0
     B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny

     do j=1,ny-1
     do i=1,nx-1
 		 s0 =B%vol(i,j)         ! 面积
		 si=0.5d0*(B%si(i,j)+B%si(i+1,j))         ! 控制体边长（两侧平均）
		 ni1=0.5d0*(B%ni1(i,j)+B%ni1(i+1,j))      ! 界面法方向（两侧平均）
		 ni2=0.5d0*(B%ni2(i,j)+B%ni2(i+1,j))
         sj=0.5d0*(B%sj(i,j)+B%sj(i,j+1))
		 nj1=0.5d0*(B%nj1(i,j)+B%nj1(i,j+1))
		 nj2=0.5d0*(B%nj2(i,j)+B%nj2(i,j+1))
         uni=uu(i,j)*ni1+v(i,j)*ni2             !法向速度
		 unj=uu(i,j)*nj1+v(i,j)*nj2

!        谱半径
		 B%Lci(i,j)=(abs(uni)+cc(i,j))*si          ! 无粘项Jocabian矩阵的谱半径 （Blazek's Book 6.1.4节）
		 B%Lcj(i,j)=(abs(unj)+cc(i,j))*sj

		 tmp1=gamma/d(i,j)*(B%Amu(i,j)/Pr+B%Amu_t(i,j)/PrT)
		 B%Lvi(i,j)=tmp1*si*si/s0         ! 粘性项Jocabian矩阵谱半径
         B%Lvj(i,j)=tmp1*sj*sj/s0
 	 enddo
	 enddo
  end subroutine comput_Lvc
