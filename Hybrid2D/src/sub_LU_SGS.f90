!  采用LU-SGS方法，计算DU=U(n+1)-U(n)
!  2D Code by Li Xinliang, 2012-4-29
!-----------------------------------------------------------------------------------------
    subroutine  du_LU_SGS_2D(nMesh,mBlock)                          ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
    use Global_Var
    use Flow_Var 
    
    implicit none
	integer:: nMesh,mBlock,Nvar1,nx,ny,nz,plane,i,j,m

    real(sp) ,dimension(:),allocatable:: alfa,dui,duj,DF   ! alfa,对角线元素值;  dui,duj,DF对流通量
	real(sp) ,parameter:: w_LU=1.d0   !  LU_SGS的松弛因子(1-2之间), 增大w_LU会提高稳定性，但会降低收敛速度
	Type (Block_TYPE),pointer:: B
    TYPE (Mesh_TYPE),pointer:: MP
     MP=> Mesh(nMesh)
	 NVAR1=MP%NVAR
	 B => MP%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny
     allocate(alfa(NVAR1),dui(NVAR1),duj(NVAR1),DF(NVAR1)) 

! LU-SGS的两次扫描
!----------------------------------
!   从i=1,j=1 到i=nx-1,j=ny-1的扫描过程  (向上扫描过程)
!   扫描 i+j=plane 的平面
!   w_LU是松弛因子（1到2之间），增大w_LU会提高稳定性，但会降低收敛速度
   do plane=6,nx+ny-6            

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m,alfa,dui,duj,DF)
   do j=3,ny-3
   i=plane-j
	 if( i .lt. 3 .or. i .gt. nx-3) cycle    ! 超出了这个平面
    
	  alfa(:)=B%vol(i,j)/B%dt(i,j)+w_LU*(B%Lci(i,j)+B%Lcj(i,j))               ! 对角元
	  if(If_viscous .eq. 1)  alfa(:)=alfa(:)+2.d0*(B%Lvi(i,j)+B%Lvj(i,j))           ! 加上粘性项

       if(NVAR1 .eq. 5) then
	     alfa(5)=alfa(5)+B%Lvi(i,j)+B%Lvj(i,j)
       elseif(NVAR1 .eq. 6) then
	     alfa(5)=alfa(5)+0.09*B%U(6,i,j)/B%U(1,i,j)
         alfa(6)=alfa(6)+2.d0*0.0828*B%U(6,i,j)/B%U(1,i,j)
	   endif

	  		   
	 if(i.ne. 3) then
!                                                      通量的差量，用来近似计算A*W (See Blazek's book, page 208)
       call comput_DFn(NVAR1,DF,B%U(:,i-1,j),B%DU(:,i-1,j),B%ni1(i,j),B%ni2(i,j),gamma)  
       dui(:)=0.5d0*(DF(:)*B%si(i,j)+w_LU*B%Lci(i-1,j)*B%dU(:,i-1,j))
       if(If_viscous .eq. 1)    dui(:)=dui(:)+B%Lvi(i-1,j)*B%dU(:,i-1,j)         ! 2012-2-29, 稳定性更好些
     else
	   dui=0.d0                             ! 左侧没有点
     endif
	 
	 if(j.ne.3) then
       call comput_DFn(NVAR1,DF,B%U(:,i,j-1),B%DU(:,i,j-1),B%nj1(i,j),B%nj2(i,j),gamma)  
       duj(:)=0.5d0*(DF(:)*B%sj(i,j)+w_LU*B%Lcj(i,j-1)*B%dU(:,i,j-1))

        if(If_viscous .eq. 1)    duj(:)=duj(:)+B%Lvj(i,j-1)*B%dU(:,i,j-1)    ! 2012-2-29

	 else
	   duj=0.d0
	 endif

	do m=1,NVAR1
    B%dU(m,i,j)=(B%Res(m,i,j)+dui(m)+duj(m))/alfa(m)
	enddo

   enddo
!$OMP END PARALLEL DO 
   enddo
!----------------------------------------------------------
!  从 (nx-1,ny-1)到(1,1)的扫描过程 （向下扫描过程）
!  plane=i+j+k
   do plane=nx+ny-6,6,-1   

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,m,alfa,dui,duj,duk,DF)
	 do j=ny-3,3,-1
	 i=plane-j
	 if( i .lt. 3 .or. i .gt. nx-3) cycle            ! 超出了这个平面
      alfa(:)=B%vol(i,j)/B%dt(i,j)+w_LU*(B%Lci(i,j)+B%Lcj(i,j))
      if(If_viscous .eq. 1)  alfa(:)=alfa(:)+2.d0*(B%Lvi(i,j)+B%Lvj(i,j) )         
       if(NVAR1 .eq. 5) then
         alfa(5)=alfa(5)+B%Lvi(i,j)+B%Lvj(i,j)
       elseif(NVAR1 .eq. 6) then
 	     alfa(5)=alfa(5)+0.09*B%U(6,i,j)/B%U(1,i,j)
         alfa(6)=alfa(6)+2.d0*0.0828*B%U(6,i,j)/B%U(1,i,j)
	   endif

	 if(i.ne. nx-3) then
!                              通量的差量，用来近似计算A*W (See Blazek's book, page 208)
       call comput_DFn(NVAR1,DF,B%U(:,i+1,j),B%DU(:,i+1,j),B%ni1(i,j),B%ni2(i,j),gamma)  
       dui(:)=-0.5d0*(DF(:)*B%si(i+1,j)-w_LU*B%Lci(i+1,j)*B%dU(:,i+1,j))
       if(If_viscous .eq. 1)    dui(:)=dui(:)+B%Lvi(i+1,j)*B%dU(:,i+1,j)
	   else
	   dui=0.d0
       endif
	 
	 if(j.ne. ny-3) then
       call comput_DFn(NVAR1,DF,B%U(:,i,j+1),B%DU(:,i,j+1),B%nj1(i,j),B%nj2(i,j),gamma)  
       duj(:)=-0.5d0*(DF(:)*B%sj(i,j+1)-w_LU*B%Lcj(i,j+1)*B%dU(:,i,j+1))
       if(If_viscous .eq. 1)    duj(:)=duj(:)+B%Lvj(i,j+1)*B%dU(:,i,j+1)
	 else
	   duj=0.d0
	 endif

	do m=1,NVAR1
      B%dU(m,i,j)=B%dU(m,i,j)+(dui(m)+duj(m))/alfa(m)
	enddo
   enddo
!$OMP END PARALLEL DO 

   enddo
   deallocate(dui,duj,DF,alfa)
   end subroutine  du_LU_SGS_2D







!  计算通量的差量 DF=F(U+DU)-F(U),  LU-SGS方法中使用，用来近似A*DU  
    subroutine comput_DFn(NV,DF,U,DU,n1,n2,gamma)
    use const_var ,only :sp,dop
	implicit none
    integer:: NV
	real(sp) ,dimension(NV):: DF,U,DU,U2
	real(sp) :: n1,n2,un1,un2,gamma,p1,p2
    U2=U+DU                                 ! 新的守恒变量
	un1=(U(2)*n1+U(3)*n2)/U(1)              !un 法向速度
	p1=(gamma-1.d0)*(U(4)-0.5d0*(U(2)**2+U(3)**2)/U(1))   ! 压力
	un2=(U2(2)*n1+U2(3)*n2)/U2(1)          !法向速度un
	p2=(gamma-1.d0)*(U2(4)-0.5d0*(U2(2)**2+U2(3)**2)/U2(1))    !压力
!  通量之差 DF=F(U+DU)-F(U)
    DF(1)=U2(1)*un2-U(1)*un1                  ! d*un
	DF(2)=(U2(2)*un2+p2*n1)-(U(2)*un1+p1*n1)
	DF(3)=(U2(3)*un2+p2*n2)-(U(3)*un1+p1*n2)
	DF(4)=(U2(4)+p2)*un2-(U(4)+p1)*un1
     if(NV .eq. 5) then
      DF(5)=U2(5)*un2-U(5)*un1
	 elseif (NV .eq. 6) then       
	  DF(5)=U2(5)*un2-U(5)*un1   ! k方程的对流通量
	  DF(6)=U2(6)*un2-U(6)*un1   ! w方程的对流通量
	 endif
	end subroutine comput_dFn