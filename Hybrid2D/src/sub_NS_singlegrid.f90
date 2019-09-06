
!---------------------------------------------------------------------------------------------
! ����LU_SGS������ʱ���ƽ�һ��ʱ�䲽 ����nMesh������ �ĵ�������
 subroutine NS_Time_advance_LU_SGS(nMesh)
   use Global_var
   use const_var ,only :sp,dop
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(sp) :: du

    call Comput_Residual_one_mesh(nMesh)     ! ���������ϼ���в�
!    if(nMesh .ne. 1) call Add_force_function(nMesh)                !  ���ǿ�Ⱥ�������������Ĵ�����ʹ��,ĿǰLU-SGS��֧�֣�

    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock)                          ! ����LU_SGS��������DU=U(n+1)-U(n)
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   ʱ���ƽ� 
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
        !call Boundary_condition_onemesh(nMesh)         ! �߽����� ���趨Ghost Cell��ֵ��
!!      call update_buffer_onemesh(nMesh)              ! ͬ������Ľ�����



 end subroutine NS_Time_advance_LU_SGS


!-------------------------------------------------------------------------------------------------------   
! ����в�����ȫ���飩
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
       MP=>Mesh(1)               ! ��������
        

!----------------------------------------------
   do mBlock=1,MP%Num_Block
      B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
      B%Res_max(:) =0.d0           ! ���в�
      B%Res_rms(:) =0.d0           ! �������в�    
     nx=B%nx; ny=B%ny
!    ������ʱ����
!     allocate(d(0:nx,0:ny),uu(0:nx,0:ny),v(0:nx,0:ny),T(0:nx,0:ny),cc(0:nx,0:ny),p(0:nx,0:ny))   ! Bug found, 2012-5-1
     allocate(d(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),uu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), &
              v(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),T(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),  &
              p(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), cc(nx-1,ny-1))
     allocate(Fluxi(4,nx,ny),Fluxj(4,nx,ny))


!------------------------------------------------------------------------------------     
	 call comput_duvtpckw(nMesh,mBlock)    ! ��������� d,u,v,T,p,cc,Kt,Wt
!---------------------------------------------------------------------------------------
!  ���N-S���̵������ģ��  
!  ����һ�������Ĳв�Ҷ�� ; ��nMesh ������ĵ�mBlock��  ������ģ��Ҳ�ڸÿ��м��㣩  

   call Residual (nMesh,mBlock)       
  
     
!--------------------------------------------------------------------------------------
!   ͳ�����в�;������в� 
    do j=3,ny-3
    do i=3,nx-3
! -------------------------------------------------------------------------------------------
!    ʱ���ƽ�
       do m=1,MP%Nvar
        Res=B%Res(m,i,j)

!--------------------------------------------------------------------------------------------------
        if(abs(Res) .gt. B%Res_max(m))then
            B%Res_max(m)=abs(Res)          ! ���в�
            B%Resmax_location_i=i
            B%Resmax_location_j=j
        end if
        
        B%Res_rms(m)=B%Res_rms(m)+Res*Res                           ! �������в�      
!--------------------------------------------------------------------------------------------------       
	   enddo
    enddo
    enddo
            

!---------------------------------------------------------------------------------------- 
     call comput_Lvc(nMesh,mBlock)
	 call comput_dt(nMesh,mBlock)        ! ����(����) ʱ�䲽��
 
 !   ɾ������ʱ����    
	 deallocate(d,uu,v,T,cc,p,Fluxi,Fluxj)
 
   enddo    
   
   B%Res_rms(:)=sqrt(B%Res_rms(:)/(1.d0*Mesh(nMesh)%Num_Cell))   !    ȫ���������ܾ������в�
  end Subroutine Comput_Residual_one_mesh



!---------------------------------------------------------
!  �����غ��������������� (d,u,v,T,p,c) 
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
     B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
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
! ���㣨���أ�ʱ�䲽��
   subroutine comput_dt(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real(sp) ,parameter:: C0=1
     Type (Block_TYPE),pointer:: B
    
     B => Mesh(nMesh)%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny

     do j=1,ny-1
     do i=1,nx-1
      if(Iflag_local_dt .eq. 0) then   ! ȫ��ʱ�䲽��
        B%dt(i,j)=dt_global   
       else                             ! ����ʱ�䲽��
         B%dt(i,j)=CFL*B%vol(i,j)/(B%Lci(i,j)+B%Lcj(i,j) +C0*(B%Lvi(i,j)+B%Lvj(i,j)) )     ! �ֲ�ʱ�䲽��
         if(B%dt(i,j) .gt. dtmax) B%dt(i,j)=dtmax
         if(B%dt(i,j) .lt. dtmin) B%dt(i,j)=dtmin
      endif
	 enddo
	 enddo
  end subroutine comput_dt

!----------------------------------------------------------
! ������ճ�Լ�ճ������װ뾶���ڼ��������������ֲ�ʱ�䲽��������ʽ���в��˳�ȣ���ʹ��
   subroutine comput_Lvc(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real(sp) :: C0,si,sj,s0,ni1,ni2,nj1,nj2,uni,unj,tmp1
     Type (Block_TYPE),pointer:: B
    
	 C0=1.d0
     B => Mesh(nMesh)%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny

     do j=1,ny-1
     do i=1,nx-1
 		 s0 =B%vol(i,j)         ! ���
		 si=0.5d0*(B%si(i,j)+B%si(i+1,j))         ! ������߳�������ƽ����
		 ni1=0.5d0*(B%ni1(i,j)+B%ni1(i+1,j))      ! ���淨��������ƽ����
		 ni2=0.5d0*(B%ni2(i,j)+B%ni2(i+1,j))
         sj=0.5d0*(B%sj(i,j)+B%sj(i,j+1))
		 nj1=0.5d0*(B%nj1(i,j)+B%nj1(i,j+1))
		 nj2=0.5d0*(B%nj2(i,j)+B%nj2(i,j+1))
         uni=uu(i,j)*ni1+v(i,j)*ni2             !�����ٶ�
		 unj=uu(i,j)*nj1+v(i,j)*nj2

!        �װ뾶
		 B%Lci(i,j)=(abs(uni)+cc(i,j))*si          ! ��ճ��Jocabian������װ뾶 ��Blazek's Book 6.1.4�ڣ�
		 B%Lcj(i,j)=(abs(unj)+cc(i,j))*sj

		 tmp1=gamma/d(i,j)*(B%Amu(i,j)/Pr+B%Amu_t(i,j)/PrT)
		 B%Lvi(i,j)=tmp1*si*si/s0         ! ճ����Jocabian�����װ뾶
         B%Lvj(i,j)=tmp1*sj*sj/s0
 	 enddo
	 enddo
  end subroutine comput_Lvc
