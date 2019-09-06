!---------------------------------------------------------------------------------------------
! ����LU_SGS������ʱ���ƽ�һ��ʱ�䲽 ����nMesh������ �ĵ�������
 subroutine GSIS_LU_SGS(nMesh)
   use Global_var
   use const_var ,only :sp,dop
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(sp) :: du

    call Comput_Residual_GSIS(nMesh)     ! ���������ϼ���в�

    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock)                          ! ����LU_SGS��������DU=U(n+1)-U(n)
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   ʱ���ƽ� 
      do j=1,ny-1
      do i=3,nx-3
      do m=1,4
        du=B%du(m,i,j)  
        B%U(m,i,j)=B%U(m,i,j)+du                      ! U(n+1)=U(n)+dU
      enddo
      enddo
      enddo
    enddo    





    end subroutine GSIS_LU_SGS

    
    !-------------------------------------------------------------------------------------------------------   
! ����в�����ȫ���飩
   Subroutine Comput_Residual_GSIS(nMesh)
   use Global_Var
   use Flow_Var 
   use const_var ,only :sp,dop
   implicit none
   real(sp) :: Res
   integer:: nMesh,mBlock,nx,ny,i,j,m,Nvar1
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B

   MP=>Mesh(1)               ! ��������
   do mBlock=1,MP%Num_Block
      B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
      B%Res_max(:) =0.d0           ! ���в�
      B%Res_rms(:) =0.d0           ! �������в�    
     nx=B%nx; ny=B%ny
!    ������ʱ����
     allocate(d(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),uu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), &
              v(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),T(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),  &
              p(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), cc(nx-1,ny-1))
     allocate(Fluxi(4,nx,ny),Fluxj(4,nx,ny))


!------------------------------------------------------------------------------------     
	 call comput_duvtpckw(nMesh,mBlock)    ! ��������� d,u,v,T,p,cc,Kt,Wt
!---------------------------------------------------------------------------------------
!  ���N-S���̵������ģ��  
!  ����һ�������Ĳв�Ҷ�� ; ��nMesh ������ĵ�mBlock��  ������ģ��Ҳ�ڸÿ��м��㣩  

   call Residual_GSIS (1,mBlock)       

!---------------------------------------------------------------------------------------- 
   !--------------------------------------------------------------------------------------
!   ͳ�����в�;������в� 
    do j=4,ny-4
    do i=4,nx-4
! -------------------------------------------------------------------------------------------
!    ʱ���ƽ�
       do m=1,MP%Nvar
        Res=B%Res(m,i,j)
!--------------------------------------------------------------------------------------------------
        if(abs(Res) .gt. B%Res_max(m))  B%Res_max(m)=abs(Res)          ! ���в�
        B%Res_rms(m)=B%Res_rms(m)+Res*Res                           ! �������в�      
!--------------------------------------------------------------------------------------------------       
	   enddo
    enddo
    enddo
   
   
     call comput_Lvc(1,mBlock)
	 call comput_dt(1,mBlock)        ! ����(����) ʱ�䲽��
 
 !   ɾ������ʱ����    
	 deallocate(d,uu,v,T,cc,p,Fluxi,Fluxj)
 
   enddo    
   
  
    end Subroutine Comput_Residual_GSIS
    
       subroutine Residual_GSIS(nMesh,mBlock)
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
   
   call get_viscous(nMesh,mBlock)         ! �������ճ��ϵ��
   B%Amu_t(:,:)=0.d0         
   nx=B%nx ; ny=B%ny

   Scheme=MP%Iflag_Scheme         ! ��ֵ��ʽ ����ͬ�����ϲ��ò�ͬ��ʽ��
   IFlux=MP%Iflag_Flux            ! ͨ������
   Reconstruction=MP%IFlag_Reconstruction  ! �ع���ʽ
  
   Fluxi=0.d0  ! ����
   Fluxj=0.d0


continue

!----- i- direction ----------------------------------------------------------------------------------

   do j=1,ny-1
   do i=3,nx-2
       
	 si=B%si(i,j)   ! ������߽糤��
	 ni1=B%ni1(i,j); ni2=B%ni2(i,j)   ! ������
!-----��ճ�� -----------Inviscous term-------------------------------------------------------------

!-----ѡ���ع�����ֵ������ֵ�� ��ʽ  ��ʹ��ԭʼ�������غ���������������ع���

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
   

!-------������ת����ֱ�ڽ��� ������-���� ����ϵ��
      QL(1)=UL(1); QL(2)=UL(2)*ni1+UL(3)*ni2 ; QL(3)= -UL(2)*ni2+UL(3)*ni1 ; QL(4)=UL(4) ! �ܶȡ�ѹ���������ٶȡ������ٶ� ����ֵ��
      QR(1)=UR(1); QR(2)=UR(2)*ni1+UR(3)*ni2 ; QR(3)= -UR(2)*ni2+UR(3)*ni1 ; QR(4)=UR(4) ! �ܶȡ�ѹ���������ٶȡ������ٶ� ����ֵ��

!-------ѡ��ͨ�����ѷ���----(Steger_Warming,Van Leer, Roe, AUSM, HLL, HLLC----
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
!   ���ͨ�� ���任��x-y����ϵ��  
        Fluxi(1,i,j)=-Flux0(1)*si                            ! ����ͨ��
        Fluxi(2,i,j)=-(Flux0(2)*ni1-Flux0(3)*ni2)*si         ! x-������ͨ��
        Fluxi(3,i,j)=-(Flux0(2)*ni2+Flux0(3)*ni1)*si         ! y-������ͨ��
        Fluxi(4,i,j)=-Flux0(4)*si                            ! ����ͨ��
        
        !write(*,*) i,Fluxi(1,i,j),Fluxi(2,i,j),Fluxi(3,i,j),Fluxi(4,i,j)
      continue
	  if (If_viscous .eq. 0) cycle   ! �����Euler���̣���ճ����������ճ�������

!--i-������ճͨ���������������i-����ճ��ͨ��
!--------------------------------------------------------------------------------------------------------- 
!----------- Viscous term---------------------------------------------------------------------------------  
! ��ɢϵ����ճ��ϵ�����ȴ���ϵ����= ��������ֵ��ƽ�� ���߽紦���õ���ֵ)    
    if( i .eq. 1) then 
     Amu1=B%Amu(i,j) + B%Amu_t(i,j)                  ! ճ��ϵ�� (����+����), �����ϵ�ֵ=����ֵ��ƽ��
     Amk1=Cp*(B%Amu(i,j)/Pr + B%Amu_t(i,j)/PrT )     ! �ȴ���ϵ��
    else if (i .eq. nx ) then
      Amu1=B%Amu(i-1,j) + B%Amu_t(i-1,j)                  ! ճ��ϵ�� (����+����), �����ϵ�ֵ=����ֵ��ƽ��
      Amk1=Cp*(B%Amu(i-1,j)/Pr + B%Amu_t(i-1,j)/PrT )     ! �ȴ���ϵ��
    else
      Amu1=(B%Amu(i-1,j)+B%Amu(i,j) + B%Amu_t(i-1,j)+B%Amu_t(i,j) )*0.5d0                 ! ճ��ϵ�� (����+����), �����ϵ�ֵ=����ֵ��ƽ��
      Amk1=Cp*((B%Amu(i-1,j)+B%Amu(i,j))/Pr + (B%Amu_t(i-1,j)+B%Amu_t(i,j))/PrT )*0.5d0   ! �ȴ���ϵ��
    endif

!----Jocabianϵ�� ����������Լ�������ĵ���, ���ڼ����������ĵ�����
    Dix=B%x1(i,j)-B%x1(i-1,j)
    Diy=B%y1(i,j)-B%y1(i-1,j)
    Djx=(B%x1(i-1,j+1)+B%x1(i,j+1)-B%x1(i-1,j-1)-B%x1(i,j-1))*0.25d0
    Djy=(B%y1(i-1,j+1)+B%y1(i,j+1)-B%y1(i-1,j-1)-B%y1(i,j-1))*0.25d0
    Ds=1.d0/(Dix*Djy-Djx*Diy)
! �������Լ�������ĵ���    
    Diu=uu(i,j)-uu(i-1,j)
    Div=v(i,j)-v(i-1,j)
    DiT=T(i,j)-T(i-1,j)
    Dju=(uu(i-1,j+1)+uu(i,j+1)-uu(i-1,j-1)-uu(i,j-1))*0.25d0
    Djv=(v(i-1,j+1)+v(i,j+1)-v(i-1,j-1)-v(i,j-1))*0.25d0
    DjT=(T(i-1,j+1)+T(i,j+1)-T(i-1,j-1)-T(i,j-1))*0.25d0
! ��������x,y����ĵ���
    ux=(Diu*Djy-Dju*Diy)*Ds
    vx=(Div*Djy-Djv*Diy)*Ds
    Tx=(DiT*Djy-DjT*Diy)*Ds
    uy=(-Diu*Djx+Dju*Dix)*Ds
    vy=(-Div*Djx+Djv*Dix)*Ds
    Ty=(-DiT*Djx+DjT*Dix)*Ds
!  ճ��Ӧ��������ͨ��
    u1=(uu(i,j)+uu(i-1,j))*0.5d0
    v1=(v(i,j)+v(i-1,j))*0.5d0
    
    t11=((4.d0/3.d0)*ux-(2.d0/3.d0)*vy)*Amu1
    t22=((4.d0/3.d0)*vy-(2.d0/3.d0)*ux)*Amu1
    t12=(uy+vx)*Amu1
    t33=Amu1*(-(2.d0/3.d0)*ux-(2.d0/3.d0)*vy)
    qx=Amk1*Tx
    qy=Amk1*Ty

     
      if((B%solver==GKUA))then
         t11=t11+B%TaoNSi_G(1,i,j)
         t22=t22+B%TaoNSi_G(2,i,j)
         t12=t12+B%TaoNSi_G(3,i,j)
         t33=t33+B%TaoNSi_G(4,i,j)
         
         qx=qx+B%Qi_G(1,i,j)
         qy=qy+B%Qi_G(2,i,j)
      end if 
        !B%TaoNSi(1,i,j)=t11
        !B%TaoNSi(2,i,j)=t22
        !B%TaoNSi(3,i,j)=t12
        !B%TaoNSi(4,i,j)=t33
        !B%Qi(1,i,j)=qx
        !B%Qi(2,i,j)=qy
      
    
    E1=u1*t11+v1*t12+qx
    E2=u1*t12+v1*t22+qy
! ���ճ��ͨ��
    Fluxi(2,i,j)=Fluxi(2,i,j)   +(t11*ni1+t12*ni2)*si   
    Fluxi(3,i,j)=Fluxi(3,i,j)   +(t12*ni1+t22*ni2)*si
    Fluxi(4,i,j)=Fluxi(4,i,j)  +(E1*ni1+ E2*ni2)*si
  
   enddo
   enddo
   
!==================================================================================================================
!                                         j-�������ճ��ճ��ͨ��    
!-----------------------------------------j- direction -------------------------------------------------------------	 
   do j=1,ny
   do i=1,nx-1
  
   ! �߳���������   
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
 
!-------ѡ��ͨ�����ѷ���----(Steger_Warming,Van Leer, Roe, AUSM, HLL, HLLC)----
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
     
!   ��ճͨ��   	   
    Fluxj(1,i,j)=-Flux0(1)*sj 
    Fluxj(2,i,j)=-(Flux0(2)*nj1-Flux0(3)*nj2)*sj
    Fluxj(3,i,j)=-(Flux0(2)*nj2+Flux0(3)*nj1)*sj
    Fluxj(4,i,j)=-Flux0(4)*sj

    if (If_viscous .eq. 0) cycle  ! �����Euler���̣���ճ����������ճ�������
 !---------Viscous term -----------------------------------------------------------------------------
 ! ��ɢϵ����ճ��ϵ�����ȴ���ϵ����= ��������ֵ��ƽ�� ���߽紦���õ���ֵ)    
    if( j .eq. 1) then
     Amu1=B%Amu(i,j)+B%Amu_t(i,j)
     Amk1=Cp*(B%Amu(i,j)/Pr +B%Amu_t(i,j)/PrT )   ! �ȴ���ϵ��
    else if (j .eq. ny) then
      Amu1=B%Amu(i,j-1)+B%Amu_t(i,j-1)
      Amk1=Cp*(B%Amu(i,j-1)/Pr +B%Amu_t(i,j-1)/PrT )   ! �ȴ���ϵ��
    else
     Amu1=(B%Amu(i,j)+B%Amu(i,j-1) +B%Amu_t(i,j)+B%Amu_t(i,j-1))*0.5d0
     Amk1=Cp*((B%Amu(i,j-1)+B%Amu(i,j))/Pr + (B%Amu_t(i,j-1)+B%Amu_t(i,j))/PrT )*0.5d0   ! �ȴ���ϵ��
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
    
  
    !if((B%solver==GKUA))then
    !     t11=t11+B%TaoNSj_G(1,i,j)
    !     t22=t22+B%TaoNSj_G(2,i,j)
    !     t12=t12+B%TaoNSj_G(3,i,j)
    !     t33=t33+B%TaoNSj_G(4,i,j)
    !     
    !     qx=qx+B%Qj_G(1,i,j)
    !     qy=qy+B%Qj_G(2,i,j)
    !end if
    !B%TaoNSj(1,i,j)=t11
    !B%TaoNSj(2,i,j)=t22
    !B%TaoNSj(3,i,j)=t12
    !B%TaoNSj(4,i,j)=t33
    !
    !B%Qj(1,i,j)=qx  
    !B%Qj(2,i,j)=qy
 
    
    u1=(uu(i,j)+uu(i,j-1))*0.5d0
    v1=(v(i,j)+v(i,j-1))*0.5d0

    E1=u1*t11+v1*t12+qx
    E2=u1*t12+v1*t22+qy

    Fluxj(2,i,j)=Fluxj(2,i,j) +(t11*nj1+t12*nj2)*sj
    Fluxj(3,i,j)=Fluxj(3,i,j) +(t12*nj1+t22*nj2)*sj
    Fluxj(4,i,j)=Fluxj(4,i,j) +(E1*nj1+ E2*nj2)*sj

 
   enddo
   enddo
  
 !-------����в� ����������----------------------------------------------------  
 !                                    ����ǿ�Ⱥ��� ����������ʱ�������ǿ�Ⱥ�����
    do j=1,ny-1
    do i=3,nx-4
    do m=1,4
    !  B%Res(m,i,j)= Fluxi(m,i+1,j)-Fluxi(m,i,j)+Fluxj(m,i,j+1)-Fluxj(m,i,j)
        B%Res(m,i,j)= Fluxi(m,i+1,j)-Fluxi(m,i,j)
       
    enddo
  !   write(*,*)i,j,B%Res(1,i,j),B%Res(2,i,j),B%Res(3,i,j),B%Res(4,i,j)
	enddo
	enddo
continue
  end  
