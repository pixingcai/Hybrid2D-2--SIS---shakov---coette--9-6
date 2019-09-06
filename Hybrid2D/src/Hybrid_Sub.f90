
    !---------------------------------------------------------------------------------------------
! Gsis main subroutine for DVM
    subroutine Hybrid_Time_advance_LU_SGS(nMesh)
    use Global_var
    use const_var 
    use Com_ctrl
    implicit none
    integer::iter_count
    integer::nMesh,mBlock,i,j,m,nx,ny,iv,jv,iv1,jv1,k,flag
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    real(sp):: res_iter,temp1,temp2
    real(sp) :: du,temprho,tempU,tempV,tempT,temprhoNS,tempUNS,tempVNS,tempTNS,crt1,vjv1,viu1,vij21,GDVijv1,vjv2,viu2,vij22,GDVijv2,crt2
    real(sp) :: txx_tmp,tyy_tmp,txy_tmp,tzz_tmp,qx_tmp,qy_tmp,Coef
    real(sp) :: txx_tmp1,tyy_tmp1,txy_tmp1,tzz_tmp1,qx_tmp1,qy_tmp1,T_tmp1,P_tmp1
    real(sp) :: txx_tmp2,tyy_tmp2,txy_tmp2,tzz_tmp2,qx_tmp2,qy_tmp2,T_tmp2,P_tmp2
    real(sp):: coe11,coe12,coe21,coe22,Key1,Key2
    real(sp):: coe11_2,coe12_2,coe21_2,coe22_2,Key12,Key22,rho_L,rhou_L,rhov_L,e_L,rho_R,rhou_R,rhov_R,e_R
    res_iter=1.0
    iter_count=0
    MP=>Mesh(1) 
   
    B => MP%Block(1)
    nx=B%nx; ny=B%ny

  
  
    do mBlock=1,MP%Num_Block
        B => MP%Block(mBlock)
        nx=B%nx; ny=B%ny
        do i=1,nx-1
        do j=1,ny-1
            B%S0_GKUA(1:10,i,j)=B%S_GKUA(1:10,i,j)         
        end do
        end do
        
    end do
    call GKUA_Phy_boundary
    call Hybrid_GKUA_Comput_Block(1)
    call DETFP(Time_Method)
    call MassCorrect
    call Convert_NS(1)
    
    do mBlock=1,MP%Num_Block
            B => MP%Block(mBlock)
            nx=B%nx; ny=B%ny
            do i=1,nx-1
            do j=1,ny-1
                B%S0_GKUA(1:10,i,j)=B%S_GKUA(1:10,i,j)-B%S0_GKUA(1:10,i,j)
                B%S00_GKUA(1:10,i,j)=B%S_GKUA(1:10,i,j)
            end do
            end do
        
    end do  
    

    

    do mBlock=1,MP%Num_Block
        B => MP%Block(mBlock)
        nx=B%nx; ny=B%ny
        do i=3,nx-3
        do j=3,ny-3
            temprho=B%S0_GKUA(1,i,j)
            tempU=B%S0_GKUA(2,i,j)/Ma*sqrt(2.0/gamma)
            tempV=B%S0_GKUA(3,i,j)/Ma*sqrt(2.0/gamma)
            tempT=B%S0_GKUA(4,i,j)
            
            temprhoNS=B%U0(1,i,j)
            tempUNS=B%U0(2,i,j)/B%U0(1,i,j)
            tempVNS=B%U0(3,i,j)/B%U0(1,i,j)
            tempTNS=(B%U0(4,i,j)-0.5*B%U0(1,i,j)*(tempUNS*tempUNS+tempVNS*tempVNS))/(Cv*B%U0(1,i,j))
            
            temprhoNS=temprhoNS+temprho
            tempUNS=tempUNS+tempU
            tempVNS=tempVNS+tempV
            tempTNS=tempTNS+tempT
            
            B%U(1,i,j)=temprhoNS
            B%U(2,i,j)=temprhoNS*tempUNS
            B%U(3,i,j)=temprhoNS*tempVNS
            B%U(4,i,j)=Cv*tempTNS*temprhoNS+0.5*temprhoNS*(tempUNS*tempUNS+tempVNS*tempVNS)
      
            B%U2(1,i,j)=B%U(1,i,j)
            B%U2(2,i,j)=B%U(2,i,j)
            B%U2(3,i,j)=B%U(3,i,j)
            B%U2(4,i,j)=B%U(4,i,j)
         
     
        end do
        end do
        
    end do
    
   

  
    do while(res_iter.gt.1E-6)
        iter_count=iter_count+1  
        call NS_Time_advance_LU_SGS(1)
        call get_res(res_iter)

    end do
    MP%iter=iter_count
    print * ,iter_count
    
    do mBlock=1,MP%Num_Block
          B => MP%Block(mBlock)
          nx=B%nx; ny=B%ny
          do i=3,nx-3
          do j=3,ny-3
            B%U0(1:4,i,j)=B%U(1:4,i,j)  
          end do
          end do
          
    end do 

    
    do mBlock=1,MP%Num_Block
        B => MP%Block(mBlock)
        nx=B%nx; ny=B%ny
        do i=3,nx-3
        do j=3,ny-3
            temprhoNS=B%U(1,i,j)
            tempUNS=(B%U(2,i,j)/B%U(1,i,j))
            tempVNS=(B%U(3,i,j)/B%U(1,i,j))
            tempTNS=(B%U(4,i,j)-0.5*B%U(1,i,j)*(tempUNS*tempUNS+tempVNS*tempVNS))/(Cv*B%U(1,i,j))
              
            temprho=B%U2(1,i,j)
            tempU=(B%U2(2,i,j)/B%U2(1,i,j))
            tempV=(B%U2(3,i,j)/B%U2(1,i,j))
            tempT=(B%U2(4,i,j)-0.5*B%U2(1,i,j)*(tempU*tempU+tempV*tempV))/(Cv*B%U2(1,i,j))
      
            B%S0_GKUA(1,i,j)=temprhoNS-temprho
            B%S0_GKUA(2,i,j)=(tempUNS-tempU)*Ma/sqrt(2.0/gamma)
            B%S0_GKUA(3,i,j)=(tempVNS-tempV)*Ma/sqrt(2.0/gamma)
            B%S0_GKUA(4,i,j)=tempTNS-tempT
              
            B%S_GKUA(1:4,i,j)=B%S_GKUA(1:4,i,j)+B%S0_GKUA(1:4,i,j)
              
            B%S_GKUA(5,i,j)= B%S_GKUA(5,i,j)+B%TaoNS_out(1,i,j)
            B%S_GKUA(6,i,j)= B%S_GKUA(6,i,j)+B%TaoNS_out(2,i,j)
            B%S_GKUA(7,i,j)= B%S_GKUA(7,i,j)+B%TaoNS_out(3,i,j)
            B%S_GKUA(10,i,j)= B%S_GKUA(10,i,j)+B%TaoNS_out(4,i,j)
            B%S_GKUA(8,i,j)= B%S_GKUA(8,i,j)+B%Q_out(1,i,j)
            B%S_GKUA(9,i,j)= B%S_GKUA(9,i,j)+B%Q_out(2,i,j)
              
            B%S0_GKUA(1:4,i,j)=B%S0_GKUA(1:4,i,j)+B%S00_GKUA(1:4,i,j)
        end do
        end do
          
    end do  
      

    Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! ʱ�� ��ʹ��ȫ��ʱ�䲽����ʱ�����壩
    Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! ���㲽��


    end subroutine Hybrid_Time_advance_LU_SGS
    

! ����в�����NS�飩
    Subroutine Hybrid_GKUA_Comput_Block(nMesh)
    use Global_Var
    use Flow_Var  
    use Com_ctrl
    use const_var ,only :sp,dop
    implicit none
    real(sp) :: Res
    integer:: nMesh,mBlock,nx,ny,i,j,m,Nvar1,jv,iv,ird
    integer :: iv1,jv1,NVIN,NVJN,index_1(4),index_2(4),index_3(4),index_4(4)
    real(sp)  :: Tij,crt,viu,vjv,vij2,GDVijv,err111(4)
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B

    MP=>Mesh(1)             

!----------------------------------------------
    do mBlock=1,MP%Num_Block
        B => MP%Block(mBlock)                 
        if(B%solver==GKUA)then
            nx=B%nx; ny=B%ny
            call GKUA_Residual (1,mBlock)       
        end if
    enddo   

    continue
    end Subroutine Hybrid_GKUA_Comput_Block
    
    Subroutine Convert_NS(nMesh)
   use Global_Var
   use Flow_Var  
   use Com_ctrl
   use const_var ,only :sp,dop
   implicit none
   real(sp) ::tempU,tempV,tempP,tempT,RONij , Uij, Vij, Ttraij,Trotij,Tvibij,Pij,toxxmid,toxymid,toyymid,tozzmid
   integer:: nMesh,mBlock,nx,ny,i,j,m,ksub,ijwall(2),ijwall_1(2),ijwall_0(2),pos_x,pos_y
   integer:: ibegin,iend,jbegin,jend,icheck,jcheck,Block_check,Block_i
   real(sp) :: s,swnx,swny,wnx,wny,taox,taoy,tao,pijt,qtratmp,qrottmp,qvibtmp,qtottmp,qtottmp1,xtmp,ytmp,qxtmp,qytmp
   real(sp) :: minmod,dx,dy,si,ni1,ni2,sj,nj1,nj2
   real(sp) :: Diu,Div,Dju,Djv,DiT,DjT,ux,vx,Tx,uy,vy,Ty,t11,t12,t22,t33,E1,E2,qx,qy
   real(sp) :: Dix,Diy,Djx,Djy,Ds,Amu1,u1,v1,Amk1,tmp
   real(sp) :: DiV1,DiV2,DiV3,DiV4,DiV5,DiV6,DiQ1,DiQ2,DiQ3
   real(sp) :: DjV1,DjV2,DjV3,DjV4,DjV5,DjV6,DjQ1,DjQ2,DjQ3
   real(sp) :: V1x,V1y,V2x,V2y,V3x,V3y,V4x,V4y,V5x,V5y,V6x,V6y
   real(sp) :: Q1x,Q1y,Q2x,Q2y,Q3x,Q3y
   real(sp) :: t11_g,t12_g,t22_g,t33_g,qx_g,qy_g,CF_g
   real(sp) :: u_g,v_g
   real(sp) :: dl,uul,vl,ppl,dr,uur,vr,ppr,pr1,temp1,temp2,temp3,temp4,temp5,temp6
   real(sp) ,pointer,dimension(:,:) :: filter_tmp
    character(len=50):: filename,filename_Dem
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   MP=>Mesh(1)
   do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
             if(B%solver==GKUA)then
                       B%S_GKUA(:,0,0)=B%S_GKUA(:,1,0)+B%S_GKUA(:,0,1)-B%S_GKUA(:,1,1)
                       B%S_GKUA(:,nx,0)=B%S_GKUA(:,nx,1)+B%S_GKUA(:,nx-1,0)-B%S_GKUA(:,nx-1,1)
                       B%S_GKUA(:,0,ny)=B%S_GKUA(:,1,ny)+B%S_GKUA(:,0,ny-1)-B%S_GKUA(:,1,ny-1)
                       B%S_GKUA(:,nx,ny)=B%S_GKUA(:,nx-1,ny)+B%S_GKUA(:,nx,ny-1)-B%S_GKUA(:,nx-1,ny-1)
                       do i=0,nx
                       do j=0,ny
                          B%U(1,i,j)=B%S_GKUA(1,i,j)
                          tempU=B%S_GKUA(2,i,j)/Ma*sqrt(2.0/gamma)
                          B%U(2,i,j)= B%U(1,i,j)*tempU
                          tempV=B%S_GKUA(3,i,j)/Ma*sqrt(2.0/gamma) 
                          B%U(3,i,j)= B%U(1,i,j)*tempV
                          tempT=B%S_GKUA(4,i,j)
                          B%U(4,i,j)=B%U(1,i,j)*Cv*tempT+0.5*B%U(1,i,j)*(tempU*tempU+tempV*tempv)      
                       end do
                       end do
                       B%U(:,0,:)=2.0*B%U(:,1,:)-B%U(:,2,:)
                       B%U(:,nx,:)=2.0*B%U(:,nx-1,:)-B%U(:,nx-2,:)
                       B%U(:,:,0)=2.0*B%U(:,:,1)-B%U(:,:,2)
                       B%U(:,:,ny)=2.0*B%U(:,:,ny-1)-B%U(:,:,ny-2)
                       
         
!    ������ʱ����
!     allocate(d(0:nx,0:ny),uu(0:nx,0:ny),v(0:nx,0:ny),T(0:nx,0:ny),cc(0:nx,0:ny),p(0:nx,0:ny))   ! Bug found, 2012-5-1
     allocate(d(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),uu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), &
              v(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),T(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),  &
              p(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), cc(nx-1,ny-1))
     allocate(Fluxi(4,nx,ny),Fluxj(4,nx,ny))


!------------------------------------------------------------------------------------     
	 call comput_duvtpckw(1,mBlock)    ! ��������� d,u,v,T,p,cc,Kt,Wt
     call get_viscous(1,mBlock)         ! �������ճ��ϵ��


!----- i- direction ----------------------------------------------------------------------------------
    B%Amu_t=0.0
   do j=1,ny-1
   do i=1,nx
	 si=B%si(i,j)   ! ������߽糤��
	 ni1=B%ni1(i,j); ni2=B%ni2(i,j)   ! ������   
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
    
    B%TaoNSi_in(1,i,j)=t11
    B%TaoNSi_in(2,i,j)=t22
    B%TaoNSi_in(3,i,j)=t12
    B%TaoNSi_in(4,i,j)=t33
    B%Qi_in(1,i,j)=qx   !!-qx
    B%Qi_in(2,i,j)=qy   !!-qy          
    if((i/=1).and.(i/=nx))then

              tmp= (B%S_GKUA(5,i,j)+B%S_GKUA(6,i,j)+B%S_GKUA(10,i,j)+B%S_GKUA(5,i-1,j)+B%S_GKUA(6,i-1,j)+B%S_GKUA(10,i-1,j))/6.0
              B%TaoNSi_G(1,i,j)= -((B%S_GKUA(5,i,j)+B%S_GKUA(5,i-1,j))/2.0-tmp)/(gamma*Ma**2)-t11
              B%TaoNSi_G(2,i,j)= -((B%S_GKUA(6,i,j)+B%S_GKUA(6,i-1,j))/2.0-tmp)/(gamma*Ma**2)-t22
              B%TaoNSi_G(3,i,j)= -(B%S_GKUA(7,i,j)+B%S_GKUA(7,i-1,j))/2.0/(gamma*Ma**2)-t12
              B%TaoNSi_G(4,i,j)=-((B%S_GKUA(10,i,j)+B%S_GKUA(10,i-1,j))/2.0-tmp)/(gamma*Ma**2)-t33
              
              B%Qi_G(1,i,j)=  -(B%S_GKUA(8,i-1,j)+B%S_GKUA(8,i,j))/2.0*sqrt(2.0)/(Ma**3*gamma**1.5)-qx
              B%Qi_G(2,i,j)=  -(B%S_GKUA(9,i-1,j)+B%S_GKUA(9,i,j))/2.0*sqrt(2.0)/(Ma**3*gamma**1.5)-qy
    else if(i==1)then
        
              tmp= (B%S_GKUA(5,i,j)+B%S_GKUA(6,i,j)+B%S_GKUA(10,i,j))/3.0
              B%TaoNSi_G(1,i,j)= -(B%S_GKUA(5,i,j)-tmp)/(gamma*Ma**2)-t11
              B%TaoNSi_G(2,i,j)= -(B%S_GKUA(6,i,j)-tmp)/(gamma*Ma**2)-t22
              B%TaoNSi_G(3,i,j)= -B%S_GKUA(7,i,j)/(gamma*Ma**2)-t12
              B%TaoNSi_G(4,i,j)= -(B%S_GKUA(10,i,j)-tmp)/(gamma*Ma**2)-t33
              
              B%Qi_G(1,i,j)=  -B%S_GKUA(8,i,j)*sqrt(2.0)/(Ma**3*gamma**1.5)-qx
              B%Qi_G(2,i,j)=  -B%S_GKUA(9,i,j)*sqrt(2.0)/(Ma**3*gamma**1.5)-qy
    else if(i==nx)then
        
              tmp= (B%S_GKUA(5,i-1,j)+B%S_GKUA(6,i-1,j)+B%S_GKUA(10,i-1,j))/3.0
              B%TaoNSi_G(1,i,j)= -(B%S_GKUA(5,i-1,j)-tmp)/(gamma*Ma**2)-t11
              B%TaoNSi_G(2,i,j)= -(B%S_GKUA(6,i-1,j)-tmp)/(gamma*Ma**2)-t22
              B%TaoNSi_G(3,i,j)= -B%S_GKUA(7,i-1,j)/(gamma*Ma**2)-t12
              B%TaoNSi_G(4,i,j)= -(B%S_GKUA(10,i-1,j)-tmp)/(gamma*Ma**2)-t33
              
              B%Qi_G(1,i,j)=  -B%S_GKUA(8,i-1,j)*sqrt(2.0)/(Ma**3*gamma**1.5)-qx
              B%Qi_G(2,i,j)=  -B%S_GKUA(9,i-1,j)*sqrt(2.0)/(Ma**3*gamma**1.5)-qy
            

        end if
        
   enddo
   enddo
!     call filter_i(1,mBlock)

!==================================================================================================================
!                                         j-�������ճ��ճ��ͨ��    
!-----------------------------------------j- direction -------------------------------------------------------------	 
   do j=1,ny
   do i=1,nx-1
  
   ! �߳���������   
	sj=B%sj(i,j)
	nj1=B%nj1(i,j); nj2=B%nj2(i,j)
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
   !!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Dix=(B%x1(i+1,j-1)+B%x1(i+1,j)-B%x1(i-1,j-1)-B%x1(i-1,j))*0.25d0
    Diy=(B%y1(i+1,j-1)+B%y1(i+1,j)-B%y1(i-1,j-1)-B%y1(i-1,j))*0.25d0
    Djx=B%x1(i,j)-B%x1(i,j-1)
    Djy=B%y1(i,j)-B%y1(i,j-1)
   
    Ds=1.d0/(Dix*Djy-Djx*Diy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Diu=(uu(i+1,j-1)+uu(i+1,j)-uu(i-1,j-1)-uu(i-1,j))*0.25d0
    Div=(v(i+1,j-1)+v(i+1,j)-v(i-1,j-1)-v(i-1,j))*0.25d0
    DiT=(T(i+1,j-1)+T(i+1,j)-T(i-1,j-1)-T(i-1,j))*0.25d0
    Dju=uu(i,j)-uu(i,j-1)
    Djv=v(i,j)-v(i,j-1)
    DjT=T(i,j)-T(i,j-1)
    
    
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
    
    B%TaoNSj_in(1,i,j)=t11
    B%TaoNSj_in(2,i,j)=t22
    B%TaoNSj_in(3,i,j)=t12
    B%TaoNSj_in(4,i,j)=t33
    
    B%Qj_in(1,i,j)=qx    !!-qx
    B%Qj_in(2,i,j)=qy    !!-qy
    

    
    if((j/=1).and.(j/=ny))then
             
              tmp= (B%S_GKUA(5,i,j)+B%S_GKUA(6,i,j)+B%S_GKUA(10,i,j)+B%S_GKUA(5,i,j-1)+B%S_GKUA(6,i,j-1)+B%S_GKUA(10,i,j-1))/6.0
              B%TaoNSj_G(1,i,j)= -((B%S_GKUA(5,i,j)+B%S_GKUA(5,i,j-1))/2.0-tmp)/(gamma*Ma**2)-t11
              B%TaoNSj_G(2,i,j)= -((B%S_GKUA(6,i,j)+B%S_GKUA(6,i,j-1))/2.0-tmp)/(gamma*Ma**2)-t22
              B%TaoNSj_G(3,i,j)= -(B%S_GKUA(7,i,j)+B%S_GKUA(7,i,j-1))/2.0/(gamma*Ma**2)-t12
              B%TaoNSj_G(4,i,j)= -((B%S_GKUA(10,i,j)+B%S_GKUA(10,i,j-1))/2.0-tmp)/(gamma*Ma**2)-t33
              
              B%Qj_G(1,i,j)=  -(B%S_GKUA(8,i,j-1)+B%S_GKUA(8,i,j))/2.0*sqrt(2.0)/(Ma**3*gamma**1.5)-qx
              B%Qj_G(2,i,j)=  -(B%S_GKUA(9,i,j-1)+B%S_GKUA(9,i,j))/2.0*sqrt(2.0)/(Ma**3*gamma**1.5)-qy
    else if(j==1)then
              tmp= (B%S_GKUA(5,i,j)+B%S_GKUA(6,i,j)+B%S_GKUA(10,i,j))/3.0
              B%TaoNSj_G(1,i,j)= -(B%S_GKUA(5,i,j)-tmp)/(gamma*Ma**2)-t11
              B%TaoNSj_G(2,i,j)= -(B%S_GKUA(6,i,j)-tmp)/(gamma*Ma**2)-t22
              B%TaoNSj_G(3,i,j)= -B%S_GKUA(7,i,j)/(gamma*Ma**2)-t12
              B%TaoNSj_G(4,i,j)= -(B%S_GKUA(10,i,j)-tmp)/(gamma*Ma**2)-t33
    
              B%Qj_G(1,i,j)=  -B%S_GKUA(8,i,j)*sqrt(2.0)/(Ma**3*gamma**1.5)-qx
              B%Qj_G(2,i,j)=  -B%S_GKUA(9,i,j)*sqrt(2.0)/(Ma**3*gamma**1.5)-qy
    else if(j==ny)then
             tmp= (B%S_GKUA(5,i,j-1)+B%S_GKUA(6,i,j-1)+B%S_GKUA(10,i,j-1))/3.0
              B%TaoNSj_G(1,i,j)= -(B%S_GKUA(5,i,j-1)-tmp)/(gamma*Ma**2)-t11
              B%TaoNSj_G(2,i,j)= -(B%S_GKUA(6,i,j-1)-tmp)/(gamma*Ma**2)-t22
              B%TaoNSj_G(3,i,j)= -B%S_GKUA(7,i,j-1)/(gamma*Ma**2)-t12
              B%TaoNSj_G(4,i,j)= -(B%S_GKUA(10,i,j-1)-tmp)/(gamma*Ma**2)-t33
    
              B%Qj_G(1,i,j)=  -B%S_GKUA(8,i,j-1)*sqrt(2.0)/(Ma**3*gamma**1.5)-qx
              B%Qj_G(2,i,j)=  -B%S_GKUA(9,i,j-1)*sqrt(2.0)/(Ma**3*gamma**1.5)-qy
    
    end if
    


   enddo
   enddo
   ! call filter_j(1,mBlock)  
    deallocate(d,uu,v,T,cc,p,Fluxi,Fluxj)
    end if   
   end do      
   

          
end Subroutine Convert_NS
    
    subroutine filter_i(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
 
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B
	 integer nMesh,mBlock,nx,ny,i,j,k,l,ibeg,ied,jbeg,jed,i_off,j_off,f_iter
     real(sp) :: tmp1(500,500),tmp2(500,500),tmp3(500,500),tmp4(500,500),tmp5(500,500)
     real(sp) :: W(-1:1,-1:1),Var(500,500)
     
     MP=> mesh(nMesh)
     B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny
  
    W(-1,-1)=1.0/16.0;W(-1,0)=2.0/16.0;W(-1,1)=1.0/16.0
    W(0,-1)=2.0/16.0;W(0,0)=4.0/16.0;W(0,1)=2.0/16.0
    W(1,-1)=1.0/16.0;W(1,0)=2.0/16.0;W(1,1)=1.0/16.0
  
    
    ibeg=1;ied=nx
    jbeg=1;jed=ny-1
    do f_iter=1,2
    tmp1=0.0;tmp2=0.0;tmp3=0.0;tmp4=0.0;tmp5=0.0
    do j=jbeg,jed
    do i=ibeg,ied
       
        i_off=0;j_off=0
            if(i==ibeg.and.j==jbeg)then
                 i_off=1;j_off=1
            else if(i==ibeg.and.j==jed)then
                i_off=1;j_off=-1
            else if(i==ied.and.j==jbeg)then
                i_off=-1;j_off=1
            else if(i==ied.and.j==jed)then
                i_off=-1;j_off=-1
            else if(i==ibeg.and.j/=jbeg.and.j/=jed)then
                i_off=1;j_off=0
            else if (i==ied.and.j/=jbeg.and.j/=jed)then
                i_off=-1;j_off=0
            else if (j==jbeg.and.i/=ibeg.and.i/=ied)then
                i_off=0;j_off=1
            else if(j==jed.and.i/=ibeg.and.i/=ied)then
                i_off=0;j_off=-1
            else 
                i_off=0;j_off=0
            endif
       
            do k=-1,1
            do l=-1,1 
                tmp1(i,j)=tmp1(i,j)+B%TaoNSi_G(1,i,j)*W(k,l)
                tmp2(i,j)=tmp2(i,j)+B%TaoNSi_G(2,i,j)*W(k,l)
                tmp3(i,j)=tmp3(i,j)+B%TaoNSi_G(3,i,j)*W(k,l)
                tmp4(i,j)=tmp4(i,j)+ B%Qi_G(1,i,j)*W(k,l)
                tmp5(i,j)=tmp5(i,j)+ B%Qi_G(2,i,j)*W(k,l)
            end do
            end do
            
            
    enddo
    enddo
    do j=jbeg,jed
    do i=ibeg,ied
                B%TaoNSi_G(1,i,j)=tmp1(i,j)
                B%TaoNSi_G(2,i,j)=tmp2(i,j)
                B%TaoNSi_G(3,i,j)= tmp3(i,j)
                B%Qi_G(1,i,j)=tmp4(i,j)
                B%Qi_G(2,i,j)=tmp5(i,j)
    end do
    end do
    end do
    
    return
    
    end subroutine filter_i
    
    
    subroutine filter_j(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
 
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B
	 integer nMesh,mBlock,nx,ny,i,j,k,l,ibeg,ied,jbeg,jed,i_off,j_off,f_iter
     real(sp) :: tmp1(500,500),tmp2(500,500),tmp3(500,500),tmp4(500,500),tmp5(500,500)
     real(sp) :: W(-1:1,-1:1),Var(500,500)
     
     MP=> mesh(nMesh)
     B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny
  
    W(-1,-1)=1.0/16.0;W(-1,0)=2.0/16.0;W(-1,1)=1.0/16.0
    W(0,-1)=2.0/16.0;W(0,0)=4.0/16.0;W(0,1)=2.0/16.0
    W(1,-1)=1.0/16.0;W(1,0)=2.0/16.0;W(1,1)=1.0/16.0
  
    
    ibeg=1;ied=nx-1
    jbeg=1;jed=ny
    do f_iter=1,2
    tmp1=0.0;tmp2=0.0;tmp3=0.0;tmp4=0.0;tmp5=0.0
    do j=jbeg,jed
    do i=ibeg,ied
       
        i_off=0;j_off=0
            if(i==ibeg.and.j==jbeg)then
                 i_off=1;j_off=1
            else if(i==ibeg.and.j==jed)then
                i_off=1;j_off=-1
            else if(i==ied.and.j==jbeg)then
                i_off=-1;j_off=1
            else if(i==ied.and.j==jed)then
                i_off=-1;j_off=-1
            else if(i==ibeg.and.j/=jbeg.and.j/=jed)then
                i_off=1;j_off=0
            else if (i==ied.and.j/=jbeg.and.j/=jed)then
                i_off=-1;j_off=0
            else if (j==jbeg.and.i/=ibeg.and.i/=ied)then
                i_off=0;j_off=1
            else if(j==jed.and.i/=ibeg.and.i/=ied)then
                i_off=0;j_off=-1
            else 
                i_off=0;j_off=0
            endif
       
            do k=-1,1
            do l=-1,1 
                tmp1(i,j)=tmp1(i,j)+B%TaoNSj_G(1,i,j)*W(k,l)
                tmp2(i,j)=tmp2(i,j)+B%TaoNSj_G(2,i,j)*W(k,l)
                tmp3(i,j)=tmp3(i,j)+B%TaoNSj_G(3,i,j)*W(k,l)
                tmp4(i,j)=tmp4(i,j)+ B%Qj_G(1,i,j)*W(k,l)
                tmp5(i,j)=tmp5(i,j)+ B%Qj_G(2,i,j)*W(k,l)
            end do
            end do
            
            
    enddo
    enddo
    do j=jbeg,jed
    do i=ibeg,ied
                B%TaoNSj_G(1,i,j)=tmp1(i,j)
                B%TaoNSj_G(2,i,j)=tmp2(i,j)
                B%TaoNSj_G(3,i,j)= tmp3(i,j)
                B%Qj_G(1,i,j)=tmp4(i,j)
                B%Qj_G(2,i,j)=tmp5(i,j)
    end do
    end do
    end do
    
    return
    
    end subroutine filter_j    
    subroutine getNSF
     use Global_Var
     use Flow_Var 
     use const_var ,only :sp,dop
     implicit none
 
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B
	 integer nMesh,mBlock,nx,ny,i,j,k,l,ibeg,ied,jbeg,jed,i_off,j_off,f_iter
     real(sp) :: tmp1(500,500),tmp2(500,500),tmp3(500,500),tmp4(500,500),tmp5(500,500)
     real(sp) :: W(-1:1,-1:1),Var(500,500)
     
     MP=> mesh(1)
     B => MP%Block(1)  
     B%TaoNS_out=0;B%Q_out=0
     nx=B%nx; ny=B%ny
     do i=1,nx-1
         do j=1,ny-1
             do k=1,4
             B%TaoNS_out(k,i,j)=0.25*(B%TaoNSi_out(k,i,j)+B%TaoNSj_out(k,i,j)+B%TaoNSi_out(k,i+1,j)+B%TaoNSj_out(k,i,j+1))&
                 & -0.25*(B%TaoNSi_in(k,i,j)+B%TaoNSj_in(k,i,j)+B%TaoNSi_in(k,i+1,j)+B%TaoNSj_in(k,i,j+1))
             B%TaoNS_out(k,i,j)=-1.0*B%TaoNS_out(k,i,j)*(gamma*Ma**2)
             end do
             do k=1,2
            B%Q_out(k,i,j)=0.25*(B%Qi_out(k,i,j)+B%Qj_out(k,i,j)+B%Qi_out(k,i+1,j)+B%Qj_out(k,i,j+1))&
                 & -0.25*(B%Qi_in(k,i,j)+B%Qj_in(k,i,j)+B%Qi_in(k,i+1,j)+B%Qj_in(k,i,j+1))
            B%Q_out(k,i,j)=-1.0*B%Q_out(k,i,j)*(Ma**3*gamma**1.5)/sqrt(2.0)
            end do
         end do
     end do
     
         
  
 
    
    return
    
    end subroutine getNSF 
