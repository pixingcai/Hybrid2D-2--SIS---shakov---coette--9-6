! modules for Boundary layer condition
! For OpenCFD-EC 2D ver 1.1
! Copyright by Li Xinliang, lixl@imech.ac.cn
! Code by Li Xinliang, 2011-3-27 
! -------modified-----------------------------------------------------
! 2011-11-23: Symmetry boundary condition is adding
! 2011-12-9:  Isothormal wall boundary condition is adding
!---------------------------------------------------------------------
! ����߽����������ڱ߽磩 ������һ������
     subroutine Boundary_condition_onemesh(nMesh)
     use Global_Var
     implicit none
     
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: nMesh,mBlock,ksub

 ! ------------------------------------------------------ 
     do mBlock=1,Mesh(nMesh)%Num_Block
        B => Mesh(nMesh)%Block(mBlock)
        do  ksub=1,B%subface
            Bc=> B%bc_msg(ksub)
           if(Bc%neighb .lt. 0 ) then               ! ���ڱ߽�
                  if( Bc%neighb .eq. BC_Wall ) then
                      !call boundary_wall(nMesh,mBlock,ksub)         ! �̱ڱ߽�
                  else if( Bc%neighb .eq. BC_Farfield ) then
                      call boundary_Farfield(nMesh,mBlock,ksub)     ! Զ���߽�
                  else if( Bc%neighb .eq. BC_Symmetry_or_slidewall ) then
                      call Symmetry_or_slidewall(nMesh,mBlock,ksub)   ! �Գƣ����ƣ��߽�
                  else if ( Bc%neighb .eq. BC_inlet )then
                      call boundary_inlet(nMesh,mBlock,ksub)
                  else if(Bc%neighb.eq.BC_outlet)then
                      call boundary_outlet(nMesh,mBlock,ksub)
                   endif
           endif
        enddo
     enddo
     
  end subroutine Boundary_condition_onemesh

!------------------------------------------------------------






!-------------------------------------------------------------------
!  �������߽�����
!  ʹ��LAP��������
    subroutine boundary_wall(nMesh,mBlock,ksub)
    use Global_Var
     
    implicit none
    Type (Block_TYPE),pointer:: B
    Type (BC_MSG_TYPE),pointer:: Bc
    Type (Mesh_TYPE),pointer:: MP
    integer::icheck,jcheck
    integer:: nMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i2,j2,i1,j1,k,n1,pos_x,pos_y
    real(sp) :: d1,p1,T1,d2,p2,T2,u1,v1,u2,v2,Tsb,Amu2,d_wall,temp(2,4),temp1(2),temp2(2),dnw,VK(2)
    real(sp) ::temp_nw(2),s_nw,dwall_tmp,uwall_tmp,vwall_tmp,twall_tmp
    real(sp) ,parameter:: beta1_SST=0.075d0
     
    MP => Mesh(nMesh)
    B  => MP%Block(mBlock)
    Bc => B%bc_msg(ksub)
     
    Bc%F_Wall=0.0
  
    ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
    n1=iend-ibegin+jend-jbegin
 
    do j=jbegin,jend
    do i=ibegin,iend
        icheck=max(ibegin,iend)
        jcheck=max(jbegin,jend)
        if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
            cycle
        end if
        if(Bc%face .eq. 1 ) then 
            pos_x=1;pos_y=j
        else if(Bc%face .eq. 2 ) then 
            pos_x=i;pos_y=1
        else if(Bc%face .eq. 3 ) then 
            pos_x=1;pos_y=j
        else if(Bc%face .eq. 4 ) then 
            pos_x=i;pos_y=1
        end if
           
        do k=1,LAP
    !      (i1,j1)  �ڵ㣻  (i2,j2) Ϊ��Ӧ��buffer���ĵ�   
            if(Bc%face .eq. 1) then
                i1=i+k-1+2; j1=j; i2=i-k+2 ; j2=j  
            else if(Bc%face .eq. 2) then
                i1=i; j1=j+k-1+2; i2=i; j2=j-k+2
            else if(Bc%face .eq. 3) then
                i1=i-k-2; j1=j;  i2=i+k-1-2; j2=j
            else
                i1=i; j1=j-k-2; i2=i; j2=j+k-1-2
            endif
            d1=B%U(1,i1,j1)    ! �ڵ㴦���ܶȡ�ѹ�����¶ȡ��ٶ�
            u1=B%U(2,i1,j1)/d1 
            v1=B%U(3,i1,j1)/d1
            p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)             
            T1=gamma*Ma*Ma*p1/d1     

            d2=-d1+2.0*dwall_tmp
            p2=d2*T2/(gamma*Ma*Ma)
            u2=-u1 +2.0* uwall_tmp            ! �޻��Ʊ�
            v2=-v1 +2.0* vwall_tmp
            
            B%U(1,i2,j2)=d2
            B%U(2,i2,j2)=d2*u2
            B%U(3,i2,j2)=d2*v2
            B%U(4,i2,j2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2)
        
       enddo
    enddo
    enddo 

    continue
       
           
    end



! ��߽�������ֻ�������/���ڣ������ֳ�/������ (�������������������������д���֤��
   subroutine boundary_Farfield(nMesh,mBlock,ksub)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP
     integer::icheck,jcheck
     integer:: nMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i2,j2,i1,j1,k,n1,pos_x,pos_y
     real(sp) :: d1,p1,T1,d2,p2,T2,u1,v1,u2,v2,Tsb,Amu2,d_wall,temp(2,4),temp1(2),temp2(2),dnw,VK(2)
      real(sp) ::temp_nw(2),s_nw,dwall_tmp,uwall_tmp,vwall_tmp,twall_tmp
     real(sp) ,parameter:: beta1_SST=0.075d0
     
     MP => Mesh(nMesh)
     B  => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)
     

     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
     n1=iend-ibegin+jend-jbegin

     do j=jbegin,jend
     do i=ibegin,iend
          icheck=max(ibegin,iend)
          jcheck=max(jbegin,jend)
          if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
          end if
          if(Bc%face .eq. 1 ) then 
             pos_x=1;pos_y=j
          else if(Bc%face .eq. 2 ) then 
             pos_x=i;pos_y=1
          else if(Bc%face .eq. 3 ) then
             pos_x=1;pos_y=j
          else if(Bc%face .eq. 4 ) then
             pos_x=i;pos_y=1
          end if
          do k=1,LAP
        !      (i1,j1)  �ڵ㣻  (i2,j2) Ϊ��Ӧ��buffer���ĵ�   
		     if(Bc%face .eq. 1) then
               i1=i+k-1; j1=j; i2=i-k ; j2=j  
             else if(Bc%face .eq. 2) then
               i1=i; j1=j+k-1; i2=i; j2=j-k
             else if(Bc%face .eq. 3) then
               i1=i-k; j1=j;  i2=i+k-1; j2=j
             else
               i1=i; j1=j-k; i2=i; j2=j+k-1
             endif

		    if(Twall .le. 0) then   ! ���ȱ�

 		       B%U(1,i2,j2)= B%U(1,i1,j1)       ! d(0)=d(1)   �Գ�
               B%U(2,i2,j2)=-B%U(2,i1,j1)       ! u(0)=-u(1)  -> d*u ���Գ�
               B%U(3,i2,j2)=-B%U(3,i1,j1)       ! v(0)=-v(1)  -> d*v ���Գ�
               B%U(4,i2,j2)= B%U(4,i1,j1)       ! E(0)=E(1)   �Գ�
        
		    else   ! ���±�
                if (B%solver==NS_Solver) then
			        d1=B%U(1,i1,j1)    ! �ڵ㴦���ܶȡ�ѹ�����¶ȡ��ٶ�
                    u1=B%U(2,i1,j1)/d1 
                    v1=B%U(3,i1,j1)/d1
		            p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)             
                    T1=gamma*Ma*Ma*p1/d1     
			 
                    T2=2.d0*(Twall+0.0279)-T1    ! ���±�  0.5*(T1+T2)=Twall
                    if(T2.LE.0.0)THEN
                       T2=Twall
                    end if
                     p2=-p1+2.0*0.51515              ! �߽����裬���洦����ѹ���ݶ�Ϊ0   
                     d2=gamma*Ma*Ma*p2/T2 
			         u2=-u1              ! �޻��Ʊ�
			         v2=-v1+2.0*0.127
			         B%U(1,i2,j2)=d2
                     B%U(2,i2,j2)=d2*u2
                     B%U(3,i2,j2)=d2*v2
  			         B%U(4,i2,j2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2)
                else if(B%solver==GKUA)then
                     dwall_tmp=BC%Swall(1,pos_x,pos_y)
                     uwall_tmp=BC%Swall(2,pos_x,pos_y)/Ma*sqrt(2.0/gamma)
                     vwall_tmp=BC%Swall(3,pos_x,pos_y)/Ma*sqrt(2.0/gamma)
                     Twall_tmp=BC%Swall(4,pos_x,pos_y)
                    
                	 d1=B%U(1,i1,j1)    ! �ڵ㴦���ܶȡ�ѹ�����¶ȡ��ٶ�
			         u1=B%U(2,i1,j1)/d1 
			         v1=B%U(3,i1,j1)/d1
		             p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)             
                     T1=gamma*Ma*Ma*p1/d1     
			 
                     T2=2.d0*(Twall_tmp)-T1    ! ���±�  0.5*(T1+T2)=Twall
                     if(T2.LE.0.0)THEN
                        T2=Twall_tmp
                     end if
                     d2=-d1+2.0*dwall_tmp
                     p2=d2*T2/(gamma*Ma*Ma)
                     u2=-u1 +2.0* uwall_tmp            ! �޻��Ʊ�
			         v2=-v1 +2.0* vwall_tmp
			     
                     B%U(1,i2,j2)=d2
                     B%U(2,i2,j2)=d2*u2
                     B%U(3,i2,j2)=d2*v2
  			         B%U(4,i2,j2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2)
              
                endif 
                
            end if 
            
          enddo
          
	   enddo
       enddo 

       continue
       
     
    end subroutine boundary_Farfield

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    !-------------------------------------------------------------------
!  ����������ר�ã������ϱ����ٶȱ߽�
    subroutine boundary_inlet(nMesh,mBlock,ksub)
     use Global_Var
     
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP
     integer::icheck,jcheck
     integer:: nMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i2,j2,i1,j1,k,n1,pos_x,pos_y,nx,ny
     real(sp) :: d1,p1,T1,d2,p2,T2,u1,v1,u2,v2,Tsb,Amu2,d_wall,temp(2,4),temp1(2),temp2(2),dnw,VK(2)
      real(sp) ::temp_nw(2),s_nw
     real(sp) ,parameter:: beta1_SST=0.075d0
     
     MP => Mesh(nMesh)
     B  => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)
     nx=B%nx
     ny=B%ny
     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      n1=iend-ibegin+jend-jbegin

       do j=jbegin,jend
       do i=ibegin,iend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
           
             do k=1,LAP
            !      (i1,j1)  �ڵ㣻  (i2,j2) Ϊ��Ӧ��buffer���ĵ�   
		         if(Bc%face .eq. 1) then
                   i1=nx-k; j1=j; i2=i-k ; j2=j  
                 else if(Bc%face .eq. 2) then
                   i1=i; j1=ny-k; i2=i; j2=j-k
                 else if(Bc%face .eq. 3) then
                   i1=k; j1=j;  i2=i+k-1; j2=j
                 else
                   i1=i; j1=k; i2=i; j2=j+k-1
                 endif
                 B%U(1,i2,j2)=B%U(1,i1,j1)
                 B%U(2,i2,j2)=B%U(2,i1,j1)
                 B%U(3,i2,j2)=B%U(3,i1,j1)
  			     B%U(4,i2,j2)=B%U(4,i1,j1)
                 
             enddo
             
	   enddo
       enddo 
  
    end subroutine boundary_inlet
    


! �ڱ߽�������������
   subroutine boundary_outlet(nMesh,mBlock,ksub)
     use Global_Var
     
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP

     integer:: NMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i1,j1,i2,j2,k
     real(sp) :: d_inf,u_inf,v_inf,p_inf,p_out,n1,n2 ,temp_T,temp_coef,temp_P,temp_Rho 
     real(sp) :: d1,u1,v1,p1,c1,d2,u2,v2,p2,pb,db,ub,vb,un1,c_inf
     real(sp) :: dx,dy,s,Ma_n,Riemann_inf,Riemann_bound,vn_bound,c_bound,T_bound,ut1,u_bound,v_bound,p_bound
     real(sp) :: entropy_bound,d_bound
     !!����������Զ���ڱ߽�����
     d_inf=P_outlet*1.0d0; u_inf=0.0; v_inf=0.0 ; p_inf=P_outlet*1.0d0/(gamma*Ma*Ma)

     MP=>Mesh(nMesh)
     B => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)

     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      
       do j=jbegin,jend
       do i=ibegin,iend
           !   ����߽� �� ������       
          if( Bc%face .eq. 1 ) then 
             dx=B%x(i,j+1)-B%x(i,j) ;    dy=B%y(i,j+1)-B%y(i,j);  s=sqrt(dx*dx+dy*dy)
             n1=-dy/s; n2= dx/s    ! �ⷨ��  
          else if(Bc%face .eq. 3) then 
             dx=B%x(i,j+1)-B%x(i,j) ;    dy=B%y(i,j+1)-B%y(i,j);  s=sqrt(dx*dx+dy*dy)
             n1=dy/s; n2= -dx/s    ! �ⷨ��  
          else if(Bc%face .eq. 2) then
             dx=B%x(i+1,j)-B%x(i,j) ;   dy=B%y(i+1,j)-B%y(i,j) ;  s=sqrt(dx*dx+dy*dy)   
             n1=dy/s; n2=-dx/s      ! �ⷨ�� 
          else
             dx=B%x(i+1,j)-B%x(i,j) ;   dy=B%y(i+1,j)-B%y(i,j) ;  s=sqrt(dx*dx+dy*dy)   
             n1=-dy/s; n2=dx/s        ! �ⷨ�� 
          endif
!-----------------------------------------------------------------
 !  (i1,j1) �ǿ����߽���ڵ㣬 (i2,j2) �Ǳ߽����Ghost Cell��    
         if(Bc%face .eq. 1) then       ! i- �棬 
           i1=i; j1=j  
         else if(Bc%face .eq. 2) then  ! j- ��
           i1=i; j1=j
         else if(Bc%face .eq. 3) then  ! i+ �� (i=ibegin=iend=nx), i1=i-1 ���ڵ�, i2=i=nx��Ghost Cell
           i1=i-1; j1=j
         else
           i1=i; j1=j-1
         endif
            d1=B%U(1,i1,j1) ;   u1=B%U(2,i1,j1)/d1 ;     v1=B%U(3,i1,j1)/d1
            p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)              ! �ڵ㴦��ֵ
            c1=sqrt(gamma*p1/d1) 
            un1=(u1*n1+v1*n2)

            
            ut1=(-1.0*u1*n2+v1*n1) 
            Riemann_Bound=un1+2.0*c1/(gamma-1.0)
            
            c_inf=sqrt(gamma*p_inf/d_inf) !!!����¶ȵ��ڳ����¶ȵ��ڲο��¶�,�������¶ȱ߽���������������
            Riemann_inf=-2.0*c_inf/(gamma-1.0)!!����Զ�ٶ�0
            
            vn_bound=(Riemann_Bound+Riemann_inf)/2.0
            c_bound=(Riemann_bound-Riemann_inf)*(gamma-1.0)/4.0
         !!!���Ի���
                entropy_bound=p1/(d1**gamma)
                d_bound=(c_bound**2/(gamma*entropy_bound))**(1.0/(gamma-1))
                p_bound=entropy_bound*d_bound**gamma
                u_bound=vn_bound*n1+ut1*(-n2)
                v_bound=vn_bound*n2+ut1*n1
    

!-------------------------------------------------------------------------------------------
          do k=1,LAP
	       if(Bc%face .eq. 1) then       ! i- �棬 
             i2=i-k ; j2=j  
           else if(Bc%face .eq. 2) then  ! j- ��
             i2=i; j2=j-k
           else if(Bc%face .eq. 3) then  ! i+ �� (i=ibegin=iend=nx), i1=i-1 ���ڵ�, i2=i=nx��Ghost Cell
             i2=i+k-1; j2=j
           else
            i2=i; j2=j+k-1
           endif
	   
	     
              !B%U(1,i2,j2)=d_bound
              !B%U(2,i2,j2)=d_bound*u_bound
              !B%U(3,i2,j2)=d_bound*v_bound
              !B%U(4,i2,j2)=p_bound/(gamma-1.d0)+0.5d0*d_bound*(u_bound*u_bound+v_bound*v_bound)
           
              B%U(1,i2,j2)=B%U(1,i1,j1)
              B%U(2,i2,j2)=B%U(2,i1,j1)
              B%U(3,i2,j2)=B%U(3,i1,j1)
              B%U(4,i2,j2)= B%U(4,i1,j1)


          
        enddo
	   enddo
      enddo 
    end subroutine boundary_outlet

!-------------------------------------------------------------------------


!-------------------------------------------------------------------
!  �ԳƱ߽�����(���ƹ̱�)
!  ʹ��LAP��������
    subroutine Symmetry_or_slidewall(nMesh,mBlock,ksub)
     use Global_Var
     
     implicit none
     
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP

     integer:: nMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i2,j2,i1,j1,k
     real(sp) :: p2,dx,dy,si,n1,n2,Vn
     
     MP=>Mesh(nMesh)
     B => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      

       do j=jbegin,jend
       do i=ibegin,iend
!      ����߽編����    
	     if(Bc%face .eq. 1 .or. Bc%face .eq. 3) then    ! i+ or i- boundary
	      dx=B%x(i,j+1)-B%x(i,j) ;      dy=B%y(i,j+1)-B%y(i,j)
	      si=sqrt(dx*dx+dy*dy)
          n1=dy/si; n2=-dx/si   ! normal vector at (i,j) or (I-1/2,J) 
      
		 else   ! j+ or j-
	      dx=B%x(i+1,j)-B%x(i,j) ;      dy=B%y(i+1,j)-B%y(i,j)
  	      si=sqrt(dx*dx+dy*dy)
          n1=-dy/si; n2=dx/si   ! normal vector at i, j+1/2
 		 endif


       do k=1,LAP
         if(Bc%face .eq. 1) then
           i1=i+k-1; j1=j; i2=i-k ; j2=j  
         else if(Bc%face .eq. 2) then
           i1=i; j1=j+k-1; i2=i; j2=j-k
         else if(Bc%face .eq. 3) then
           i1=i-k; j1=j;  i2=i+k-1; j2=j
         else
           i1=i; j1=j-k; i2=i; j2=j+k-1
         endif
           Vn=B%U(2,i1,j1)*n1+B%U(3,i1,j1)*n2   ! ������

         B%U(1,i2,j2)= B%U(1,i1,j1)       ! d(0)=d(1)   �Գ�
         B%U(2,i2,j2)= B%U(2,i1,j1)-2.d0*Vn*n1       ! �������෴������������
         B%U(3,i2,j2)= B%U(3,i1,j1)-2.d0*Vn*n2       ! �������෴������������
         B%U(4,i2,j2)= B%U(4,i1,j1)       ! E(0)=E(1)   �Գ�
	   
	    if(MP%Nvar .eq. 5) then
	      B%U(5,i2,j2)=B%U(5,i1,j1)     ! �������Գ�
	    elseif(MP%Nvar .eq. 6) then         ! k,w
	      B%U(5,i2,j2)=B%U(5,i1,j1)     ! �������Գ�
	      B%U(6,i2,j2)=B%U(6,i1,j1)
	    endif
	    
	   enddo
	   enddo
       enddo 
      
    end


!---------------------------------------------------------------------------
! inner boundary condition  �ڱ߽�
! �����������ӹ�ϵ�����»�����������������Ϣ (����һ������)
! ������ΪLAP������ (Ŀǰ�汾�趨LAP=2)

     subroutine update_buffer_onemesh(nMesh)
     use Global_Var
     
     implicit none
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc,Bc1
     real(sp) ,allocatable:: Ua(:,:,:)
     integer::flag_Hybrid
     integer:: nMesh,Nvar1,i,j,k,k1,m,mBlock,ksub,n1,m_neighbour,msub,orient,Kflag_initial,nx1,ny1
     integer:: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2

!----------------------------------------------------------
 MP=>Mesh(nMesh)   ! ���� ��nMesh=3,2,1���� �֡��С������� 
 Nvar1=MP%Nvar     ! ������Ŀ ��4��6������������ֻ��4�������� �����������4��6������������SSTģ�͵�����������
 do mBlock=1,MP%Num_Block
  B => MP%Block(mBlock)
  
  if( (B%solver==NS_Solver))then  !1-25
      do  ksub=1,B%subface
          Bc=> B%bc_msg(ksub)
    if( Bc%neighb .ge. 0 ) then   ! inner boundary
         ibegin2=Bc%ist; iend2=Bc%iend; jbegin2=Bc%jst; jend2=Bc%jend       ! write
         n1=iend2-ibegin2+jend2-jbegin2     ! +1 (number of cells)
         !  �����ӵ����Ϣ������ʱ����Ua
         allocate(Ua(Nvar1,n1,LAP))             !�洢������Ϣ����ʱ����

           m_neighbour=Bc%neighb;  msub= Bc%subface   
           B1 =>  Mesh(nMesh)%Block(m_neighbour)
           Bc1 => B1%bc_msg(msub)
       if( (B1%solver==NS_Solver))then   !1-25
           ibegin1=Bc1%ist; iend1=Bc1%iend; jbegin1=Bc1%jst; jend1=Bc1%jend   ! read
 
           if(Bc1%face .eq. 1 ) then                     !  boundary  i-
             do j=jbegin1,jend1-1
		     do i=1,LAP   
		       Ua(:,j-jbegin1+1,i)=B1%U(:,ibegin1+i-1,j) 
             enddo
	         enddo

	       else if ( Bc1%face .eq. 2) then  !  boundary  j-
             do  j=1,LAP
		     do  i=ibegin1,iend1-1  
		       Ua(:,i-ibegin1+1,j)=B1%U(:,i,jbegin1+j-1)
             enddo
		     enddo

	       else if( Bc1%face .eq. 3) then    !  boundary  i+
             do j=jbegin1, jend1-1
		     do i=1,LAP  
               Ua(:,j-jbegin1+1,i)=B1%U(:,iend1-i,j)  
             enddo
		     enddo

	       else                                !  boundary  j+
             do j=1,LAP
	         do i=ibegin1,iend1-1
               Ua(:,i-ibegin1+1,j)=B1%U(:,i,jend1-j)
             enddo
		     enddo

	       endif

!                                    ����ʱ����Ua�е���Ϣд�뻺����
       orient=Bc%orient              ! orient==2 ˳��;   ����ֵΪ ����
       do k=1,n1
        if(orient .eq. 2) then
          k1=k
        else
          k1=n1+1-k
        endif
       if(Bc%face .eq. 1 ) then             !  boundary  i-
         do i=1,LAP  
	            B%U(:,ibegin2-i,jbegin2+k-1)=Ua(:,k1,i)
         enddo

	   else if(Bc%face .eq. 2 ) then        !  boundary  j-
         do j=1,LAP   
	          B%U(:,ibegin2+k-1,jbegin2-j)=Ua(:,k1,j) 
         enddo
	   else if (Bc%face .eq. 3 ) then       !  boundary  i+
         do i=1,LAP
	       B%U(:,iend2+i-1,jbegin2+k-1)=Ua(:,k1,i)
         enddo
	   else                                          !  boundary  j+
         do j=1,LAP
	       B%U(:,ibegin2+k-1,jend2+j-1)=Ua(:,k1,j)
         enddo
      endif

      enddo
      
        end if  !1-25
        deallocate(Ua)
     endif
    enddo
 
  !  ���������Ľǵ㣬�ò�ֵ�ķ�����ֵ     
        nx1=B%nx; ny1=B%ny
        call U_average_conner(Nvar1,B%U(:,1,0),B%U(:,1,1),B%U(:,0,1),B%U(:,0,0),Cv)
        call U_average_conner(Nvar1,B%U(:,1,ny1),B%U(:,1,ny1-1),B%U(:,0,ny1-1),B%U(:,0,ny1),Cv)
        call U_average_conner(Nvar1,B%U(:,nx1-1,0),B%U(:,nx1-1,1),B%U(:,nx1,1),B%U(:,nx1,0),Cv)
        call U_average_conner(Nvar1,B%U(:,nx1-1,ny1),B%U(:,nx1-1,ny1-1),B%U(:,nx1,ny1-1),B%U(:,nx1,ny1),Cv)

   end if        !1-25
   
    enddo


   end  subroutine update_buffer_onemesh




!-------------------------------
! ���㻺�����ǵ㴦��ֵ ������U1(0,0)���ֵ���� ����巽������
  subroutine U_average_conner(Nvar1,U1,U2,U3,U4,Cv)
  use const_var ,only :sp,dop
  implicit none
  
  integer::Nvar1
  real(sp) ,dimension(Nvar1):: U1,U2,U3,U4
  real(sp) :: Cv,d1,uu1,v1,T1,d2,uu2,v2,T2,d3,uu3,v3,T3,d4,uu4,v4,T4
    d1=U1(1); uu1=U1(2)/d1; v1=U1(3)/d1; T1=(U1(4)-(uu1*U1(2)+v1*U1(3))*0.5d0)/(d1*Cv)  ! density, velocity, Temperature 
    d2=U2(1); uu2=U2(2)/d2; v2=U2(3)/d2; T2=(U2(4)-(uu2*U2(2)+v2*U2(3))*0.5d0)/(d2*Cv)   
    d3=U3(1); uu3=U3(2)/d3; v3=U3(3)/d3; T3=(U3(4)-(uu3*U3(2)+v3*U3(3))*0.5d0)/(d3*Cv) 
    d4=d1+d3-d2; uu4=uu1+uu3-uu2; v4=v1+v3-v2; T4=T1+T3-T2
    U4(1)=d4; U4(2)=d4*uu4; U4(3)=d4*v4; U4(4)=d4*(Cv*T4+(uu4*uu4+v4*v4)*0.5d0)
    if(Nvar1 .eq. 6) then
      U4(5)=U1(5)+U3(5)-U2(5)
      U4(6)=U1(6)+U3(6)-U2(6)
    endif
  end subroutine U_average_conner
