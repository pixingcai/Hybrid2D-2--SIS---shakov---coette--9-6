    
    Subroutine GKUA_Phy_boundary
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      integer:: iv,jv,iv1,jv1,i,j
      integer:: nx,ny
      integer:: NMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i1,j1,i2,j2,k
      real(sp) ::FRDtmp,TIJtmp,viu,vjv,vij2,GDVijv
      real(sp) :: d1,u1,v1,p1,c1,d2,u2,v2,p2,pb,db,ub,vb
      real(sp) :: dx,dy,s,Ma_n
       real(sp) ::n1,n2
      
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
     
       MP=>Mesh(1)   
       do mBlock=1,MP%Num_Block
           B => MP%Block(mBlock)
          if(B%solver==GKUA)then
           do  ksub=1,B%subface
               Bc=> B%bc_msg(ksub)
               if(Bc%neighb.lt.0)then  !!外边界
                   if(Bc%neighb .eq. BC_Wall)then
                      call GKUA_boundary_wall(1,mBlock,ksub)  
                   else if( Bc%neighb .eq. BC_Farfield ) then
                      call GKUA_boundary_farfield(1,mBlock,ksub)
                   else if( Bc%neighb .eq. BC_inlet)then
                      call GKUA_boundary_inlet(1,mBlock,ksub)
                   else if (Bc%neighb .eq. BC_outlet) then
                      call GKUA_boundary_outlet(1,mBlock,ksub)
                      
                   end if
               end if
               
               
             
           end do
           end if
       end do
       

    
    end subroutine GKUA_Phy_boundary
    
     
    Subroutine GKUA_boundary_wall(nMesh,mBlock,ksub)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Integer :: nMesh,mBlock,ksub
      integer:: iv,jv,iv1,jv1,i,j,nx,ny,ird,ibegin ,iend,jbegin,jend
      integer:: pos_x,pos_y,pos_x1,pos_y1,ijwall(2),icheck,jcheck
      real(sp) ::FRDtmp,wall_u,wall_v
    
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(ksub)
       wall_u=0.0
       wall_v=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
       do i=ibegin,iend    
       do j=jbegin,jend
                 icheck=max(ibegin,iend)
                 jcheck=max(jbegin,jend)
                 if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                     cycle
                 end if
                if(Bc%face .eq. 1 ) then 
                    pos_x=1;pos_y=j;pos_x1=2;pos_y1=j;ijwall(1)=-3;ijwall(2)=j
                else if(Bc%face .eq. 2 ) then 
                    pos_x=i;pos_y=1;pos_x1=i;pos_y1=2;ijwall(1)=i;ijwall(2)=-3
                else if(Bc%face .eq. 3 ) then 
                    pos_x=nx-1;pos_y=j;pos_x1=nx-2;pos_y1=j;ijwall(1)=nx+3;ijwall(2)=j
                else if(Bc%face .eq. 4 ) then 
                    pos_x=i;pos_y=ny-1;pos_x1=i;pos_y1=ny-2;ijwall(1)=i;ijwall(2)=ny+3
                end if
         !   write(*,*)ijwall(1),ijwall(2),mBlock,ksub
                  DO jv=GKUA_JST,GKUA_JEND
                  DO iv=GKUA_IST,GKUA_IEND
                     do ird=1,NRD
                        FRDtmp=(-B%FRD(ird,pos_x1,pos_y1,iv,jv)+ 3.*B%FRD(ird,pos_x,pos_y,iv,jv))*0.5
                        if(FRDtmp.LE.0.)FRDtmp=B%FRD(ird,pos_x,pos_y,iv,jv)
                        if(FRDtmp.GE.1.5*B%FRD(ird,pos_x,pos_y,iv,jv))FRDtmp=1.5*B%FRD(ird,pos_x,pos_y,iv,jv)
                        B%FRD(ird,ijwall(1),ijwall(2),iv,jv)=FRDtmp
                     enddo
                  ENDDO
                  ENDDO
       enddo
       enddo
       call GKUA_wall_boundary (mBlock,ksub,wall_u,wall_v)
            
    end Subroutine GKUA_boundary_wall   
    
    !-------------------------------------------------------------------------
    Subroutine GKUA_boundary_farfield(nMesh,mBlock,ksub)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Integer :: nMesh,mBlock,ksub
      integer:: iv,jv,iv1,jv1,i,j,nx,ny,ird,ibegin ,iend,jbegin,jend
      integer:: pos_x,pos_y,pos_x1,pos_y1,ijwall(2),icheck,jcheck
      real(sp) ::FRDtmp,wall_u,wall_v
    
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
      
       MP=>Mesh(1)   
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(ksub)
       wall_u=0.14828
       wall_v=0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
       do i=ibegin,iend    
       do j=jbegin,jend
                 icheck=max(ibegin,iend)
                 jcheck=max(jbegin,jend)
                 if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                     cycle
                 end if
                if(Bc%face .eq. 1 ) then 
                    pos_x=1;pos_y=j;pos_x1=2;pos_y1=j;ijwall(1)=-3;ijwall(2)=j
                else if(Bc%face .eq. 2 ) then 
                    pos_x=i;pos_y=1;pos_x1=i;pos_y1=2;ijwall(1)=i;ijwall(2)=-3
                else if(Bc%face .eq. 3 ) then 
                    pos_x=nx-1;pos_y=j;pos_x1=nx-2;pos_y1=j;ijwall(1)=nx+3;ijwall(2)=j
                else if(Bc%face .eq. 4 ) then 
                    pos_x=i;pos_y=ny-1;pos_x1=i;pos_y1=ny-2;ijwall(1)=i;ijwall(2)=ny+3
                end if
         !   write(*,*)ijwall(1),ijwall(2),mBlock,ksub
                  DO jv=GKUA_JST,GKUA_JEND
                  DO iv=GKUA_IST,GKUA_IEND
                     do ird=1,NRD
                        FRDtmp=(-B%FRD(ird,pos_x1,pos_y1,iv,jv)+ 3.*B%FRD(ird,pos_x,pos_y,iv,jv))*0.5
                        if(FRDtmp.LE.0.)FRDtmp=B%FRD(ird,pos_x,pos_y,iv,jv)
                        if(FRDtmp.GE.1.5*B%FRD(ird,pos_x,pos_y,iv,jv))FRDtmp=1.5*B%FRD(ird,pos_x,pos_y,iv,jv)
                        B%FRD(ird,ijwall(1),ijwall(2),iv,jv)=FRDtmp
                     enddo
                  ENDDO
                  ENDDO
       enddo
       enddo
       call GKUA_wall_boundary (mBlock,ksub,wall_u,wall_v)
      ! call DETFPwall(mBlock,ksub)
             
                     
        end subroutine GKUA_boundary_farfield
!-------------------------------------------------------------------------


   subroutine GKUA_boundary_inlet(nMesh,mBlock,ksub)
      use Global_Var
      use Com_ctrl
      use const_var ,only :sp,dop
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP
     integer:: iv,jv,iv1,jv1,nx,ny
     integer:: NMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i1,j1,i2,j2,k
     real(sp) :: d_inf,u_inf,v_inf,p_tmp,p_out,n1,n2,temp_T,temp_coef,temp_P,temp_Rho 
     real(sp) :: tp_ron,tp_u,tp_v
     real(sp) :: d1,u1,v1,p1,c1,d2,u2,v2,p2,pb,db,ub,vb
     real(sp) :: dx,dy,s,Ma_n
     real(sp) ::FRDtmp,TIJtmp,viu,vjv,vij2,GDVijv
     real(sp) :: un1,ut1,Riemann_bound,c_inf,Riemann_inf,vn_bound,c_bound,T_bound,U_bound,V_bound,P_bound
     real(sp)  :: entropy_inf,d_Bound
     MP=>Mesh(1)
     B => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)
     nx=B%nx
     ny=B%ny
     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      
       do j=jbegin,jend
       do i=ibegin,iend              
!-------------------------------------------------------------------------------------------
          do k=1,LAP
               !      (i1,j1)  内点；  (i2,j2) 为对应的buffer区的点   
		     if(Bc%face .eq. 1) then
               i1=nx-k; j1=j; i2=i-k ; j2=j  
             else if(Bc%face .eq. 2) then
               i1=i; j1=ny-k; i2=i; j2=j-k
             else if(Bc%face .eq. 3) then
               i1=k; j1=j;  i2=i+k-1; j2=j
             else
               i1=i; j1=k; i2=i; j2=j+k-1
             endif
	         B%S_GKUA(:,i2,j2)=B%S_GKUA(:,i1,j1)
             DO jv=GKUA_JST,GKUA_JEND          
                jv1=jv+GKUA_JPT
             DO iv=GKUA_IST,GKUA_IEND          
                iv1=iv+GKUA_IPT
                         B%FRD(1,i2,j2,iv,jv)= B%FRD(1,i1,j1,iv,jv)
                         B%FRD(2,i2,j2,iv,jv)= B%FRD(2,i1,j1,iv,jv)
             ENDDO
             ENDDO
             
        
           
        enddo
	   enddo
      enddo 

    end subroutine GKUA_boundary_inlet

!-------------------------------------------------------------------------


   subroutine GKUA_boundary_outlet(nMesh,mBlock,ksub)
      use Global_Var
      use Com_ctrl
      use const_var ,only :sp,dop
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP
     integer:: iv,jv,iv1,jv1
     integer:: NMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i1,j1,i2,j2,k
     real(sp) :: d_inf,u_inf,v_inf,p_tmp,p_out,n1,n2,temp_T,temp_coef,temp_P,temp_Rho 
     real(sp) :: tp_ron,tp_u,tp_v
     real(sp) :: d1,u1,v1,p1,c1,d2,u2,v2,p2,pb,db,ub,vb
     real(sp) :: dx,dy,s,Ma_n
     real(sp) ::FRDtmp,TIJtmp,viu,vjv,vij2,GDVijv
     
     real(sp) :: Riemann_inf,Riemann_bound,vn_bound,c_bound,T_bound,ut1,u_bound,v_bound,p_bound
     real(sp) :: entropy_bound,d_bound,un1,c_inf
!  本软件目前用来计算外流，给定无穷远条件
     d_inf=P_outlet*1.0d0; u_inf=0.0; v_inf=0.0 ; p_inf=P_outlet*1.0d0/(gamma*Ma*Ma)

     MP=>Mesh(nMesh)
     B => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)

     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      
       do j=jbegin,jend
       do i=ibegin,iend
           !   计算边界 外 法方向       
          if( Bc%face .eq. 1 ) then 
             dx=B%x(i,j+1)-B%x(i,j) ;    dy=B%y(i,j+1)-B%y(i,j);  s=sqrt(dx*dx+dy*dy)
             n1=-dy/s; n2= dx/s    ! 外法线  
          else if(Bc%face .eq. 3) then 
             dx=B%x(i,j+1)-B%x(i,j) ;    dy=B%y(i,j+1)-B%y(i,j);  s=sqrt(dx*dx+dy*dy)
             n1=dy/s; n2= -dx/s    ! 外法线  
          else if(Bc%face .eq. 2) then
             dx=B%x(i+1,j)-B%x(i,j) ;   dy=B%y(i+1,j)-B%y(i,j) ;  s=sqrt(dx*dx+dy*dy)   
             n1=dy/s; n2=-dx/s      ! 外法线 
          else
             dx=B%x(i+1,j)-B%x(i,j) ;   dy=B%y(i+1,j)-B%y(i,j) ;  s=sqrt(dx*dx+dy*dy)   
             n1=-dy/s; n2=dx/s        ! 外法线 
          endif
!-----------------------------------------------------------------
 !  (i1,j1) 是靠近边界的内点， (i2,j2) 是边界外的Ghost Cell点    
         if(Bc%face .eq. 1) then       ! i- 面， 
           i1=i; j1=j  
         else if(Bc%face .eq. 2) then  ! j- 面
           i1=i; j1=j
         else if(Bc%face .eq. 3) then  ! i+ 面 (i=ibegin=iend=nx), i1=i-1 是内点, i2=i=nx是Ghost Cell
           i1=i-1; j1=j
         else
           i1=i; j1=j-1
         endif
            d1=B%U(1,i1,j1) ;   u1=B%U(2,i1,j1)/d1 ;     v1=B%U(3,i1,j1)/d1
            p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)              ! 内点处的值
            c1=sqrt(gamma*p1/d1) 
            un1=(u1*n1+v1*n2)
            ut1=(-1.0*u1*n2+v1*n1)
            Riemann_Bound=un1+2.0*c1/(gamma-1.0)
            
            c_inf=sqrt(gamma*p_inf/d_inf) !!!入口温度等于出口温度等于参考温度,即出口温度边界条件等于来流。
            Riemann_inf=-2.0*c_inf/(gamma-1.0)!!无穷远速度0
            
            vn_bound=(Riemann_Bound+Riemann_inf)/2.0
            c_bound=(Riemann_bound-Riemann_inf)*(gamma-1.0)/4.0
            !!不考虑回流对速度边界的影响
                entropy_bound=p1/(d1**gamma)
                d_bound=(c_bound**2/(gamma*entropy_bound))**(1.0/(gamma-1))
                p_bound=entropy_bound*d_bound**gamma
                u_bound=vn_bound*n1+ut1*(-n2)
                v_bound=vn_bound*n2+ut1*n1
    

!-------------------------------------------------------------------------------------------
          do k=1,LAP
	       if(Bc%face .eq. 1) then       ! i- 面， 
             i2=i-k ; j2=j  
           else if(Bc%face .eq. 2) then  ! j- 面
             i2=i; j2=j-k
           else if(Bc%face .eq. 3) then  ! i+ 面 (i=ibegin=iend=nx), i1=i-1 是内点, i2=i=nx是Ghost Cell
             i2=i+k-1; j2=j
           else
            i2=i; j2=j+k-1
          endif
	   
	     
              B%U(1,i2,j2)=d_bound
              B%U(2,i2,j2)=d_bound*u_bound
              B%U(3,i2,j2)=d_bound*v_bound
              B%U(4,i2,j2)=p_bound/(gamma-1.d0)+0.5d0*d_bound*(u_bound*u_bound+v_bound*v_bound)
              tp_ron=B%U(1,i2,j2)
              tp_u=B%U(2,i2,j2)/B%U(1,i2,j2)
              tp_v=B%U(3,i2,j2)/B%U(1,i2,j2)
              
               B%S_GKUA(1,i2,j2)= B%U(1,i2,j2)
               B%S_GKUA(2,i2,j2)= (B%U(2,i2,j2)/B%U(1,i2,j2))*Ma/sqrt(2.0/gamma)
               B%S_GKUA(3,i2,j2)=(B%U(3,i2,j2)/B%U(1,i2,j2))*Ma/sqrt(2.0/gamma)
               B%S_GKUA(4,i2,j2)=(B%U(4,i2,j2)-0.5*B%U(1,i2,j2)*(tp_u*tp_u+tp_v*tp_v))/(Cv*tp_ron)
               B%S_GKUA(5,i2,j2)=B%S_GKUA(4,i2,j2)
                 DO jv=GKUA_JST,GKUA_JEND          
                   jv1=jv+GKUA_JPT
                   vjv=B%vj(jv1)-B%S_GKUA(3,i2,j2)
                 DO iv=GKUA_IST,GKUA_IEND          
                        iv1=iv+GKUA_IPT
                        viu=B%vi(iv1)-B%S_GKUA(2,i2,j2)
                        vij2=(viu**2+vjv**2)/B%S_GKUA(4,i2,j2)
                        GDVijv=B%S_GKUA(1,i,j)/(PI*B%S_GKUA(4,i2,j2)) *EXP(-vij2)
                         B%FRD(1,i2,j2,iv,jv)=GDVijv
                         B%FRD(2,i2,j2,iv,jv)=B%S_GKUA( 4,i2,j2)*GDVijv/2.
                         B%FRD(3,i2,j2,iv,jv)=Rdof*B%S_GKUA( 5,i2,j2)*GDVijv/2.
                 ENDDO
                 ENDDO
          
        enddo
	   enddo
      enddo 

    end subroutine GKUA_boundary_outlet
!

    
  Subroutine GKUA_inner_mesco_boundary(nMesh)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Integer :: nMesh,mBlock
      integer::iv,jv,iv1,jv1,i,j,nx,ny,ird,ksub,k,k1
      integer :: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2
      integer :: N1,m_neighbour,msub,orient
      real(sp) ,allocatable:: Ua(:,:,:)
      Type (Block_TYPE),pointer:: B,B1
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc,BC1
       MP=>Mesh(1)   
      DO jv=GKUA_JST,GKUA_JEND
         jv1=jv+GKUA_JPT
      DO iv=GKUA_IST,GKUA_IEND
         iv1=iv+GKUA_IPT
       
         do mBlock=1,MP%Num_Block
            B => MP%Block(mBlock)
            nx=B%nx; ny=B%ny
            if(B%solver==GKUA)then
            do  ksub=1,B%subface
                Bc=> B%bc_msg(ksub)
             if( Bc%neighb .ge. 0 ) then   ! inner boundary
                  ibegin2=Bc%ist; iend2=Bc%iend; jbegin2=Bc%jst; jend2=Bc%jend       ! write
                  n1=iend2-ibegin2+jend2-jbegin2     ! +1 (number of cells)
    !              把连接点的信息读入临时数组Ua
                  allocate(Ua(3,n1,LAP))             !存储交换信息的临时数组
                  m_neighbour=Bc%neighb;  msub= Bc%subface   
                  B1 =>  Mesh(nMesh)%Block(m_neighbour)
                  if(B1%solver==GKUA)then
                  Bc1 => B1%bc_msg(msub)
                 ibegin1=Bc1%ist; iend1=Bc1%iend; jbegin1=Bc1%jst; jend1=Bc1%jend   ! read

                
                   if(Bc1%face .eq. 1 ) then                     !  boundary  i-
                     do j=jbegin1,jend1-1
		             do i=1,LAP   
                          Ua(:,j-jbegin1+1,i)=B1%FRD(:,ibegin1+i-1,j,iv,jv) 
                          
                     enddo
	                 enddo

	               else if ( Bc1%face .eq. 2) then  !  boundary  j-
                     do  j=1,LAP
		             do  i=ibegin1,iend1-1  
		                 Ua(:,i-ibegin1+1,j)=B1%FRD(:,i,jbegin1+j-1,iv,jv)
                     enddo
		             enddo

	               else if( Bc1%face .eq. 3) then    !  boundary  i+
                     do j=jbegin1, jend1-1
		             do i=1,LAP  
                        Ua(:,j-jbegin1+1,i)=B1%FRD(:,iend1-i,j,iv,jv)  
                     enddo
		             enddo

	               else                                !  boundary  j+
                     do j=1,LAP
	                 do i=ibegin1,iend1-1
                       Ua(:,i-ibegin1+1,j)=B1%FRD(:,i,jend1-j,iv,jv)
                     enddo
		             enddo

	               endif

!                                    把临时数组Ua中的信息写入缓冲区
       orient=Bc%orient              ! orient==2 顺序;   其他值为 逆序
       do k=1,n1
        if(orient .eq. 2) then
          k1=k
        else
          k1=n1+1-k
        endif
       if(Bc%face .eq. 1 ) then             !  boundary  i-
         do i=1,LAP  
            B%FRD(:,ibegin2-i,jbegin2+k-1,iv,jv)=Ua(:,k1,i)
         enddo
	   else if(Bc%face .eq. 2 ) then        !  boundary  j-
         do j=1,LAP   
          B%FRD(:,ibegin2+k-1,jbegin2-j,iv,jv)=Ua(:,k1,j) 
         enddo
	   else if (Bc%face .eq. 3 ) then       !  boundary  i+
         do i=1,LAP
	       B%FRD(:,iend2+i-1,jbegin2+k-1,iv,jv)=Ua(:,k1,i)
         enddo
	   else                                          !  boundary  j+
         do j=1,LAP
	       B%FRD(:,ibegin2+k-1,jend2+j-1,iv,jv)=Ua(:,k1,j)
         enddo
      endif

       enddo
       
       end if
      deallocate(Ua)
     endif
    enddo

            end if
            
            
    enddo
       
       
      END DO
      END DO
      
       
       
    end subroutine GKUA_inner_mesco_boundary
    
      Subroutine GKUA_inner_macro_boundary(nMesh) !!7-2:19=>10
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Integer :: nMesh,mBlock
      integer::iv,jv,iv1,jv1,i,j,nx,ny,ird,ksub,k,k1,l
      integer :: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2
      integer :: N1,m_neighbour,msub,orient
      real(sp) ,allocatable:: Ua(:,:,:)
      Type (Block_TYPE),pointer:: B,B1
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc,BC1
       MP=>Mesh(1)   

       
         do mBlock=1,MP%Num_Block
            B => MP%Block(mBlock)
            nx=B%nx; ny=B%ny
            if(B%solver==GKUA)then
            do  ksub=1,B%subface
                Bc=> B%bc_msg(ksub)
             if( Bc%neighb .ge. 0 ) then   ! inner boundary
                  ibegin2=Bc%ist; iend2=Bc%iend; jbegin2=Bc%jst; jend2=Bc%jend       ! write
                  n1=iend2-ibegin2+jend2-jbegin2     ! +1 (number of cells)
    !              把连接点的信息读入临时数组Ua
                  allocate(Ua(19,n1,LAP))             !存储交换信息的临时数组
                  m_neighbour=Bc%neighb;  msub= Bc%subface   
                  B1 =>  Mesh(1)%Block(m_neighbour)
                  if(B1%solver==GKUA)then
                  Bc1 => B1%bc_msg(msub)
                 ibegin1=Bc1%ist; iend1=Bc1%iend; jbegin1=Bc1%jst; jend1=Bc1%jend   ! read
                   if(Bc1%face .eq. 1 ) then                     !  boundary  i-
                     do j=jbegin1,jend1-1
		             do i=1,LAP   
                        Ua(1:10,j-jbegin1+1,i)=B1%S_GKUA(1:10,ibegin1+i-1,j) 
                         
                     enddo
	                 enddo

                   else if ( Bc1%face .eq. 2) then  !  boundary  j-
                      do  i=ibegin1,iend1-1    
                     do  j=1,LAP
		                 Ua(1:10,i-ibegin1+1,j)=B1%S_GKUA(1:10,i,jbegin1+j-1)
                               
                     enddo
		             enddo

	               else if( Bc1%face .eq. 3) then    !  boundary  i+
                     do j=jbegin1, jend1-1
		             do i=1,LAP  
                         Ua(1:10,j-jbegin1+1,i)=B1%S_GKUA(1:10,iend1-i,j)
                 
                     enddo
		             enddo

                   else    
                      do i=ibegin1,iend1-1       !  boundary  j+
                      do j=1,LAP
	                      Ua(1:10,i-ibegin1+1,j)=B1%S_GKUA(1:10,i,jend1-j)
                     
                     enddo
		             enddo

	               endif

!                                    把临时数组Ua中的信息写入缓冲区
       orient=Bc%orient              ! orient==2 顺序;   其他值为 逆序
       do k=1,n1
        if(orient .eq. 2) then
          k1=k
        else
          k1=n1+1-k
        endif
       if(Bc%face .eq. 1 ) then             !  boundary  i-
         do i=1,LAP  
        
               B%S_GKUA(1:10,ibegin2-i,jbegin2+k-1)=Ua(1:10,k1,i)
         
            
         enddo
	   else if(Bc%face .eq. 2 ) then        !  boundary  j-
         do j=1,LAP   
             
                B%S_GKUA(1:10,ibegin2+k-1,jbegin2-j)=Ua(1:10,k1,j) 
            
         enddo
	   else if (Bc%face .eq. 3 ) then       !  boundary  i+
         do i=1,LAP
             
	            B%S_GKUA(1:10,iend2+i-1,jbegin2+k-1)=Ua(1:10,k1,i)
             
         enddo
	   else                                          !  boundary  j+
         do j=1,LAP
           
	          B%S_GKUA(1:10,ibegin2+k-1,jend2+j-1)=Ua(1:10,k1,j)
         
              
         enddo
      endif

       enddo
       
       end if
      deallocate(Ua)
     endif
    enddo

  
            end if
            
    enddo
       
       

      
       
       
    end subroutine GKUA_inner_macro_boundary
    

    
