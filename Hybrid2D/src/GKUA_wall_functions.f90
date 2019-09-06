      SUBROUTINE GKUA_wall_boundary(nm,nn,wall_u,wall_v)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      integer ::nm,nn
      real(sp) :: wall_u,wall_v
      integer :: ij(ndim),ijws(ndim),ijwall(ndim),steptmp,stedtmp(2,ndim)
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y
      integer :: iv,jv,i,j,k,iv1,jv1,drct,posi,nx,ny,icheck,jcheck
      real(sp)  :: sumnw0(1000),sumnw(1000)
      real(sp) :: SIJ,swnx,swny,cvuv,gdvij,vnw,cuvw,TAT,TAR
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
       steptmp=0.0
       sumnw=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
      call DETNW(nm,nn,sumnw,wall_u,wall_v)

     
     continue
      do i=ibegin,iend    
        do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if(Bc%face .eq. 1 ) then 
                    ijwall(1)=-3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=1.0;k=j 
             else if(Bc%face .eq. 2 ) then 
                    ijwall(1)=i;ijwall(2)=-3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=1.0;k=i     
             else if(Bc%face .eq. 3 ) then 
                    ijwall(1)=nx+3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=-1.0;k=j
             else if(Bc%face .eq. 4 ) then 
                    ijwall(1)=i;ijwall(2)=ny+3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=-1.0;k=i
             end if
           
         DO jv=GKUA_JST,GKUA_JEND
            jv1=jv+GKUA_JPT
         DO iv=GKUA_IST,GKUA_IEND
            iv1=iv+GKUA_IPT
  
            vnw= ( (B%vi(iv1)-wall_u) *swnx  +   (B%vj(jv1)-wall_v) *swny )*steptmp
            if (vnw .ge. 0.0) then
               cvuv=(B%vi(iv1)-wall_u)**2+(B%vj(jv1)-wall_v)**2
               gdvij=sumnw(k)*exp(-cvuv/Twall)/(pi*Twall)
         
               if(sumnw(k).LT.0.)sumnw(k)=0.
               B%FRD(1,ijwall(1),ijwall(2),iv,jv)=gdvij
               B%FRD(2,ijwall(1),ijwall(2),iv,jv)=(Twall*gdvij/2.)
              ! B%FRD(3,ijwall(1),ijwall(2),iv,jv)=(Rdof*Twall*gdvij/2.)  !7-2

            endif
         ENDDO
         ENDDO
        end do
      end do
      continue

  
      return
      END SUBROUTINE GKUA_wall_boundary
    
      SUBROUTINE DETNW(nm,nn,sumnw,wall_u,wall_v)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   ::wall_u,wall_v
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy
      real(sp) :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny
      real(sp) :: starttime,endtime,steptimebegin
      real(sp)  :: sumnw0(1000),sumnw(1000),sumnw1(1000),sumnw01(1000)
      integer :: ij(ndim),ijws(ndim),ijwall(ndim),steptmp,stedtmp(2,ndim)
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
       steptmp=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
       sumnw =0.;sumnw0=0.;sumnw1 =0.;sumnw01 =0.
    IF (KGH.EQ.1) THEN
          call NWGaussHerimte(nm,nn,sumnw,wall_u,wall_v)
      ELSEIF (KGL.EQ.1) THEN
          call NWGaussLegendre(nm,nn,sumnw,wall_u,wall_v)
      ELSE
     !     call NWNewtonCotes(nm,nn,sumnw,wall_u)
      ENDIF
          
  
      RETURN
      END SUBROUTINE DETNW

      SUBROUTINE NWGaussLegendre(nm,nn,sumnw,wall_u,wall_v)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy,wall_u,wall_v
      real(sp) :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny,jvp,ivp
      real(sp) :: Sjifenxs,vnn,FBEI
      real(sp)  :: sumnw0(1000),sumnw(1000),sumnw1(1000),sumnw01(1000)
      integer :: ij(ndim),ijws(ndim),ijwall(ndim),steptmp,stedtmp(2,ndim)
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y,icheck,jcheck,mpi_begin
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
      
     
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
       steptmp=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
      
      sumnw =0.;sumnw0=0.;sumnw1 =0.;sumnw01 =0.
        do i=ibegin,iend    
        do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if(Bc%face .eq. 1 ) then 
                 k=j; ijwall(1)=-3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=1.0;mpi_begin=jbegin
             else if(Bc%face .eq. 2 ) then 
                 k=i;  ijwall(1)=i;ijwall(2)=-3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=1.0;mpi_begin=ibegin
             else if(Bc%face .eq. 3 ) then 
                 k=j;   ijwall(1)=nx+3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=-1.0;mpi_begin=jbegin
             else if(Bc%face .eq. 4 ) then 
                 k=i;  ijwall(1)=i;ijwall(2)=ny+3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=-1.0;mpi_begin=ibegin
             end if

         DO jv=GKUA_JST,GKUA_JEND
            jv1=jv+GKUA_JPT
            jvp=(jv1-1)/KGLN+1
            vy=B%vj(jv1)
            vs1 = B%FBAVY(jvp)*B%wj(jv1)
         DO iv=GKUA_IST,GKUA_IEND
            iv1=iv+GKUA_IPT
            ivp=(iv1-1)/KGLN+1
            vx=B%vi(iv1)
            Sjifenxs=B%FBAVX(ivp)*B%wi(iv1)*vs1 

            vnw= ( (B%vi(iv1) -wall_u)*swnx  +   (B%vj(jv1)-wall_v) *swny )*steptmp
            vnn=(vnw-abs(vnw))/2.
            FBEI=vnn*B%FRD(1,ijwall(1),ijwall(2),iv,jv)
 
            sumnw(k)=sumnw(k)+Sjifenxs*FBEI
            FBEI=(vnw+abs(vnw))/2.*exp(-((B%vi(iv1)-wall_u) **2+(B%vj(jv1)-wall_v) **2)/Twall)/(pi*Twall)
            sumnw1(k)=sumnw1(k)+Sjifenxs*FBEI 
         ENDDO
         continue
         ENDDO
         
        enddo
        enddo
         
 
     icount = iend-ibegin+jend-jbegin+1
      call MPI_allreduce(sumnw(mpi_begin),sumnw0(mpi_begin),icount,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      sumnw=sumnw0
      call MPI_allreduce(sumnw1(mpi_begin),sumnw01(mpi_begin),icount,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      sumnw1=sumnw01
      

         do i=ibegin,iend    
         do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if((Bc%face .eq. 1).or.(Bc%face .eq. 3) ) then 
                    k=j
             else if((Bc%face .eq. 2 ).or.(Bc%face .eq. 4)) then 
                    k=i
             end if
             sumnw(k)=-sumnw(k)/sumnw1(k)
   
         end do
         end do
    continue
      RETURN
      END SUBROUTINE NWGaussLegendre    
    
      SUBROUTINE NWGaussHerimte(nm,nn,sumnw,wall_u,wall_v)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy,wall_u,wall_v
      real(sp) :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny,jvp,ivp
      real(sp) :: Sjifenxs,vnn,FBEI
      real(sp)  :: sumnw0(1000),sumnw(1000),sumnw1(1000),sumnw01(1000)
      integer :: ij(ndim),ijws(ndim),ijwall(ndim),steptmp,stedtmp(2,ndim)
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y,icheck,jcheck,mpi_begin
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
      
     
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
       steptmp=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
      
      sumnw =0.;sumnw0=0.;sumnw1 =0.;sumnw01 =0.
        do i=ibegin,iend    
        do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if(Bc%face .eq. 1 ) then 
                 k=j; ijwall(1)=-3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=1.0;mpi_begin=jbegin
             else if(Bc%face .eq. 2 ) then 
                 k=i;  ijwall(1)=i;ijwall(2)=-3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=1.0;mpi_begin=ibegin
             else if(Bc%face .eq. 3 ) then 
                 k=j;   ijwall(1)=nx+3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=-1.0;mpi_begin=jbegin
             else if(Bc%face .eq. 4 ) then 
                 k=i;  ijwall(1)=i;ijwall(2)=ny+3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=-1.0;mpi_begin=ibegin
             end if

         DO jv=GKUA_JST,GKUA_JEND
            jv1=jv+GKUA_JPT
            vy=B%vj(jv1)
            WDVY=B%wj(jv1)*EXP(vy*vy)
         DO iv=GKUA_IST,GKUA_IEND
            iv1=iv+GKUA_IPT
            vx=B%vi(iv1)
            WDVX=B%wi(iv1)*EXP(vx*vx)
            WDVXY=WDVX*WDVY
            vnw= ( (B%vi(iv1)-wall_u )*swnx  +   (B%vj(jv1)-wall_v) *swny )*steptmp
            vnn=(vnw-abs(vnw))/2.
            FBEI=vnn*B%FRD(1,ijwall(1),ijwall(2),iv,jv)
 
            sumnw(k)=sumnw(k)+WDVXY*FBEI
            FBEI=(vnw+abs(vnw))/2.*exp(-((B%vi(iv1)-wall_u) **2+(B%vj(jv1)-wall_v) **2)/Twall)/(pi*Twall)
            sumnw1(k)=sumnw1(k)+WDVXY*FBEI 
         ENDDO
         ENDDO    
        enddo
        enddo
         
 
     icount = iend-ibegin+jend-jbegin+1
      call MPI_allreduce(sumnw(mpi_begin),sumnw0(mpi_begin),icount,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      sumnw=sumnw0
      call MPI_allreduce(sumnw1(mpi_begin),sumnw01(mpi_begin),icount,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      sumnw1=sumnw01
      

         do i=ibegin,iend    
         do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if((Bc%face .eq. 1).or.(Bc%face .eq. 3) ) then 
                    k=j
             else if((Bc%face .eq. 2 ).or.(Bc%face .eq. 4)) then 
                    k=i
             end if
             sumnw(k)=-sumnw(k)/sumnw1(k)
     
            
         end do
         end do
    continue
      RETURN
      END SUBROUTINE NWGaussHerimte 
    
    
!c-----Use the equal-space Newton-Cotes repeated-Simpson integration formula to solve the 
!c       moments from micro-kinetics to macro-properties on the discrete velocity space. 

