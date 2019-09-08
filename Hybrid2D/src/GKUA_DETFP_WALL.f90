       Subroutine GKUA_DET_wall(nMesh)
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
      do mBlock=1,MP%Num_Block
         B => MP%Block(mBlock)
         if(B%solver==GKUA)then
         do  ksub=1,B%subface
             Bc=> B%bc_msg(ksub)
             if(Bc%neighb .eq. BC_Wall)then
                 call DETFPwall(mBlock,ksub)
             else if  ( Bc%neighb .eq. BC_Farfield ) then
                 call DETFPwall(mBlock,ksub)
             end if
             
         end do 
         end if
      end do
      
      return      
    end Subroutine GKUA_DET_wall   
    
    
    SUBROUTINE DETFPwall(nm,nn)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: RONIJ,Uij,Vij,Ttij,Trij,CF0ij,CFij
      real(sp)   :: Tij,Ptij,Pij,qxij,qyij
      real(sp)   :: CONQ,EPR
      real(sp) :: TOXXK,TOYYK,TOZZK,TOXYK,TOYZK,TOXZK
      integer:: i,j,kps,ilp,iv,jv,iv1,jv1
      integer :: ibd,idim,njpblk,njpbd,njpdim,ijktmp1,ijktmp2,ijktmp3,insum,ijp
      real(sp),external::fi 
      integer :: nm,nn,drct,nega,posi
      integer :: ij(ndim),ijwall(ndim),ijwall_1(ndim),ijwall_0(ndim)
      real(sp) :: sum_uv,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,Trotij,TEMP
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y,nx,ny,icheck,jcheck
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
    
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
      
    !   call DETQnwall(nm,nn)


      do i=ibegin,iend    
        do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             
             if((Bc%face .eq. 1)  ) then 
                 pos_x=1;pos_y=j;ijwall_1(1)=ibegin-1;ijwall_1(2)=j;ijwall_0(1)=ibegin;ijwall_0(2)=j;ijwall(1)=ibegin+1;ijwall(2)=j
             else if((Bc%face .eq. 2) ) then
                 pos_x=i;pos_y=1;ijwall_1(1)=i;ijwall_1(2)=jbegin-1;ijwall_0(1)=i;ijwall_0(2)=jbegin;ijwall(1)=i;ijwall(2)=jbegin+1
             else if(Bc%face .eq.3)then
                pos_x=1;pos_y=j;ijwall_1(1)=iend;ijwall_1(2)=j;ijwall_0(1)=iend-1;ijwall_0(2)=j;ijwall(1)=iend-2;ijwall(2)=j
             else if((Bc%face .eq. 4))then
                pos_x=i;pos_y=1;ijwall_1(1)=i;ijwall_1(2)=jend;ijwall_0(1)=i;ijwall_0(2)=jend-1;ijwall(1)=i;ijwall(2)=jend-2
             end if
             DO kps=1,10  !7-2
               
               BC%Swall(kps,pos_x,pos_y)=(-B%S_GKUA(kps,ijwall(1),ijwall(2))+3.0*B%S_GKUA(kps,ijwall_0(1),ijwall_0(2)))*0.5
               B%S_GKUA(kps,ijwall_1(1),ijwall_1(2))=2.0*BC%Swall(kps,pos_x,pos_y)-B%S_GKUA(kps,ijwall_0(1),ijwall_0(2))
             
             ENDDO
             
           

       enddo
     enddo
     ! 
     !
     ! continue
      RETURN
    END SUBROUTINE DETFPwall
    
    
    
       
    SUBROUTINE GaussLegendrewall(nm,nn)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy,viiv2,vjjv2,viiv3,vjjv3,viiv,vjjv
      real(sp)     :: HDVIJ,GDVIJ,WDGVXY,WDHVXY
      real(sp)   :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      real(sp)   :: FBEI(Nsum)
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny,kps
      integer:: jvp,ivp,mpi_begin
      real(sp) :: Sjifenxs,GDV1,HDV1,LDV1,VX1GDV,VY1GDV,VX2GDV,VY2GDV,VX1HDV,VY1HDV
      real(sp)  :: sumnw0(maxnode),sumnw(maxnode),sumnw1(maxnode),sumnw01(maxnode)
      integer :: ij(ndim),ijws(ndim),ijwall(ndim),steptmp,stedtmp(2,ndim)
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y,icheck,jcheck
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
       steptmp=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
      
      BC%Swall =0.
     ! BC%Swall0=0.
      FBEI=0.
      
      
      DO jv=GKUA_JST,GKUA_JEND
         jv1=jv+GKUA_JPT
         VY=B%vj(jv1)
         jvp=(jv1-1)/KGLN+1
         vs1 = B%FBAVY(jvp)*B%wj(jv1)
      DO iv=GKUA_IST,GKUA_IEND
         iv1=iv+GKUA_IPT
         VX=B%vi(iv1)
         ivp=(iv1-1)/KGLN+1
         Sjifenxs=B%FBAVX(ivp)*B%wi(iv1)*B%FBAVY(jvp)*B%wj(jv1) 

        do i=ibegin,iend    
        do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             
             if(Bc%face .eq. 1 ) then 
                    pos_x=1;pos_y=j;ijwall(1)=-3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=1.0;mpi_begin=jbegin
             else if(Bc%face .eq. 2 ) then 
                    pos_x=i;pos_y=1;ijwall(1)=i;ijwall(2)=-3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=1.0;mpi_begin=ibegin
             else if(Bc%face .eq. 3 ) then 
                    pos_x=1;pos_y=j;ijwall(1)=nx+3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=-1.0;mpi_begin=jbegin
             else if(Bc%face .eq. 4 ) then 
                   pos_x=i;pos_y=1;ijwall(1)=i;ijwall(2)=ny+3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=-1.0;mpi_begin=ibegin
             end if
            GDV1=B%FRD(1,ijwall(1),ijwall(2),iv,jv)
            HDV1=B%FRD(2,ijwall(1),ijwall(2),iv,jv)
            LDV1=B%FRD(3,ijwall(1),ijwall(2),iv,jv)
               
            VX1GDV=VX*GDV1
            VY1GDV=VY*GDV1
 
            VX2GDV=VX*VX1GDV
            VY2GDV=VY*VY1GDV

            FBEI( 1)=GDV1    !n
            FBEI( 2)=VX1GDV  !U'
            FBEI( 3)=VY1GDV  !V'
            FBEI( 4)=HDV1+VX2GDV+VY2GDV  !Ttra'
            FBEI( 5)=LDV1 !Trot'
            FBEI( 7)=VX2GDV  !Taoxx'
            FBEI( 8)=VX*VY*GDV1 !Taoxy'
            FBEI( 9)=VY2GDV !Taoyy'
             FBEI(10)=HDV1 !Taozz'
            FBEI( 8)=VX*(HDV1+VX2GDV+VY2GDV)
            FBEI( 9)=VY*(HDV1+VX2GDV+VY2GDV)
          
     
            vnw= ( B%vi(iv1) *swnx  +   B%vj(jv1) *swny )*steptmp
            FBEI(11)=vnw*GDV1

            DO kps=1,11
              BC%Swall(kps,pos_x,pos_y)=BC%Swall(kps,pos_x,pos_y)+Sjifenxs*FBEI(kps)
            ENDDO
        enddo
        enddo 
      ENDDO
      ENDDO
      icount = iend-ibegin+jend-jbegin
      
      DO kps=1,11
      
        do i=ibegin,iend    
        do j=jbegin,jend
             icheck=max(ibegin,iend)
             jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if((Bc%face .eq. 1).or.(Bc%face .eq. 3) ) then 
                    pos_x=1;pos_y=j;k=j
             else if((Bc%face .eq. 2).or.(Bc%face .eq. 4) ) then 
                    pos_x=i;pos_y=1;k=i
             end if
             sumnw0(k)= BC%Swall(kps,pos_x,pos_y)
        ENDDO
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        call MPI_allreduce(sumnw0(mpi_begin),sumnw(mpi_begin),icount,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        do i=ibegin,iend    
        do j=jbegin,jend
             icheck=max(ibegin,iend)
             jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if((Bc%face .eq. 1).or.(Bc%face .eq. 3) ) then 
                    pos_x=1;pos_y=j;k=j
             else if((Bc%face .eq. 2).or.(Bc%face .eq. 4) ) then 
                    pos_x=i;pos_y=1;k=i
             end if
             BC%Swall(kps,pos_x,pos_y)=sumnw(k)
       
        ENDDO
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      enddo 
    
     RETURN
    END SUBROUTINE GaussLegendrewall 
    
    
    !c
    SUBROUTINE DETQnwall(nm,nn)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy
      real(sp) :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny
      real(sp) :: starttime,endtime,steptimebegin
      real(sp)  :: sumnw0(maxnode),sumnw(maxnode),sumnw1(maxnode),sumnw01(maxnode)
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
       BC%Qnwall0=0.;BC%Qnwall =0.

         call QnwallGaussLegendre(nm,nn)
   
      RETURN
    END SUBROUTINE DETQnwall

    SUBROUTINE QnwallGaussLegendre(nm,nn)  !!waiting for further debug
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy
      real(sp) :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny,jvp,ivp,kps,mpi_begin
      real(sp) :: starttime,endtime,steptimebegin,Sjifenxs
      real(sp)  :: sumnw0(maxnode),sumnw(maxnode),sumnw1(maxnode),sumnw01(maxnode)
      integer :: ij(ndim),ijws(ndim),ijwall(ndim),steptmp,stedtmp(2,ndim)
      integer :: ibegin,iend,jbegin,jend,pos_x,pos_y,icheck,jcheck
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
       Bc=> B%bc_msg(nn)
       steptmp=0.0
       ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
      
      BC%Qnwall0=0.;BC%Qnwall =0.
      
      do i=ibegin,iend    
      do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
              if(Bc%face .eq. 1 ) then 
                    ijwall(1)=-3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=1.0
                    ijws(1)=1;ijws(2)=j;mpi_begin=jbegin
             else if(Bc%face .eq. 2 ) then 
                    ijwall(1)=i;ijwall(2)=-3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=1.0
                    ijws(1)=i;ijws(2)=1;mpi_begin=ibegin
             else if(Bc%face .eq. 3 ) then 
                    ijwall(1)=nx+3;ijwall(2)=j;swnx=B%ni1(i,j);swny=B%ni2(i,j);steptmp=-1.0
                     ijws(1)=1;ijws(2)=j;mpi_begin=jbegin
             else if(Bc%face .eq. 4 ) then 
                   ijwall(1)=i;ijwall(2)=ny+3;swnx=B%nj1(i,j);swny=B%nj2(i,j);steptmp=-1.0
                     ijws(1)=i;ijws(2)=1;mpi_begin=ibegin
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
            vnw=( B%vi(iv1)*swnx+B%vj(jv1)*swny )*steptmp
            VXYUVtmp= (B%vi(iv1)**2+B%vj(jv1)**2)*B%FRD(1,ijwall(1),ijwall(2),iv,jv)+B%FRD(2,ijwall(1),ijwall(2),iv,jv)
            VXYUVtmp= vnw*VXYUVtmp
            BC%Qnwall(1,ijws(1),ijws(2))=  BC%Qnwall(1,ijws(1),ijws(2))+Sjifenxs*VXYUVtmp
            VXYUVtmp= vnw*B%FRD(3,ijwall(1),ijwall(2),iv,jv)
            BC%Qnwall(2,ijws(1),ijws(2))=BC%Qnwall(2,ijws(1),ijws(2))+Sjifenxs*VXYUVtmp
         ENDDO
         ENDDO
      enddo
      enddo
      icount=iend-ibegin+jend-jbegin
     DO kps=1,2
      
        do i=ibegin,iend    
        do j=jbegin,jend
 
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if((Bc%face .eq. 1).or.(Bc%face .eq. 3) ) then 
                    pos_x=1;pos_y=j;k=j
             else if((Bc%face .eq. 2).or.(Bc%face .eq. 4) ) then 
                    pos_x=i;pos_y=1;k=i
             end if
             sumnw0(k)= BC%Qnwall(kps,pos_x,pos_y)
        ENDDO
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        call MPI_allreduce(sumnw0(mpi_begin),sumnw(mpi_begin),icount,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        
        do i=ibegin,iend    
        do j=jbegin,jend
            icheck=max(ibegin,iend)
            jcheck=max(jbegin,jend)
             if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                 cycle
             end if
             if((Bc%face .eq. 1).or.(Bc%face .eq. 3) ) then 
                    pos_x=1;pos_y=j;k=j
             else if((Bc%face .eq. 2).or.(Bc%face .eq. 4) ) then 
                    pos_x=i;pos_y=1;k=i
             end if
             BC%Qnwall(kps,pos_x,pos_y)=sumnw(k)
       
        ENDDO
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      enddo 
      
      RETURN
    END SUBROUTINE QnwallGaussLegendre    


    
    
    SUBROUTINE IntegralGaussHermitewall(nm,nn)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp)   :: WDVX,WDVY,vnw,vnnFDV,vx,vy,viiv2,vjjv2,viiv3,vjjv3,viiv,vjjv
      real(sp)     :: HDVIJ,GDVIJ,WDGVXY,WDHVXY
      real(sp) :: SIJ,swnx,swny,WDVXY,VNNGDV,us1,vs1,VXYUVtmp
      integer:: nm,nn,drct,posi,i,j,k,iv,jv,iv1,jv1,icount,nx,ny
      real(sp) :: starttime,endtime,steptimebegin
      real(sp)  :: sumnw0(maxnode),sumnw(maxnode),sumnw1(maxnode),sumnw01(maxnode)
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
      
      BC%Swall =0.
      BC%Swall0=0.
    
      DO jv=GKUA_JST,GKUA_JEND
         jv1 =jv+GKUA_JPT
         vjjv=B%vj(jv1)
         vjjv2=vjjv*vjjv
         vjjv3=vjjv2*vjjv
         WDVY=B%wj(jv1)*EXP(vjjv2)
      DO iv=GKUA_IST,GKUA_IEND
         iv1=iv+GKUA_IPT
         viiv=B%vi(iv1)
         viiv2=viiv*viiv
         viiv3=viiv2*viiv
         WDVX=B%wi(iv1)*EXP(viiv2)
    
        do i=ibegin,iend    
        do j=jbegin,jend
             if(Bc%face .eq. 1 ) then 
                    pos_x=1;pos_y=j;ijwall(1)=-3;ijwall(2)=j
             else if(Bc%face .eq. 2 ) then 
                    pos_x=i;pos_y=1;ijwall(1)=i;ijwall(2)=-3
             else if(Bc%face .eq. 3 ) then 
                    pos_x=nx-1;pos_y=j;ijwall(1)=nx+3;ijwall(2)=j
             else if(Bc%face .eq. 4 ) then 
                   pos_x=i;pos_y=ny-1;ijwall(1)=i;ijwall(2)=ny+3
             end if
         
            GDVIJ=B%FRD(1,ijwall(1),ijwall(2),iv,jv)
            HDVIJ=B%FRD(2,ijwall(1),ijwall(2),iv,jv)
            
            WDGVXY=WDVX*WDVY*GDVIJ
            WDHVXY=WDVX*WDVY*HDVIJ
          
            BC%Swall( 1,ijwall(1),ijwall(2))=BC%Swall( 1,ijwall(1),ijwall(2))+WDGVXY        !n
            BC%Swall( 2,ijwall(1),ijwall(2))=BC%Swall( 2,ijwall(1),ijwall(2))+WDGVXY*viiv   !U
            BC%Swall( 3,ijwall(1),ijwall(2))=BC%Swall( 3,ijwall(1),ijwall(2))+WDGVXY*vjjv   !V
            BC%Swall( 4,ijwall(1),ijwall(2))=BC%Swall( 4,ijwall(1),ijwall(2))+WDGVXY*(viiv2+vjjv2)+WDHVXY  !t''
            BC%Swall( 5,ijwall(1),ijwall(2))=BC%Swall( 5,ijwall(1),ijwall(2))+WDGVXY*viiv2   !txx''
            BC%Swall( 6,ijwall(1),ijwall(2))=BC%Swall( 6,ijwall(1),ijwall(2))+WDGVXY*vjjv2   !tyy''
            BC%Swall( 7,ijwall(1),ijwall(2))=BC%Swall( 7,ijwall(1),ijwall(2))+WDGVXY*viiv*vjjv  !txy''
            BC%Swall( 8,ijwall(1),ijwall(2))=BC%Swall( 8,ijwall(1),ijwall(2))+viiv*(WDGVXY*(viiv2+vjjv2)+WDHVXY)  !qx''
            BC%Swall( 9,ijwall(1),ijwall(2))=BC%Swall( 9,ijwall(1),ijwall(2))+vjjv*(WDGVXY*(viiv2+vjjv2)+WDHVXY)  !qy''
            BC%Swall(10,ijwall(1),ijwall(2))=BC%Swall(10,ijwall(1),ijwall(2))+WDHVXY  !tzz''
        enddo
        enddo
        
      ENDDO
      ENDDO
    
      !icount = Nsum*iend-ibegin+jend-jbegin 
      !
      !call MPI_allreduce(BC%Swall,BC%Swall0,icount,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
      !
      !BC%Swall=BC%Swall0
    
      RETURN
    END SUBROUTINE IntegralGaussHermitewall    
    
    
    