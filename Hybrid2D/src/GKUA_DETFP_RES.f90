
!C-1     This routine "DETFP" solve the flow properties using the reduced velocity
!C         distributions and dicrete velocity coordinate method.
    SUBROUTINE DETFP_res
    use Global_var
    use Com_ctrl
    use const_var ,only :sp,dop
    implicit none
    integer :: i,j,m,nn,idim,drct,posi,iv,jv,iv1,jv1,ITIME,NIstep,nx,ny,kps
    integer :: steptmp(ndim),stedtmp(2,ndim)
    integer :: steptmp2(ndim),ij(ndim),ijws(ndim) !,ijp(ndim)
    character*80 fname
    real(sp)  :: IAMatrix(6,6),Amatrix(6,6),res
    real(sp) ,external::Omiga22
    real(sp)  :: S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,ERRm,EPR,CF0IJ,BRHSTMP
    real(sp)  :: RONIJ,Uij,Vij,us1,vs1,sum_uv,Ttraij,Trotij,Toxxij,Toyyij,Toxyij,Qtra_x,Qtra_y,Qrot_x,Qrot_y,Pij,CFij,CFij_tmp
    real(sp)  :: taoxx,taoyy,taoxy,taozz,coetmp
    real(sp)  :: S1tmp,S2tmp,S3tmp,S4tmp,S5tmp,RONintmp,TEMP1
    integer :: ijp,njpblk,njpbd,njpdim,ijtmp1,IJTMP2
    real(sp)  :: DRONtmp,DTtrtmp,DUtmp,DVtmp,DQXTMP,DQYTMP
    real(sp)  :: Sentropyinf,cinf,cpt,FVXDV,TIJtmp,viu,vjv,Machlocal,Miu_tmp
    real(sp)  :: SBNX1,SBNY1,SBNX2,SBNY2,SBNXY,VTBRY1,VTBRY2,VNBRYACU,CBRYACU
    real(sp)  :: vij2,GDVijv,SIJ,sbnx,sbny,vninf,vtinf1,vtinf2,vnbry,vnbry1,vnbry2,cbry,riemann1,riemann2,unriemann,Sentropybry
    integer stats(MPI_STATUS_SIZE)
      
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
    MP=>Mesh(1)   
  
    do m=1,MP%Num_Block
        B => MP%Block(m)
        B%Res_max=0.0;B%Res_rms=0.0
        nx=B%nx; ny=B%ny
        if(B%solver==GKUA)then
            if(Quadrature_Mode==GKUA_GH)then
                call  GaussHermite_rs(m)
            elseif(Quadrature_Mode==GKUA_GL)then
                call GaussLegendre_rs(m)
            else if(Quadrature_Mode==GKUA_NC)then 
                write(*,*)" waiting for further coding"
                stop
            end if
  
        !   统计最大残差和均方根残差 
            do j=1,ny-1
            do i=1,nx-1
                do kps=1,4
                    Res= B%S00_GKUA(kps,i,j)
                    if(abs(Res) .gt. B%Res_max(kps))  B%Res_max(kps)=abs(Res)          ! 最大残差
                        B%Res_rms(kps)=B%Res_rms(kps)+Res*Res                           ! 均方根残差      
                enddo
            enddo
            enddo
            B%Res_rms(:)=sqrt(B%Res_rms(:)/(1.d0*Mesh(1)%Num_Cell))   !    全部网格点的总均方根残差
        endif
        
    enddo   
    continue
    return
    end subroutine DETFP_res

!c-----Use 4,9,10-points" Gauss-Legendre integration formula to solve the moments from 
!c       micro-kinetics to macro-properties in each discrete velocity son-space. 
      SUBROUTINE GaussLegendre_rs(nm)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
    
      integer :: iv,jv,iv1,jv1,i,j,icount,kpp,nm,ivp,jvp
      real(sp)  :: VX,VY,Sjifenxs,GDV1,HDV1,LDV1,VX1GDV,VY1GDV,VX1HDV,VY1HDV,VX2GDV,VY2GDV
      real(sp) :: SY(Nsum),FBEI(Nsum),temp
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      
      
       MP=>Mesh(1)   
       B => MP%Block(nm)
       B%S00_GKUA=0.0
       
      DO jv=GKUA_JST,GKUA_JEND  !1:nyab*KGLN
         jv1=jv+GKUA_JPT
         VY=B%vj(jv1)
         jvp=(jv1-1)/KGLN+1         
      DO iv=GKUA_IST,GKUA_IEND  !1:nxab*KGLN
         iv1=iv+GKUA_IPT
         VX=B%vi(iv1)
         ivp=(iv1-1)/KGLN+1

         Sjifenxs=B%FBAVX(ivp)*B%wi(iv1)*B%FBAVY(jvp)*B%wj(jv1)
         DO j=1,B%ny-1
         DO i=1,B%nx-1
            GDV1=B%FRD00(1,i,j,iv,jv)
            HDV1=B%FRD00(2,i,j,iv,jv)
         !   LDV1=B%FRD00(3,i,j,iv,jv)
            
            VX1GDV=VX*GDV1
            VY1GDV=VY*GDV1
            
            
            VX2GDV=VX*VX1GDV
            VY2GDV=VY*VY1GDV
            
            FBEI(1)=GDV1 !n
            FBEI(2)=VX1GDV !U'
            FBEI(3)=VY1GDV !V'
            FBEI(4)=HDV1+VX2GDV+VY2GDV !Ttra'
        !    FBEI(5)=LDV1
            
            
            
            DO kpp=1,4
               B%S00_GKUA(kpp,i,j)=B%S00_GKUA(kpp,i,j)+Sjifenxs*FBEI(kpp)
            ENDDO
 
            
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      
      DO j=1,B%ny-1
      DO i=1,B%nx-1
         DO kpp=1,4
             call MPI_allreduce(B%S00_GKUA(kpp,i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      
              B%S00_GKUA(kpp,i,j)=temp
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         END DO
         
      END DO
      END DO
      continue
      RETURN
    END subroutine GaussLegendre_rs
  
    !c-----Use 4,9,10-points" Gauss-Legendre integration formula to solve the moments from 
!c       micro-kinetics to macro-properties in each discrete velocity son-space. 
      SUBROUTINE GaussHermite_rs(nm)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      IMPLICIT NONE
      integer :: iv,jv,iv1,jv1,i,j,icount,kpp,nm
      real(sp) :: viiv,vjjv,viiv2,vjjv2,WDVX,WDVY,WDVXY,GDVIJ,HDVIJ,LDVIJ,WDGVXY,WDHVXY,WDLVXY,temp
      real(sp),allocatable::S0(:,:,:)
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
       MP=>Mesh(1)   
       B => MP%Block(nm)
       B%S00_GKUA=0.0
    
      
      DO jv=GKUA_JST,GKUA_JEND    
         jv1=jv+GKUA_JPT
      DO iv=GKUA_IST,GKUA_IEND   
         iv1=iv+GKUA_IPT
         viiv=B%vi(iv1)
         vjjv=B%vj(jv1)
         viiv2=viiv*viiv
         vjjv2=vjjv*vjjv
      
         WDVX=B%wi(iv1)*EXP(viiv2)
         WDVY=B%wj(jv1)*EXP(vjjv2)
         WDVXY=WDVX*WDVY
         DO j=0,B%ny
         DO i=0,B%nx
            GDVIJ=B%FRD00(1,i,j,iv,jv)
            HDVIJ=B%FRD00(2,i,j,iv,jv)

            
            WDGVXY=WDVX*WDVY*GDVIJ
            WDHVXY=WDVX*WDVY*HDVIJ
 
            
            B%S00_GKUA(1,i,j)=B%S00_GKUA(1,i,j)+WDGVXY        !n
            B%S00_GKUA(2,i,j)=B%S00_GKUA(2,i,j)+WDGVXY*viiv   !U
            B%S00_GKUA(3,i,j)=B%S00_GKUA(3,i,j)+WDGVXY*vjjv   !V
            B%S00_GKUA(4,i,j)=B%S00_GKUA(4,i,j)+WDGVXY*(viiv2+vjjv2)+WDHVXY  !t''
  
            

         ENDDO
         ENDDO
      ENDDO
      ENDDO
      
      DO j=1,B%ny-1
      DO i=1,B%nx-1
         DO kpp=1,4
             call MPI_allreduce(B%S00_GKUA(kpp,i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      
              B%S00_GKUA(kpp,i,j)=temp
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         END DO
         B%S00_GKUA(4,i,j)= B%S00_GKUA(4,i,j)+ B%S00_GKUA(5,i,j)
         
      END DO
      END DO
      END subroutine GaussHermite_rs
