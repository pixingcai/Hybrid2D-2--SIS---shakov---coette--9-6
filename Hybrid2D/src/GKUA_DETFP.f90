
!     This routine "DETFP" solve the flow properties using the reduced velocity
!        distributions and dicrete velocity coordinate method.
    SUBROUTINE DETFP(iTime)
    use Global_var
    use Com_ctrl
    use const_var ,only :sp,dop
    implicit none
    integer :: i,j,m,nn,idim,drct,posi,iv,jv,iv1,jv1,ITIME,NIstep
    integer :: steptmp(ndim),stedtmp(2,ndim)
    integer :: steptmp2(ndim),ij(ndim),ijws(ndim) !,ijp(ndim)
    character*80 fname
    real(sp)  :: IAMatrix(6,6),Amatrix(6,6)
    real(sp) ,external::Omiga22
    real(sp)  :: S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,ERRm,EPR,CF0IJ,BRHSTMP,tmp
    real(sp)  :: RONIJ,Uij,Vij,us1,vs1,sum_uv,Ttraij,Trotij,Toxxij,Toyyij,Toxyij,Qtra_x,Qtra_y,Qrot_x,Qrot_y,Pij,CFij,CFij_tmp
    real(sp)  :: taoxx,taoyy,taoxy,taozz,coetmp,qxij,qyij,tij
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
    ERR=0.0
    ERRm=-1.E30
    Amatrix=0.;IAMatrix=0.
    TEMP1=0.0
    do m=1,MP%Num_Block
        B => MP%Block(m)
        if(B%solver==GKUA)then
            if(Quadrature_Mode==GKUA_GH)then
                call IntegralGaussHermite(m)
            elseif(Quadrature_Mode==GKUA_GL)then  
                call GaussLegendre(m)
            else if(Quadrature_Mode==GKUA_NC)then   
                  !    call NewtonCotes(m)
                write(*,*)" waiting for further coding"
                stop
            else
                stop
            endif
            
        end if

    enddo
    
    out_mass=0.0
    do m=1,MP%Num_Block
        B => MP%Block(m)
        if(B%solver==GKUA)then
            ERR=0.0
            DO j=0,B%ny
            DO i=0,B%nx
                S1=B%S_GKUA(1,i,j); S2=B%S_GKUA(2,i,j); S3=B%S_GKUA(3,i,j); S4=B%S_GKUA(4,i,j)
                S5=B%S_GKUA(5,i,j); S6=B%S_GKUA(6,i,j); S7=B%S_GKUA(7,i,j); S8=B%S_GKUA(8,i,j)
                S9=B%S_GKUA(9,i,j); S10=B%S_GKUA(10,i,j)
                RONIJ=S1
                if(RONij.LE.0.)RONij=R00
                Uij=B%S_GKUA(2,i,j)/RONIJ
                B%S_GKUA(2,i,j)=Uij
                Vij=B%S_GKUA(3,i,j)/RONIJ
                B%S_GKUA(3,i,j)=Vij
                us1 = Uij * Uij
                vs1 = Vij * Vij
                sum_uv = us1+vs1
                Ttraij=(S4-RONIJ*sum_uv)/(1.5*RONIJ)
                if(Ttraij.LE.0.)Ttraij=T00
                B%S_GKUA(4,i,j)=Ttraij
                Toxxij=(B%S_GKUA(5,i,j)-RONIJ*us1)*2.
                B%S_GKUA(5,i,j)=Toxxij
                Toyyij=(B%S_GKUA(6,i,j)-RONIJ*vs1)*2.
                B%S_GKUA(6,i,j)=Toyyij
                Toxyij=(B%S_GKUA(7,i,j)-RONIJ*Uij*Vij)*2.
                B%S_GKUA(7,i,j)=Toxyij
                qxij=B%S_GKUA(8,i,j)-2.*Uij*S5-2.*Vij*S7+RONIJ*Uij*(sum_uv-1.5*Ttraij)
                B%S_GKUA(8,i,j)=qxij
                qyij=B%S_GKUA(9,i,j)-2.*Vij*S6-2.*Uij*S7+RONIJ*Vij*(sum_uv-1.5*Ttraij)
                B%S_GKUA(9,i,j)=qyij
                B%S_GKUA(10,i,j)=B%S_GKUA(10,i,j) *2. !taozz
       
                B%CF0=B%CF
                if(KHS==VHS_Mode)then
                    if(x1s .le. 1.e-8) then
                        B%CF(i,j)=cfc*B%S_GKUA(1,i,j)
                    else
                        B%CF(i,j)=cfc*B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j) **x1s
                    endif
                else if(KHS==Sutherland_Mode)then  
                    Miu_tmp=(vis_ref/vis_inf)*(B%S_GKUA(4,i,j)*T_inf/Tref)**(3.0/2.0)*(Sgas+Tref)/(Sgas+B%S_GKUA(4,i,j)*T_inf)
                    B%CF(i,j)=cfc*B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j)/Miu_tmp
                else if (KHS==LJ_Mode)then
                    tmp=Omiga22(B%S_GKUA(4,i,j)*T_inf/TSLJ)/Omiga22(T_inf/TSLJ)
                    Miu_tmp=sqrt(B%S_GKUA(4,i,j))/tmp
                    B%CF(i,j)=cfc*B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j)/Miu_tmp
                else  
                    if(myid==0)then
                        write(*,*)"Viscosity Model error"
                        stop
                    end if
                end if
            

     
                if(iTime.EQ.Time_Method) then
                    if (B%CF0(i,j) .GT. 1.E-12) then
                        EPR=((B%CF(i,j)-B%CF0(i,j))/B%CF0(i,j))*((B%CF(i,j)-B%CF0(i,j))/B%CF0(i,j))
                    else
                        EPR=(B%CF(i,j)-B%CF0(i,j))*(B%CF(i,j)-B%CF0(i,j))
                    endif
                    write(*,*)i,j,EPR,B%CF(i,j),B%CF0(i,j)
                    ERR=ERR+EPR
       
                endif  
                

            end do
            end do
            continue
            if(iTime.EQ.Time_Method) then
                if (nf .ge. 1) then
                    ERR=SQRT(ERR/(dt_global *dt_global ))/(1.0*(B%nx+3)*(B%ny+3))
                else
                    ERR=SQRT(ERR)/(1.0*(B%nx+3)*(B%ny+3))
                endif
                ERRm=max(ERRm,ERR)
         endif
         
         
         do j=1,B%ny-1
	     do i=1,B%nx-1
                out_mass=out_mass+B%S_GKUA(1,i,j)*B%vol(i,j)
         end do
         end do
         end if  
      end do

      
       ERR=ERRm 
       

       
      call GKUA_DET_wall(1)       
      continue
      
      RETURN
      END
      
    subroutine MassCorrect
    use Global_var
    use Com_ctrl
    use const_var ,only :sp,dop
    implicit none
    integer :: i,j,m
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
    MP=>Mesh(1)
    do m=1,MP%Num_Block
        B => MP%Block(m)
        if(B%solver==GKUA)then
            do j=1,B%ny-1
            do i=1,B%nx-1
                B%S_GKUA(1,i,j)=B%S_GKUA(1,i,j)*Init_mass/out_mass
            end do
            end do
        endif
    end do
    end subroutine MassCorrect

      SUBROUTINE IntegralGaussHermite(nm)
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
     ! allocate(S0(Nsum,0:B%nx,0:B%ny))
      B%S_GKUA=0.0
      B%HOT_GKUA=0.0
      B%HOQ_GKUA=0.0
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
            GDVIJ=B%FRD(1,i,j,iv,jv)
            HDVIJ=B%FRD(2,i,j,iv,jv)
          !  LDVIJ=B%FRD(3,i,j,iv,jv)
            
            WDGVXY=WDVX*WDVY*GDVIJ
            WDHVXY=WDVX*WDVY*HDVIJ
       !     WDLVXY=WDVX*WDVY*LDVIJ
            
            B%S_GKUA(1,i,j)=B%S_GKUA(1,i,j)+WDGVXY        !n
            B%S_GKUA(2,i,j)=B%S_GKUA(2,i,j)+WDGVXY*viiv   !U
            B%S_GKUA(3,i,j)=B%S_GKUA(3,i,j)+WDGVXY*vjjv   !V
            B%S_GKUA(4,i,j)=B%S_GKUA(4,i,j)+WDGVXY*(viiv2+vjjv2)+WDHVXY  !t''
  
            
            B%S_GKUA(5,i,j)=B%S_GKUA(5,i,j)+WDGVXY*viiv2   !txx''
            B%S_GKUA(7,i,j)=B%S_GKUA(7,i,j)+WDGVXY*viiv*vjjv  !txy''
            B%S_GKUA(6,i,j)=B%S_GKUA(6,i,j)+WDGVXY*vjjv2   !tyy''
            
            B%S_GKUA(8,i,j)=B%S_GKUA(8,i,j)+viiv*(WDGVXY*(viiv2+vjjv2)+WDHVXY)  !qx''
            B%S_GKUA(9,i,j)=B%S_GKUA(9,i,j)+vjjv*(WDGVXY*(viiv2+vjjv2)+WDHVXY)  !qy''
            B%S_GKUA(10,i,j)=B%S_GKUA(10,i,j)+WDHVXY  !tzz''
            
            B%HOT_GKUA(1,i,j)=B%HOT_GKUA(1,i,j)+WDGVXY*viiv2*viiv   !!VX**3   
            B%HOT_GKUA(2,i,j)=B%HOT_GKUA(2,i,j)+WDGVXY*viiv2*vjjv    !!VX**2*VY
            B%HOT_GKUA(3,i,j)=B%HOT_GKUA(3,i,j)+WDGVXY*viiv*vjjv2   !!VX*VY**2
            B%HOT_GKUA(4,i,j)=B%HOT_GKUA(4,i,j)+WDGVXY*vjjv*vjjv2   !!VY*VY**2
            B%HOT_GKUA(5,i,j)=B%HOT_GKUA(5,i,j)+WDHVXY*viiv   !!VXV*VZ**2
            B%HOT_GKUA(6,i,j)=B%HOT_GKUA(6,i,j)+WDHVXY*vjjv   !!VXV*VZ**2
            
            B%HOQ_GKUA(1,i,j)=B%HOQ_GKUA(1,i,j)+(WDGVXY*(viiv2+vjjv2)+WDHVXY)*viiv2        !!!vx**2*V2
            B%HOQ_GKUA(2,i,j)=B%HOQ_GKUA(2,i,j)+(WDGVXY*(viiv2+vjjv2)+WDHVXY)*viiv*vjjv     !!vx*vy*V2
            B%HOQ_GKUA(3,i,j)=B%HOQ_GKUA(3,i,j)+(WDGVXY*(viiv2+vjjv2)+WDHVXY)*vjjv2       !!vy**2*V2
            
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      DO j=0,B%ny
      DO i=0,B%nx
      
         DO kpp=1,Nsum
             
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
             call MPI_allreduce(B%S_GKUA(kpp,i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      
              B%S_GKUA(kpp,i,j)=temp
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         END DO
      END DO
      END DO
      
      DO j=0,B%ny
      DO i=0,B%nx
      
         DO kpp=1,6
             
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
             call MPI_allreduce(B%HOT_GKUA(kpp,i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      
              B%HOT_GKUA(kpp,i,j)=temp
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         END DO
      END DO
      END DO
      
      DO j=0,B%ny
      DO i=0,B%nx
      
         DO kpp=1,3
             
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
             call MPI_allreduce(B%HOQ_GKUA(kpp,i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      
              B%HOQ_GKUA(kpp,i,j)=temp
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         END DO
      END DO
      END DO
      
      RETURN
      END SUBROUTINE
!c-----Use 4,9,10-points" Gauss-Legendre integration formula to solve the moments from 
!c       micro-kinetics to macro-properties in each discrete velocity son-space. 
      SUBROUTINE GaussLegendre(nm)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp) ,allocatable::S0(:,:,:)
      integer :: iv,jv,iv1,jv1,i,j,icount,kpp,nm,ivp,jvp
      real(sp)  :: VX,VY,Sjifenxs,GDV1,HDV1,LDV1,VX1GDV,VY1GDV,VX1HDV,VY1HDV,VX2GDV,VY2GDV
      real(sp) :: SY(Nsum),FBEI(Nsum),temp
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      
      
       MP=>Mesh(1)   
       B => MP%Block(nm)
       B%S_GKUA=0.0
       
      DO jv=GKUA_JST,GKUA_JEND  !1:nyab*KGLN
         jv1=jv+GKUA_JPT
         VY=B%vj(jv1)
         jvp=(jv1-1)/KGLN+1         
      DO iv=GKUA_IST,GKUA_IEND  !1:nxab*KGLN
         iv1=iv+GKUA_IPT
         VX=B%vi(iv1)
         ivp=(iv1-1)/KGLN+1

         Sjifenxs=B%FBAVX(ivp)*B%wi(iv1)*B%FBAVY(jvp)*B%wj(jv1)
         DO j=0,B%ny
         DO i=0,B%nx
            GDV1=B%FRD(1,i,j,iv,jv)
            HDV1=B%FRD(2,i,j,iv,jv)
            
            VX1GDV=VX*GDV1
            VY1GDV=VY*GDV1


            VX2GDV=VX*VX1GDV
            VY2GDV=VY*VY1GDV

            FBEI(1)=GDV1 !n
            FBEI(2)=VX1GDV !U'
            FBEI(3)=VY1GDV !V'
            FBEI(4)=HDV1+VX2GDV+VY2GDV !Ttra'
            FBEI(5)=VX2GDV
            FBEI(6)=VY2GDV
            FBEI(7)=VX*VY*GDV1
            FBEI(8)=VX*FBEI(4) !(HDV1+VX2GDV+VY2GDV)
            FBEI(9)=VY*FBEI(4) !(HDV1+VX2GDV+VY2GDV)
            FBEI(10)=HDV1
            DO kpp=1,Nsum
               B%S_GKUA(kpp,i,j)=B%S_GKUA(kpp,i,j)+Sjifenxs*FBEI(kpp)
            ENDDO
 
            
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      
      DO j=0,B%ny
      DO i=0,B%nx
     
         DO kpp=1,Nsum
             
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
             call MPI_allreduce(B%S_GKUA(kpp,i,j),temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 
              B%S_GKUA(kpp,i,j)=temp
             call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         END DO
      END DO
      END DO
      continue
      RETURN
      END subroutine GaussLegendre
    
!c-----Use the equal-space Newton-Cotes repeated-Simpson integration formula to solve the 
!c       moments from micro-kinetics to macro-properties on the discrete velocity space. 
!c

      SUBROUTINE NewtonCotes(nm)   
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp),allocatable::S0(:,:,:)
      real(sp) QSY(Nsum),QF(Nsum)
      integer :: iv,jv,iv1,jv1,i,j,icount,kps,nm
      real(sp) :: VX,VY,GDV1,HDV1,VX1GDV,VY1GDV,VX1HDV,VY1HDV,VX2GDV,VY2GDV
      real(sp) :: xsxjifen,xsyjifen
      
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
       MP=>Mesh(1)   
       B => MP%Block(nm)
      allocate(S0(Nsum,0:B%nx,0:B%ny))
      B%S_GKUA=0.0;S0=0.
      xsxjifen=DVX/3.0;xsyjifen=DVY/3.0
   
      DO j=0,B%ny
      DO i=0,B%nx
         DO iv=GKUA_IST,GKUA_IEND
            iv1=iv+GKUA_IPT
            VX=B%vi(iv1)           
            DO kps=1,Nsum
               QSY(kps)=0.0
            ENDDO
         DO jv=GKUA_JST,GKUA_JEND    !!! DO jv=1,nvj
            jv1=jv+GKUA_JPT
            VY=B%vj(jv1)
            GDV1=B%FRD(1,i,j,iv,jv)
            HDV1=B%FRD(2,i,j,iv,jv)
            VX1GDV=VX*GDV1
            VY1GDV=VY*GDV1
            VX1HDV=VX*HDV1
            VY1HDV=VY*HDV1
            VX2GDV=VX*VX1GDV
            VY2GDV=VY*VY1GDV

            QF(1)=GDV1
            QF(2)=VX1GDV
            QF(3)=VY1GDV
            QF(4)=HDV1+VX2GDV+VY2GDV
            QF(5)=VX2GDV
            QF(6)=VY2GDV
            QF(7)=VX*VY*GDV1
            QF(8)=VX*QF(4)
            QF(9)=VY*QF(4)
            QF(10)=HDV1
            DO kps=1,Nsum
               if (jv1.eq.1) then           !!! if (jv.eq.1) then
                  QSY(kps)=QSY(kps)+xsyjifen*QF(kps)
               else if (jv1.eq.nvjt) then    !!! else if (jv.eq.nvj) then
                  QSY(kps)=QSY(kps)+xsyjifen*QF(kps)
               else
                  if (mod(jv1,2).eq.0) then   !!! if (mod(jv,2).eq.0) then
                    QSY(kps)=QSY(kps)+xsyjifen*4.0*QF(kps)
                  else
                    QSY(kps)=QSY(kps)+xsyjifen*2.0*QF(kps)
                  endif
               endif
            ENDDO
         ENDDO

         DO kps=1,Nsum
            if (iv1.eq.1) then             !!!  if (iv.eq.1) then
               B%S_GKUA(kps,i,j)=B%S_GKUA(kps,i,j)+xsxjifen*QSY(kps)
            else if (iv1.eq.nvit) then      !!!  else if (iv.eq.nvi) then 
               B%S_GKUA(kps,i,j)=B%S_GKUA(kps,i,j)+xsxjifen*QSY(kps)
            else
               if (mod(iv1,2).eq.0) then    !!! if (mod(iv,2).eq.0) then
                  B%S_GKUA(kps,i,j)=B%S_GKUA(kps,i,j)+xsxjifen*4.*QSY(kps)
               else
                  B%S_GKUA(kps,i,j)=B%S_GKUA(kps,i,j)+xsxjifen*2.*QSY(kps)
               end if
            end if
         ENDDO
         ENDDO
      ENDDO
      ENDDO



      RETURN
      END    
