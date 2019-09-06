!C
!C---- This routine "lisansudu" set discrete velocity points using the discrete velocity
!C       coordinate method.
      SUBROUTINE GKUA_lisansudu
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none 
      real(sp)  ::vd(nvt),wd(nvt),fbav(nvt)
      integer :: iv,jv,NV0,m,NV0X,NV0Y,iv0
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
       MP=>Mesh(1)   
       
!C---- Check dimension size on X-direction discrete velocity points whether to overflow: 
      do iv=1,nvt
         vd(iv)=0.
         wd(iv)=0.
      enddo
      IF (NVI.GT.NVIT) THEN
         if (myid == 0) then
            WRITE(*,*)'Why X-direction(NVI>NVIT)? & 
                   $   To continue RUN please input"enter"'
         endif
      ENDIF
!---- Check dimension size on Y-direction discrete velocity points whether or not to overflow: 
      IF (NVJ.GT.NVJT) THEN
         if (myid == 0) then
            WRITE(*,*)'Why Y-direction(NVJ>NVJT)? &
                   &   To continue RUN please input"enter"'
         endif
      ENDIF

      call lisanv(NVI,dvx,vxdown,vxup,VD,WD,FBAV,NV0)
      do m=1,MP%Num_Block
           B => MP%Block(m)
           if(B%solver==GKUA)then 
           do iv=1,NVI
              B%vi(iv)=vd(iv)
              B%wi(iv)=wd(iv)
           enddo
           end if
      end do
      if (kgl.eq.1) then
          
         NV0X=NV0
        do m=1,MP%Num_Block
            B => MP%Block(m)
            if(B%solver==GKUA)then 
               do iv0=1,nv0x
                  B%fbavx(iv0)=fbav(iv0)
               enddo
            end if
        end do
        WRITE(*,*)NV0
        CONTINUE
      endif
      call lisanv(NVJ,dvy,vydown,vyup,VD,WD,FBAV,NV0)
     do m=1,MP%Num_Block
        B => MP%Block(m)
          if(B%solver==GKUA)then 
            do jv=1,NVJ
               B%vj(jv)=vd(jv)
               B%wj(jv)=wd(jv)
            enddo
          end if
     end do
     
      if (kgl.eq.1) then
         NV0Y=NV0
         do m=1,MP%Num_Block
             B => MP%Block(m)
             if(B%solver==GKUA)then 
                do iv0=1,nv0y
                  B%fbavy(iv0)=fbav(iv0)
                enddo
             end if
         end do
        
      endif
      IF (myid == 0) then
        WRITE(*,*) 'X-direction-velocity point:NVIT,NVI=',NVIT,NVI
        WRITE(*,*) 'Y-direction-velocity point:NVJT,NVJ=',NVJT,NVJ
      endif
      IF (NVJ.GT.NVJT .or. NVI.GT.NVIT) THEN
         IF (myid == 0) then
         WRITE(*,*)'Why(NVI>NVIT,or NVJ>NVJT)-CONTINUE RUN,input"enter"'
         ENDIF
      END IF
      RETURN
      END   SUBROUTINE GKUA_lisansudu          

      SUBROUTINE lisanv(nv,dv,vdown,vup,VD,WD,FBAV,NV0)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      integer:: nv ,NV0,IV
      real(sp)  ::dv ,vdown,vup
      real(sp) :: vd(nvt),wd(nvt),VD0(nvt),WD0(nvt),FBAV(nvt),DVMAX,DVMIN,dviv1
      IF (KGH.EQ.1) THEN
        IF (NV.EQ.7 .OR. NV.EQ.8 .OR. NV.EQ.16) THEN
           !if (myid == 0 .and. KP.EQ.0 ) then
           !   WRITE(21,*)'Using Gauss-Hermite table to acquire & 
           !            &  discrete velocity ordinate points and weights'
           !endif
          call GaussHermite(VD0,WD0,NV)
        ENDIF
        DO IV=1,NV
           VD(IV)=-VD0(NV+1-IV)
           WD(IV)=WD0(NV+1-IV)
           VD(NV+IV)=VD0(IV)
           WD(NV+IV)=WD0(IV)
        ENDDO
        NV=2*NV
        if (myid == 0 .and. KPtmp.EQ.0 ) then
           !WRITE(21,*)'DISCRETE VELOCITY POINTS &
           !       &    AND WEIGHTS:VD(IV),WD(IV)'
           DO IV=1,NV
              !WRITE(21,*) iv,VD(IV),WD(IV)
           END DO
        endif
!c        vspace=vd(nv)
        vdown=-vd(nv)
        vup=vd(nv)
      ELSE
         DO IV=1,NV+1
            VD(IV)=vdown+REAL(IV-1)*DV
            VD0(IV)=VD(IV)
         ENDDO
         if (myid == 0) then
            WRITE(*,*)'As NVT>8,roots and weights of GAUSS-HERMITE"s &
                    &  integral are uneasily determined.'
         endif
         IF (KGL.EQ.1) THEN
            if (myid == 0 .and. KPtmp.EQ.0 ) then
               !WRITE(21,*)'Use Newton-Legendre integration fomula with &
               !        &   non-equal nodes and weights!'
               !WRITE(21,*)'use space nonequal nodes:KGLN=3,4,5,6,9,10'
            endif
            call VNODEQUAN(VD0,NV,VD,WD,FBAV,NV0)
            DVMAX=0.0
            DVMIN=10.
            DO IV=2,NV
               dviv1=VD(IV)-VD(iv-1)
               IF (DViv1.GT.DVMAX) DVMAX=dviv1
               IF (DVIV1.LT.DVMIN) DVMIN=DVIV1
            ENDDO
            DO IV=1,NV
               if (myid == 0 .and. KPtmp.EQ.0 ) then
                  !WRITE(21,*) iv,VD(IV),WD(IV)
               endif
            ENDDO
            if (myid == 0 .and. KPtmp.EQ.0 ) then
               WRITE( *,*)'Use KGLN=3,4,5,6,9,10-points" Gauss-Legendre&
          & integral"', ' in each DV son-space:DVMAX,DVMIN=',DVMAX,DVMIN
          !     WRITE(21,*)'Use KGLN=3,4,5,6,9,10-points" Gauss-Legendre&
          !& integral"', ' in each DV son-space:DVMAX,DVMIN=',DVMAX,DVMIN
            endif
         ELSE
           NV=NV+1
           if (myid == 0 .and. KPtmp.EQ.0 ) then
             WRITE(*,*)'   --Use Newton-cotes fomula with equal space--'
             WRITE(*,*) 'VD(IV)='
             !WRITE(21,*)'  --Use Newton-cotes fomula with equal space--'
             !WRITE(21,*) 'VD(IV)='
             DO IV=1,NV
                WRITE( *,*) iv,VD(IV)
                !WRITE(21,*) iv,VD(IV)
             END DO
             WRITE(*,*)'discrete vel. cellsize:DV=',DV
           endif
         END IF
      END IF
      RETURN
      END     

        
!c---- Discrete (Vx,Vy) velocity space by Gauss-Legendre integrating formula.
      SUBROUTINE VNODEQUAN(VD0,NV,VD,WD,FBAV,NV0)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp) :: vd0(nvt),vd(nvt),wd(nvt),FBAV(nvt)
      real(sp) :: VNODE(NGLT),WGAUL(NGLT)
      real(sp) :: AGL(NGLT),HGL(NGLT),VBMA
      integer:: IIV,NV0,NV,IV ,IVI
      call GAULEG(AGL,HGL)
!c---- Discrete Vx-velocity comparement:
      IIV=0

      NV0=NV
      DO iv=1,NV0
         call DISVEL(vd0(iv),vd0(iv+1),AGL,HGL,VNODE,WGAUL,VBMA)
         FBAV(iv)=VBMA
         DO IVI=1,KGLN
!c      write(*,*)IIV+IVI,IVI,NV0,KGLN
            VD(IIV+IVI)=VNODE(IVI)
            WD(IIV+IVI)=WGAUL(IVI)
         END DO
         IIV=IIV+KGLN
      END DO
      NV=IIV
      RETURN
      END
!C--The nodes and weights of Gauss-Legendre integration in space:-1 to +1.
      SUBROUTINE DISVEL(VA,VB,AGL,HGL,VNODE,WGAUL,VBMA)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
!      PARAMETER(NGLT=10)
      real(sp) :: AGL(NGLT),HGL(NGLT),VNODE(NGLT),WGAUL(NGLT),VA,VB,VAPB,VBMA

      VAPB=(VA+VB)/2.
      VBMA=(VB-VA)/2.

      IF (KGLN.EQ.3) THEN
         VNODE(1)=VAPB+VBMA*(-AGL(2))
         WGAUL(1)=HGL(2)
         VNODE(2)=VAPB
         WGAUL(2)=HGL(1)
         VNODE(3)=VAPB+VBMA*AGL(2)
         WGAUL(3)=HGL(2)

      ELSE IF (KGLN.EQ.4) THEN
         VNODE(1)=VAPB+VBMA*(-AGL(2))
         WGAUL(1)=HGL(2)
         VNODE(2)=VAPB+VBMA*(-AGL(1))
         WGAUL(2)=HGL(1)
         VNODE(3)=VAPB+VBMA*(AGL(1))
         WGAUL(3)=HGL(1)
         VNODE(4)=VAPB+VBMA*AGL(2)
         WGAUL(4)=HGL(2)

      ELSE IF (KGLN.EQ.5) THEN
         VNODE(3)=VAPB
         WGAUL(3)=HGL(1)
         VNODE(2)=VAPB+VBMA*(-AGL(2))
         WGAUL(2)=HGL(2)
         VNODE(4)=VAPB+VBMA*AGL(2)
         WGAUL(4)=HGL(2)
         VNODE(1)=VAPB+VBMA*(-AGL(3))
         WGAUL(1)=HGL(3)
         VNODE(5)=VAPB+VBMA*AGL(3)
         WGAUL(5)=HGL(3)

      ELSE IF (KGLN.EQ.6) THEN
         VNODE(1)=VAPB+VBMA*(-AGL(3))
         WGAUL(1)=HGL(3)
         VNODE(2)=VAPB+VBMA*(-AGL(2))
         WGAUL(2)=HGL(2)
         VNODE(3)=VAPB+VBMA*(-AGL(1))
         WGAUL(3)=HGL(1)
         VNODE(4)=VAPB+VBMA*AGL(1)
         WGAUL(4)=HGL(1)
         VNODE(5)=VAPB+VBMA*AGL(2)
         WGAUL(5)=HGL(2)
         VNODE(6)=VAPB+VBMA*AGL(3)
         WGAUL(6)=HGL(3)

      ELSE IF (KGLN.EQ.9) THEN
         VNODE(5)=VAPB
         WGAUL(5)=HGL(1)
         VNODE(4)=VAPB+VBMA*(-AGL(2))
         WGAUL(4)=HGL(2)
         VNODE(6)=VAPB+VBMA*AGL(2)
         WGAUL(6)=HGL(2)
         VNODE(3)=VAPB+VBMA*(-AGL(3))
         WGAUL(3)=HGL(3)
         VNODE(7)=VAPB+VBMA*AGL(3)
         WGAUL(7)=HGL(3)
         VNODE(2)=VAPB+VBMA*(-AGL(4))
         WGAUL(2)=HGL(4)
         VNODE(8)=VAPB+VBMA*AGL(4)
         WGAUL(8)=HGL(4)
         VNODE(1)=VAPB+VBMA*(-AGL(5))
         WGAUL(1)=HGL(5)
         VNODE(9)=VAPB+VBMA*AGL(5)
         WGAUL(9)=HGL(5)

      ELSE IF (KGLN.EQ.10) THEN
         VNODE(1)=VAPB+VBMA*(-AGL(5))
         WGAUL(1)=HGL(5)
         VNODE(2)=VAPB+VBMA*(-AGL(4))
         WGAUL(2)=HGL(4)
         VNODE(3)=VAPB+VBMA*(-AGL(3))
         WGAUL(3)=HGL(3)
         VNODE(4)=VAPB+VBMA*(-AGL(2))
         WGAUL(4)=HGL(2)
         VNODE(5)=VAPB+VBMA*(-AGL(1))
         WGAUL(5)=HGL(1)
         VNODE(6)=VAPB+VBMA*AGL(1)
         WGAUL(6)=HGL(1)
         VNODE(7)=VAPB+VBMA*AGL(2)
         WGAUL(7)=HGL(2)
         VNODE(8)=VAPB+VBMA*AGL(3)
         WGAUL(8)=HGL(3)
         VNODE(9)=VAPB+VBMA*AGL(4)
         WGAUL(9)=HGL(4)
         VNODE(10)=VAPB+VBMA*AGL(5)
         WGAUL(10)=HGL(5)
      END IF
      RETURN
      END

      SUBROUTINE GAULEG(AGL,HGL)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      real(sp) :: AGL(NGLT),HGL(NGLT)

      IF (KGLN.EQ.3) THEN
         AGL(1)=0.0
         HGL(1)=0.8888888889
         AGL(2)=0.7745966692
         HGL(2)=0.5555555556
      ELSE IF (KGLN.EQ.4) THEN
         AGL(1)=0.3399810436
         HGL(1)=0.6521451549
         AGL(2)=0.8611363116
         HGL(2)=0.3478548451
      ELSE IF (KGLN.EQ.5) THEN
         AGL(1)=0.0
         HGL(1)=0.5688888889
         AGL(2)=0.5384693101
         HGL(2)=0.4786286705
         AGL(3)=0.9061798459
         HGL(3)=0.2369268851
      ELSE IF (KGLN.EQ.6) THEN
         AGL(1)=0.2386191861
         HGL(1)=0.4679139346
         AGL(2)=0.6612093865
         HGL(2)=0.3607615730
         AGL(3)=0.9324695142
         HGL(3)=0.1713244924
      ELSE IF (KGLN.EQ.9) THEN
         AGL(1)=0.0
         HGL(1)=0.3302393550
         AGL(2)=0.3242534234
         HGL(2)=0.3123470770
         AGL(3)=0.6133714327
         HGL(3)=0.2606106964
         AGL(4)=0.8360311073
         HGL(4)=0.1806481607
         AGL(5)=0.9681602395
         HGL(5)=0.0812743884
      ELSE IF (KGLN.EQ.10) THEN
         AGL(1)=0.1488743390
         HGL(1)=0.2955242247
         AGL(2)=0.4333953941
         HGL(2)=0.2692667193
         AGL(3)=0.6794095683
         HGL(3)=0.2190863625
         AGL(4)=0.8650633667
         HGL(4)=0.1494513492
         AGL(5)=0.9739065285
         HGL(5)=0.0666713443
      END IF
      RETURN
      END






!C---  A new table for a modified (half-range) Gauss-Hermite duadrature
!C       with an evaluation of the integral.
!c           Roots(a(j)--VXD(IV)) and Weights(H(j)--WD(IV)) of modified
!c         Gauss-Hermite quadrature are stored in two arrays:
      SUBROUTINE GaussHermite(AGH,HGH,NV)
       use const_var ,only :sp,dop
      implicit none
      real(sp)  ::AGH(NV),HGH(NV)
      integer ::NV
      IF (NV.EQ.7) THEN
         AGH(1)=0.637163079E-1
         HGH(1)=0.160614533
         AGH(2)=0.318193583
         HGH(2)=0.306303720
         AGH(3)=0.724183347
         HGH(3)=0.275545914
         AGH(4)=1.238085494
         HGH(4)=0.120624443
         AGH(5)=1.838440401
         HGH(5)=0.218895217E-1
         AGH(6)=2.531553434
         HGH(6)=0.123780543E-2
         AGH(7)=3.373437976
         HGH(7)=0.109889525E-4 
      ELSE IF (NV.EQ.8) THEN
         AGH(1)=0.5297864393185113E-1
         HGH(1)=0.1341091884533595
         AGH(2)=0.2673983721677653
         HGH(2)=0.2683307544726388
         AGH(3)=0.6163028841823999
         HGH(3)=0.2759533979884218
         AGH(4)=1.064246312116224
         HGH(4)=0.1574482826187903
         AGH(5)=1.588855862270055
         HGH(5)=0.4481410991746290E-1
         AGH(6)=2.183921153095858
         HGH(6)=0.5367935756025333E-2
         AGH(7)=2.863133883708075
         HGH(7)=0.2020636491324107E-3
         AGH(8)=3.686007162724397
         HGH(8)=0.1192596926595344E-5
      ELSE IF (NV.EQ.16) THEN
         AGH(1)=0.1975365846007727E-1
         HGH(1)=0.5052463202137790E-1
         AGH(2)=0.1028022452379175
         HGH(2)=0.113608556894151
         AGH(3)=0.2473976694524551
         HGH(3)=0.1629212923145150
         AGH(4)=0.4466962259616832
         HGH(4)=0.1835628011162462
         AGH(5)=0.6930737203019995
         HGH(5)=0.1654386377556098
         AGH(6)=0.9794041703307299
         HGH(6)=0.1165724905535033
         AGH(7)=1.299789321277036
         HGH(7)=0.6199969609915657E-1
         AGH(8)=1.649854240397434
         HGH(8)=0.2391970961868355E-1
         AGH(9)=2.026808152168867
         HGH(9)=0.6409914424050132E-2
         AGH(10)=2.429450491602143
         HGH(10)=0.1135695310688778E-2
         AGH(11)=2.858266528543266
         HGH(11)=0.1252862213295624E-3
         AGH(12)=3.315769275038698
         HGH(12)=0.7950495719622457E-5
         AGH(13)=3.807377116755898
         HGH(13)=0.2590007619415064E-6
         AGH(14)=4.343606345470172
         HGH(14)=0.3611549139742782E-8
         AGH(15)=4.946377204048386
         HGH(15)=0.1537677916189839E-10
         AGH(16)=5.675017934041922
         HGH(16)=0.8674204452494624E-14
      END IF
      RETURN
      END 
