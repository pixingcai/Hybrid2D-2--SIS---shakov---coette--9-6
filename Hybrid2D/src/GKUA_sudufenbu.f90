!---- This routine "sudufenbu" set discrete velocity points using the discrete velocity
!      coordinate method,and check reliability of the the discrete velocity points
!       by the integrating-obtained macro-properties using Maxwellian velocity distribution.
!     Notice:in the following arrays:vd(nvt),wd(nvt),nvt=max(NVIt,NVJt).
    SUBROUTINE GKUA_sudufenbu
    use Global_var
    use Com_ctrl
    use const_var ,only :sp,dop
    implicit none
    integer :: m,iv,iv1,jv,jv1,i,j,NVIN,NVJN,index_1(4),index_2(4),index_3(4),index_4(4)
    integer :: ivp
    real(sp)  :: Tij,crt,viu,vjv,vij2,GDVijv,err111(4)
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
    MP=>Mesh(1)   
!---- Set discrete velocity points using the discrete velocity coordinate method:
1 call GKUA_lisansudu

!!---- To check reliability of the discrete velocity points,
!!---- Set Maxwellian velocity distribution and Store intial macro-properties:
       
    do m=1,MP%Num_Block
        B => MP%Block(m)
        if(B%solver==GKUA)then 
            do j=-1,B%ny+1
            do i=-1,B%nx+1      
                crt=B%S_GKUA(1,i,j)/(PI*B%S_GKUA(4,i,j))  
                do jv=GKUA_JST,GKUA_JEND          !!!  (do jv=1,nvj)
                   jv1=jv+GKUA_JPT
                   vjv=B%vj(jv1)-B%S_GKUA(3,i,j)
                do iv=GKUA_IST,GKUA_IEND          !!!  (do iv=1,nvi)
                    iv1=iv+GKUA_IPT
                    viu=B%vi(iv1)-B%S_GKUA(2,i,j)
                    vij2=(viu**2+vjv**2)/B%S_GKUA(4,i,j)
                    GDVijv=crt*EXP(-vij2)
                    B%FRD(1,i,j,iv,jv)=GDVijv
                    B%FRD(2,i,j,iv,jv)=B%S_GKUA( 4,i,j)*GDVijv/2.
                      !   B%FRD(3,i,j,iv,jv)=Rdof*B%S_GKUA( 5,i,j)*GDVijv/2.
                enddo
                enddo
            enddo
            enddo
            call GKUA_Phy_boundary
        end if
        continue
    end do


!!---- Decide the macro-properties related to the Maxwellian velocity distribution:
    if(myid==0) then
        write(*,*)'By initial "sudufenbu", cal. macro-properties-"detfp"'
    endif
    call DETFP(Time_Method)
    if(myid==0) then
        write(*,*) 'Initial vel. distribution ERR--<eps-"enter"=',err,eps
    endif
      CONTINUE
      IF (ERR .GT. EPS) THEN
         if(myid==0) then
            WRITE(*,*)'To integrate-solve macro properties by the Discrete velocity coordinate method,'
            WRITE(*,*) 'THE DISCRETE VELOCITY ZONE IS TOO SMALL!'
            WRITE(*,*)'Please increase discrete zone or decrease DVX,DVY!'
         endif
         NVIN=(VXup-vxdown+1.)/DVX
         NVJN=(VYup-vydown+1.)/DVY
         IF (NVIN.GT.NVIT .and. NVJN.LT.NVJT) THEN
            if(myid==0) then
               write(*,*)'NVI>=NVIT:stop increasing X-direction vel. zone!'
            endif
            NVJ=NVJN
            Vyup=Vyup+0.5
            Vydown=Vydown-0.5
            GOTO 1
         ELSE IF (NVJN.GT.NVJT .and. NVIN.LT.NVIT) THEN
            if(myid==0) then
               write(*,*)'NVJ>=NVJT:pause increase Y-direction vel. zone!'
            endif
            NVI=NVIN
            VXup=VXup+0.5
            VXdown=VXdown-0.5
            GOTO 1
         ELSE IF (NVIN.GT.NVIT .and. NVJN.GT.NVJT) THEN
            if(myid==0) then
               write(*,*)'Although the precision of integrating-solving macro properties isn"t satisfied,' 
               write(*,*)'But NVI>=NVIT,NVJ>=NVJT:had to stop increasing X- and Y- direction discrete velocity zone!'
            endif
         ELSE
            if(myid==0) then
               WRITE(*,*)'We use Vxspace=Vxspace+0.5,Vyspace=Vyspace+0.5 as one iteration-increasing.'
            endif
            NVI=NVIN
            VXup=VXup+0.5
            VXdown=VXdown-0.5
            NVJ=NVJN
            Vyup=Vyup+0.5
            Vydown=Vydown-0.5
            GOTO 1
         END IF 
      ELSE
         if(myid==0) then
            write(*,*) 'By now,the precision of the discrete velocity points and weights has been satisfied!'
         endif
      END IF  
      if(myid==0) then
         WRITE(*,*)'At last,NVI,VXup,VXdown',NVI,VXup,VXdown
         WRITE(*,*)'        NVJ,VYup,VYdown',NVJ,VYup,VYdown
      endif

      continue
   
      RETURN
      END  SUBROUTINE GKUA_sudufenbu       

