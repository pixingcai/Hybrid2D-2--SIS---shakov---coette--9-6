
      
      
      SUBROUTINE GKUA_Residual(nMesh,mBlock)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Integer :: nMesh,mBlock
      real(sp)  :: RHSERROR,DFRDI,DFRDJ
      real(sp)  :: AFRD(NumABFRD*2),BFRD(NumABFRD*2)
      real(sp)  :: alfa,ATSVXP,ATSVYP,ATSVXN,ATSVYN,ATSVXP1,ATSVYP1,ATSVXN1,&
          & ATSVYN1,FRDtmp,crt,vjv,viu,vij2,GDVijv
      integer::plane,iv,jv,iv1,jv1,i,j,nx,ny,ird,m
    
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
       
      AFRD=0.;BFRD=0.
      RHSERROR=-1.E20

      call Get_Feq (mBlock)          
      dt_global_GKUA=dt_global
      DO jv=GKUA_JST,GKUA_JEND
         jv1=jv+GKUA_JPT
      DO iv=GKUA_IST,GKUA_IEND
         iv1=iv+GKUA_IPT
         call RHSFDVsub(mBlock,iv,jv,iv1,jv1)
      
          do plane=2,nx-1+ny-1
            do j=1,ny-1
               i=plane-j
               if( i .LT. 1 .or. i .GT. nx-1) cycle    ! ���������ƽ��
               !write(*,*)"plane,i",plane,i
               ATSVXP=(B%ATSVX(i+1,j)+abs(B%ATSVX(i+1,j)))*0.5
               ATSVYP=(B%ATSVY(i,j+1)+abs(B%ATSVY(i,j+1)))*0.5
               ATSVXN=(B%ATSVX(i,j)-abs(B%ATSVX(i,j)))*0.5
               ATSVYN=(B%ATSVY(i,j)-abs(B%ATSVY(i,j)))*0.5
               alfa=B%vol(i,j)/dt_global_GKUA +(ATSVXP-ATSVXN+ATSVYP-ATSVYN)+B%vol(i,j)*B%CF(i,j)      ! �Խ�����  
                  
               ATSVXP1=(B%ATSVX(i,j)+abs(B%ATSVX(i,j)))*0.5  !ATSVX�����Ѿ�������dS(I-1/2)
               ATSVYP1=(B%ATSVY(i,j)+abs(B%ATSVY(i,j)))*0.5
               ATSVXN1=(B%ATSVX(i+1,j)-abs(B%ATSVX(i+1,j)))*0.5
               ATSVYN1=(B%ATSVY(i,j+1)-abs(B%ATSVY(i,j+1)))*0.5
               do ird=1,NRD
            
                     DFRDI=B%FRD00(ird,i-1,j,iv,jv) 
                     DFRDI=ATSVXP1*DFRDI
          
                     DFRDJ=B%FRD00(ird,i,j-1,iv,jv) 
                     DFRDJ=ATSVYP1*DFRDJ
          
                  B%FRD00(ird,i,j,iv,jv)=(B%vol(i,j)*B%RHSFRD(ird,i,j)+DFRDI+DFRDJ)/alfa !*
               
          
               enddo
            enddo
            enddo
            
            do plane=nx-1+ny-1,2,-1
            do j=ny-1,1,-1
               i=plane-j
               if( i .LT. 1 .or. i .GT. nx-1) cycle    ! ���������ƽ��
               !write(*,*)"plane,i",plane,i
               ATSVXP=(B%ATSVX(i+1,j)+abs(B%ATSVX(i+1,j)))*0.5
               ATSVYP=(B%ATSVY(i,j+1)+abs(B%ATSVY(i,j+1)))*0.5
               ATSVXN=(B%ATSVX(i,j)-abs(B%ATSVX(i,j)))*0.5
               ATSVYN=(B%ATSVY(i,j)-abs(B%ATSVY(i,j)))*0.5
               alfa=B%vol(i,j)/dt_global_GKUA +(ATSVXP-ATSVXN+ATSVYP-ATSVYN)+B%vol(i,j)*B%CF(i,j)      ! �Խ�����  
                  
               ATSVXP1=(B%ATSVX(i,j)+abs(B%ATSVX(i,j)))*0.5  !ATSVX�����Ѿ�������dS(I-1/2)
               ATSVYP1=(B%ATSVY(i,j)+abs(B%ATSVY(i,j)))*0.5
               ATSVXN1=(B%ATSVX(i+1,j)-abs(B%ATSVX(i+1,j)))*0.5
               ATSVYN1=(B%ATSVY(i,j+1)-abs(B%ATSVY(i,j+1)))*0.5
               do ird=1,NRD
          
               
                     DFRDI=-ATSVXN1*B%FRD00(ird,i+1,j,iv,jv)
          
                     DFRDJ=-ATSVYN1*B%FRD00(ird,i,j+1,iv,jv)
          
                  B%DUFRD(ird,i,j)=B%FRD00(ird,i,j,iv,jv)+(DFRDI+DFRDJ)/alfa
                
               enddo
            enddo
            enddo
            do j=1,ny-1
            do i=1,nx-1
               do ird=1,NRD
                  B%FRD(ird,i,j,iv,jv)=B%FRD(ird,i,j,iv,jv)+B%DUFRD(ird,i,j)
                 !  B%FRD(ird,i,j,iv,jv)=B%FRD(ird,i,j,iv,jv)+B%RHSFRD(ird,i,j)*dt_global_GKUA
                  B%FRD00(ird,i,j,iv,jv)= B%RHSFRD(ird,i,j)
                  
               enddo
              end do
            enddo
            do j=-1,ny+1
            do i=-1,nx+1
               do ird=1,NRD
                  IF(B%FRD(ird,i,j,iv,jv) .lt. 1.E-30) THEN
                     B%FRD(ird,i,j,iv,jv)=1.E-30
                  ENDIF
               enddo
            enddo
            enddo
           
      ENDDO
      ENDDO

    END SUBROUTINE GKUA_Residual
    

      SUBROUTINE RHSFDVsub(nm,iv,jv,iv1,jv1)
      use Global_var
      use Com_ctrl
      use const_var 
      IMPLICIT NONE
      integer:: nx,ny,ibegin,iend,jbegin,jend
      real(sp) ,external::FMINMD,FVALMTer
      real(sp) ,external :: GDVN 
      real(sp) :: FRDNES(3),RPitmp,HDVN,taotmp,feq1,heq1,vij2,GDVijv
      
      real(sp)  :: ASP,ASN,uvrb00,RHSERROR,RHSERROR0
      real(sp)  :: cpt,tmp,crt,viu,vjv,vkw,vijk2,cfijk,FDVM,FVQ,FDVNijk
      integer :: i,j,iv,iv1,jv,jv1,ITIME,ntyp,icount,ird
      integer :: nm,nn,drct,nega,posi,num,srtag,ibd,iad,bdnum
      real(sp)  :: DFRD(NRD,-nad:nmt),FluxFRD(NRD,-nad:nmt,-nad:nmt)
      real(sp)  :: SMUSCL,ULG,URG,ULH,URH,FVALMT,EPSMUSCL
      real(sp)  :: EPSWENO,Cweno1,Cweno2,IS1,IS2,IS3,AWENO1,AWENO2,AWENO,OWENO1,OWENO2
      real(sp)  :: GDVtmp,HDVtmp
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       nx=B%nx; ny=B%ny
 
              
      EPSMUSCL=1.E-30;EPSWENO=1.E-15
      Cweno1= 1./3.; Cweno2= 2./3.
      do j=1,ny-1
      do i=0,nx+1 
         do ird=1,NRD
            DFRD(ird,i)=B%FRD(ird,i,j,iv,jv)-B%FRD(ird,i-1,j,iv,jv)
         enddo
      enddo
      do i=1,nx
         B%ATSVX(i,j)=(B%vi(iv1)*B%ni1(i,j)+B%vj(jv1)*B%ni2(i,j))*B%si(i,j)
         ASP=(B%ATSVX(i,j)+ABS(B%ATSVX(i,j)))*0.5
         ASN=(B%ATSVX(i,j)-ABS(B%ATSVX(i,j)))*0.5
         IF(Iflag_scheme.EQ.GKUA_Scheme_UP1)THEN
            do ird=1,NRD
               FluxFRD(ird,i,j)=ASP*B%FRD(ird,i-1,j,iv,jv)+ASN*B%FRD(ird,i,j,iv,jv)
            enddo
         ENDIF
         IF(Iflag_scheme.EQ.GKUA_Scheme_NND4A2)THEN
            do ird=1,NRD
               ULG=B%FRD(ird,i-1,j,iv,jv)+0.5*FMINMD(DFRD(ird,i),DFRD(ird,i-1))
               URG=B%FRD(ird,i  ,j,iv,jv)-0.5*FMINMD(DFRD(ird,i),DFRD(ird,i+1))
               FluxFRD(ird,i,j)=ASP*ULG+ASN*URG
            enddo
         ENDIF
         IF(Iflag_scheme.EQ.GKUA_Scheme_WENO3)THEN
            do ird=1,NRD
              IS1= DFRD(ird,i-1)*DFRD(ird,i-1); IS2= DFRD(ird,i)*DFRD(ird,i)
              IS3= DFRD(ird,i+1)*DFRD(ird,i+1) 
              !Cweno1= 1./3.; Cweno2= 2./3.
              Aweno1= Cweno1/(IS1+EPSWENO)**2;Aweno2= Cweno2/(IS2+EPSWENO)**2
              Aweno = Aweno1+Aweno2
              Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
              ULH=Oweno1*(-0.5*B%FRD(ird,i-2,j,iv,jv)+1.5*B%FRD(ird,i-1,j,iv,jv))&
                  &+Oweno2*(0.5*B%FRD(ird,i-1,j,iv,jv)+0.5*B%FRD(ird,i,j,iv,jv))
              !Cweno1= 2./3.; Cweno2= 1./3.
              Aweno1= Cweno2/(IS2+EPSWENO)**2;Aweno2= Cweno1/(IS3+EPSWENO)**2
              Aweno = Aweno1+Aweno2
              Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
              URH=Oweno1*(0.5*B%FRD(ird,i-1,j,iv,jv)+0.5*B%FRD(ird,i,j,iv,jv))+&
                  &Oweno2*(1.5*B%FRD(ird,i,j,iv,jv)-0.5*B%FRD(ird,i+1,j,iv,jv))
              FluxFRD(ird,i,j)=ASP*ULH+ASN*URH
            enddo
         ENDIF
      enddo
      enddo
      
     do nn=1,B%subface
            Bc=> B%bc_msg(nn)
            if(Bc%neighb .lt. 0 ) then               ! ���ڱ߽�
               if( ((Bc%neighb .eq. BC_Wall).or.(Bc%neighb .eq. BC_Farfield)).and.(Bc%face .eq. 1) ) then
                    ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend
                    do j=jbegin,jend-1
                 
                    ASP=(B%ATSVX(2,j)+ABS(B%ATSVX(2,j)))*0.5
                    ASN=(B%ATSVX(2,j)-ABS(B%ATSVX(2,j)))*0.5
                  if(Iflag_scheme.EQ.GKUA_Scheme_UP1)then
                     do ird=1,NRD
                        FluxFRD(ird,2,j)=ASP*B%FRD(ird,1,j,iv,jv)+ASN*B%FRD(ird,2,j,iv,jv)
                     enddo
                  elseif(Iflag_scheme.EQ.GKUA_Scheme_WENO3)then
                     do ird=1,NRD
                       IS2= (B%FRD(ird,2,j,iv,jv)-B%FRD(ird,1,j,iv,jv))**2 
                       IS3= (B%FRD(ird,3,j,iv,jv)-B%FRD(ird,2,j,iv,jv))**2 
                       ULH=0.5*(B%FRD(ird,1,j,iv,jv)+B%FRD(ird,2,j,iv,jv))
                       Aweno1= Cweno2/(IS2+EPSWENO)**2;Aweno2= Cweno1/(IS3+EPSWENO)**2
                       Aweno = Aweno1+Aweno2
                       Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
                       URH=Oweno1*(0.5*B%FRD(ird,1,j,iv,jv)+0.5*B%FRD(ird,2,j,iv,jv))+&
                           &Oweno2*(1.5*B%FRD(ird,2,j,iv,jv)-0.5*B%FRD(ird,3,j,iv,jv))
                       FluxFRD(ird,2,j)=ASP*ULH+ASN*URH
                     enddo
                  else
                     do ird=1,NRD
                        ULG=B%FRD(ird,1,j,iv,jv)+0.5*(B%FRD(ird,2,j,iv,jv)-B%FRD(ird,1,j,iv,jv))
                        URG=B%FRD(ird,2,j,iv,jv)-0.5*FMINMD(B%FRD(ird,2,j,iv,jv)-&
                            &B%FRD(ird,1,j,iv,jv),B%FRD(ird,3,j,iv,jv)-B%FRD(ird,2,j,iv,jv))
                        FluxFRD(ird,2,j)=ASP*ULG+ASN*URG
                     enddo
                  endif
                  do ird=1,NRD
                     FluxFRD(ird,1,j)=B%ATSVX(1,j)*B%FRD(ird,-3,j,iv,jv)
                  enddo
                  end do
                    
               endif
            end if
               if(Bc%neighb .lt. 0 ) then 
               if( ((Bc%neighb .eq. BC_Wall).or.(Bc%neighb .eq. BC_Farfield)).and.(Bc%face .eq. 3) ) then
                    ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend
                    do j=jbegin,jend-1
                           ASP=(B%ATSVX(nx-1,j)+ABS(B%ATSVX(nx-1,j)))*0.5
                           ASN=(B%ATSVX(nx-1,j)-ABS(B%ATSVX(nx-1,j)))*0.5
                          if(Iflag_scheme.EQ.GKUA_Scheme_UP1)then
                             do ird=1,NRD
                                 FluxFRD(ird,nx-1,j)=ASP*B%FRD(ird,nx-2,j,iv,jv)+ASN*B%FRD(ird,nx-1,j,iv,jv)
                             enddo
                          elseif(Iflag_scheme.EQ.GKUA_Scheme_WENO3)then
                             do ird=1,NRD
                                   IS1= (B%FRD(ird,nx-2,j,iv,jv)-B%FRD(ird,nx-3,j,iv,jv))**2
                                   IS2= (B%FRD(ird,nx-1,j,iv,jv)-B%FRD(ird,nx-2,j,iv,jv))**2 
                                   Aweno1= Cweno1/(IS1+EPSWENO)**2;Aweno2= Cweno2/(IS2+EPSWENO)**2
                                   Aweno = Aweno1+Aweno2
                                   Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
                                   ULH=Oweno1*(-0.5*B%FRD(ird,nx-3,j,iv,jv)+1.5*B%FRD(ird,nx-2,j,iv,jv))+&
                                       & Oweno2*(0.5*B%FRD(ird,nx-2,j,iv,jv)+0.5*B%FRD(ird,nx-1,j,iv,jv))
                                   URH=0.5*B%FRD(ird,nx-2,j,iv,jv)+0.5*B%FRD(ird,nx-1,j,iv,jv)
                                   FluxFRD(ird,nx-1,j)=ASP*ULH+ASN*URH
                             enddo
                          else
                              do ird=1,NRD
                                    ULG=B%FRD(ird,nx-2,j,iv,jv)+0.5*FMINMD(B%FRD(ird,nx-1,j,iv,jv)-&
                                        & B%FRD(ird,nx-2,j,iv,jv),B%FRD(ird,nx-2,j,iv,jv)-B%FRD(ird,nx-3,j,iv,jv))
                                    URG=B%FRD(ird,nx-1,j,iv,jv)-0.5*(B%FRD(ird,nx-1,j,iv,jv)-B%FRD(ird,nx-2,j,iv,jv))
                                    FluxFRD(ird,nx-1,j)=ASP*ULG+ASN*URG
                             enddo
                          endif
                          
                          do ird=1,NRD
                             FluxFRD(ird,nx,j)=B%ATSVX(nx,j)*B%FRD(ird,nx+3,j,iv,jv)
                          enddo
                    end do
               endif
               endif
     enddo   
     
  

      do j=1,ny-1
      do i=1,nx-1
         do ird=1,NRD
           B%RHSFRD(ird,i,j)=FluxFRD(ird,i+1,j)-FluxFRD(ird,i,j)
         enddo
      enddo
      enddo
      
      do i=1,nx-1
      do j=0,ny+1 
         do ird=1,NRD
            DFRD(ird,j)=B%FRD(ird,i,j,iv,jv)-B%FRD(ird,i,j-1,iv,jv)
         enddo
      enddo
      do j=1,ny
         B%ATSVY(i,j)=(B%vi(iv1)*B%nj1(i,j)+B%vj(jv1)*B%nj2(i,j))*B%sj(i,j)
         ASP=(B%ATSVY(i,j)+ABS(B%ATSVY(i,j)))*0.5
         ASN=(B%ATSVY(i,j)-ABS(B%ATSVY(i,j)))*0.5
         IF(Iflag_scheme.EQ.GKUA_Scheme_UP1)THEN
            do ird=1,NRD
               FluxFRD(ird,i,j)=ASP*B%FRD(ird,i,j-1,iv,jv)+ASN*B%FRD(ird,i,j,iv,jv)
            enddo
         ENDIF
         IF(Iflag_scheme.EQ.GKUA_Scheme_NND4A2)THEN
            do ird=1,NRD
               ULG=B%FRD(ird,i,j-1,iv,jv)+0.5*FMINMD(DFRD(ird,j),DFRD(ird,j-1))
               URG=B%FRD(ird,i,j  ,iv,jv)-0.5*FMINMD(DFRD(ird,j),DFRD(ird,j+1))
               FluxFRD(ird,i,j)=ASP*ULG+ASN*URG
            enddo
         ENDIF
         IF(Iflag_scheme.EQ.GKUA_Scheme_WENO3)THEN
            do ird=1,NRD
               IS1= DFRD(ird,j-1)*DFRD(ird,j-1); IS2= DFRD(ird,j)*DFRD(ird,j) 
               IS3= DFRD(ird,j+1)*DFRD(ird,j+1) 
               !Cweno1= 1./3.; Cweno2= 2./3.
               Aweno1= Cweno1/(IS1+EPSWENO)**2;Aweno2= Cweno2/(IS2+EPSWENO)**2
               Aweno = Aweno1+Aweno2
               Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
               ULH=Oweno1*(-0.5*B%FRD(ird,i,j-2,iv,jv)+1.5*B%FRD(ird,i,j-1,iv,jv))+&
                   & Oweno2*(0.5*B%FRD(ird,i,j-1,iv,jv)+0.5*B%FRD(ird,i,j,iv,jv))
               !Cweno1= 2./3.; Cweno2= 1./3.
               Aweno1= Cweno2/(IS2+EPSWENO)**2;Aweno2= Cweno1/(IS3+EPSWENO)**2
               Aweno = Aweno1+Aweno2
               Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
               URH=Oweno1*(0.5*B%FRD(ird,i,j-1,iv,jv)+0.5*B%FRD(ird,i,j,iv,jv))+&
                   & Oweno2*(1.5*B%FRD(ird,i,j,iv,jv)-0.5*B%FRD(ird,i,j+1,iv,jv))
               FluxFRD(ird,i,j)=ASP*ULH+ASN*URH
            enddo
         ENDIF
      enddo
      enddo
      do nn=1,B%subface
            Bc=> B%bc_msg(nn)
            if(Bc%neighb .lt. 0 ) then            
            if( ((Bc%neighb .eq. BC_Wall).or.(Bc%neighb .eq. BC_Farfield)).and.(Bc%face .eq. 2) ) then
                ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend
                  do i=ibegin,iend-1
                       ASP=(B%ATSVY(i,2)+ABS(B%ATSVY(i,2)))*0.5
                       ASN=(B%ATSVY(i,2)-ABS(B%ATSVY(i,2)))*0.5 
                  IF(Iflag_scheme.EQ.GKUA_Scheme_UP1)THEN
                     do ird=1,NRD
                        FluxFRD(ird,i,2)=ASP*B%FRD(ird,i,1,iv,jv)+ASN*B%FRD(ird,i,2,iv,jv)
                     enddo
                  ELSEIF(Iflag_scheme.EQ.GKUA_Scheme_WENO3)THEN
                     do ird=1,NRD
                        IS2= (B%FRD(ird,i,2,iv,jv)-B%FRD(ird,i,1,iv,jv))**2 
                        IS3= (B%FRD(ird,i,3,iv,jv)-B%FRD(ird,i,2,iv,jv))**2 
                        ULH=0.5*B%FRD(ird,i,1,iv,jv)+0.5*B%FRD(ird,i,2,iv,jv)
                        Aweno1= Cweno2/(IS2+EPSWENO)**2;Aweno2= Cweno1/(IS3+EPSWENO)**2
                        Aweno = Aweno1+Aweno2
                        Oweno1= Aweno1/Aweno;Oweno2=Aweno2/Aweno
                        URH=Oweno1*(0.5*B%FRD(ird,i,1,iv,jv)+0.5*B%FRD(ird,i,2,iv,jv))+&
                            & Oweno2*(1.5*B%FRD(ird,i,2,iv,jv)-0.5*B%FRD(ird,i,3,iv,jv))
                        FluxFRD(ird,i,2)=ASP*ULH+ASN*URH
                     enddo
                  ELSE
                     do ird=1,NRD
                        ULG=B%FRD(ird,i,1,iv,jv)+0.5*(B%FRD(ird,i,2,iv,jv)-B%FRD(ird,i,1,iv,jv))
                        URG=B%FRD(ird,i,2,iv,jv)-0.5*FMINMD(B%FRD(ird,i,2,iv,jv)-B%FRD(ird,i,1,iv,jv),&
                            & B%FRD(ird,i,3,iv,jv)-B%FRD(ird,i,2,iv,jv))
                        FluxFRD(ird,i,2)=ASP*ULG+ASN*URG
                     enddo
                  ENDIF
                  do ird=1,NRD
                     FluxFRD(ird,i,1)=B%ATSVY(i,1)*B%FRD(ird,i,-3,iv,jv)
                  enddo
                  end do
                  
            end if
            if( ((Bc%neighb .eq. BC_Wall).or.(Bc%neighb .eq. BC_Farfield)).and.(Bc%face .eq. 4) ) then
                ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend
                  do i=ibegin,iend-1
                     ASP=(B%ATSVY(i,ny-1)+ABS(B%ATSVY(i,ny-1)))*0.5
                     ASN=(B%ATSVY(i,ny-1)-ABS(B%ATSVY(i,ny-1)))*0.5
      
                  IF(Iflag_scheme.EQ.GKUA_Scheme_UP1)THEN
                     do ird=1,NRD
                        FluxFRD(ird,i,ny-1)=ASP*B%FRD(ird,i,ny-2,iv,jv)+ASN*B%FRD(ird,i,ny-1,iv,jv)
                     enddo
                  ELSEIF(Iflag_scheme.EQ.GKUA_Scheme_WENO3)THEN
                     do ird=1,NRD
                        IS1= (B%FRD(ird,i,ny-2,iv,jv)-B%FRD(ird,i,ny-3,iv,jv))**2
                        IS2= (B%FRD(ird,i,ny-1,iv,jv)-B%FRD(ird,i,ny-2,iv,jv))**2 
                        Aweno1= Cweno1/(IS1+EPSWENO)**2;Aweno2= Cweno2/(IS2+EPSWENO)**2
                        Aweno = Aweno1+Aweno2
                        Oweno1=Aweno1/Aweno;Oweno2=Aweno2/Aweno
                        ULH=Oweno1*(-0.5*B%FRD(ird,i,ny-3,iv,jv)+1.5*B%FRD(ird,i,ny-2,iv,jv))+&
                            & Oweno2*(0.5*B%FRD(ird,i,ny-2,iv,jv)+0.5*B%FRD(ird,i,ny-1,iv,jv))
                        URH=0.5*B%FRD(ird,i,ny-2,iv,jv)+0.5*B%FRD(ird,i,ny-1,iv,jv)
                        FluxFRD(ird,i,ny-1)=ASP*ULH+ASN*URH
                     enddo
                  ELSE
                     do ird=1,NRD
                        ULG=B%FRD(ird,i,ny-2,iv,jv)+0.5*FMINMD(B%FRD(ird,i,ny-1,iv,jv)-B%FRD(ird,i,ny-2,iv,jv),&
                            & B%FRD(ird,i,ny-2,iv,jv)-B%FRD(ird,i,ny-3,iv,jv))
                        URG=B%FRD(ird,i,ny-1,iv,jv)-0.5*(B%FRD(ird,i,ny-1,iv,jv)-B%FRD(ird,i,ny-2,iv,jv))
                        FluxFRD(ird,i,ny-1)=ASP*ULG+ASN*URG
                     enddo
                  ENDIF
                  do ird=1,NRD
                     FluxFRD(ird,i,ny)=B%ATSVY(i,ny)*B%FRD(ird,i,ny+3,iv,jv)
                  enddo
                  end do
                  
            end if
            end if
       enddo    
         
      do i=1,nx-1
      do j=1,ny-1
         do ird=1,NRD
            B%RHSFRD(ird,i,j)=B%RHSFRD(ird,i,j)+FluxFRD(ird,i,j+1)-FluxFRD(ird,i,j)
         enddo
      enddo
      enddo
      !
      do j=1,ny-1
      do i=1,nx-1
          
        B%RHSFRD(1,i,j)=-B%RHSFRD(1,i,j)/B%vol(i,j)+B%CF(i,j)*(B%FRD0(1,i,j,iv,jv)-B%FRD(1,i,j,iv,jv)) 
        B%RHSFRD(2,i,j)=-B%RHSFRD(2,i,j)/B%vol(i,j)+B%CF(i,j)*(B%FRD0(2,i,j,iv,jv)-B%FRD(2,i,j,iv,jv))


      enddo
      enddo
      
     

    END SUBROUTINE RHSFDVsub
    
    subroutine Get_Feq(mBlock)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Integer :: nMesh,mBlock
      real(sp)  :: RHSERROR,DFRDI,DFRDJ,cpt,viuq,vq,FVXDV
      real(sp)  :: AFRD(NumABFRD*2),BFRD(NumABFRD*2)
      real(sp)  :: alfa,ATSVXP,ATSVYP,ATSVXN,ATSVYN,ATSVXP1,ATSVYP1,ATSVXN1,ATSVYN1,FRDtmp,crt,vjv,viu,vij2,GDVijv
      integer::plane,iv,jv,iv1,jv1,i,j,nx,ny,ird,m
    
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
    
          do m=1,MP%Num_Block
           B => MP%Block(m)
           if(B%solver==GKUA)then 
           do j=1,B%ny-1
              do i=1,B%nx-1      
                 DO jv=GKUA_JST,GKUA_JEND          !!!  (do jv=1,nvj)
                   jv1=jv+GKUA_JPT
                   vjv=B%vj(jv1)-B%S_GKUA(3,i,j)
                     DO iv=GKUA_IST,GKUA_IEND          !!!  (do iv=1,nvi)
                        iv1=iv+GKUA_IPT
                        viu=B%vi(iv1)-B%S_GKUA(2,i,j)
                        crt=B%S_GKUA(1,i,j)/(PI*B%S_GKUA(4,i,j)) 
                        cpt=cpr/(B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j)*B%S_GKUA(4,i,j))
                        viuq=viu*B%S_GKUA(8,i,j)
                        vq=viuq+vjv*B%S_GKUA(9,i,j)
                        
                        vij2=(viu**2+vjv**2)/B%S_GKUA(4,i,j)
                         FVXDV=cpt*Vq
                         GDVijv=crt*EXP(-vij2)
                         B%FRD0(1,i,j,iv,jv)=GDVijv*(1.0+FVXDV*(2.*vij2-4.))
                         B%FRD0(2,i,j,iv,jv)=(B%S_GKUA( 4,i,j)*GDVijv/2.)*(1.0+FVXDV*(2.*vij2-2.))
                        ENDDO
                       ENDDO

                     ENDDO 
           ENDDO 
         
           end if
           
          continue
      end do
    end  subroutine Get_Feq
    
