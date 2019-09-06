!    subroutine write_FRD
!       use Global_Var
!       use Com_ctrl
!       implicit none
!       TYPE (Mesh_TYPE),pointer:: MP
!       Type (Block_TYPE),pointer:: B
!       Type (BC_MSG_TYPE),pointer:: Bc
!       integer stats(MPI_STATUS_SIZE)
!       character*4  fno,fnum
!       character(len=50):: fname,filename_Dem
!       integer :: i,j, m,i_x,i_y,iv,jv,iv1,jv1,ids,k,num_coord,ICOUNT
!       real(sp) ::  coord_FRD(100,2),FRD_OUT(200,200,3),GRD_OUT(200,200,3),LRD_OUT(200,200,3),TMP_OUT(200,200,3)
!       real(sp) ::  send_real_buffer(5000),recv_real_buffer(5000),distance
!       integer coord_FRD_INDEX(100,3)
!       INTEGER mpi_err,send_int_buffer(20),recv_int_buffer(20),num_count(4,1000)
!        MP=>Mesh(1)
!       fname="input/"//"Coordinates.dat"
!        if(myid==0)then
!           write(*,*)"read coordinates"
!       end if
!       OPEN(UNIT=115,FILE=trim(fname),STATUS='OLD',share="denynone")
!       read(115,*)num_coord
!       do i=1,num_coord
!          read(115,*)coord_FRD(i,1),coord_FRD(i,2)
!      
!          distance=1E10
!          do m=1,MP%Num_Block
!           B => Mesh(1)%Block(m)
!            do i_x=1,B%nx
!                do i_y=1,B%ny
!                    
!                    if(sqrt((coord_FRD(i,1)-B%x(i_x,i_y))**2+(coord_FRD(i,2)-B%y(i_x,i_y))**2).LT.distance)then
!                        distance=sqrt((coord_FRD(i,1)-B%x(i_x,i_y))**2+(coord_FRD(i,2)-B%y(i_x,i_y))**2)
!                        coord_FRD_INDEX(I,1)=m
!                        coord_FRD_INDEX(I,2)=i_x
!                        coord_FRD_INDEX(I,3)=i_y
!                        
!                    end if
!                end do
!            end do
!          end do
!          
!          
!          B => Mesh(1)%Block(coord_FRD_INDEX(I,1))
!          if(B%solver==GKUA)then
!          DO jv=GKUA_JST,GKUA_JEND
!                    jv1=jv+GKUA_JPT
!                    DO iv=GKUA_IST,GKUA_IEND
!                    iv1=iv+GKUA_IPT
!                    FRD_OUT(iv1,jv1,1)=B%FRD(1,coord_FRD_INDEX(I,2),coord_FRD_INDEX(I,3),iv,jv)
!                    FRD_OUT(iv1,jv1,2)=B%FRD(2,coord_FRD_INDEX(I,2),coord_FRD_INDEX(I,3),iv,jv)  
!                    FRD_OUT(iv1,jv1,3)=B%FRD(3,coord_FRD_INDEX(I,2),coord_FRD_INDEX(I,3),iv,jv) 
!                    end do
!          end do
!           call MPI_Barrier(mpi_comm_world,ierr)
!            send_int_buffer(1)=(GKUA_IEND-GKUA_IST+1)
!            send_int_buffer(2)=(GKUA_JEND-GKUA_JST+1)
!            send_int_buffer(3)=GKUA_IPT
!            send_int_buffer(4)=GKUA_JPT
!      !      WRITE(*,*),IST,IEND,JST,JEND,IPT,JPT
!        if(myid .NE.0)then
!         call MPI_Send(send_int_buffer,4,mpi_integer,0,0,MPI_COMM_WORLD,mpi_err)
!        endif
!        if(myid ==0)then
!            do ids=1,nprocs-1
!                call MPI_Recv(recv_int_buffer,4,mpi_integer,ids,0,MPI_COMM_WORLD,stats,mpi_err)
!                 !   write(*,*) "id 0 finish recv",recv_int_buffer(1:4)
!                    num_count(1,ids)=recv_int_buffer(1)
!                    num_count(2,ids)=recv_int_buffer(2)
!                    num_count(3,ids)=recv_int_buffer(3)
!                    num_count(4,ids)=recv_int_buffer(4)
!            end do
!           ! WRITE(*,*)num_count(1:4,1:10)
!        endif
!        
!        call MPI_Barrier(mpi_comm_world,ierr)
!        DO k=1,3
!        if(myid .NE.0)then
!           DO jv=1,GKUA_JEND-GKUA_JST+1
!              jv1=jv+GKUA_JPT
!              DO iv=1,GKUA_IEND-GKUA_IST+1
!                 iv1=iv+GKUA_IPT
!                 send_real_buffer((jv-1)*(GKUA_IEND-GKUA_IST+1)+iv)=B%FRD(k,coord_FRD_INDEX(I,2),coord_FRD_INDEX(I,3),iv,jv)
!              end do 
!           end do
!           call MPI_Send(send_real_buffer,(GKUA_IEND-GKUA_IST+1)*(GKUA_JEND-GKUA_JST+1),mpi_real,0,myid,MPI_COMM_WORLD,mpi_err)
!        end if
!        if(myid ==0)then
!            do ids=1,nprocs-1
!                call MPI_Recv(recv_real_buffer,num_count(1,ids)*num_count(2,ids),mpi_real,ids,ids,MPI_COMM_WORLD,stats,mpi_err)
!                do jv=1,num_count(2,ids)
!                      jv1=jv+num_count(4,ids)
!                    do iv=1,num_count(1,ids)
!                        iv1=iv+num_count(3,ids)
!                      FRD_OUT(iv1,jv1,k)=recv_real_buffer((jv-1)*num_count(1,ids)+iv)
!                    end do
!                end do
!                
!            end do
!          
!        endif
!         END DO 
!        call MPI_Barrier(mpi_comm_world,ierr)
!        GRD_OUT=0;LRD_OUT=0;
!        call GET_CE_FRD(coord_FRD_INDEX(I,1),coord_FRD_INDEX(I,2),coord_FRD_INDEX(I,3),GRD_OUT,LRD_OUT)
!        icount = 200*200*3
!        call MPI_allreduce(GRD_OUT,TMP_OUT,icount,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
!        GRD_OUT=TMP_OUT
!         call MPI_Barrier(mpi_comm_world,ierr)
!        call MPI_allreduce(LRD_OUT,TMP_OUT,icount,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
!        LRD_OUT=TMP_OUT
!         call MPI_Barrier(mpi_comm_world,ierr)       
!    !     write(*,*)"i is :",i
!         WRITE(fnum,'(I3)'),i
!         IF(myid==0)THEN
!       !     write(*,*) "fnum is ",fnum
!            fname="output/"//"FRD_COORD"//fnum//".dat" 
!          !  write(*,*) "file name is ",fname
!            open(UNIT=112,file=trim(fname),FORM='FORMATTED',ACTION='READWRITE',STATUS='REPLACE')
!            WRITE(112,1121)
!            write(112,1122)  '"0"',NVIT,NVJT
!            DO iv=1,NVIT
!                DO jv=1,NVJT
!                     WRITE(112,1124)B%vi(iv),B%vj(jv),FRD_out(Iv,Jv,1),FRD_out(Iv,Jv,2),FRD_out(Iv,Jv,3),GRD_OUT(iv,jv,1),GRD_OUT(iv,jv,2),GRD_OUT(iv,jv,3) &
!                     &,LRD_OUT(iv,jv,1),LRD_OUT(iv,jv,2),LRD_OUT(iv,jv,3)
!                end do
!            end do
!                 CLOSE(112)
!              END IF
!
!      end if    
!      end do
!      
!      close(115)
!     
!1121  format(' VARIABLES="VX","VY","FRD1","FRD2","FRD3","GRD1","GRD2","GRD3","LRD1","LRD2","LRD3" ')
!1122  format('ZONE T=',A10,' ,I=',I6,', J=',I6,',   F=POINT')
!1124  format(11f15.7)   
!    end subroutine write_FRD
!    
    subroutine Get_CE_FRD(nBlock,i,j,FRD_OUT,FFRD_OUT)
           use Global_Var
           use Com_ctrl
           implicit none
           TYPE (Mesh_TYPE),pointer:: MP
           Type (Block_TYPE),pointer:: B
           integer :: i,j ,iv ,jv,iv1,jv1,nBlock
           real(sp) ::Ttraij,Ronij,Coetmp,Taoxx,Taoxy,Taoyy,Taozz,Qx,Qy,Crt,RPitmp,vjv,viu
           real(sp) :: corrtmp11,corrtmp12,corrtmp21,corrtmp22,U,V,VIJ2,GDVIJV,corrtmp31,corrtmp32
           real(sp) :: FRD_OUT(200,200,3),FFRD_OUT(200,200,3),taotmp
               MP=>Mesh(1)
               B => MP%Block(nBlock)
               FRD_OUT(200,200,3)=0.0;FFRD_OUT(200,200,3)=0.0
               !B%S_GKUA(1,tempi,tempj)=Ua(1,k1,i)    !!RHO
               !B%S_GKUA(2,tempi,tempj)= Ua(3,k1,i)   !!U
               !B%S_GKUA(3,tempi,tempj)=Ua(4,k1,i)    !!V
               !B%S_GKUA(4,tempi,tempj)=Ua(2,k1,i)    !!Ttran
               !B%S_GKUA(5,tempi,tempj)=Ua(2,k1,i)
                                    
               Ttraij=B%S_GKUA(4,i,j)
               RONIJ =B%S_GKUA(1,i,j)
               coetmp=(1.-1./Pr)/(1.-ZrotD) 
               taoxx=B%S_GKUA(7,i,j)-(B%S_GKUA(7,i,j)+B%S_GKUA(9,i,j)+B%S_GKUA(10,i,j))/3.0
               taoyy=B%S_GKUA(9,i,j)-(B%S_GKUA(7,i,j)+B%S_GKUA(9,i,j)+B%S_GKUA(10,i,j))/3.0
               taoxy=B%S_GKUA(8,i,j)
               taozz=B%S_GKUA(10,i,j)-(B%S_GKUA(7,i,j)+B%S_GKUA(9,i,j)+B%S_GKUA(10,i,j))/3.0
				U=B%S_GKUA(2,i,j)
                V=B%S_GKUA(3,i,j) 
                QX=B%S_GKUA(11,i,j)+B%S_GKUA(13,i,j)
                QY=B%S_GKUA(12,i,j)+B%S_GKUA(14,i,j)
				
     
               crt=B%S_GKUA(1,i,j)/(PI*B%S_GKUA(4,i,j)) 
                                    
                     
               DO jv=GKUA_JST,GKUA_JEND         
                   jv1=jv+GKUA_JPT
               DO iv=GKUA_IST,GKUA_IEND   
                   iv1=iv+GKUA_IPT
                   RPitmp= RONIJ/pi
                   vjv=B%vj(jv1)-V 
                   viu=B%vi(iv1)-U
                   Corrtmp11=2.0*(gamma-1.0)/(gamma*RONij*Ttraij**2)*(viu*Qx+vjv*Qy)*((viu**2+vjv**2)/Ttraij-2.0)
                   Corrtmp12=(taoxx*viu**2+taoyy*vjv**2+2.0*taoxy*vjv*viu+taozz*Ttraij/2.0)/(RONij*Ttraij**2)
                   Corrtmp11=(1.0+Corrtmp11+Corrtmp12)
                   Corrtmp21=2.0*(gamma-1.0)/(gamma*RONij*Ttraij**2)*(viu*Qx+vjv*Qy)*((viu**2+vjv**2)/Ttraij-1.0)
                   Corrtmp22=(taoxx*viu**2+taoyy*vjv**2+2.0*taoxy*vjv*viu+3.0*taozz*Ttraij/2.0)/(RONij*Ttraij**2)
                   Corrtmp21=(1.0+Corrtmp21+Corrtmp22)
                   
                   Corrtmp31=2.0*(gamma-1.0)/(gamma*RONij*Ttraij**2)*(viu*Qx+vjv*Qy)*((viu**2+vjv**2)/Ttraij-1.0)
                   Corrtmp32=(taoxx*viu**2+taoyy*vjv**2+2.0*taoxy*vjv*viu+taozz*Ttraij/2.0)/(RONij*Ttraij**2)
                   Corrtmp31=(1.0+Corrtmp31+Corrtmp32)
                                              
                   FRD_OUT(iv1,jv1,1)=RPitmp/Ttraij*exp(-(vjv**2+viu**2)/Ttraij)*Corrtmp11
                   FRD_OUT(iv1,jv1,2)=Ttraij/2.0*RPitmp/Ttraij*exp(-(vjv**2+viu**2)/Ttraij)*Corrtmp21
                   FRD_OUT(iv1,jv1,3)=Rdof*Ttraij/2.*RPitmp/Ttraij*exp(-(vjv**2+viu**2)/Ttraij)*Corrtmp31
                   
                   !  vij2=(viu**2+vjv**2)/B%S_GKUA(4,i,j)
                   !  GDVijv=crt*EXP(-vij2)
                   !FFRD_OUT(iv1,jv1,1)=GDVijv
                   !FFRD_OUT(iv1,jv1,2)=B%S_GKUA( 4,i,j)*GDVijv/2.
                   !FFRD_OUT(iv1,jv1,3)=Rdof*B%S_GKUA( 5,i,j)*GDVijv/2.  
                   
                   RPitmp= B%S_GKUA(1,i,j)/pi
                   viu=B%vi(iv1)- U;vjv=B%vj(jv1)- V 
                   taotmp=B%TaoXXES(i,j)*B%TaoYYES(i,j)-B%TaoXYES(i,j)*B%TaoXYES(i,j)
                   FFRD_OUT(iv1,jv1,1)=RPitmp/sqrt(taotmp)*exp(-(viu*viu*B%TaoYYES(i,j)-&
                       &2.*viu*vjv*B%TaoXYES(i,j)+vjv*vjv*B%TaoXXES(i,j))/taotmp)
                   FFRD_OUT(iv1,jv1,2)=B%TaoZZES(i,j)*0.5*FFRD_OUT(iv1,jv1,1)
                   FFRD_OUT(iv1,jv1,3)=Rdof*B%S_GKUA(18,i,j)/2.*FFRD_OUT(iv1,jv1,1)
                                                
               ENDDO
               ENDDO
    end subroutine Get_CE_FRD      
        
    
!---------------------------------------
!         打印残差（最大残差和均方根残差）
    subroutine output_Res_GKUA(nMesh)
    use Global_var
     use Com_ctrl
     use Flow_Var 
	implicit none
    integer::time_label(8)
	integer:: nMesh,i,block_check,i_max,j_max
    real(sp) ::Res_moutg(6),Res_routg(6),Res_mouts(6),Res_routs(6),time_tmp
    TYPE (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    MP=>Mesh(1) 
    Res_moutg=0;Res_routg=0;Res_mouts=0;Res_routs=0
    do i=1,MP%Num_Block
         B => MP%Block(i)
         Res_moutg(:)=max(Res_moutg(:),B%Res_max(:))
         Res_routg(:)=max(Res_routg(:),B%Res_rms(:))
    end do
    if(isnan(Res_routg(1)))then
        if(myid==0)then
            write(*,*)"NaN detected ...."
        end if
        call MPI_GROUP_FREE(new_group,IERR)
        call MPI_COMM_FREE(new_comm,IERR)
        call MPI_FINALIZE(IERR)
        stop
    end if
   
    if(myid==0)then
      print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
      print*, "  The R.M.S Residuals of GKUA are "
      write(*, "(4E20.10)") Res_routg(1:4)

      
      call DATE_AND_TIME(VALUES=time_label)
      time_tmp=((time_label(5)-initialtime(5))*3600+(time_label(6)-initialtime(6))*60+(time_label(7)-initialtime(7)))/3600.0
      open(99,file="Residual_gkua.dat",position="append")
      write(99,"(I9,5E24.5,I4)") Mesh(nMesh)%Kstep, time_tmp ,Res_routg(1:4),MP%iter
      close(99) 
    end if
    end subroutine output_Res_GKUA

    !---------------------------------------
!         打印残差（最大残差和均方根残差）
    subroutine output_Res_NS(nMesh)
    use Global_var
     use Com_ctrl
     use Flow_Var 
	implicit none
	integer:: nMesh,i,block_check,i_max,j_max
    real(sp) ::Res_moutg(6),Res_routg(6),Res_mouts(6),Res_routs(6)
    TYPE (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    MP=>Mesh(1) 
    Res_mouts=0;Res_routs=0
    do i=1,MP%Num_Block
         B => MP%Block(i)
         Res_mouts(:)=max(Res_mouts(:),B%Res_max(:))
         Res_routs(:)=max(Res_routs(:),B%Res_rms(:))  
 
    end do
   
    if(myid==0)then
      print*, "  The Max Residuals of NS are "
      write(*, "(4E20.10)") Res_mouts(1:4)
    end if
    end subroutine output_Res_NS
    
    subroutine get_res(res_iter)
    use Global_var
     use Com_ctrl
     use Flow_Var 
	implicit none
	integer:: nMesh,i,block_check,i_max,j_max
    real(sp) ::Res_moutg(6),Res_routg(6),Res_mouts(6),Res_routs(6),res_iter
    TYPE (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    MP=>Mesh(1) 
    Res_mouts=0;Res_routs=0
    do i=1,MP%Num_Block
         B => MP%Block(i)
         Res_routs(:)=max(Res_routs(:),B%Res_rms(:))  
         Res_mouts(:)=max(Res_mouts(:),B%Res_max(:))
    end do
    res_iter=Res_mouts(4)
    end subroutine get_res
!----------------------------------------------------------------------
!  输出几何及物理量 （tecplot格式）, 最细网格flow2d.dat; 粗网格 flow2d-2.dat ; 最粗网格 flow2d-3.dat
   
      subroutine output (nMesh)
       use Global_Var
       use Com_ctrl
       implicit none
       TYPE (Mesh_TYPE),pointer:: MP
       Type (Block_TYPE),pointer:: B
       Type (BC_MSG_TYPE),pointer:: Bc
       integer  nMesh,i,j,m,ksub,n1
        integer:: ibegin,iend,jbegin,jend,icheck,jcheck,pos_x,pos_y
       real(sp) :: d1,u1,v1,T1,txx,tyy,txy,tzz,qx,qy,errors,temp,V_tmp(9)
       character(len=50):: filename,filename_Dem
       
	   if(nMesh .eq. 1) then
          filename="flow2d"
          filename_Dem="flow2d_Dem"
	   else
          write(filename,"('flow2d-'I1.1'')") nMesh    
       endif

	   print*, "write data file ...", filename
       
	   MP=>Mesh(1)
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,d,u,v,T,p,Ma"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx-5, " j= ", B%ny-5
  	     do j=3,B%ny-3
	     do i=3,B%nx-3
            d1=B%U(1,i,j)
            u1=B%U(2,i,j)
            v1=B%U(3,i,j)
            temp=B%U(4,i,j)
            T1=(temp-0.5d0*d1*(u1*u1+v1*v1))/(Cv*d1)
      
             write(99,"(8f20.10)") B%x1(i,j),B%y1(i,j), d1,u1,v1,T1,d1*T1/(gamma*Ma*Ma),Ma*sqrt(u1**2+v1**2)/sqrt(T1)
         enddo
         enddo


        enddo
         do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         if((MODE==Hybrid_Mode).and.(B%solver==GKUA))then
             write(99,*)  " zone ", "i= ", B%nx, " j= ", B%ny
  	         do j=1,B%ny
	         do i=1,B%nx
               d1=(B%S_GKUA(1,i,j)+B%S_GKUA(1,i-1,j)+B%S_GKUA(1,i,j-1)+B%S_GKUA(1,i-1,j-1))*0.25d0
               u1=(B%S_GKUA(2,i,j)+B%S_GKUA(2,i-1,j)+B%S_GKUA(2,i,j-1)+B%S_GKUA(2,i-1,j-1))*0.25d0/Ma*sqrt(2.0/gamma)
               v1=(B%S_GKUA(3,i,j)+B%S_GKUA(3,i-1,j)+B%S_GKUA(3,i,j-1)+B%S_GKUA(3,i-1,j-1))*0.25d0/Ma*sqrt(2.0/gamma)
               temp=(B%S_GKUA(4,i,j)+B%S_GKUA(4,i-1,j)+B%S_GKUA(4,i,j-1)+B%S_GKUA(4,i-1,j-1))*0.25d0
               T1=temp
      
             write(99,"(8f20.10)") B%x(i,j),B%y(i,j), d1,u1,v1,T1,d1*T1/(gamma*Ma*Ma),Ma*sqrt(u1**2+v1**2)/sqrt(T1)
             enddo
             enddo
         end if


        enddo
        
       close(99)
       filename="flow2d-GKUA"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,d,u,v,T,p,txx,tyy,txy,qx,qy,tzz"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         if(B%solver==GKUA)then
         write(99,*)  " zone ", "i= ", B%nx-1, " j= ", B%ny-1
  	     do j=1,B%ny-1
	     do i=1,B%nx-1
           d1=B%S_GKUA(1,i,j);u1=B%S_GKUA(2,i,j);v1=B%S_GKUA(3,i,j)
           T1=B%S_GKUA(4,i,j);txx=B%S_GKUA(5,i,j)-d1*T1;tyy=B%S_GKUA(6,i,j)-d1*T1
           txy=B%S_GKUA(7,i,j);qx=B%S_GKUA(8,i,j);qy=B%S_GKUA(9,i,j)
           tzz=B%S_GKUA(10,i,j)-d1*T1
           temp=-(qx*Pr+v1*txy)*B%CF(i,j)
         write(99,"(13f15.10)") B%x1(i,j),B%y1(i,j), d1,u1,v1,T1,d1*T1,txx,tyy,txy,qx,qy,tzz
    
         enddo
         enddo
         end if
        enddo
       close(99)
       
       filename="flow2d-S0"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,d,u,v,T,p,txx,tyy,txy,qx,qy,tzz"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         if(B%solver==GKUA)then
         write(99,*)  " zone ", "i= ", B%nx-5, " j= ", B%ny-5
  	     do j=3,B%ny-3
	     do i=3,B%nx-3
           d1=B%S0_GKUA(1,i,j);u1=B%S0_GKUA(2,i,j);v1=B%S0_GKUA(3,i,j)
           T1=B%S0_GKUA(4,i,j);txx=B%S0_GKUA(5,i,j);tyy=B%S0_GKUA(6,i,j)
           txy=B%S0_GKUA(7,i,j);qx=B%S0_GKUA(8,i,j);qy=B%S0_GKUA(9,i,j)
           tzz=B%S0_GKUA(10,i,j)
         write(99,"(13f15.10)") B%x1(i,j),B%y1(i,j), d1,u1,v1,T1,d1*T1,txx,tyy,txy,qx,qy,tzz
    
         enddo
         enddo
         end if
        enddo
       close(99)
       
       filename="flow2d-S00"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,d,u,v,T,p,txx,tyy,txy,qx,qy,tzz"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         if(B%solver==GKUA)then
         write(99,*)  " zone ", "i= ", B%nx-5, " j= ", B%ny-5
  	     do j=3,B%ny-3
	     do i=3,B%nx-3
           d1=B%S00_GKUA(1,i,j);u1=B%S00_GKUA(2,i,j);v1=B%S00_GKUA(3,i,j)
           T1=B%S00_GKUA(4,i,j);txx=B%S00_GKUA(5,i,j);tyy=B%S00_GKUA(6,i,j)
           txy=B%S00_GKUA(7,i,j);qx=B%S00_GKUA(8,i,j);qy=B%S00_GKUA(9,i,j)
           tzz=B%S00_GKUA(10,i,j)
         write(99,"(13f15.10)") B%x1(i,j),B%y1(i,j), d1,u1,v1,T1,d1*T1,txx,tyy,txy,qx,qy,tzz
    
         enddo
         enddo
         end if
        enddo
       close(99)
       
       filename="HOT_Q"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,I11,I12,I13,I14,I15,I16,I17,I18,I19,CF"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         if(B%solver==GKUA)then
         write(99,*)  " zone ", "i= ", B%nx-1, " j= ", B%ny-1
  	     do j=1,B%ny-1
	     do i=1,B%nx-1
            V_tmp(1:6)=B%HOT_GKUA(1:6,i,j)
            V_tmp(7)=B%HOQ_GKUA(1,i,j)
            V_tmp(8)=B%HOQ_GKUA(2,i,j)
            V_tmp(9)=B%HOQ_GKUA(3,i,j)
         write(99,"(11f15.10)") B%x1(i,j),B%y1(i,j), V_tmp(1:9),B%CF(i,j)
    
         enddo
         enddo
         end if
        enddo
       close(99)
       
       filename="Tao_g_i"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,T1_i,T1_G_i,T2_i,T2_G_i,T12_i,T12_G_i,T3_i,T3_G_i"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         
         write(99,*)  " zone ", "i= ", B%nx, " j= ", B%ny-1
  	     do j=1,B%ny-1
	     do i=1,B%nx

           
         write(99,"(10f20.10)") (B%x1(i,j)+B%x1(i-1,j))/2.0,B%y1(i,j),B%TaoNSi_in(1,i,j),B%TaoNSi_G(1,i,j),&
             &B%TaoNSi_in(2,i,j),B%TaoNSi_G(2,i,j),B%TaoNSi_in(3,i,j),B%TaoNSi_G(3,i,j),B%TaoNSi_in(4,i,j),B%TaoNSi_G(4,i,j)

    
         enddo
         enddo
         
        enddo
       close(99)
       filename="Tao_g_j"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,T1_j,T1_G_j,T2_j,T2_G_j,T12_j,T12_G_j,T3_j,T3_G_j"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
        
         write(99,*)  " zone ", "i= ", B%nx-1, " j= ", B%ny
  	     do j=1,B%ny
	     do i=1,B%nx-1

           
         write(99,"(10f20.10)") B%x1(i,j),(B%y1(i,j)+B%y1(i,j-1))/2.0,B%TaoNSj_in(1,i,j),B%TaoNSj_G(1,i,j),&
             &B%TaoNSj_in(2,i,j),B%TaoNSj_G(2,i,j),B%TaoNSj_in(3,i,j),B%TaoNSj_G(3,i,j),B%TaoNSj_in(4,i,j),B%TaoNSj_G(4,i,j)

    
         enddo
         enddo
        
        enddo
       close(99)
       
       filename="q_g_i"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,qx_i,qx_G_i,qy_i,qy_G_i"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         
         write(99,*)  " zone ", "i= ", B%nx, " j= ", B%ny-1
  	     do j=1,B%ny-1
	     do i=1,B%nx

           
         write(99,"(6f20.10)") (B%x1(i,j)+B%x1(i-1,j))/2.0,B%y1(i,j),B%Qi_in(1,i,j),B%Qi_G(1,i,j),B%Qi_in(2,i,j) ,B%Qi_G(2,i,j)
    
         enddo
         enddo
        
        enddo
       close(99)
       
       filename="q_g_j"
       open(99,file="output/"//trim(filename)//".dat")
       write(99,*) "variables=x,y,qx_j,qx_G_j,qy_j,qy_G_j"
        do m=1,MP%Num_Block
         B => Mesh(1)%Block(m)
         
         write(99,*)  " zone ", "i= ", B%nx-1, " j= ", B%ny
  	     do j=1,B%ny
	     do i=1,B%nx-1

           
         write(99,"(5f20.10)") B%x1(i,j),(B%y1(i,j)+B%y1(i,j-1))/2.0,B%Qj_in(1,i,j),B%Qj_G(1,i,j),B%Qj_in(2,i,j),B%Qj_G(2,i,j)
    
         enddo
         enddo
         
        enddo
       close(99)
       
       open(99,file="output/"//trim(filename_Dem)//".dat")
       write(99,*) "variables=x,y,d,u,v,T,p,Amut"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx, " j= ", B%ny
  	     do j=1,B%ny
	     do i=1,B%nx
           d1=(B%U(1,i,j)+B%U(1,i-1,j)+B%U(1,i,j-1)+B%U(1,i-1,j-1))*0.25d0
           u1=(B%U(2,i,j)+B%U(2,i-1,j)+B%U(2,i,j-1)+B%U(2,i-1,j-1))*0.25d0/d1
           v1=(B%U(3,i,j)+B%U(3,i-1,j)+B%U(3,i,j-1)+B%U(3,i-1,j-1))*0.25d0/d1
           temp=(B%U(4,i,j)+B%U(4,i-1,j)+B%U(4,i,j-1)+B%U(4,i-1,j-1))*0.25d0
           T1=(temp-0.5d0*d1*(u1*u1+v1*v1))/(Cv*d1)
         write(99,"(8f20.10)") B%x(i,j)*Lref,B%y(i,j)*Lref, d1*ND_inf*Matom*Mole_Mass,u1*Ma*sqrt(gamma*RCON*T_inf),&
       &  v1*Ma*sqrt(gamma*RCON*T_inf),T1*T_inf,d1*T1*ND_inf*Matom*Mole_Mass*RCON*T_inf,B%Amu_t(i,j)*Re
         enddo
         enddo
        enddo
       close(99)
       

       

       
      open(99,file="output/"//"VAR_wall.dat")
      write(99,*) "variables=x,y,rho_wall,U_wall,V_wall,Ttra_wall,P_wall,P_t_wall,tau_wall,q_wall,qwall2"
      do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         
        do ksub=1, B%subface
         Bc => B%bc_msg(ksub)
         ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend  
         if( Bc%neighb== BC_Wall)then
            write(99,*)'ZONE I=',abs(iend-ibegin+jend-jbegin),', J=',1,', F=point'
            do i=ibegin,iend    
            do j=jbegin,jend
                
               icheck=max(ibegin,iend)
               jcheck=max(jbegin,jend)
                if((jbegin==jend .and.i==icheck).or.(ibegin==iend.and.j==jcheck))then
                    cycle
                end if
               
                if(Bc%face .eq. 1) then 
                        pos_x=1;pos_y=j
                else if(Bc%face .eq. 2 ) then 
                        pos_x=i;pos_y=1                        
                else if(Bc%face .eq. 3)then
                        pos_x=1;pos_y=j                       
                else if(Bc%face .eq. 4)then
                        pos_x=i;pos_y=1                      
                end if
                
                temp=ND_inf*Matom*Mole_Mass
               write(99,"(11f20.10)") BC%Swall0(1,pos_x,pos_y)*Lref,BC%Swall0(2,pos_x,pos_y)*Lref,BC%Swall0(3,pos_x,pos_y)*temp,&
            &   BC%Swall0(4,pos_x,pos_y)*Ma*sqrt(gamma*RCON*T_inf),BC%Swall0(5,pos_x,pos_y)*Ma*sqrt(gamma*RCON*T_inf),&
                   &BC%Swall0(6,pos_x,pos_y)*T_inf, BC%Swall0(9,pos_x,pos_y)*temp*gamma*RCON*T_inf*Ma**2,BC%Swall0(10,&
                   &pos_x,pos_y)*temp*gamma*RCON*T_inf*Ma**2,BC%Swall0(11,pos_x,pos_y)*temp*gamma*RCON*T_inf*Ma**2,&
             &  BC%Swall0(12,pos_x,pos_y)*temp*(Ma*sqrt(gamma*RCON*T_inf))**3,&
                   &BC%Swall0(13,pos_x,pos_y)*temp*(Ma*sqrt(gamma*RCON*T_inf))**3
               

            end do
            end do
         end if
        end do
      end do
      close(99)
      
      open(99,file="output/"//"Force_wall.dat")
     
      do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
        do ksub=1, B%subface
         Bc => B%bc_msg(ksub)
         if( Bc%neighb== BC_Wall)then
               temp=ND_inf*Matom*Mole_Mass
               write(99,"(4f20.10)") BC%F_Wall(1,1)*temp*gamma*RCON*T_inf*Ma**2,BC%F_Wall(2,1)*temp*gamma*RCON*T_inf*Ma**2,&
                   &BC%F_Wall(1,2)*temp*gamma*RCON*T_inf*Ma**2,BC%F_Wall(2,2)*temp*gamma*RCON*T_inf*Ma**2
         
            
         end if
        end do
      end do
      close(99)
       
      if(MP%Nvar .eq. 5) then
      open(99,file="output/"//"SA2d.dat")
      write(99,*) "variables=x,y,vt"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx+1, " j= ", B%ny+1
  	     do j=0,B%ny
	     do i=0,B%nx
         write(99,"(3f20.10)") B%x1(i,j),B%y1(i,j),B%U(5,i,j) 
         enddo
         enddo
        enddo
       close(99)
      endif

      if(MP%Nvar .eq. 6) then
      open(99,file="output/"//"SST2d.dat")
       write(99,*) "variables=x,y,Kt,Wt"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx+1, " j= ", B%ny+1
  	     do j=0,B%ny
	     do i=0,B%nx
         write(99,"(4f20.10)") B%x1(i,j),B%y1(i,j),B%U(5,i,j),B%U(6,i,j) 
         enddo
         enddo
        enddo
       close(99)
      endif
      
      
      ! open(99,file="output/"//trim(filename)//"tao_ns"//".dat")
      ! write(99,*) "variables=x,y,tao11,tao22,tao12,tao33,qx,qy"
      !  do m=1,MP%Num_Block
      !   B => Mesh(nMesh)%Block(m)
      !   write(99,*)  " zone ", "i= ", B%nx, " j= ", B%ny
  	   !  do j=1,B%ny
	     !do i=1,B%nx
      !     d1=(B%TaoNSi_in(1,i,j)+B%TaoNSj_in(1,i,j)+B%TaoNSi_in(1,i,j-1)+B%TaoNSj(1,i-1,j))*0.25d0
      !     u1=(B%TaoNSi(2,i,j)+B%TaoNSj(2,i,j)+B%TaoNSi(2,i,j-1)+B%TaoNSj(2,i-1,j))*0.25d0
      !     v1=(B%TaoNSi(3,i,j)+B%TaoNSj(3,i,j)+B%TaoNSi(3,i,j-1)+B%TaoNSj(3,i-1,j))*0.25d0
      !     t1=(B%TaoNSi(4,i,j)+B%TaoNSj(4,i,j)+B%TaoNSi(4,i,j-1)+B%TaoNSj(4,i-1,j))*0.25d0
      !     qx=(B%Qi(1,i,j)+B%Qj(1,i,j)+B%Qi(1,i,j-1)+B%Qj(1,i-1,j))*0.25d0
      !     qy=(B%Qi(2,i,j)+B%Qj(2,i,j)+B%Qi(2,i,j-1)+B%Qj(2,i-1,j))*0.25d0
      !     write(99,"(8f20.10)") B%x(i,j),B%y(i,j), d1,u1,v1,t1,qx,qy
      !   enddo
      !   enddo
      !  enddo
      ! close(99)
    end subroutine output
    
       subroutine output2 (nMesh)
       use Global_Var
       
       implicit none
       TYPE (Mesh_TYPE),pointer:: MP
       Type (Block_TYPE),pointer:: B
       Type (BC_MSG_TYPE),pointer:: Bc
       integer  nMesh,i,j,m,ksub,n1
       real(sp) :: d1,u1,v1,T1,temp
       character(len=50):: filename
       
	   if(nMesh .eq. 1) then
          filename="flow2d"
	   else
          write(filename,"('flow2d-'I1.1'')") nMesh    ! flow2d-2.dat ; flow2d-3.dat 
       endif

	   print*, "write data file ...", filename
    
	   MP=>Mesh(nMesh)
       open(99,file="output/"//trim(filename)//"_gkua"//".dat")
       write(99,*) "variables=x,y,d,u,v,T"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx, " j= ", B%ny
  	     do j=1,B%ny
	     do i=1,B%nx
         write(99,"(6f20.10)") B%x(i,j),B%y(i,j), B%GKUA_RONpt(i,j),B%GKUA_Upt(i,j),B%GKUA_Vpt(i,j),B%GKUA_Tpt(i,j)
         enddo
         enddo
        enddo
       close(99)
       
       

    end subroutine output2

