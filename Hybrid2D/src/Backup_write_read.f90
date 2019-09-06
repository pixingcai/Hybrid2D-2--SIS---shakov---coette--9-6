
!-----------------------------------------------
!  read flowfield information for restart
    subroutine init_flow_read
    use Global_var
    use Com_ctrl
    real(sp)  ,external::Omiga22
    real(sp) :: d1,u1,v1,T1,temp,qx,qy,errors
    character(len=50):: fname,fno2,fname2,fname3
    Type (Mesh_TYPE),pointer:: MP
    integer:: i,j,m,step,nMesh,flag,Kflag
    Type (Block_TYPE),pointer:: B
    if(myid==0)then           
        print*, "Init from 'flow2d.dat' ......"
    end if

    Twall=Twall/T_inf
    fno2 = ''
    write(fno2,'(I2.2)'), myid
    fname2="restart2"//"/flag"//".dat"
     
    open(110,file=trim(fname2),action='read')
    read(110,*)  flag
    read(110,*)  Mesh(1)%tt
    read(110,*)  Mesh(1)%Kstep
    close(110)
    if(flag==1)then  
        fname="restart2"//"/macro"//trim(fno2)//".dat"
        fname3="restart2"//"/macro"//".dat"
        Kflag=2
    else 
        fname2="restart1"//"/flag"//".dat"
        open(110,file=trim(fname2),action='read')
        read(110,*)  flag
        read(110,*)  Mesh(1)%tt
        read(110,*)  Mesh(1)%Kstep
        close(110)
        if(flag==1)then
            fname="restart1"//"/macro"//trim(fno2)//".dat"  
            fname3="restart1"//"/macro"//".dat"  
            Kflag=1
        else
            write(*,*)"No backup data is available"
            stop
        end if
            
    end if
    write(*,*)"Reading flowfield data from :",fname
    call backup_Mesco_read(Kflag) 
    call Boundary_condition_onemesh(1)                   ! �����߽����� 
    call update_buffer_onemesh(1)                        ! �ڱ߽�����

     
   end  subroutine init_flow_read
    
    SUBROUTINE backup_Mesco_read(flag)
    use Global_var
    use Com_ctrl
    use const_var ,only :sp,dop
    implicit none
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    Integer :: mBlock,flag
    character(len=50):: fname
    integer::sour_id,alloc_err,count
    integer::dis_list,i,j,iv,jv,iv1,jv1
    integer::nums_file_parall_io,ids

    integer stats(MPI_STATUS_SIZE)
    character*30 start,end
    character*4 fno
    character*6 fno2
    real(sp) beginsave,endsave
    integer :: NUMPROC,iv0,jv0
      
    MP=>Mesh(1)               ! ��������
      
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        
    if (myid == 0)print *,"restart-------kpbak=",kpbak
        fno = ''
        fno2 = ''
        write(fno,'(I2.2)'), myid/1024
        write(fno2,'(I5.5)'), myid
        if(flag==1) then
            fname="restart1/star"//trim(fno)//"/restart"//trim(fno2)//".dat"
        else if(flag==2)then
            fname="restart2/star"//trim(fno)//"/restart"//trim(fno2)//".dat"
        else 
            write(*,*) "Restart reading error"
            stop
    end if
    open(198,file=trim(fname),form='unformatted')        
    do mBlock=1,MP%Num_Block
        B => MP%Block(mBlock) 
        if(B%solver==GKUA)then
            read (198) B%FRD
        end if
    enddo
          
    close(198)
    call GKUA_lisansudu
    call DETFP(Time_Method)
    Init_mass = out_mass!!For mass correct
    call GET_SYM_NUM
    call GKUA_inner_mesco_boundary(1)
    call GKUA_inner_macro_boundary(1)  !������
    call GKUA_Phy_boundary
    call output (1)

         
    end subroutine backup_Mesco_read
    
   subroutine backup_Macro_write (flag)
   use Global_Var
   use Com_ctrl
   implicit none
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   integer  flag,i,j,m,ksub,n1,temp
   real(sp) :: d1,u1,v1,T1,qx,qy,errors
   character(len=50):: fname,fno2,fname2,fname3
   MP=>Mesh(1)
   B => MP%Block(1)
   if(myid==0)then
      write(*,*)"Backing up ......"
   end if
       
   call backup_Mesco_write(flag)

   if(flag==1)then
      fname2="restart1"//"/flag"//".dat"
   else if(flag==2)then
      fname2="restart2"//"/flag"//".dat"
   else 
      write(*,*)"flag error in data back up"
   end if
       
   if(myid==0)then
      open(110,file=trim(fname2),status='replace')
      write(110,*)  1
      write(110,*) Mesh(1)%tt
      write(110,*) Mesh(1)%Kstep
   end if
       
   end subroutine backup_Macro_write
    
    !c
   SUBROUTINE backup_Mesco_write(flag)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      Type (Mesh_TYPE),pointer:: MP
      Type (Block_TYPE),pointer:: B
      Integer :: mBlock,ioflag,i,j,iv,jv,iv1,jv1
      integer::sour_id,alloc_err,count
      integer::dis_list,flag
      integer::nums_file_parall_io,ids

      integer stats(MPI_STATUS_SIZE)

      character*80 fname
      character*30 start,end
      character*4 fno
      character*6 fno2
      real(sp):: beginsave,endsave
      integer :: NUMPROC,iv0,jv0
      
      MP=>Mesh(1)               ! ��������
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      
      nums_file_parall_io = 4
      DO ids=0, nums_file_parall_io-1
         IF(mod(myid,nums_file_parall_io) .eq. ids) THEN 
            fno = ''
            fno2 = ''
            write(fno,'(I2.2)'), myid/1024
            write(fno2,'(I5.5)'), myid
            if(flag==1) then
               fname="restart1/star"//trim(fno)//"/restart"//trim(fno2)//".dat"
            else if(flag==2) then
               fname="restart2/star"//trim(fno)//"/restart"//trim(fno2)//".dat"
            else 
                write(*,*)"flag error in backuping"
                stop
            end if
            open(198,file=trim(fname),form='unformatted')
            do mBlock=1,MP%Num_Block
                   B => MP%Block(mBlock) 
                   if(B%solver==GKUA)then
                      write (198) B%FRD
                   end if            
            enddo    
           close(198)
         end if
      end do
      
      
         call MPI_BARRIER(MPI_COMM_WORLD,IERR)
         if(myid.eq.0) then 
            print *,'the data save file!'
            call flush(6)
         endif
     

      RETURN
      END SUBROUTINE backup_Mesco_write

