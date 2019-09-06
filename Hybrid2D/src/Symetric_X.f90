  !!!Get the geometry infomation on X_symetric boundary    
      
  Subroutine GET_SYM_NUM
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      integer :: j,ird,jv ,iv,jv1,iv1,ksub
      integer ::st,ed,nm
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       
      NumABFRD=0
      DO jv=GKUA_JST,GKUA_JEND
         jv1=jv+GKUA_JPT
      DO iv=GKUA_IST,GKUA_IEND
         iv1=iv+GKUA_IPT
         do nm=1,MP%Num_Block
            B => MP%Block(nm)
            do  ksub=1,B%subface
                Bc=> B%bc_msg(ksub)
                
                if( Bc%neighb .eq. BC_Symmetry_or_slidewall ) then
             		  if((Bc%face .eq. 1).or.(Bc%face .eq. 3)) then
                         st=Bc%jst;ed=Bc%jend-1
                      else if((Bc%face .eq. 2).or.(Bc%face .eq. 2)) then
                         st=Bc%ist;ed=Bc%iend-1
                      else
                            write(*,*)"error in SUB:GET_SYM_NUM "
                            stop
                      endif
                  do j=st,ed
                      do ird=1,NRD
                         NumABFRD=NumABFRD+1
                      enddo
                  enddo
               endif
            enddo
         enddo
      end do 
      end do
      
      if(myid==0)write(*,*)"on symmetry boundary:NumABFRD=",NumABFRD
 
      end subroutine GET_SYM_NUM