      FUNCTION Omiga22(Ts)
      use const_var ,only :sp,dop 
       implicit none
      real(sp)  :: Omiga22,Ts
     
      Omiga22=1.16145/Ts**0.14874+0.52487/exp(0.77320*Ts)+2.16178/exp(2.43787*Ts)
      RETURN
    END
    
      FUNCTION FMINMD(X,Y)
      use const_var ,only :sp,dop 
      implicit none
      real(sp)  :: FMINMD,X,Y
      
         IF ((X.GT.0. .and.  Y.GT.0.).or.(X.LT.0. .and.  Y.LT.0.) ) THEN
            IF (ABS(X) .LE. ABS(Y)) THEN 
               FMINMD=X
            ELSE
               FMINMD=Y
            END IF
         ELSE
            FMINMD=0.
         END IF
      RETURN
    END
    
     FUNCTION FVALMTer(X,Y)
      use const_var ,only :sp,dop
      implicit none
      real(sp)  :: FVALMTer,X,Y,EPSVA
         EPSVA=1.E-30
         FVALMTer=(y*y+0.5*EPSVA)/(x*x+y*y+EPSVA)
      RETURN
    END
    

      
      FUNCTION GDVN(nm,i,j,iv1,jv1,HDVN)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      implicit none
      integer ::nm,i,j,iv1,jv1
      real(sp) ::GDVN,HDVN,tij,crt,vjv,viu,viuq,vij2,GDVM1,vq,FVXDV,HDVM1,&
          & cpt,Ttij,Cij,pij, GDVijv
      Type (Block_TYPE),pointer:: B
      Type (Mesh_TYPE),pointer:: MP
      Type (BC_MSG_TYPE),pointer:: Bc
       MP=>Mesh(1)   
       B => MP%Block(nm)
       
          pij=B%S_GKUA(1,i,j)*B%S_GKUA(4,i,j)
            
            cpt=cpr/(pij*B%S_GKUA(4,i,j))
            crt=B%S_GKUA(1,i,j)/(PI*B%S_GKUA(4,i,j))
            vjv=B%vj(jv1)-B%S_GKUA(3,i,j)
            viu=B%vi(iv1)-B%S_GKUA(2,i,j)
      
            viuq=viu*B%S_GKUA(8,i,j)!qx*Cx
            vij2=(viu*viu+vjv*vjv)/B%S_GKUA(4,i,j)
            GDVM1=crt*EXP(-vij2)
            vq=viuq+vjv*B%S_GKUA(9,i,j) !qx*Cx+qy*Cy

            
            FVXDV=cpt*Vq
            GDVN=GDVM1*(1.+FVXDV*(2.*vij2-4.))
            GDVN=GDVM1
            HDVM1=B%S_GKUA( 4,i,j)*GDVM1/2.
            HDVN=HDVM1*(1.+FVXDV*(2.*vij2-2.))
            HDVN=HDVM1
         
                 crt=B%S_GKUA(1,i,j)/(PI*B%S_GKUA(4,i,j))  
                    vjv=B%vj(jv1)-B%S_GKUA(3,i,j)
                    viu=B%vi(iv1)-B%S_GKUA(2,i,j)
                        vij2=(viu**2+vjv**2)/B%S_GKUA(4,i,j)
                        GDVijv=crt*EXP(-vij2)
                         GDVN=GDVijv
                         HDVN=B%S_GKUA( 4,i,j)*GDVijv/2.
                
          

     

            
      RETURN
    END FUNCTION GDVN
    
        
    FUNCTION Vdof_f(Ttmp)
      use Global_var
      use Com_ctrl
      use const_var ,only :sp,dop
      real(sp)  :: Vdof_f,Ttmp
      Vdof_f=0. 
    END