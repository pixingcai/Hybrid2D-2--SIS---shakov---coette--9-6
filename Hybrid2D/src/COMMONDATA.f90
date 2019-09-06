  module const_var
  
  implicit none
  integer,parameter:: ndim=2,nad=3,NGLT=10,Nsum=11,NRD=2  !7-2
  integer,parameter::  VHS_Mode=0, Sutherland_Mode=1,LJ_Mode=2
  integer, parameter :: sp =selected_real_kind(8)
  integer, parameter :: dop =selected_real_kind(8)
  real(sp) ,parameter::  PI=3.1415926535897932d0, P_rato_limit=2.d0
  real(sp) ,parameter::  Boltz=1.38065E-23,Matom=1.66053873E-27
  real(sp) ::            RCON
  integer,parameter:: Scheme_UD1=0,Scheme_NND2=1, Scheme_UD3=2,Scheme_WENO3=3,Scheme_MUSCL2=4,Scheme_MUSCL3=5, Scheme_OMUSCL2=6
  integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_VanLeer=5,Flux_Ausm=6
  integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2
  integer,parameter:: BC_Wall=-2, BC_Farfield=-4, BC_Symmetry_or_slidewall=-3,BC_inlet=-5,BC_outlet=-6 !GKUA֧��x��Ϊ�Գ���
  integer,parameter:: Hybrid_Mode=1,NS_Mode=2,ES_Mode=3
  integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0
  integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2,Turbulence_SST=3
   integer,parameter:: GKUA_Scheme_UP1=0,GKUA_Scheme_NND4A2=1,GKUA_Scheme_WENO3=3
  integer,parameter:: LAP=2                   ! ��Ϳ�֮��Ľ�����(overlap)��� ��LAP=2���֧��4�ף�LAP=3���֧��6�ף�LAP=4���֧��8�׾��ȣ�
  integer,parameter:: NS_solver=1,GKUA=2 
  real(sp) ,parameter:: PrT=0.9d0                ! PrT ����Plandtl��
  real(sp) :: Init_mass,out_mass
    end module const_var
 
! ȫ�ֱ�����������������������������
!========================================================================================
  module Global_Var    
  use const_var
 
  implicit none
  integer ::MODE
  TYPE BC_MSG_TYPE              ! �߽�������Ϣ
   integer::  f_no, face, ist, iend, jst, jend,  neighb, subface, orient
   integer:: ist_neighb,iend_neighb,jst_neighb,jend_neighb
   real(sp) ,pointer,dimension(:,:) :: Var_Wall !@CARDC
   real(sp) ,pointer,dimension(:,:) :: F_Wall   !@CARDC
   real(sp) ,pointer,dimension(:,:,:)::Swall,Swall0   !@CARDC
   real(sp) ,pointer,dimension(:,:,:)::Qnwall0,Qnwall   !@CARDC
  END TYPE BC_MSG_TYPE

!------------------------------------�����--------------------------------------
   TYPE Block_TYPE                               ! ���ݽṹ������� ���������α����������������Ϣ 
     integer :: Block_no,nx,ny,subface,solver           ! ��ţ�������nx,ny��������
     integer :: Resmax_location_i,Resmax_location_j
     real(sp) ,pointer,dimension(:,:):: x,y,x1,y1   ! (x,y) : coordinate of vortex; (x1,y1): coordinate of cell center  
	 real(sp) ,pointer,dimension(:,:,:) :: U,Un,U1,U0,U2,U3     ! �غ���� (��ʱ�䲽����һʱ�䲽��ֵ), conversation variables 
     real(sp) ,pointer,dimension(:,:,:,:,:) :: FRD,FRD0,FRD00,FRD000 !@CARDC
     real(sp) ,pointer,dimension(:,:) :: CF ,CF0,CF00  !@CARDC
     real(sp) ,pointer,dimension(:,:) ::ATSVX,ATSVY !@CARDC
     real(sp) ,pointer,dimension(:) :: vi,vj,wi,wj,FBAVX,FBAVY  !@CARDC
     real(sp) ,pointer,dimension(:,:,:) :: DUFRD,RHSFRD,S_GKUA,S0_GKUA,S00_GKUA,HOT_GKUA,HOQ_GKUA ,S1_GKUA!@CARDC
     real(sp) ,pointer,dimension(:,:) :: Taoxx,Taoyy,Taoxy,Qx,Qy
     real(sp) ,pointer,dimension(:,:) ::QXES,QYES
     real(sp) ,pointer,dimension(:,:) ::TaoXXES,TaoYYES,TaoZZES,TaoXYES,TaoXXES0,TaoYYES0,TaoZZES0,TaoXYES0
     real(sp) ,pointer,dimension(:,:) ::TLaXX,TLaYY,TLaZZ,TLaXY  !@CARDC
     real(sp) ,pointer,dimension(:,:) :: GKUA_P,GKUA_RON,GKUA_T,GKUA_U,GKUA_V  !@CARDC
      real(sp) ,pointer,dimension(:,:) :: GKUA_Ppt,GKUA_RONpt,GKUA_Tpt,GKUA_Upt,GKUA_Vpt  !@CARDC
   
     real(sp) ,pointer,dimension(:,:,:):: BNIRHS!@CARDC
     real(sp) ,pointer,dimension(:,:,:)::TaoNSi_in,TaoNSj_in,Qi_in,Qj_in,TaoNSi_G,TaoNSj_G,Qi_G,Qj_G !NSӦ������
     real(sp) ,pointer,dimension(:,:,:)::TaoNSi_out,TaoNSj_out,Qi_out,Qj_out
     real(sp) ,pointer,dimension(:,:,:)::TaoNS_out,Q_out
     real(sp) ,pointer,dimension(:,:,:) :: Res ,Res_i     ! �в� ����ͨ����
     real(sp) ,pointer,dimension(:,:)::    dt       ! (�ֲ�)ʱ�䲽��
     real(sp) ,pointer,dimension(:,:,:) :: QF       ! ǿ�Ⱥ��� (���������д�����ʹ��)
     real(sp) ,pointer,dimension(:,:,:) :: deltU    ! �غ�����Ĳ�ֵ, dU=U(n+1)-U(n)  ��������ʹ��
	 real(sp) ,pointer,dimension(:,:,:)::  dU       ! �غ�����Ĳ�ֵdU=U(n+1)-U(n), LU-SGS������ʹ��
	 real(sp) ,pointer,dimension(:,:):: Amu,Amu_t ! ����ճ��ϵ��������ճ��ϵ��
	        ! ������: vol���������; si,sj i,j-���������߽糤��; (ni1,ni2) i-���������߽編����; (nj1,nj2) j-���������߽編����;
     real(sp) ,pointer,dimension(:,:):: vol,si,sj,ni1,ni2,nj1,nj2        
     real(sp) ,pointer,dimension(:,:):: Lci,Lcj,Lvi,Lvj       ! �װ뾶 (Lci,Lcj i-,j-������ճ���װ뾶; Lvi,Lvj i-, j-����ճ�����װ뾶)
     real(sp) ,pointer,dimension(:,:):: dw                    ! ������ľ��� ��k-w SSTģ����ʹ�ã�
     real(sp) ::  Res_max(6),Res_rms(6)                 ! ���в�������в�
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg    ! �߽�������Ϣ
   End TYPE Block_TYPE  

!---------------------------���� -------------------------------------------------------- 
                                      ! (�絥������ֻ��1�ף���������񣬿����ж���) 
   TYPE Mesh_TYPE                     ! ���ݽṹ�����񡱣� �������α��������������Ϣ
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep,Iter          ! ������ (1��Ϊ��ϸ����2��Ϊ������ 3��Ϊ��������...)�����������������Ŀ,ʱ�䲽 
     integer:: Nvar  ! ���������̣���Ŀ����ʹ��k-wģ�ͣ�����6������ ��ϡ����ʹ����ģ�ͣ������������4����
     real(sp) ::  tt                 !  �ƽ���ʱ��
	 TYPE (Block_TYPE),pointer,dimension(:):: Block     ! ������顱  �������ڡ����񡱣�

!                                                       ���Ʋ��������ڿ�����ֵ������ͨ������������ģ�͵�    
!             ��Щ���Ʋ��������ڡ����񡱣���ͬ�����񡱿��Բ��ò�ͬ�ļ��㷽��������ģ�͵ȡ�	 �����磬�������õ;��ȷ�����������ʹ������ģ��,...��
	integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction
   End TYPE Mesh_TYPE
!---------------------------------------------------------------------------------------------


! global variables                                       ���ӳ�����ɼ���ȫ�ֱ���
!----------------------------------------------------------------------------
   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh       ! ������ ������
   integer,save:: Num_Mesh          ! ���������  
   integer,save:: Nvar   ! ���̣�����������Ŀ�� ��ʹ��BLģ��Nvar=4;  ʹ��SAģ�� Nvar=5, ��ʹ��K-W SSTģ�� Nvar=6 (4����������+k����+w���̣�                                                        
 ! ���Ʊ���  If_viscous=0 Euler���̣�1 N_S���̣� Iflag_turbulence_model ����ģ�ͣ�Ŀǰ�汾ֻ֧��BL��
 ! Iflag_Scheme ��ֵ��ʽ�� Iflag_flux ͨ�������� Iflag_local_dt �Ƿ���þֲ�ʱ�䲽�� ; Num_Threads OpenMP���õ��߳���
   integer,save::  Kstep_save,If_viscous, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method,Kstep_show,Num_Threads


 ! global parameter (for all Meshes )                     ��������, ��ȫ�塰���񡱶�����
 ! Ma: Mach��;  Re: Reynolds��; Pr: Prandtl��; Cp,Cv: ��ѹ�����ݱ���;
 ! t_end: End time ���������ʱ��; P_outlet: ����ѹ������������������������,����������;  Twall: ���£�������������λK��
 ! T_inf: �����¶� (������ֵ����λK),��surthland��ʽ��ʹ��; 
 ! Kt_inf, Wt_inf: �����Ķ��ܡ����ܱȺ�ɢ��  (Amut_inf=Kt_inf/Wt_inf)
   
     real(sp) ,save:: Ma,Re,gamma,Pr,AoA,Cp,Cv,t_end,p_outlet,T_inf,Twall,vt_inf,Kt_inf,Wt_inf
 
 !                                         ȫ�ֿ��Ʋ�����������ֵ������ͨ������������ģ�͵� ����Щֻ����ϸ������Ч��
    real(sp) ,save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin ,dt_global_GKUA   ! RK�����еĳ�������ʱ�䲽���йص���
    integer,save:: Pre_Step_Mesh(3)                                                 ! ������ֵʱ��������Ԥ��������

  end module Global_Var
   
!----------------------------------------------------------------------------
! ���������� ������ÿ��ʱ�����ڴ棬�ÿ����������ͷţ�������ʱ������ 
! ��Щ��������ֵ���豣��;  
module Flow_Var
 use const_var ,only :sp,dop
 implicit none 
   real(sp) , save,pointer,dimension(:,:)::  d,uu,v,T,p,cc  ! �ܶȡ�x-�ٶȡ�y-�ٶȡ�ѹ�������٣�
   real(sp) , save,pointer,dimension(:,:,:):: Fluxi,Fluxj         ! i- ��j-�����ͨ��
    end module Flow_Var
    
    

    

      module Com_ctrl
          use const_var ,only :sp,dop
         use mpi
         implicit none
         real(sp)  :: KN,X1S,cfc,cpr,cnw,cbnw,xs,OMGgas,ALPgas,cfc_Su,Sgas
         integer orig_group,new_group,new_rank
         integer :: maxnode,initialtime(8),finaltime(8)
          !maxnode    : ���п顢����ά�������������
         integer  ::  maxnbd
         integer,allocatable::rank1(:,:),rank2(:)
         integer  :: nvi,nvj,nmt,nvt !,imt,jmt
         integer  :: nvit,nvjt,kbody,KMAX_SAVE,kpbak,KP,KPFLOW
         integer  :: XDIM, YDIM, I1MAX, J1MAX
         integer  :: GKUA_IST, GKUA_IEND, GKUA_JST, GKUA_JEND, GKUA_IST1,GKUA_IEND1,GKUA_JST1, &
                    & GKUA_JEND1, GKUA_IPT, GKUA_JPT, GKUA_IPT1, GKUA_JPT1
         integer  :: myid,nprocs,ierr,prt_id,prt_x,prt_y,status(MPI_STATUS_SIZE)
         integer  :: myid_x,myid_y,new_comm
         integer  :: KGL,KGLN,KGH,KHW,NFM,NFM1,nf,nft,KCDOM,NIstepmax ,KPtmp
         real(sp)  :: NIerrormax
         
         real(sp)  :: dvx,vxdown,vxup,dvy,vydown,vyup
         real(sp)  :: LMD00,MUREF,ROREF,TREF,Lref,vis_ref,vis_inf,p_ref
         real(sp)  :: P00,R00,T00,u00,v00,M00,CRT00,S00,TwT00
         real(sp)  :: Zrot,Zvib,Tcv,Rdof,ZrotD,ZvibD,ZtraD
         real(sp)  :: kfly,ND_inf,Tw_inf,p_inf,DSLJ,TSLJ
         real(sp)  :: MFP_REF,Kn_REF,MFP_inf,KN_inf,Cgasmodel,Dgas_REF,Dgas_inf
         real(sp)  :: time,timet
         real(sp)  :: CFLst,CFLed,EPS,ERR,DTcol,PSt,RSt,TSt,Tomigaref
         real(sp) :: ND_DSMC,T00_DSMC,Mole_Mass,mach,fLamda !,EPSMUSCL
         integer  :: KHS,KDSMC,KShakhov,KES,khalf,kall,KChangeCFL
         LOGICAL  :: LP
         integer  :: Iflag_Scheme0,NFCFLst,NFCFLed,NFUP1,DNFCFL !NRK
         integer  :: NFCFLst0,NFCFLed0,NumABFRD
         integer  :: tnumnbd,numnbd !numnbd:�ܹ����ٱ�  
         integer  :: next,ntyptmp,CPUst,CPUed,nextwall  !Ϊ�Գ���׼����ǰ���̶�Ӧ��CPU
         real(sp)  :: characterlength

      end module Com_ctrl
