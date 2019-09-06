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
  integer,parameter:: BC_Wall=-2, BC_Farfield=-4, BC_Symmetry_or_slidewall=-3,BC_inlet=-5,BC_outlet=-6 !GKUA支持x轴为对称轴
  integer,parameter:: Hybrid_Mode=1,NS_Mode=2,ES_Mode=3
  integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0
  integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2,Turbulence_SST=3
   integer,parameter:: GKUA_Scheme_UP1=0,GKUA_Scheme_NND4A2=1,GKUA_Scheme_WENO3=3
  integer,parameter:: LAP=2                   ! 块和块之间的交叠区(overlap)宽度 （LAP=2最高支持4阶，LAP=3最高支持6阶，LAP=4最高支持8阶精度）
  integer,parameter:: NS_solver=1,GKUA=2 
  real(sp) ,parameter:: PrT=0.9d0                ! PrT 湍流Plandtl数
  real(sp) :: Init_mass,out_mass
    end module const_var
 
! 全局变量，包括参数、几何量及物理量
!========================================================================================
  module Global_Var    
  use const_var
 
  implicit none
  integer ::MODE
  TYPE BC_MSG_TYPE              ! 边界链接信息
   integer::  f_no, face, ist, iend, jst, jend,  neighb, subface, orient
   integer:: ist_neighb,iend_neighb,jst_neighb,jend_neighb
   real(sp) ,pointer,dimension(:,:) :: Var_Wall !@CARDC
   real(sp) ,pointer,dimension(:,:) :: F_Wall   !@CARDC
   real(sp) ,pointer,dimension(:,:,:)::Swall,Swall0   !@CARDC
   real(sp) ,pointer,dimension(:,:,:)::Qnwall0,Qnwall   !@CARDC
  END TYPE BC_MSG_TYPE

!------------------------------------网格块--------------------------------------
   TYPE Block_TYPE                               ! 数据结构：网格块 ；包含几何变量及物理变量的信息 
     integer :: Block_no,nx,ny,subface,solver           ! 块号；网格数nx,ny；子面数
     integer :: Resmax_location_i,Resmax_location_j
     real(sp) ,pointer,dimension(:,:):: x,y,x1,y1   ! (x,y) : coordinate of vortex; (x1,y1): coordinate of cell center  
	 real(sp) ,pointer,dimension(:,:,:) :: U,Un,U1,U0,U2,U3     ! 守恒变量 (本时间步及上一时间步的值), conversation variables 
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
     real(sp) ,pointer,dimension(:,:,:)::TaoNSi_in,TaoNSj_in,Qi_in,Qj_in,TaoNSi_G,TaoNSj_G,Qi_G,Qj_G !NS应力张量
     real(sp) ,pointer,dimension(:,:,:)::TaoNSi_out,TaoNSj_out,Qi_out,Qj_out
     real(sp) ,pointer,dimension(:,:,:)::TaoNS_out,Q_out
     real(sp) ,pointer,dimension(:,:,:) :: Res ,Res_i     ! 残差 （净通量）
     real(sp) ,pointer,dimension(:,:)::    dt       ! (局部)时间步长
     real(sp) ,pointer,dimension(:,:,:) :: QF       ! 强迫函数 (多重网格法中粗网格使用)
     real(sp) ,pointer,dimension(:,:,:) :: deltU    ! 守恒变量的差值, dU=U(n+1)-U(n)  多重网格使用
	 real(sp) ,pointer,dimension(:,:,:)::  dU       ! 守恒变量的插值dU=U(n+1)-U(n), LU-SGS方法中使用
	 real(sp) ,pointer,dimension(:,:):: Amu,Amu_t ! 层流粘性系数、湍流粘性系数
	        ! 几何量: vol控制体面积; si,sj i,j-方向控制体边界长度; (ni1,ni2) i-方向控制体边界法方向; (nj1,nj2) j-方向控制体边界法方向;
     real(sp) ,pointer,dimension(:,:):: vol,si,sj,ni1,ni2,nj1,nj2        
     real(sp) ,pointer,dimension(:,:):: Lci,Lcj,Lvi,Lvj       ! 谱半径 (Lci,Lcj i-,j-方向无粘项谱半径; Lvi,Lvj i-, j-方向粘性项谱半径)
     real(sp) ,pointer,dimension(:,:):: dw                    ! 到壁面的距离 （k-w SST模型中使用）
     real(sp) ::  Res_max(6),Res_rms(6)                 ! 最大残差，均方根残差
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg    ! 边界链接信息
   End TYPE Block_TYPE  

!---------------------------网格 -------------------------------------------------------- 
                                      ! (如单重网格，只有1套；如多重网格，可以有多套) 
   TYPE Mesh_TYPE                     ! 数据结构“网格”； 包含几何变量及物理变量信息
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep,Iter          ! 网格编号 (1号为最细网格，2号为粗网格， 3号为更粗网格...)，网格块数，网格数目,时间步 
     integer:: Nvar  ! 变量（方程）数目，如使用k-w模型，则有6个变量 （稀网格不使用湍模型，因而变量仍是4个）
     real(sp) ::  tt                 !  推进的时间
	 TYPE (Block_TYPE),pointer,dimension(:):: Block     ! “网格块”  （从属于“网格”）

!                                                       控制参数，用于控制数值方法、通量技术、湍流模型等    
!             这些控制参数从属于“网格”，不同“网格”可以采用不同的计算方法、湍流模型等。	 （例如，粗网格用低精度方法，粗网格不使用湍流模型,...）
	integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction
   End TYPE Mesh_TYPE
!---------------------------------------------------------------------------------------------


! global variables                                       各子程序均可见的全局变量
!----------------------------------------------------------------------------
   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh       ! 主数据 “网格”
   integer,save:: Num_Mesh          ! 网格的套数  
   integer,save:: Nvar   ! 方程（变量）的数目， 如使用BL模型Nvar=4;  使用SA模型 Nvar=5, 如使用K-W SST模型 Nvar=6 (4个基本方程+k方程+w方程）                                                        
 ! 控制变量  If_viscous=0 Euler方程，1 N_S方程； Iflag_turbulence_model 湍流模型（目前版本只支持BL）
 ! Iflag_Scheme 数值格式； Iflag_flux 通量技术； Iflag_local_dt 是否采用局部时间步长 ; Num_Threads OpenMP采用的线程数
   integer,save::  Kstep_save,If_viscous, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method,Kstep_show,Num_Threads


 ! global parameter (for all Meshes )                     流动参数, 对全体“网格”都适用
 ! Ma: Mach数;  Re: Reynolds数; Pr: Prandtl数; Cp,Cv: 定压、定容比热;
 ! t_end: End time 计算结束的时间; P_outlet: 出口压力（亚声速内流计算必须给定,无量纲量）;  Twall: 壁温（有量纲量，单位K）
 ! T_inf: 来流温度 (有量纲值，单位K),在surthland公式中使用; 
 ! Kt_inf, Wt_inf: 来流湍动能、湍能比耗散率  (Amut_inf=Kt_inf/Wt_inf)
   
     real(sp) ,save:: Ma,Re,gamma,Pr,AoA,Cp,Cv,t_end,p_outlet,T_inf,Twall,vt_inf,Kt_inf,Wt_inf
 
 !                                         全局控制参数，控制数值方法、通量技术及湍流模型等 （有些只对最细网格有效）
    real(sp) ,save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin ,dt_global_GKUA   ! RK方法中的常量，与时间步长有关的量
    integer,save:: Pre_Step_Mesh(3)                                                 ! 构建初值时，粗网格预迭代步数

  end module Global_Var
   
!----------------------------------------------------------------------------
! 流场物理量 （计算每块时申请内存，该块计算结束后释放；属于临时变量） 
! 这些物理量的值无需保留;  
module Flow_Var
 use const_var ,only :sp,dop
 implicit none 
   real(sp) , save,pointer,dimension(:,:)::  d,uu,v,T,p,cc  ! 密度、x-速度、y-速度、压力、声速；
   real(sp) , save,pointer,dimension(:,:,:):: Fluxi,Fluxj         ! i- 及j-方向的通量
    end module Flow_Var
    
    

    

      module Com_ctrl
          use const_var ,only :sp,dop
         use mpi
         implicit none
         real(sp)  :: KN,X1S,cfc,cpr,cnw,cbnw,xs,OMGgas,ALPgas,cfc_Su,Sgas
         integer orig_group,new_group,new_rank
         integer :: maxnode,initialtime(8),finaltime(8)
          !maxnode    : 所有块、所有维数中最大网格数
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
         integer  :: tnumnbd,numnbd !numnbd:总共多少边  
         integer  :: next,ntyptmp,CPUst,CPUed,nextwall  !为对称轴准备当前进程对应的CPU
         real(sp)  :: characterlength

      end module Com_ctrl
