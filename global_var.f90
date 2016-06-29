module global_var
implicit none

!	constants
	! Arnett Bowers: G = 6.6732D-8; M_Solar = 1.987D33; c = 2.9979D10
	real(8), parameter :: pi=3.14159265358979323846D0
	real(8), parameter :: Grav_Const =6.67428D-8
	real(8), parameter :: M_Solar =	1.98855D33
	real(8), parameter :: c =2.99792458D10
	real(8), parameter :: e = 4.8032042510D-10
	real(8), parameter :: k_b = 1.3806488D-16
	real(8), parameter ::  h_bar = 1.0545718d-27

!	input parameters
	real(8),parameter :: delr_0 = 1.d0						! starting dr for adjustable grid of static profile
	real(8),parameter :: delr_adj = 0.0115 !0.05d0
	real(8),parameter :: delr_ref = delr_adj/1.d2

	character (len=30) :: eos_tab =	"input/eos_APR_ILQ.txt" !"input/eos_ArtXVIII.txt"		!"input/eos_ArtVIII.txt"
	integer, parameter :: eos_tab_Nmax = 1000

	character (len=30) :: comp_tab = "input/comp_DH01_HP.txt"
	integer, parameter :: comp_tab_Nmax = 1000

	character (len=30) :: she_tab = "input/Rshear_Steiner_SLy4.txt"
	integer, parameter :: she_tab_Nmax = 9000
	integer, parameter :: she_format = 2

	integer, parameter :: NIface = 2	! no of interface within the crust; including CC, OC interface
	integer, parameter :: sp_Nmax = 200000
	integer, parameter :: n_poly_tail = 1000

	integer, parameter :: pg_Ni = 5000
	integer, parameter :: pg_N1 = 8000 !10000
	integer, parameter :: pg_N2 = 3000 !3000 !15000
	integer, parameter :: pg_N3 = 1500
	integer, parameter :: pg_NOu = 50000

	real(8), parameter :: pg_ri_frac = 1.d-4 !1.d-4	!4.d-2	!1.d-4
	real(8), parameter :: pg_rs_frac = 0.d0 !!!!!!!!!!!!******************** Causes a weird mode to Arise in LD formalism (Conversion from H0 K to Z Z`; a singular point is involved)

	integer,parameter :: muller_n_max = 400
	integer :: NMat = 4
	integer, parameter :: NEqt = 6

	character (len=30) :: pef_yA = "data/y_Co.txt"	! output files of pulsation equation solutions
	character (len=30) :: pef_yB = "data/y_Oc.txt"
	character (len=30) :: pef_zC = "data/z3_Cr.txt" 
	character (len=30) :: pef_zD = "data/z4_Cr.txt"

	character (len=30) :: pef_FN_Co1 = "data/y_Co1.txt" 
	character (len=30) :: pef_FN_Co2 = "data/y_Co2.txt"
	character (len=30) :: pef_FN_Cr1 = "data/y_Cr1.txt" 
	character (len=30) :: pef_FN_Cr2 = "data/y_Cr2.txt"
	character (len=30) :: pef_FN_Cr3 = "data/y_Cr3.txt" 
	character (len=30) :: pef_FN_Cr4 = "data/y_Cr4.txt"
	character (len=30) :: pef_FN_Cr5 = "data/y_Cr5.txt"
	character (len=30) :: pef_FN_Oc1 = "data/y_Oc1.txt" 
	character (len=30) :: pef_FN_Oc2 = "data/y_Oc2.txt"

	character (len=30) :: pef_LD_Co1 = "data/y_Co1.txt"	! for LD formalism
	character (len=30) :: pef_LD_Co2 = "data/y_Co2.txt"
	character (len=30) :: pef_LD_Co3 = "data/y_Co3.txt"
	character (len=30) :: pef_LD_Cr1 = "data/y_Cr1.txt"
	character (len=30) :: pef_LD_Cr2 = "data/y_Cr2.txt"
	character (len=30) :: pef_LD_Cr3 = "data/y_Cr3.txt"
	character (len=30) :: pef_LD_Cr4 = "data/y_Cr4.txt"
	character (len=30) :: pef_LD_Cr5 = "data/y_Cr5.txt"
	character (len=30) :: pef_LD_Oc1 = "data/y_Oc1.txt"
	character (len=30) :: pef_LD_Oc2 = "data/y_Oc2.txt"
	character (len=30) :: pef_LD_Oc3 = "data/y_Oc3.txt"
	character (len=30) :: pef_LD_Ou = "data/y_Ou.txt"
	
	character (len=30) :: pef_y_FM_Cr = "data/y_Cr.txt" 

!	setting variables
	!	Allowed to be changed in setting page
	!	Mainly contains the options variables in the program
	!	To introduce a new setting:
	!	(1) define the variable in this module
	!	(2) change "settings" module (UI, verifying statement)
	!	(3) change "main" module or any other module corresponding to the setting variable

	real(8) :: l_0 = 2.d0
	real(8) :: rho_0 = 2.d15 !3.103d15 !0.988d15
	real(8) :: P_t = 1.7d34   !1.7d34  !4.90424d32 !5.37190965d32 !4.54453d32 !5.24168d32	!3.847d32
	real(8) :: P_g =  1.d29 !1.151d21
	real(8) :: P_min = 5.266d21	!1.d29 age 100yr SLy4	!5.d22	!5.266d21 !5.82d15

	integer :: sp_opt = 2, pg_opt = 2, ei_opt = 2, pe_opt = 9, pes_opt = 5	 ! remember to switch NMat if pe_opt is changed
	integer :: eos_opt = 3, comp_opt = 6, shear_fcn_opt = 2
	integer :: newt_V_opt = 1	! switch V_tilde between dPdr and rho g; opt 2 - V=P'/(-dP/dr)/r; opt 3 - V=P'/rho/g/r
	real(8) :: poly_K = 180.d0	! K in geometrized unit 
	real(8) :: poly_n = 1.d0	! polytropic index
	real(8) :: drho = 0.d0		! density jump (fraction) in polytropic model
	real(8) :: rho_t = 9.d14	! transition density in polytropic model	!6.15d13 !1.d12	check shear and drho effects
	real(8) :: poly_cs = 2.d8 !2.d8
	real(8) :: QS_a4 = 0.7d0	! parameters of quark star models; in MeV unit; ref: PRD 89 103014 (2014)
	real(8) :: QS_a2 = 0.d0
	real(8) :: B_eff = 145.d0 ** 4
	real(8) :: QS_Gap = 5.d0
	real(8) :: QS_mu0 = 2.47d0	! unit MeV/fm^3; ref: PRD 89 103014 (2014)
	real(8) :: Ki_poly = 0.015	! shear modulus proportional contant Andersson 2011
	real(8) :: mu_poly = 1.d-14	! shear modulus contant Andersson 2011 (in km^-2)

	real(8) :: mu_red = 1.d0	! artificially reduce shear modulus
	real(8) :: rel = 1.d0	! adjust relativistic correction for relat cowling, range 0 < rel <= 1

	real(8) :: R_inf_f = 25.d0		! Zerilli function final r x omega (theoretically should be infinity)

	real(8) :: R_t= 0.9d0				!	use with sp_opt 5; fixing Rc/R0; added after v10
	real(8) :: ROc_t = 0.99d0
!	global variables


	real(8), dimension(0:sp_Nmax) :: nb, rho, P, m, mu, sp_r, delr, nu
	real(8) :: M0, R0, R_mid, R_g
	integer :: eos_n, comp_n, shear_n
	integer :: sp_N1, sp_N2, sp_Nt
	logical :: sp_eos_fcall

	character (len=30), dimension(1:NIface-1) :: pef_ziC, pef_ziD
	character (len=30), dimension(1:NIface-1) :: pef_LD_Cri1, pef_LD_Cri2, pef_LD_Cri3, pef_LD_Cri4, pef_LD_Cri5
	real(8), dimension(1:NIface-1) :: P_i
	integer :: sp_Ni(1:NIface-1)
	real(8), dimension(1:NIface-1, 0:pg_Ni*2) :: r_Cri, x_Cri, P_Cri, rho_Cri, m_Cri, mu_Cri, nu_Cri
	real(8), dimension(1:NIface-1) :: dx_Cri, XCri_i, XCri_f

	real(8), dimension(0:pg_N1*2) :: r_Co, x_Co, P_Co, rho_Co, m_Co, mu_Co, nu_Co
	real(8), dimension(0:pg_N2*2) :: r_Cr, x_Cr, P_Cr, rho_Cr, m_Cr, mu_Cr, nu_Cr
	real(8), dimension(0:pg_N3*2) :: r_Oc, x_Oc, P_Oc, rho_Oc, m_Oc, nu_Oc
	real(8) :: P_Ou, rho_Ou, m_Ou, nu_Ou
	real(8) :: dx_Co, dx_Cr, dx_Oc, XCo_i, XCo_f, XCr_i, XCr_f, XOc_i, XOc_f
	real(8) :: dr_Ou, ROu_i, ROu_f
	integer :: pg_ji, pg_js
	
	real(8) :: afreqH, afreqL, afreq, Omega_sq
	complex(8) :: OmeC_sq
	logical :: pe_fcall

	real(8), allocatable :: Coe(:)
	complex(8), allocatable :: CoeC(:)
	integer :: gr_opt

	complex(8) :: k2
	logical :: zero_freq

endmodule global_var