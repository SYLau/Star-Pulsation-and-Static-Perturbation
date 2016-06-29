module settings_opt01
use global_var
contains

	subroutine se01_Ctrl
	implicit none
	integer :: ans
	logical :: pass

14		continue
		pass = .false.
		
write(*,*) "Note: module - subroutine ei01_o_soln_C1L1R1 - pes_opt == 1 .and. newt_V_opt == 2 modified!!"
write(*,*) "Note: module - subroutine pes03_E_Fluid - freq = 0 V = 0 not sure"
		if (eos_opt == 1) then
			write(*,*) "EOS: ", eos_tab
		elseif (eos_opt == 2) then
			write(*,*) "EOS: polytrope"
		elseif (eos_opt == 3) then
			write(*,*) "EOS: MIT Bag Model"
		elseif (eos_opt == 4) then
			write(*,*) "EOS: Hybrid Star: MIT + NS"
		elseif (eos_opt == 5) then
			write(*,*) "EOS: Hybrid Star: Poly + NS"
		elseif (eos_opt == 6) then
			write(*,*) "EOS: Hybrid Star: MIT + NS"
		else
			write(*,*) "err: eos_tab"
			pause
		endif

		write(*,*) "=========================================================="
		write(*,*) "Options:"
		write(*,*) "1. Begin"
		write(*,*) "2. Display settings"
		write(*,*) "3. Change settings"
		write(*,*) "4. Display NS parameters"
		write(*,*) "5. Change NS parameters"
		write(*,*) "=========================================================="
		read(*,*) ans
		write(*,*) "=========================================================="
		if (ans == 1) then
			call se01_01(pass)
		elseif (ans == 2) then
			call se01_02
		elseif (ans == 3) then
			call se01_03
		elseif (ans == 4) then
			call se01_04
		elseif (ans == 5) then
			call se01_05
		else
			write(*,*) "invalid input"
		endif

		if (pass == .false.) then
			write(*,*) "=========================================================="
			goto 14
		else
			return
		endif
	endsubroutine se01_Ctrl

	subroutine se01_01(pass)
	! Allow the code to run after briefly verifying the settings
	implicit none
	logical :: pass
		!	Verifying statement
		if (sp_opt > 6 .or. sp_opt < 1 .or. pg_opt > 4 .or. pg_opt < 1 &
			.or. ei_opt > 2 .or. ei_opt < 1 &
			.or. pe_opt > 9 .or. pe_opt < 1  .or. pes_opt > 5 & 
			.or. pes_opt < 1 .or. rel > 1.d0 .or. rel < 0.d0 &
			.or. comp_opt < 1 .or. comp_opt > 7 .or. shear_fcn_opt < 1 &
			.or. shear_fcn_opt > 2) then
				write(*,*) "invalid settings"
				pass = .false.
			else
				pass = .true.
				write(*,*) "=========================================================="
				write(*,*) "		Program Begins"
				write(*,*) "=========================================================="
				return
		endif
	endsubroutine se01_01

	subroutine se01_02
	implicit none
			write(*,*) "Static Profile Options:"
			write(*,*) "	(1) 2 comp TOV model, EOS table, standard"
			write(*,*) "	(2) 2 comp TOV model, EOS table, including compute nu(r)"
			write(*,*) "	(3) 3 comp TOV model, EOS table, including compute nu(r)"
			write(*,*) "	(4) 3 comp TOV model, EOS table, includes nu(r), multi-interfaces"
			write(*,*) "	(5) 2 comp TOV model, EOS table, including compute nu(r); fixed Rc/R0 or Mc/M0 (Adjust in module)"
			write(*,*) "	Choice: ", sp_opt
			write(*,*)
			write(*,*) "Pulsation Grid Options:"
			write(*,*) "	(1) 2 comp model, standard"
			write(*,*) "	(2) 2 comp model, interpolate nu(r) for Relativistic Cowling"
			write(*,*) "	(3) 3 comp model, interpolate nu(r) for Relativistic Cowling"
			write(*,*) "	(4) 3 comp model, interpolate nu(r), multi-interfaces"
			write(*,*) "	Choice: ", pg_opt
			write(*,*) 
			write(*,*) "EigenValue Problem Solver Options:"
			write(*,*) "	(1) Newtonian or Relativisitic CA"
			write(*,*) "	(2) LD formalism"
			write(*,*) "	Choice: ", ei_opt
			write(*,*) 
			write(*,*) "Pulsation Equation Iterations Options:"
			write(*,*) "	(1) 2 comp model; shooting"
			write(*,*) "	(2) 3 comp model; shooting"
			write(*,*) "	(3) 1 comp model; LD formalisms; shooting"
			write(*,*) "	(4) 1 comp model"
			write(*,*) "	(5) 2 comp model; LD & Andersson; shooting"
			write(*,*) "	(6) 3 comp model; LD & Andersson; shooting"
			write(*,*) "	(7) 2 comp model; Full Newtonian; shooting"
			write(*,*) "	(8) 1 comp model; Full Newtonian; shooting"
			write(*,*) "	Choice: ", pe_opt
			write(*,*) 
			write(*,*) "Pulsation Equation Set Options:"
			write(*,*) "	(1) Newtonian Cowling"
			write(*,*) "	(2) Relativistic Cowling"
			write(*,*) "	(3) Full Relativistic (LD)"
			write(*,*) "	(4) Full Newtonian"
			write(*,*) "	Choice: ", pes_opt
			write(*,*) 
			write(*,*) "EOS Options:"
			write(*,*) "	(1) Using EOS Table"
			write(*,*) "	(2) Using Polytropic EOS"
			write(*,*) "	(3) Using Quark Star EOS"
			write(*,*) "	(4) Using Hybrid Star EOS: Quark Inner Core + NS For the Rest"
			write(*,*) "	(5) Using Hybrid Star EOS: Poly Inner Core + NS For the Rest"
			write(*,*) "	(6) Ref (5); input rho_t"
			write(*,*) "	Choice: ", eos_opt
			write(*,*) 
			write(*,*) "Shear Modulus Input Options:"
			write(*,*) "	(1) Using Composition Table"
			write(*,*) "	(2) Using Shear Modulus Table"
			write(*,*) "	(3) 0 Shear Modulus"
			write(*,*) "	(4) Polytrope: Constant Shear Speed"
			write(*,*) "	(5) Quark Star: Formula"
			write(*,*) "	(6) Solid Quark Star: Formula"
			write(*,*) "	Choice: ", comp_opt
			write(*,*) 
			write(*,*) "Shear Modulus Formula Options:"
			write(*,*) "	(1) McDermott (1988)"
			write(*,*) "	(2) D Tsang (2012)"
			write(*,*) "	Choice: ", shear_fcn_opt
			write(*,*) 
			write(*,*) "Number of matching conditions:"
			write(*,*) "	NMat= ", NMat
			write(*,*) 
			write(*,*) "Shear reduction factor:"
			write(*,*) "	mu_red= ", mu_red
			write(*,*) 
			write(*,*) "Relativistic Corrrection factor:"
			write(*,*) "	Rel= ", rel
	endsubroutine se01_02

	subroutine se01_03
	implicit none
	integer :: ans
15			continue
			write(*,*) "Inputs:"
			write(*,*) "	(1) Static Profile Options"
			write(*,*) "	(2) Pulsation Grid Options"
			write(*,*) "	(3) EigenValue Problem Solver Options"
			write(*,*) "	(4) Pulsation Equation Iterations Options"
			write(*,*) "	(5) Pulsation Equation Set Options"
			write(*,*) "	(6) EOS Options"
			write(*,*) "	(7) Shear Modulus Input Options"
			write(*,*) "	(8) Shear Modulus Formula Options"
			write(*,*) "	(9) Shear reduction factor"
			write(*,*) "	(10) Relativistic Corrrection factor"
			write(*,*) "	(0) Back"
			write(*,*) "=========================================================="
			read(*,*) ans
			write(*,*) "=========================================================="
			if (ans == 1) then
				write(*,*) "Static Profile Options:"
				write(*,*) "-Previous value: ", sp_opt
				read(*,*) sp_opt
			elseif (ans == 2) then
				write(*,*) "Pulsation Grid Options:"
				write(*,*) "-Previous value: ", pg_opt
				read(*,*) pg_opt
			elseif (ans == 3) then
				write(*,*) "EigenValue Problem Solver Options:"
				write(*,*) "-Previous value: ", ei_opt
				read(*,*) ei_opt			
			elseif (ans == 4) then
				write(*,*) "Pulsation Equation Iterations Options:"
				write(*,*) "-Previous value: ", pe_opt
				read(*,*) pe_opt
			elseif (ans == 5) then
				write(*,*) "Pulsation Equation Set Options:"
				write(*,*) "-Previous value: ", pes_opt
				read(*,*) pes_opt
			elseif (ans == 6) then
				write(*,*) "EOS Options:"
				write(*,*) "-Previous value: ", eos_opt
				read(*,*) eos_opt
			elseif (ans == 7) then
				write(*,*) "Shear Modulus Input Options:"
				write(*,*) "-Previous value: ", comp_opt
				read(*,*) comp_opt
			elseif (ans == 8) then
				write(*,*) "Shear Modulus Formula Options:"
				write(*,*) "-Previous value: ", shear_fcn_opt
				read(*,*) shear_fcn_opt
			elseif (ans == 9) then
				write(*,*) "Shear reduction factor:"
				write(*,*) "-Previous value: ", mu_red
				read(*,*) mu_red
			elseif (ans == 10) then
				write(*,*) "Relativistic Corrrection factor:"
				write(*,*) "-Previous value: ", rel
				read(*,*) rel
			elseif (ans == 0) then
				return
			else
				write(*,*) "invalid input"
			endif
			goto 15
		endsubroutine se01_03

		subroutine se01_04
		implicit none
			write(*,*) "Angular momentum number:"
			write(*,*) l_0
			write(*,*) "Central Density:"
			write(*,*) rho_0
			write(*,*) "Crust Core Transition Pressure:"
			write(*,*) P_t
			write(*,*) "Ocean Crust Transition Pressure:"
			write(*,*) P_g
			write(*,*) "Cutoff Pressure:"
			write(*,*) P_min
			write(*,*) "Polytrope: Density Jump Fraction"
			write(*,*) drho
			write(*,*) "Polytrope: Crust Core Transition Density"
			write(*,*) rho_t
			write(*,*) "Polytrope: Shear Speed throughout the crust"
			write(*,*) poly_cs
			write(*,*) "Quark Star: a4:"
			write(*,*) QS_a4
			write(*,*) "Quark Star: a2:"
			write(*,*) QS_a2
			write(*,*) "Quark Star: B_eff:"
			write(*,*) B_eff
			write(*,*) "Quark Star: Gap Parameter:"
			write(*,*) QS_Gap
		endsubroutine se01_04

		subroutine se01_05
		implicit none
		integer :: ans
16			continue
			write(*,*) "Inputs:"
			write(*,*) "	(1) Angular momentum number"
			write(*,*) "	(2) Central Density"
			write(*,*) "	(3) Crust Core Transition Pressure"
			write(*,*) "	(4) Ocean Crust Transition Pressure"
			write(*,*) "	(5) Cutoff Pressure"
			write(*,*) "	(6) Polytrope: Density Jump Fraction"
			write(*,*) "	(7) Polytrope: Crust Core Transition Density"
			write(*,*) "	(8) Polytrope: Shear Speed throughout the crust"
			write(*,*) "	(9) Quark Star: a4"
			write(*,*) "	(10) Quark Star: a2"
			write(*,*) "	(11) Quark Star: B_eff"
			write(*,*) "	(12) Quark Star: Gap Parameter"
			write(*,*) "	(0) Back"
			write(*,*) "=========================================================="
			read(*,*) ans
			write(*,*) "=========================================================="
			if (ans == 1) then
				write(*,*) "Angular momentum number:"
				write(*,*) "-Previous value: ", l_0
				read(*,*) l_0
			elseif (ans == 2) then
				write(*,*) "Central Density:"
				write(*,*) "-Previous value: ", rho_0
				read(*,*) rho_0
			elseif (ans == 3) then
				write(*,*) "Crust Core Transition Pressure:"
				write(*,*) "-Previous value: ", P_t
				read(*,*) P_t
			elseif (ans == 4) then
				write(*,*) "Ocean Crust Transition Pressure:"
				write(*,*) "-Previous value: ", P_g
				read(*,*) P_g
			elseif (ans == 5) then
				write(*,*) "Cutoff Pressure:"
				write(*,*) "-Previous value: ", P_min
				read(*,*) P_min
			elseif (ans == 6) then
				write(*,*) "Polytrope: Density Jump Fraction:"
				write(*,*) "-Previous value: ", drho
				read(*,*) drho
			elseif (ans == 7) then
				write(*,*) "Polytrope: Crust Core Transition Density:"
				write(*,*) "-Previous value: ", rho_t
				read(*,*) rho_t
			elseif (ans == 8) then
				write(*,*) "Polytrope: Shear Speed throughout the crust:"
				write(*,*) "-Previous value: ", poly_cs
				read(*,*) poly_cs
			elseif (ans == 9) then
				write(*,*) "Quark Star: a4:"
				write(*,*) "-Previous value: ", QS_a4
				read(*,*) QS_a4
			elseif (ans == 10) then
				write(*,*) "Quark Star: a2:"
				write(*,*) "-Previous value: ", QS_a2
				read(*,*) QS_a2
			elseif (ans == 11) then
				write(*,*) "Quark Star: B_eff:"
				write(*,*) "-Previous value: ", B_eff
				read(*,*) B_eff
			elseif (ans == 12) then
				write(*,*) "Quark Star: Gap Parameter:"
				write(*,*) "-Previous value: ", QS_Gap
				read(*,*) QS_Gap
			elseif (ans == 0) then
				return
			else
				write(*,*) "invalid input"
			endif
			goto 16
		endsubroutine se01_05

endmodule settings_opt01