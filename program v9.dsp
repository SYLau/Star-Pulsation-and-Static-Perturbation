# Microsoft Developer Studio Project File - Name="program v9" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=program v9 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "program v9.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "program v9.mak" CFG="program v9 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "program v9 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "program v9 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "program v9 - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0xc04 /d "NDEBUG"
# ADD RSC /l 0xc04 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "program v9 - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0xc04 /d "_DEBUG"
# ADD RSC /l 0xc04 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib imsl.lib imsls_err.lib imslmpistub.lib imsl.lib imsls_err.lib imslmpistub.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /pdb:none

!ENDIF 

# Begin Target

# Name "program v9 - Win32 Release"
# Name "program v9 - Win32 Debug"
# Begin Group "Library"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\library\Derivatives.f90
# End Source File
# Begin Source File

SOURCE=.\library\Eigenvalues_Eigenvectors.f90
# End Source File
# Begin Source File

SOURCE=.\Fitting.f90
DEP_F90_FITTI=\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\library\Format_IO.f90
# End Source File
# Begin Source File

SOURCE=.\library\FRead.f90
# End Source File
# Begin Source File

SOURCE=.\library\FWrite.f90
# End Source File
# Begin Source File

SOURCE=".\library\Hydrostatic Functions.f90"
# End Source File
# Begin Source File

SOURCE=.\library\Interpolation.f90
# End Source File
# Begin Source File

SOURCE=.\library\metric_var.f90
DEP_F90_METRI=\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\library\Minimum_Maximum.f90
# End Source File
# Begin Source File

SOURCE=.\library\RK4_Set.F90
# End Source File
# Begin Source File

SOURCE=.\library\Root_Finding.f90
# End Source File
# Begin Source File

SOURCE=.\library\Verify.f90
# End Source File
# End Group
# Begin Group "Notes"

# PROP Default_Filter ""
# Begin Source File

SOURCE=".\Development Notes.txt"
# End Source File
# Begin Source File

SOURCE=".\Input notes.txt"
# End Source File
# Begin Source File

SOURCE=.\introduction.txt
# End Source File
# End Group
# Begin Group "00 settings"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\settings_opt01.f90
DEP_F90_SETTI=\
	".\Debug\global_var.mod"\
	
# End Source File
# End Group
# Begin Group "01 stat_profile"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\EOS.F90
DEP_F90_EOS_F=\
	".\Debug\Format_IO.mod"\
	".\Debug\FRead.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Interpolation.mod"\
	".\Debug\verify.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_hyd.f90
DEP_F90_STAT_=\
	".\Debug\EOS.MOD"\
	".\Debug\Hydrostatic.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_opt01.f90
DEP_F90_STAT_P=\
	".\Debug\derivatives.mod"\
	".\Debug\EOS.MOD"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\sp_hyd.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_opt02.f90
DEP_F90_STAT_PR=\
	".\Debug\derivatives.mod"\
	".\Debug\EOS.MOD"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\sp_hyd.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_opt03.f90
DEP_F90_STAT_PRO=\
	".\Debug\derivatives.mod"\
	".\Debug\EOS.MOD"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\sp_hyd.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_opt04.f90
DEP_F90_STAT_PROF=\
	".\Debug\derivatives.mod"\
	".\Debug\EOS.MOD"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\sp_hyd.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_opt05.f90
DEP_F90_STAT_PROFI=\
	".\Debug\derivatives.mod"\
	".\Debug\EOS.MOD"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\sp_hyd.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\stat_profile_opt06.f90
DEP_F90_STAT_PROFIL=\
	".\Debug\derivatives.mod"\
	".\Debug\EOS.MOD"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\sp_hyd.mod"\
	
# End Source File
# End Group
# Begin Group "02 puls_grid"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\puls_grid_opt01.f90
DEP_F90_PULS_=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Interpolation.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_grid_opt02.f90
DEP_F90_PULS_G=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Interpolation.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_grid_opt03.f90
DEP_F90_PULS_GR=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Interpolation.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_grid_opt04.f90
DEP_F90_PULS_GRI=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Interpolation.mod"\
	
# End Source File
# End Group
# Begin Group "03 eigen_freq"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\eigen_freq_opt01.f90
DEP_F90_EIGEN=\
	".\Debug\eigen_freq_opt01_ex1.mod"\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_fcn.mod"\
	".\Debug\global_var.mod"\
	".\Debug\moI_opt01.mod"\
	".\Debug\puls_eqt_solve_opt01.mod"\
	".\Debug\puls_eqt_solve_opt02.mod"\
	".\Debug\puls_eqt_solve_opt04.mod"\
	".\Debug\puls_eqt_solve_opt07.mod"\
	".\Debug\puls_eqt_solve_opt08.mod"\
	".\Debug\Root_Finding.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\eigen_freq_opt01_ex1.f90
DEP_F90_EIGEN_=\
	".\Debug\Format_IO.mod"\
	".\Debug\FRead.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\eigen_freq_opt01_ex2.f90
DEP_F90_EIGEN_F=\
	".\Debug\Format_IO.mod"\
	".\Debug\FRead.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\eigen_freq_opt02.f90
DEP_F90_EIGEN_FR=\
	".\Debug\eigen_freq_opt02_ex1.mod"\
	".\Debug\Fitting.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\lo_eqt.mod"\
	".\Debug\love_opt01_P05.mod"\
	".\Debug\love_opt01_P08.mod"\
	".\Debug\Minimum_Maximum.mod"\
	".\Debug\moI_opt01.mod"\
	".\Debug\puls_eqt_solve_opt03.mod"\
	".\Debug\puls_eqt_solve_opt05.mod"\
	".\Debug\puls_eqt_solve_opt06.mod"\
	".\Debug\puls_eqt_solve_opt09.mod"\
	".\Debug\Root_Finding.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\eigen_freq_opt02_ex1.f90
DEP_F90_EIGEN_FRE=\
	".\Debug\Format_IO.mod"\
	".\Debug\FRead.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	
# End Source File
# End Group
# Begin Group "04 puls_eqt_solve"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\puls_eqt_solve_opt01.f90
DEP_F90_PULS_E=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt01.mod"\
	".\Debug\puls_eqt_set_opt02.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt02.f90
DEP_F90_PULS_EQ=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt01.mod"\
	".\Debug\puls_eqt_set_opt02.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt03.f90
DEP_F90_PULS_EQT=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt03.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt04.f90
DEP_F90_PULS_EQT_=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt01.mod"\
	".\Debug\puls_eqt_set_opt02.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt05.f90
DEP_F90_PULS_EQT_S=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt03.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt06.f90
DEP_F90_PULS_EQT_SO=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt03.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt07.f90
DEP_F90_PULS_EQT_SOL=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt04.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt08.f90
DEP_F90_PULS_EQT_SOLV=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt04.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_solve_opt09.f90
DEP_F90_PULS_EQT_SOLVE=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt05.mod"\
	".\Debug\puls_eqt_set_opt05_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# End Group
# Begin Group "05 puls_eqt_set"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\puls_eqt_set_opt01.f90
DEP_F90_PULS_EQT_SE=\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_set_opt02.f90
DEP_F90_PULS_EQT_SET=\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_set_opt03.f90
DEP_F90_PULS_EQT_SET_=\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_set_opt03_ex1.f90
DEP_F90_PULS_EQT_SET_O=\
	".\Debug\Eigenvalues_Eigenvectors.mod"\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_set_opt04.f90
DEP_F90_PULS_EQT_SET_OP=\
	".\Debug\global_var.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_set_opt05.f90
DEP_F90_PULS_EQT_SET_OPT=\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt05_ex1.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\puls_eqt_set_opt05_ex1.f90
DEP_F90_PULS_EQT_SET_OPT0=\
	".\Debug\global_var.mod"\
	
# End Source File
# End Group
# Begin Group "06 o_puls_soln"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\o_puls_soln_opt01.f90
DEP_F90_O_PUL=\
	".\Debug\global_var.mod"\
	
# End Source File
# End Group
# Begin Group "07 veri_soln"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\veri_soln_opt01.f90
DEP_F90_VERI_=\
	".\Debug\Format_IO.mod"\
	".\Debug\FRead.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Interpolation.mod"\
	".\Debug\puls_eqt_set_opt01.mod"\
	".\Debug\puls_eqt_set_opt02.mod"\
	
# End Source File
# End Group
# Begin Group "08 I_Love_Q"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\love_eqt.f90
DEP_F90_LOVE_=\
	".\Debug\EOS.MOD"\
	".\Debug\global_var.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\love_match.f90
DEP_F90_LOVE_M=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\lo_eqt.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\Root_Finding.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\love_match_P02.f90
DEP_F90_LOVE_MA=\
	".\Debug\global_var.mod"\
	".\Debug\lo_eqt.mod"\
	".\Debug\lo_match.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\love_opt01_P01.f90
DEP_F90_LOVE_O=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\lo_match.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\love_opt01_P05.f90
DEP_F90_LOVE_OP=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\lo_eqt.mod"\
	".\Debug\lo_match.mod"\
	".\Debug\lo_match_P02.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	".\Debug\Root_Finding.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\love_opt01_P07.f90
DEP_F90_LOVE_OPT=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\lo_match.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	{$(INCLUDE)}"lin_sol_gen_int.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\love_opt01_P08.f90
DEP_F90_LOVE_OPT0=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\lo_eqt.mod"\
	".\Debug\lo_match.mod"\
	".\Debug\love_opt01_P05.mod"\
	".\Debug\puls_eqt_set_opt03_ex1.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\moI_opt01.f90
DEP_F90_MOI_O=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\RK4_Set.mod"\
	
# End Source File
# End Group
# Begin Source File

SOURCE=.\global_fcn.f90
DEP_F90_GLOBA=\
	".\Debug\Format_IO.mod"\
	".\Debug\FWrite.mod"\
	".\Debug\global_var.mod"\
	".\Debug\Hydrostatic.mod"\
	".\Debug\puls_eqt_set_opt02.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\global_var.f90
# End Source File
# Begin Source File

SOURCE=.\main.f90
DEP_F90_MAIN_=\
	".\Debug\eigen_freq_opt01.mod"\
	".\Debug\eigen_freq_opt02.mod"\
	".\Debug\global_var.mod"\
	".\Debug\puls_grid_opt01.mod"\
	".\Debug\puls_grid_opt02.mod"\
	".\Debug\puls_grid_opt03.mod"\
	".\Debug\puls_grid_opt04.mod"\
	".\Debug\RK4_Set.mod"\
	".\Debug\settings_opt01.mod"\
	".\Debug\stat_profile_opt01.mod"\
	".\Debug\stat_profile_opt02.mod"\
	".\Debug\stat_profile_opt03.mod"\
	".\Debug\stat_profile_opt04.mod"\
	".\Debug\stat_profile_opt05.mod"\
	".\Debug\stat_profile_opt06.mod"\
	".\Debug\veri_soln_opt01.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\program v9.f90"
DEP_F90_PROGR=\
	".\Debug\main.mod"\
	
# End Source File
# End Target
# End Project
