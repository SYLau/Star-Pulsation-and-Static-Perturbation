module vf_plots
use dfwin

contains

subroutine graphicsmode( )
 use dflib
 logical modestatus
 integer(2) maxx, maxy
 type (windowconfig) myscreen
 common maxx, maxy
 ! Set highest resolution graphics mode.
 myscreen.numxpixels=-1
 myscreen.numypixels=-1
 myscreen.numtextcols=-1
 myscreen.numtextrows=-1
 myscreen.numcolors=-1
 myscreen.fontsize=-1
 myscreen.title = " "C ! blank
 modestatus=SETWINDOWCONFIG(myscreen)
 ! Determine the maximum dimensions.
 modestatus=GETWINDOWCONFIG(myscreen)
 maxx=myscreen.numxpixels - 1
 maxy=myscreen.numypixels - 1
endsubroutine graphicsmode




endmodule vf_plots