module derivatives
contains 

	function dydx_2p(dx,dy)
	REAL(8)	:: dx,dy
	REAL(8)	:: dydx_2p
		dydx_2p = dy/dx
	endfunction dydx_2p

endmodule derivatives