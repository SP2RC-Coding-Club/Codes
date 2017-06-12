FUNCTION  complex_dbeselk,X,n,p,A,B
	;http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html Eq. (5)
	; This functions solves for intergral using the simpson method for real and complex
	; X is for the complex input
	; n is the integer order
	; p for number of steps
	; A define lower limit of integration.
	; B define upper limit of integration.
	; Note for defining A and B in the function make sure the limilts end with d eg(4.0d), otherwise delta_x = 0.
	; Define lower limit of integration:
	checker = p mod 2
	if checker eq 1 then begin
		; n must be even for the Simpson rule
		print, 'input n must be even'
	endif
	if checker eq 0 then begin
		delta_x = (B-A)/p ; Step size.
		f_re = make_array(p+1,/DOUBLE) ; Makes array to fill.
		f_Im = make_array(p+1,/DOUBLE) ; Make array to fill.
		complex_f = make_array(p+1,/DOUBLE) ; Make array to fill.
		for k= 0,p do begin
			step = A+k*delta_x ; step is the quanty which you wish to intergrate with respect to
			;function to intergrate is K_v(z) = ((gamma(v+0.5)*(2*z)^v)/sqrt(!PI))*int(cos(t)/((t^2+z^2)^(v+0.5)))
			f_c = 1/((step^2+X^2)^(n+0.5d))
			; This is the function which we wish to integrate and will split into real and imaginary parts.
			f_re_Im = cos(step)*CONJ(f_c)/(((step^2+X^2)^(n+0.5d))*CONJ(f_c))
			f_re[k] = real_part(f_re_Im) ; Real part.
			f_Im[k] = imaginary(f_re_Im) ; Imaginary part.
		endfor

		int_f_re = 0
		int_f_Im = 0
		; This part is where we apply the composite Simpson rule for both the real and imaginary values.
		for j = 2, p,2 do begin
			int_f_re = int_f_re+ (f_re[j-2]+4*f_re[j-1]+f_re[j])
			int_f_Im = int_f_Im+ (f_Im[j-2]+4*f_Im[j-1]+f_Im[j])
		endfor
	endif

	full_int_f_re = (delta_x/3.0d)*int_f_re
	full_int_f_Im = (delta_x/3.0d)*int_f_Im
	; Adding the real and imaginary parts together.
	full_int_f_re_Im = dcomplex(full_int_f_re,full_int_f_Im)
	; Multiplying by the common factor outside the intergral to obtain the bessel K function.
	k_n = ((gamma(n+0.5d)*(2.0d*X)^n)/sqrt(!PI))*full_int_f_re_Im
	RETURN, k_n
end