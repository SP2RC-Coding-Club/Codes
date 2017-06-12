FUNCTION complex_dbeselj_deriv,X,n,p
  if n eq 0 then return, - complex_dbeselj(X,1,p,0.0,!PI)
  if n gt 0 then return, 0.5d*(complex_dbeselj(X,n-1,p,0.0,!PI)-complex_dbeselj(X,n+1,p,0.0,!PI))
; http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/introductions/Bessels/05/
; Checked by DY 26 Jan 2015
end