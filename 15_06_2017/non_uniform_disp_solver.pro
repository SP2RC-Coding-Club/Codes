;+
; :Author: Fionnlagh Mackenzie Dover
;-
;-------- plasma parameters -----------
; Converstions and contants
Mega=1.0d6 ; conversion to mega
kilo=1.0d3 ; conversion to kilo
Kb=1.380658d-23 ;m2 kg s-2 K-1, Boltzman constant
mp=1.67262178d-27; kg, proton mass
mu0=4.0d-7*!DPI ;[H·m−1] or [N·A−2]
gama=5.d/3.d

;internal
a=4.0d*Mega;; 4.0d*Mega; loop radius  (Mm)
;Lz=100.d*Mega; loop length (Mm)
Lz=25.d*a; loop length (Mm)
Ltr=1.d*a ;; Length of transition region (Mm)
Te0=0.8*Mega; internal plasma temperature [MK]
ne0=1.5d15 ; m-3
rho0=ne0*mp;internal plasma density [kg m-3]
B0=0.0015d; internal magnetic field [Tesla]
Cs0=sqrt(2*gama*kb/mp*Te0) ; m s-1
Va0=B0/sqrt(mu0*rho0) ; m s-1
Ct0=Cs0*Va0/sqrt(Va0^2+Cs0^2)  ; m s-1
beta0=2*Cs0^2/(gama*Va0^2)
p0=B0^2/(2*mu0)*beta0 ; internal gas pressure [Pa]

;external
delrho=5.d; rho0/rho1, density ratio
delTe=1.5d ; Te0/Te1
;; delTe=Te0/Te1=rho1/rho0*p0/p1=beta0/beta1*delB^2
beta1=(delTe*delrho*(1+1/beta0)-1)^(-1) ;
delB=sqrt((1+beta1)/(1+beta0)) ; B0/B1
delbeta=beta0/beta1
Te1=Te0/delTe ;; external plasma temperature [MK]
B1=B0/delB  ; external magnetic field [Tesla]
ne1=ne0/delrho
rho1=ne1*mp ; external plasma density [kg m-3]
p1=B1^2/(2*mu0)*beta1 ; external gas pressure [Pa]
Cs1=sqrt(2*gama*kb/mp*Te1) ; external sound speed [m s-1]
Va1=B1/sqrt(mu0*rho1)
Ct1=Cs1*Va1/sqrt(Cs1^2+Va1^2)
Ck=sqrt((rho0*Va0^2+rho1*Va1^2)/(rho0+rho1))

;; wave mode setup
;; cos(mi*theta)*cos(ni*pi*z/L)*cos(2*pi*t/P)
;; mi: azimuthal mode number
;; ni: longitudinal mode number
;;
;; solve mi=1, ni=1
;; fundamental fast kink mode Va0<Vp<Ck;

mi=1 ; cos (m*theta) azimuthal mode number
ni=1 ; longitudinal mode number
Kz=ni*!dpi/Lz ; m-1
Ka=Kz*a; m-1*m
nvp=1000000
Vpmin=(Ck-Va0)*0.8+Va0
Vp=Vpmin+(1+dindgen(nvp))*(Ck-Vpmin)/(nvp+1)
k0a=Ka*sqrt(-(Cs0^2-Vp^2)*(Va0^2-Vp^2)/((Cs0^2+Va0^2)*(Ct0^2-Vp^2))) ; k0a^2>0 unit [ka]; note the minus sign
K1a=Ka*sqrt((Cs1^2-Vp^2)*(Va1^2-Vp^2)/((Cs1^2+Va1^2)*(Ct1^2-Vp^2))) ; K1a^2>0 unit [ka]
;;;Homogenous dispersion relation
disp=delrho*((Va0/Va0)^2-(Vp/Va0)^2)*k1a*dbeselk_deriv(k1a,mi)/dbeselk(k1a,mi)$
  -((Va1/Va0)^2-(Vp/Va0)^2)*k0a*dbeselj_deriv(k0a,mi)/dbeselj(k0a,mi)

disp_min=min(abs(disp),sol)
vp0 = vp(sol) ; phase speed w0/kz
k0a=k0a[sol]
k0=k0a/a ; kz internal
k1a=k1a[sol]
k1=k1a/a ; kz external
vp0=vp0[0]; m/s
w0=vp0*kz ; rad/s
wc0=Cs0*kz; acoustic frequency
wc1=Cs1*kz
wa0=Va0*kz; alfven frequency
wa1=Va1*kz
wt0=Ct0*kz; cusp frequency
wt1=Ct1*kz

;To make unitless quantities 
unit_length = a ;[m]
Unit_density = rho1 ;[kg m-3]
unit_velocity = va0 ;[m s-1]
unit_B = B1 ;[Tesla]
unit_Temp = Te0 ;[K]
unit_time = unit_length/unit_velocity ; [s]
unit_pressure = unit_velocity^2*Unit_density ;[pa] 

Ltr_norm = ltr/unit_length
rho1_norm = rho1/Unit_density
rho0_norm = rho0/Unit_density
Va0_norm = Va0/unit_velocity
Va1_norm = Va1/unit_velocity
a_norm = a/unit_length
B0_norm = B0/unit_B 
kz_norm = kz*unit_length
w_norm = 1/unit_time;va0_norm*kz_norm
   
;rho_k_tot_list= DCOMPLEXARR(N, 1)

wk = Ck*kz ; kink frequency

i = dcomplex(0,1) ; gives complex number
N=52.0d ; Upper limmit tke for sumtion terms. This currently is the highest limmit that can be used before certian quanities in the code become infinte
; In Soler et al. (2013) they use N = 51, as this gives convergence for W0

; matrix = MAKE_ARRAY(rows,columns,precision)
rho_k =  DCOMPLEXARR(N, 1)
rho_k_test =  DCOMPLEXARR(N, 1)
rho_k_tot_check = DCOMPLEXARR(N,1)


w0_1 = dcomplex(0.164967079252836d,-3.263644161202017E-002)*w_norm;*wk ;[s-1] ;once work this will be iterated over the Re and Im part
w0_1_norm = w0_1/w_norm;w0_1*unit_time 

ra_norm=(Ltr_norm/!DPI)*asin((1+rho1_norm/rho0_norm-2*(Va0_norm*kz_norm/w0_1_norm)^2)/(1-rho1_norm/rho0_norm))+a_norm ;works :)
ra=((Ltr/!DPI)*asin((1+rho1/rho0-2*(Va0*kz/w0_1)^2)/(1-rho1/rho0))+a);/a
arg = ((rho0+rho1)/(rho0-rho1))-(2.0d/(mu0*(rho0-rho1)))*(kz*B0/w0_1)^2
arg_norm = ((rho0_norm+rho1_norm)/(rho0_norm-rho1_norm))-(2.0d/((rho0_norm-rho1_norm)))*(kz_norm/w0_1_norm)^2*((B0^2/mu0)/unit_pressure) 
ra_1= (Ltr/!DPI)*asin(arg)+a
ra_1_norm= (Ltr_norm/!DPI)*asin(arg_norm)+a_norm
;ra_norm = 1.068-0.0005*i

rho_k[0] = ((1/mu0)*(kz*(B0)/w0_1)^2)/Unit_density ; [kg/m^3] = unit_density 
;rho_k_0_check = ((1/mu0)*(kz*(B0)/w0_1)^2)/Unit_density
; part of the taylor series expansion. See equation (21) on Soler et al (2013) and the density profile from eq (63).
FOR k=1,N-1 DO BEGIN
  modular = k MOD 2
  if modular EQ 0 then begin
    rho_k[k] = (1/FACTORIAL(k))*(rho0_norm/2.0d)*(1-rho1_norm/rho0_norm)*((!DPI/Ltr_norm)^k)*(i^(k+2))*sin(!DPI/Ltr_norm*(ra_norm-a_norm)) ; This is for even solutions and no dim
  ENDIF
  if modular EQ 1 then begin
    rho_k[k] = (1/FACTORIAL(k))*(rho0_norm/2.0d)*(1-rho1_norm/rho0_norm)*(i^(k+1))*(!DPI/Ltr_norm)^k*cos(!DPI/Ltr_norm*(ra_norm-a_norm)) ; This is for odd solutions and no dim
  ENDIF
  ;print, k, rho_k[k],  rho_k_test[k]
ENDFOR

;========================================================================================
;DENSITY CHECK
;========================================================================================
;a_i_plot = (a_norm-ltr_norm/2)*((dindgen(N+1))/N)
;a_tr_plot =  (a_norm-ltr_norm/2)+(Ltr_norm)*((dindgen(N))/N)
;a_e_plot = a_norm+ltr_norm/2+(a_norm)*((dindgen(N+1))/N)
;a_e_jump = a+a_i_plot ; for step function
;a_ra =  MAKE_ARRAY(1, VALUE = ra_norm, /DOUBLE)
;rho_k_tot = 0
;;; This part is equation (20) in Soler et al (2013)
;for jj=0, N-1 do begin
;  rho_k_tot = 0
;  pos = a_tr_plot[jj]
;  ;print, pos
;  FOR k=0,N-1 DO BEGIN
;    rho_k_tot= rho_k_tot + rho_k[k]*(pos-ra_norm)^k
;  ENDFOR
;  ;print, pos, rho_k_tot*Unit_density
;  rho_k_tot_check[jj] = rho_k_tot
;endfor

;rho_int = MAKE_ARRAY(N, VALUE = rho0_norm, /DOUBLE)
;rho_ext = MAKE_ARRAY(N, VALUE = rho1_norm, /DOUBLE)
;rho_tr = (rho0_norm/2)*((1+rho1_norm/rho0_norm)-(1-rho1_norm/rho0_norm)*sin(!PI/Ltr_norm*(a_tr_plot-a_norm)))

;rho_transitional = [rho_int*Unit_density,rho_tr*Unit_density,rho_ext*Unit_density]
;a_transitional = [a_i_plot,a_tr_plot,a_e_plot]

;cgPlot, a_transitional, rho_transitional, Color='red7',thick = 4,XRange=[0.25,1.75], YRange=[min(rho_transitional)-min(rho_transitional)/10,max(rho_transitional)+max(rho_transitional)/10], YStyle=1, $
;  XTitle='r/R', YTitle='Density [kg m^-3]',CHARSIZE = 2.5, CHARTHICK = 1, /Window
;-----------------------------------------
;------------------------------------------
;Plots just transiot reg
;------------------------------------------
;cgPlot, a_tr_plot, rho_tr*Unit_density, Color='red7',thick = 4,XRange=[min(a_tr_plot),max(a_tr_plot)], YRange=[min(rho_tr*Unit_density)-min(rho_tr*Unit_density)/10,max(rho_tr*Unit_density)+max(rho_tr*Unit_density)/10], YStyle=1, $
;    XTitle='r/R', YTitle='Density [kg m^-3]',CHARSIZE = 2.5, CHARTHICK = 1, /Window
;;cgPlot, a_tr_plot, rho_k_tot_check*Unit_density, Color='blue',thick =4, YStyle=1, CHARSIZE=2.5, /Overplot, /Window
;;-----------------------------------------
;;To check density profile
;;delta_rho = (rho_k_tot_check*Unit_density-rho_tr*Unit_density)*100/rho_tr*Unit_density ; Delta rho is very small therefore the profile is calculated correctly 
;;print, delta_rho
;;========================================================================================
;END OF Density CHECK
;========================================================================================

;;========================================================================================
;CALC OF A_K AND S_K COEFFICIENTS 
;========================================================================================
a_k = DCOMPLEXARR(N-1, 1)
s_k = DCOMPLEXARR(N-1, 1)
s_k_2 = DCOMPLEXARR(N-1, 1)
;-----------------------------------------------------------------------------
;Expansion coefficient a_k and s_k given in the apendix of Soler et al. (2013).
; Eq. (A1) in Soler et al. (2013).
a_k[0] = 1.0d
; Eq. (A2) in Soler et al. (2013).
a_k[1] = -((2.0d*rho_k[1]-2.0d*rho_k[2]*ra_norm)/(3.0d*rho_k[1]*ra_norm))*a_k[0]
; Eq. (A3) in Soler et al. (2013).
a_k[2] = -(9.0d*rho_k[1]*ra_norm*a_k[1]+(2.0d*rho_k[1]-2.0d*rho_k[2]*ra_norm-4.0d*rho_k[3]*ra_norm^2-mi^2*rho_k[1])*a_k[0])/(8.0d*rho_k[1]*ra_norm^2.0d)
; Eq. (A4) in Soler et al. (2013).
; the 2nd last term should be ra.
; NOTE:(w0_1^2.0d*ra^2.0d/(B0^2/mu0)) = [s^-2 m^2/(T^2/(H m^-1))] = [m^3/kg] = Unit_density^-1
a_k[3] = (-1.0d/(15.0d*rho_k[1]*ra_norm^2))*((4.0d*rho_k[2]*ra_norm^2+20.0d*rho_k[1]*ra_norm)*a_k[2]+(-3.0d*rho_k[3]*ra_norm^2+3.0d*rho_k[2]*ra_norm+6.0d*rho_k[1]-mi^2*rho_k[1])*a_k[1] $
  +(-6.0d*rho_k[4]*ra_norm^2-6.0d*rho_k[3]*ra_norm+(w0_1_norm^2*ra_norm^2/((B0^2/mu0)/unit_pressure))*rho_k[1]^2-mi^2*rho_k[2])*a_k[0])
;Calculation of the sumation terms in Eq. (A5) in Soler et al. (2013).
ak4_p1 = 0.0d
ak4_p1_2 = 0.0d
ak4_p2 = 0.0d
ak4_p3 = 0.0d
FOR j=0.0d,3.0d DO BEGIN
  ;print,'j', j
  ak4_p1 = ak4_p1+(j+2.0d)*(2.0d*j-4.0d)*ra_norm^2*rho_k[5-j]*a_k[j] 
  ak4_p1_2 = ak4_p1_2+(j+2.0d)*(4.0d*j-5.0d)*ra_norm*rho_k[4-j]*a_k[j]
  if j le 2 then begin
    ak4_p2 = ak4_p2+((j+2.0d)*(2.0d*j-1.0d)-mi^2)*rho_k[3-j]*a_k[j]
  ENDIF
  if j le 1 then begin
    FOR l= 0,1-j DO BEGIN
      ;print,'l', l
      ak4_p3 = ak4_p3 + (w0_1_norm^2*ra_norm^2/((B0^2/mu0)/unit_pressure))*rho_k[l+1]*rho_k[2-j-l]*a_k[j]+2.0d*(w0_1_norm^2*ra_norm/((B0^2/mu0)/unit_pressure))*rho_k[1]^2*a_k[0]
    ENDFOR
  ENDIF
ENDFOR
; Eq. (A5) in Soler et al. (2013)
a_k[4] = -(1.0d/(24.0d*rho_k[1]*ra_norm^2))*(ak4_p1+ak4_p1_2+ak4_p2+ak4_p3)
;;==========
;;test of the summations. Produces results above
;;==========
;a_4_1 = 0.0d
;for j=0.0d,3.0d do begin
;  a_4_1 = a_4_1 + (j+2.0d)*(2.0d*j-4.0d)*ra_norm^2*rho_k[5-j]*a_k[j]
;endfor
;a_4_2=0.0d
;for j=0.0d,3.0d do begin
;  a_4_2 = a_4_2 + (j+2.0d)*(4.0d*j-5.0d)*ra_norm*rho_k[4-j]*a_k[j]
;endfor
;a_4_3=0.0d
;for j=0.0d,2.0d do begin
;  a_4_3= a_4_3 + ((j+2.0d)*(2.d0*j-1.0d)-mi^2)*rho_k[3-j]*a_k[j] 
;endfor
;a_4_4=0.0d
;for j=0.0d,1.0d do begin
;  for l=0.0d,1-j do begin
;    a_4_4= a_4_4 + (w0_1_norm^2*ra_norm^2/((B0^2/mu0)/unit_pressure))*rho_k[l+1]*rho_k[2-j-l]*a_k[j]+2.0d*(w0_1_norm^2*ra_norm/((B0^2/mu0)/unit_pressure))*rho_k[1]^2*a_k[0]
;  endfor
;endfor
;a_4 = -(1.0d/(24.0d*rho_k[1]*ra_norm^2))*(a_4_1+a_4_2+a_4_3+a_4_4)
;;==========
;;end
;;==========

for k=5.0d,N-2 do begin
  akN_p1 = 0.0d
  akN_p2 = 0.0d
  akN_p3 = 0.0d
  akN_p4 = 0.0d
  akN_p5 = 0.0d
  akN_p6 = 0.0d
  for j=0.0d, k-1.0d do begin
    akN_p1 = akN_p1 + (j+2.0d)*(2.d0*j-k)*ra_norm^2*rho_k[k-j+1]*a_k[j]
    akN_p2 = akN_p2 + (j+2.d0)*(4.d0*j-2.d0*k+3.d0)*ra_norm*rho_k[k-j]*a_k[j]
    if j le k-2 then begin
      akN_p3 = akN_p3+((j+2.d0)*(2.d0*j-k+3.d0)-mi^2)*rho_k[k-j-1]*a_k[j]
    endif
    if j le k-3 then begin
      for l=0, k-j-3 do begin
        akN_p4 = akN_p4+((w0_1_norm^2*ra_norm^2)/((B0^2/mu0)/unit_pressure))*(rho_k[l+1]*rho_k[k-j-l-2]*a_k[j])
      endfor
    endif 
    if j le k-4.0d then begin
      for l=0, k-j-4.0d do begin
        akN_p5 = akN_p5 + 2.0d*(w0_1_norm^2*ra_norm/((B0^2/mu0)/unit_pressure))*rho_k[l+1]*rho_k[k-j-l-3]*a_k[j]
      endfor
    endif
    if j le k-5.0d then begin
      for l=0,k-j-5 do begin
        akN_p6= akN_p6 + (w0_1_norm^2/((B0^2/mu0)/unit_pressure))*(rho_k[l+1]*rho_k[k-j-l-4]*a_k[j])
      endfor
    endif
  endfor
  a_k[k] = -(1.0d/(k*(k+2.0d)*rho_k[1]*ra_norm^2))*(akN_p1+akN_p2+akN_p3+akN_p4+akN_p5+akN_p6)
endfor

C = mi^2/(2.0d*ra_norm^2) ; Coupling constant given by Eq. (26) Soler et al. (2013).
;Eq. (A7) in Soler et al. (2013).
s_k[0] = 1.0d
;Eq. (A8) in Soler et al. (2013).
s_k[1] = 0.0d
;Eq. (A9) in Soler et al. (2013).
s_k[2] = 0.0d
;Eq. (A10) in Soler et al. (2013).
s_k[3] = (1.0d/(3.0d*rho_k[1]*ra_norm^2))*((mi^2*rho_k[2]-w0_1_norm^2*ra_norm^2*rho_k[1]^2/((B0^2/mu0)/unit_pressure))*s_k[0]-C*(4.0d*ra_norm^2*rho_k[1]*a_k[1]+(ra_norm^2*rho_k[2]+5.0d*ra_norm*rho_k[1])*a_k[0]))
s_k[4]=-(1.0d/(8.0d*rho_k[1]*ra_norm^2))*(9.0d*rho_k[1]*ra_norm*s_k[3]+(2.0d*w0_1_norm^2*ra_norm*rho_k[1]^2/((B0^2/mu0)/unit_pressure)+2.0d*w0_1_norm^2*ra_norm^2*rho_k[1]*rho_k[2]/((B0^2.0d/mu0)/unit_pressure)-mi^2*rho_k[3])*s_k[0] $
  +C*(6.0d*ra_norm^2*rho_k[1]*a_k[2]+3.0d*ra_norm^2*rho_k[2]*a_k[1]+9.0d*ra_norm*rho_k[1]*a_k[1]+3.0d*(ra_norm*rho_k[2]+rho_k[1])*a_k[0]))
;Calculation of the sumation terms in Eq. (A12) in Soler et al. (2013)

for k=5.0d,N-2 do begin
  skN_p1 = 0.0d
  skN_p2 = 0.0d
  skN_p3 = 0.0d
  skN_p4 = 0.0d
  skN_p5 = 0.0d
  skN_p6 = 0.0d
  for j=0.0d,k-1.0d do begin
    skN_p1 = skN_p1 + j*(2.0d*j-k-2.0d)*ra_norm^2*rho_k[k-j+1]*s_k[j]
    skN_p2 = skN_p2 + j*(4.0d*j-2.0d*k-1.0d)*ra_norm*rho_k[k-j]*s_k[j]
    if j le k-2.0d then begin
      skN_p3 = skN_p3 + (j*(2.0d*j-k+1.0d)-mi^2)*rho_k[k-j-1]*s_k[j]+C*(3.0d*j-k+4.0d)*ra_norm^2*rho_k[k-j-1]*a_k[j]
    endif
    if j le k-3.0d then begin
      for l=0.0d,k-j-3.0d do begin
        skN_p4 = skN_p4 + ((w0_1_norm^2*ra_norm^2)/((B0^2/mu0)/unit_pressure))*rho_k[l+1]*rho_k[k-j-l-2]*s_k[j]+C*(6.0d*j-2.0d*k+11.0d)*ra_norm*rho_k[k-j-2]*a_k[j]
      endfor
    endif
    if j le k-4 then begin
      for l=0.0d, k-j-4.0d do begin
        skN_p5 = skN_p5 + 2.0d*((w0_1_norm^2*ra_norm)/((B0^2/mu0)/unit_pressure))*rho_k[l+1]*rho_k[k-j-l-3]*s_k[j]+C*(3.0d*j-k+7.0d)*rho_k[k-j-3]*a_k[j]
      endfor
    endif
    if j le k-5 then begin
      for l=0.0d,k-j-5.0d do begin
        skN_p6 = skN_p6 + ((w0_1_norm^2)/((B0^2/mu0)/unit_pressure))*rho_k[l+1]*rho_k[k-j-l-4]*s_k[j]
      endfor
    endif
  endfor
  s_k[k] = -(1.0d/(k*(k-2.0d)*rho_k[1]*ra_norm^2))*(skN_p1+skN_p2+skN_p3+skN_p4+skN_p5+skN_p6)
endfor


;Eq. (A7) in Soler et al. (2013).
s_k_2[0] = 1.0d
;Eq. (A8) in Soler et al. (2013).
s_k_2[1] = 0.0d
;Eq. (A9) in Soler et al. (2013).
s_k_2[2] = 0.0d
;Eq. (A10) in Soler et al. (2013).
s_k_2[3] = (1.0d/(3.0d*rho_k[1]*ra_norm^2))*((mi^2*rho_k[2]-w0_1_norm^2*ra_norm^2*rho_k[1]^2/((B0^2/mu0)/unit_pressure))*s_k_2[0]-C*(4.0d*ra_norm^2*rho_k[1]*a_k[1]+(ra_norm^2*rho_k[2]+5.0d*ra_norm*rho_k[1])*a_k[0]))

s_k_2[4]=-(1.0d/(8.0d*rho_k[1]*ra_norm^2))*(9.0d*rho_k[1]*ra_norm*s_k_2[3]+(2.0d*w0_1_norm^2*ra_norm*rho_k[1]^2/((B0^2/mu0)/unit_pressure)+2.0d*w0_1_norm^2*ra_norm^2*rho_k[1]*rho_k[2]/((B0^2.0d/mu0)/unit_pressure)-mi^2*rho_k[3])*s_k_2[0] $
  +C*(6.0d*ra_norm^2*rho_k[1]*a_k[2]+3.0d*ra_norm^2*rho_k[2]*a_k[1]+9.0d*ra_norm*rho_k[1]*a_k[1]+3.0d*(ra_norm*rho_k[2]+rho_k[1])*a_k[0]))
;Calculation of the sumation terms in Eq. (A12) in Soler et al. (2013)
mag_p_norm = ((B0^2/mu0)/unit_pressure)
for k=5.0d,N-2 do begin
  sk_p1 = 0
  sk_p2 = 0
  sk_p3 = 0
  sk_p4 = 0
  sk_p5 = 0
  sk_p6 = 0
  s_k_2[k] = 0
  for j=0.0d, k-1.0d do begin
    sk_p1 = sk_p1 + j*(2.0d*j-k-2.0d)*ra_norm^2*rho_K[k-j+1]*s_k_2[j]
    sk_p2 = sk_p2 + j*(4.0d*j-2.0d*k-1.0d)*ra_norm*rho_k[k-j]*s_k_2[j]
    if j le k-2.0d then begin
      sk_p3 = sk_p3 + (j*(2.0d*j-k+1.0d)-mi^2)*rho_k[k-j-1]*s_k_2[j]+C*(3.0d*j-k+4.0d)*ra_norm^2*rho_k[k-j-1]*a_k[j]
    endif
    if j le k-3.0d then begin
      for l=0.0d,k-j-3 do begin
        sk_p4 = sk_p4 + ((w0_1_norm*ra_norm)^2/(mag_p_norm))*rho_k[l+1]*rho_k[k-j-l-2]*s_k_2[j]+C*(6.0d*j-2.0d*k+11.0d)*ra_norm*rho_k[k-j-2]*a_k[j]
      endfor
    endif
    if j le k-4 then begin
      for l=0.0d, k-j-4.0d do begin
        sk_p5 = sk_p5 + 2.0d*(((w0_1_norm^2)*ra_norm)/mag_p_norm)*rho_k[l+1]*rho_k[k-j-l-3]*s_k_2[j]+C*(3.0d*j-k+7.0d)*rho_k[k-j-3]*a_k[j]
      endfor
    endif
    if j le k-5.0d then begin
      for l=0.0d, k-j-5.0d do begin
        sk_p6 = sk_p6 + (w0_1_norm^2/mag_p_norm)*rho_k[l+1]*rho_k[k-j-l-4]*s_k_2[j]
      endfor
    endif
  endfor
  s_k_2[k] = -(1/(k*(k-2.0d)*rho_k[1]*ra_norm^2))*(sk_p1+sk_p2+sk_p3+sk_p4+sk_p5+sk_p6)
endfor

;comment this out if your interested in running the code
;====================================================
;A_K AND S_K CHECKER!!!!!!!!!
;====================================================
lun1 = 122
lun2 = 100
openr, lun1, '/home/fionnlagh/Msc_project/Default/ak.dat'
header = strarr(106)
readf,lun1,header
close,lun1

openr, lun2, '/home/fionnlagh/Msc_project/Default/sk.dat'
header2 = strarr(106)
readf,lun2,header2
close,lun2

a_k_soler = DCOMPLEXARR(51,1)
s_k_soler = DCOMPLEXARR(51,1)
for p=0,101 do begin
  modular = p MOD 2
  if modular EQ 0 then begin
    a_k_soler[p/2] = dcomplex(header[p],header[p+1])
    s_k_soler[p/2] = dcomplex(header2[p],header2[p+1])
  ENDIF
endfor

rel_diff_a_k = 100*(a_k-a_k_soler)/a_k_soler
rel_diff_s_k = 100*(s_k-s_k_soler)/s_k_soler
rel_diff_s_k[1:2] = 0
index_plot = INDGEN(51)

;;;;Color list:http://www.idlcoyote.com/idldoc/cg/cgcolor.html
;;;;
;;;----------------------------
;;Relative diff % comparison
;cgPlot, index_plot, (rel_diff_a_k), PSym=-15, Color='cadet blue', thick = 2, YRange=[0,30], XRange=[0,52], $
;  LineStyle=0, YStyle=1, XTitle='k', YTitle='Relative difference %', CHARSIZE = 2.5, CHARTHICK = 1, /Window
;cgPlot, index_plot, imaginary(rel_diff_a_k), PSym=-16,Color='Dark red', thick = 2, $
;  LineStyle=1, /Overplot, /AddCmd
;cgLegend, Title=['Re((DIY_a_k- Soler_a_k)/Soler_a_k))', 'Im((DIY_a_k- Soler_a_k)/Soler_a_k))'], CHARSIZE = 1.5, CHARTHICK = 1, PSym=[-15,-16], $
;  LineStyle=[0,1], Color=['cadet blue','Dark red'], Location=[2,29], $
;  /Data, /AddCmd
;;  
;cgPlot, index_plot, (rel_diff_s_k), PSym=-15, Color='cadet blue', thick = 2, YRange=[-7500,500], XRange=[0,52], $
;    LineStyle=0, YStyle=1, XTitle='k', YTitle='Relative difference %', CHARSIZE = 2.5, CHARTHICK = 1, /Window
;cgPlot, index_plot, imaginary(rel_diff_s_k), PSym=-16,Color='Dark red', thick = 2, $
;    LineStyle=1, /Overplot, /AddCmd
;cgLegend, Title=['Re((DIY_s_k- Soler_s_k)/Soler_s_k))', 'Im((DIY_s_k- Soler_s_k)/Soler_s_k))'], CHARSIZE = 1.5, CHARTHICK = 1, PSym=[-15,-16], $
;    LineStyle=[0,1], Color=['cadet blue','Dark red'], Location=[2,-6000], $
;    /Data, /AddCmd
;;;  
;;;----------------------------
;;Re(A_k) comparison
;cgPlot, index_plot, (a_k_soler), PSym=-15, Color='cadet blue', thick = 2, YRange=[-0.8,1], XRange=[0,52], $
;  LineStyle=0, YStyle=1, XTitle='k', YTitle='Re(A_K)', CHARSIZE = 2.5, CHARTHICK = 1, /Window
;cgPlot, index_plot, (a_k), PSym=-16,Color='Dark Red', thick = 2, $
;  LineStyle=2, /Overplot, /AddCmd
;cgLegend, Title=['Re(Soler data)', 'Re(DIY data)' ], CHARSIZE = 1.5, CHARTHICK = 1, PSym=[-15,-16], $
;  LineStyle=[0,2], Color=['cadet blue','Dark red'], Location=[2,0.85], $
;  /Data, /AddCmd
;
;;----------------------------
;;Im(A_k) comparison
;cgPlot, index_plot, imaginary(a_k_soler), PSym=-15, Color='cadet blue', thick = 2, YRange=[-0.2,0.5], XRange=[0,52], $
;  LineStyle=0, YStyle=1, XTitle='k', YTitle='Im(A_K)', CHARSIZE = 2.5, CHARTHICK = 1, /Window
;cgPlot, index_plot, imaginary(a_k), PSym=-16,Color='Dark Red', thick = 2, $
;  LineStyle=2, /Overplot, /AddCmd
;cgLegend, Title=['Im(Soler data)', 'Im(DIY data)' ], CHARSIZE = 1.5, CHARTHICK = 1, PSym=[-15,-16], $
;  LineStyle=[0,2], Color=['cadet blue','Dark red'], Location=[2,0.45], $
;  /Data, /AddCmd
;
;
;;----------------------------
;;Re(S_k) comparison
;cgPlot, index_plot, (s_k_soler), PSym=-15, Color='cadet blue', thick = 2, YRange=[-1,1], XRange=[0,52], $
;    LineStyle=0, YStyle=1, XTitle='k', YTitle='Re(S_K)', CHARSIZE = 2.5, CHARTHICK = 1, /Window
;cgPlot, index_plot, (s_k_2), PSym=-16,Color='black', thick = 2, $
;    LineStyle=1, /Overplot, /AddCmd
;cgPlot, index_plot, (s_k), PSym=-16,Color='Dark Red', thick = 2, $
;   LineStyle=2, /Overplot, /AddCmd
;cgLegend, Title=['Re(Soler data)','test', 'Re(DIY data)' ], CHARSIZE = 1.5, CHARTHICK = 1, PSym=[-15,-16,-16], $
;      LineStyle=[0,1,2], Color=['cadet blue', 'black' ,'Dark red'], Location=[2,0.85], $
;      /Data, /AddCmd
;
;;;----------------------------
;;Im(S_k) comparison
;cgPlot, index_plot, imaginary(s_k_soler), PSym=-15, Color='cadet blue', thick = 2, YRange=[-0.8,0.8], XRange=[0,52], $
;    LineStyle=0, YStyle=1, XTitle='k', YTitle='Im(S_K)', CHARSIZE = 2.5, CHARTHICK = 1, /Window
;cgPlot, index_plot, imaginary(s_k_2), PSym=-16,Color='black', thick = 2, $
;    LineStyle=1, /Overplot, /AddCmd
;cgPlot, index_plot, imaginary(s_k), PSym=-16,Color='Dark Red', thick = 2, $
;   LineStyle=2, /Overplot, /AddCmd
;cgLegend, Title=['Im(Soler data)', 'test', 'Im(DIY data)' ], CHARSIZE = 1.5, CHARTHICK = 1, PSym=[-15,-16,-16], $
;      LineStyle=[0,1,2], Color=['cadet blue', 'black','Dark red'], Location=[2,0.7], $
;      /Data, /AddCmd

;====================================================
;END OF A_K AND S_K CHECKER!!!!!!!!!
;====================================================

;========================================================================================
;END OF COEFFICIENTS
;========================================================================================

;========================================================================================
;START OF SOLVING DISPERSION RELATION
;========================================================================================
;I may up using SVD as I done previously. This is for simple checks
;Checks for zeta
;zeta_i = a - (ltr/2.0d) -ra
;zeta_e =  a + ltr/2.0d -ra

;new radial coordinate
zeta_i_norm = a_norm - (ltr_norm/2.0d) -ra_norm ; eq 13
zeta_e_norm =  a_norm + ltr_norm/2.0d -ra_norm ; eq 14

E_R_frac_i = 0
E_R_frac_e = 0
for k=0,N-2 do begin
  E_R_frac_i = E_R_frac_i + rho_k[k+1]*zeta_i_norm^k
  E_R_frac_e = E_R_frac_e + rho_k[k+1]*zeta_e_norm^k
endfor
E_R_frac_i = 1.0d/(w0_1_norm^2*E_R_frac_i)
E_R_frac_e = 1.0d/(w0_1_norm^2*E_R_frac_i)

; start of the sumation
G_i =0 
G_e =0
F_i = 0
F_e = 0
E_i = 0
E_e = 0
R_i = 0
R_e = 0 


for k=0,N-2 do begin
  G_i = G_i + a_k[k]*zeta_i_norm^(k+2) ; eq 28
  G_e =  G_e + a_k[k]*zeta_e_norm^(k+2) 
  
  ;alog(z) = alog(abs(z))+i*atan(z,/phase) where z is complex  
  F_i = F_i + s_k[k]*zeta_i_norm^(k)+(mi^2/(2*ra_norm^2))*(alog(abs(zeta_i_norm))+i*atan(zeta_i_norm,/phase))*a_k[k]*zeta_i_norm^(k+2) ;eq 29
  F_e = F_e + s_k[k]*zeta_e_norm^(k)+(mi^2/(2*ra_norm^2))*(alog(abs(zeta_e_norm))+i*atan(zeta_e_norm,/phase))*a_k[k]*zeta_e_norm^(k+2)
  
  E_i = E_i + (k+2.0d)*a_k[k]*zeta_i_norm^2 ; eq 30
  E_e = E_e + (k+2.0d)*a_k[k]*zeta_e_norm^2
  
  R_i = R_i + k*s_k[k]+(mi^2/(2.0d*ra_norm^2))*a_k[k]*zeta_i_norm^(k)+(mi^2/(2.0d*ra_norm^2))*(alog(abs(zeta_i_norm))+i*atan(zeta_i_norm,/phase))*(k+2.0d)*a_k[k]*zeta_i_norm^(k)
  R_e = R_e + k*s_k[k]+(mi^2/(2.0d*ra_norm^2))*a_k[k]*zeta_e_norm^(k)+(mi^2/(2.0d*ra_norm^2))*(alog(abs(zeta_e_norm))+i*atan(zeta_e_norm,/phase))*(k+2.0d)*a_k[k]*zeta_e_norm^(k)
endfor
E_i = E_R_frac_i*E_i 
E_e = E_R_frac_e*E_e
R_i = E_R_frac_i*R_i
R_e = E_R_frac_e*R_e 

k_pi = sqrt((w0_1^(2.0d)-wa0^(2.0d))/Va0^(2.0d)) ; Transverse wave number given by Eq. (9) Soler et al. (2013)
k_pe = sqrt(-(w0_1^(2.0d)-wa1^(2.0d))/Va1^(2.0d)) ; Transverse wave number given by Eq. (11) Soler et al. (2013)

k_pi_norm = k_pi*unit_length ; Transverse wave number given by Eq. (9) Soler et al. (2013)
k_pe_norm = k_pe*unit_length ; Transverse wave number given by Eq. (11) Soler et al. (2013)
;;; Making it dimensionless this way introduces a slight difference to above. This could possibly be an accumulating error that is been seen earilier in the code
;;k_pi_norm = sqrt((w0_1_norm^(2.0d)-(wa0/w_norm)^(2.0d))/(Va0/unit_velocity)^(2.0d)) ; Transverse wave number given by Eq. (9) Soler et al. (2013)
;;k_pe_norm = sqrt(-(w0_1_norm^(2.0d)-(wa1/w_norm)^(2.0d))/(Va1/unit_velocity)^(2.0d)) ; Transverse wave number given by Eq. (11) Soler et al. (2013)


 
;complex_dbeselk_deriv(k_pe*(a+Ltr/2.0d),mi,1000.0d,0.0d,40.0d))/(rho1*(w0_1-wa1^2))
dis_p1 = (k_pe_norm/rho1_norm*(w0_1_norm-kz_norm*Va1_norm))*(complex_dbeselk_deriv(k_pe_norm*((a+Ltr/2.0d)/a),mi,1000.0d,0.0d,N)/complex_dbeselk(k_pe_norm*((a+Ltr/2.0d)/a),mi,1000.0d,0.0d,N))
dis_p2 = (k_pe_norm/rho0_norm*(w0_1_norm-kz_norm*Va0_norm))*(complex_dbeselj_deriv(k_pi_norm*((a-Ltr/2.0d)/a),mi,1000.0d)/complex_dbeselj(k_pi_norm*((a-Ltr/2.0d)/a),mi,1000.0d,0.0d,!PI))

disp_e = (dis_p1*G_e-E_e)/(dis_p1*F_e-R_e)
disp_i = (dis_p2*G_i-E_i)/(dis_p2*F_i-R_i)
disp_final = disp_e-disp_i
;disp_i = 
;========================================================================================
;END OF SOLVING DISPERSION RELATION
;========================================================================================
end

