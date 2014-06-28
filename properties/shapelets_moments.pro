pro shapelets_moments,decomp,dmom,printit=printit  

;$Id: shapelets_moments.pro, v2$
;
; Copyright © 2005 Richard Massey and Alexandre Refregier.
;
; This file is a part of the Shapelets analysis code.
; www.astro.caltech.edu/~rjm/shapelets/
;
; The Shapelets code is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public Licence as published
; by the Free Software Foundation; either version 2 of the Licence, or
; (at your option) any later version.
;
; The Shapelets code is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
; GNU General Public Licence for more details.
;
; You should have received a copy of the GNU General Public Licence
; along with the Shapelets code; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
;
;+
; NAME:
;      SHAPELETS_MOMENTS
;
; PURPOSE:
;      Compute the zeroth (flux) and first (centroid) moment for a basis 
;      decomposition. Also computes the characteristic order n_f and scale
;      parameter beta_f.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      decomp - decomposition produced by shapelets_decomp.pro
;
; OUTPUTS:
;      An IDL structure containing the flux and centroid.
;
; MODIFICATION HISTORY:
;      Mar 2005 - Obsoleted by Richard Massey. Now use shapelets_quadrupole.pro.
;      Aug 2003 - Modified by AR to include quadrupole moments by
;                 merging quadmom.pro routine from RJM and Tzu-Ching Chang
;      Aug 2003 - Temporarily changed to a function by Richard Massey, then
;                 changed back: beware!
;      Jul 2000 - Written by A. Refregier
;      Mar 2001 - Modified by AR to compute n_f and beta_f and r
;-

compile_opt OBSOLETE

; compute nf, the characteristic order
nf=total( (decomp.n1+decomp.n2)*decomp.coeffs^2 ) / total( decomp.coeffs^2 )

; compute the flux
f=0.
x1=0. & x2=0.
r=0.
q11=0. & q22=0. & q12=0.
for i=0,decomp.n_coeffs-1 do begin
  n1=decomp.n1(i) & n2=decomp.n2(i)
  if n1 mod 2 eq 0 and n2 mod 2 eq 0 then begin
    f=f+sqrt(!pi)*decomp.beta*2.^(.5*(2.-n1-n2))*$
       sqrt(factorial(n1)*factorial(n2))/$
       factorial(n1/2)/factorial(n2/2)*decomp.coeffs(i)
    r=r+sqrt(!pi)*decomp.beta^3*2.^(.5*(4.-n1-n2))*$
       (1.+n1+n2)*$ 
       sqrt(factorial(n1)*factorial(n2))/$
       factorial(n1/2)/factorial(n2/2)*decomp.coeffs(i)   
    q11=q11+sqrt(!pi)*decomp.beta^(3.)*2.^(.5*(2.-n1-n2))*$
       (1.+2.*n1)*$
       sqrt(factorial(n1)*factorial(n2))/$
       factorial(n1/2)/factorial(n2/2)*$
       decomp.coeffs(i)
    q22=q22+sqrt(!pi)*decomp.beta^(3.)*2.^(.5*(2.-n1-n2))*$
       (1.+2.*n2)*$
       sqrt(factorial(n1)*factorial(n2))/$
       factorial(n1/2)/factorial(n2/2)*$
       decomp.coeffs(i)
  endif
  if n1 mod 2 eq 1 and n2 mod 2 eq 0 then begin
    x1=x1+sqrt(!pi)*decomp.beta^2*sqrt(n1+1.)*2.^(.5*(2.-n1-n2))*$
       sqrt(factorial(n1+1)*factorial(n2))/$
       factorial((n1+1)/2)/factorial(n2/2)*decomp.coeffs(i) 
  endif
  if n1 mod 2 eq 0 and n2 mod 2 eq 1 then begin
    x2=x2+sqrt(!pi)*decomp.beta^2*sqrt(n2+1.)*2.^(.5*(2.-n1-n2))*$
       sqrt(factorial(n1)*factorial(n2+1))/$
       factorial(n1/2)/factorial((n2+1)/2)*decomp.coeffs(i)
  endif
  if n1 mod 2 eq 1 and n2 mod 2 eq 1 then begin
    q12=q12+sqrt(!pi)*decomp.beta^(3.)*2.^(.5*(2.-n1-n2))*$
       sqrt(n1+1.)*sqrt(n2+1.)*$
       sqrt(factorial(n1+1)*factorial(n2+1))/$
       factorial((n1+1)/2)/factorial((n2+1)/2)*$
       decomp.coeffs(i)
  endif
endfor
x1=x1/f             ; centroid [pixels]
x2=x2/f
r=sqrt(r/f)         ; rms radius [pixels]
q11=q11/f & q22=q22/f & q12=q12/f

; compute ellipticity and other quantities
e1=(q11-q22)/(q11+q22)   ; ellipticity   [1]
e2=2.*q12/(q11+q22)
a=sqrt(q11+q22+sqrt((q11-q22)^2+4.*q12^2))  ; principal axes [pixels]
b=sqrt(q11+q22-sqrt((q11-q22)^2+4.*q12^2))
alpha=.5*atan(2.*q12,q11-q22)   ; position angle, counter-clockwise from x-axis
                                ;  [rad]
; print out results
if keyword_set(printit) then begin
  print,'warning: untested'
  print,'flux:',f
  print,'centroid:',x1+decomp.x(0),x2+decomp.x(1)
  print,'rms radius r:',r
  print,'Q11, Q22, Q12:',q11,q22,q12
  print,'e1,e2:',e1,e2
  print,'a,b:',a,b
  print,'alpha (deg, c-clock from x):',alpha*!radeg
endif

; save in structure
dmom={f:f,xc:[x1+decomp.x(0),x2+decomp.x(1)],r:r,q11:q11,q22:q22,q12:q12,$
      e1:e1,e2:e2,a:a,b:b,alpha:alpha}

end
