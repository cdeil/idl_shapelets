pro shapelets_make_nvec, n_max, n1, n2, n_coeffs, POLAR=polar

;$Id: shapelets_make_nvec.pro, v2$
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
;      SHAPELETS_MAKE_NVEC
;
; PURPOSE:
;      Set up look-up tables to specify which number in a vector of shapelet
;      coefficients corresponds to which n1 and n2, or which nr and nl.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      n_max - Highest n to go up to.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      n1,n2    - Coefficient number vectors.
;      n_coeffs - Optionally returned containing total number of coefficients.
;
; EXAMPLE USE:
;      A vector of Cartesian shapelet coefficients can be indexed using:
;        shapelets_make_nvec, n_max, n_1, n_2 [,n_coeffs]
;
;      A vector of polar shapelet coefficients can be indexed using:
;        shapelets_make_nvec, n_max, n, m, /POLAR [,n_coeffs]
;      or
;        shapelets_make_nvec, n_max, n_r, n_l [,n_coeffs]
;        n=n_r+n_l
;        m=n_r-n_l
;
; MODIFICATION HISTORY:
;      Jul 05 - POLAR option added by RM.
;      Nov 01 - Written by Richard Massey.
;-

COMPILE_OPT idl2, HIDDEN

n_coeffs=(fix(n_max)+1L)*(fix(n_max)+2L)/2L
n1=lonarr(n_coeffs,/nozero)
n2=lonarr(n_coeffs,/nozero)
i=0L
for n=0L, n_max do begin
  for n1i=0L, n do begin
    n1[i]=n1i
    n2[i]=n-n1i
    i+=1
  endfor
endfor

if keyword_set(polar) then begin
  n1+=n2     ; n
  n2=n1-2*n2 ; m
endif

end
