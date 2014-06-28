function shapelets_polar_expand,polar_coeffs

;$Id: shapelets_polar_expand.pro, v2$
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
; Apr 02 - RM Removed a degeneracy in phase coefficients
; Feb 02 - Written by Richard Massey
;
; Manufactures and reinserts the degenerate complex parts of polar shapelet
; coefficients so they can be transformed back into cartesian shapelets.
; See shapelets_polar_reduce.pro for extraction technique.
;-

COMPILE_OPT OBSOLETE

n_coeffs = n_elements(polar_coeffs)     ; Background work, figuring out how much there is to do
n_max=fix(0.5*(sqrt(8.*n_coeffs+1.)-1.))-1

 ; Polar basis subscripts ( n_r=n1 and n_l=n2 )
shapelets_make_nvec, n_max, n1, n2
n=n1+n2
m=n1-n2

; Number of entries in reduced array which correspond to mags and phases
;n_mag   = n_elements(where(m ge 0))
;n_phase = n_elements(where(m gt 0))

; Recover full information, extracing and expanding from the reduced form
phase=fltarr(n_coeffs)      ; Arrays to receive final answer in {r,theta} form
magnitude=fltarr(n_coeffs)  ;
coef=0L                     ; Incremental counter variable
for i=0,n_coeffs-1 do begin ;
  if m(i) le 0 then begin   ; Recover magnitudes in bottom half of {n,m} plot
    magnitude(i)=polar_coeffs[coef]
    coef=coef+1
  endif
endfor

; Recover magnitudes in top half of {n,m} plot by symmetry
for i=0,n_coeffs-1 do begin
  if m(i) gt 0 then begin   ; Top half of {n,m} plot
    magnitude(i)=magnitude(where(m eq -m(i) and n eq n(i)))
  endif
endfor

; Recover phases
for i=0,n_coeffs-1 do begin
  if m(i) gt 0 then begin   ; Recover phases in top half of {n,m} plot
    phase(i)=polar_coeffs(coef)
    ; Recover degeneracy from shapelets_polar_reduce.pro
    ;phase(i)=phase(i)/2.
    coef=coef+1
  endif
endfor

;; Unrelativise phases in top half of {n,m} plot explicitly
;for i=0,n_coeffs-1 do begin
;  if m(i) gt 0 then begin      ; Top half of {n,m} plot
;    ;phase(i)=phase(i)*m(i)
;    phase(i)=phase(i)*2.
;    if m(i) lt n(i) then begin ; In body of {n,m} plot, cf with 2 to the left
;      phase(i)=phase(i)+phase(where(m eq m(i) and n eq n(i)-2))
;    endif else begin           ; On top diagonal, work outwards from vertex
;     phase(i)=phase(i)+phase(where(m eq m(i)-1 and n eq n(i)-1))
;    endelse
;  endif
;endfor

; Recover phases in bottom half of {n,m} plot by symmetry
for i=0,n_coeffs-1 do begin
  if m(i) lt 0 then begin      ; Bottom half of {n,m} plot
    phase(i)=-phase(where(m eq -m(i) and n eq n(i)))
  endif
endfor


; Convert back to {Real,Imaginary} form
expanded_polar_coeffs=magnitude*exp(complex(0,1)*phase)

return,expanded_polar_coeffs

end
