function shapelets_polar_reduce,polar_coeffs

;$Id: shapelets_polar_reduce.pro, v2$
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
; Apr 02 - RM removed a degeneracy in phase coefficients
; Feb 02 - Written by Richard Massey
;
; Removes the duplicated/degenerate coefficients when cartesian shapelets have
; been converted into polar form and leaves only the minimum number of
; independent parameters. This will speed up (and make possible) parametisation
; of all HDF galaxies in shapelet (probability) space. Put them back with
; shapelets_polar_expand.pro.
;
; This is achieved because, while polar shapelet coefficients are in general
; complex, diagonal (n=0) elements are wholly real and entries mirrored in the
; n=0 diagonal (ie m=-m) are complex conjugates of each other. Both of these
; constraints are necessary and sufficient requirements for the resulting object
; to be real itself.
;
; The first version of this code seperated them into Real & Imaginary parts - it
; now splits into r & theta, which I think is more physically meaningful.
;-

COMPILE_OPT OBSOLETE

n_coeffs   = n_elements(polar_coeffs)    ; Background work, figuring out how much there is to do
n_max = fix(0.5*(sqrt(8.*n_coeffs+1.)-1.))-1


; Polar basis subscripts ( n_r=n1 and n_l=n2 )
shapelets_make_nvec, n_max, n1, n2
n=n1+n2
m=n1-n2



; Convert to {mag,phase} form
magnitude = abs(polar_coeffs)
phase     = atan(imaginary(polar_coeffs),float(polar_coeffs))

;; Make phases relative
;for i=n_coeffs-1,0,-1 do begin
;  if m(i) gt 0 then begin      ; only bother with top half of {n,m} plot
;    if m(i) lt n(i) then begin ; in body of {n,m} plot, cf with 2 to the left
;      phase(i)=phase(i)-phase(where(m eq m(i) and n eq n(i)-2))
;    endif else begin           ; on top diagonal, work backwards towards vertex
;      phase(i)=phase(i)-phase(where(m eq m(i)-1 and n eq n(i)-1))
;    endelse
;    ; Not using because only comparing like m to like m at the moment
;    ;phase(i)=phase(i)/m(i)
;
;    ; Think this gets rid of a degeneracy
;    phase(i)=phase(i)/2.
;    phase(i)=(phase(i)+!pi) mod (!pi)
;    phase(i)=phase(i)*2.
;
;    ; Old one that I know works, (I think it has a degeneracy, though)
;    ;phase(i)=phase(i)/2.
;    ;phase(i)=phase(i) mod (!dpi)
;  endif
;endfor

; Keep only what is necessary
reduced_polar_coeffs=fltarr(n_coeffs)      ; Array to receive final answer
coef=0L                   ; Incremental counter variable
for i=0,n_coeffs-1 do begin
  if m[i] lt 0 then begin ; Get magnitudes from bottom half of {n,m} plot
    reduced_polar_coeffs(coef)=magnitude[i]
    coef=coef+1
  endif else if m[i] eq 0 then begin 
    reduced_polar_coeffs(coef)=polar_coeffs[i]
    coef=coef+1
  endif
endfor
for i=0,n_coeffs-1 do begin
  if m[i] gt 0 then begin ; Get phases from top half of {n,m} plot
    reduced_polar_coeffs(coef)=phase[i]
    coef=coef+1
  endif
endfor


return,reduced_polar_coeffs

end
