pro shapelets_keep_it_real, structure,           $
                            NOHISTORY=nohistory

;$Id: shapelets_keep_it_real.pro, v2$
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
;      SHAPELETS_KEEP_IT_REAL
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Makes sure a polar shapelet model is wholly real by discarding
;      any imaginary part. Has no effect on a Cartesian shapelet model.
;
; INPUTS:
;      STRUCTURE - A shapelet decomp or shapecat structure.
;
; OPTIONAL INPUTS:
;      NOHISTORY - Prevents the operation from being recorded in the
;                  object's history tag.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      STRUCTURE - the structure is returned, real.

; MODIFICATION HISTORY:
;      Dec 05 - Written by Richard Massey
;-

COMPILE_OPT idl2

; Maintain backwards compatibility
if not shapelets_structure_type(structure,message=message) then message,message

; Keep it real
if structure.polar then begin
  ; Extract shapelet coefficients from the structure
  if structure.type eq "shapecat" then begin
    coeffs=transpose(structure.coeffs)
    n_max=structure.maxn_max
  endif else if structure.type eq "decomp" then begin
    coeffs=structure.coeffs
    n_max=structure.n_max
  endif else message,"Structure type not recognised!"
  shapelets_make_nvec,n_max,nn,mm,/POLAR
  ; Loop over each coefficient pairs, and calculate their real parts
  if n_max ge 1 then begin
    for n=1,n_max do begin
      for m=n mod 2,n,2 do begin
        positive=where(nn eq n and mm eq m, n_positive)
        negative=where(nn eq n and mm eq -m,n_negative)
        if n_negative ne 1 or n_positive ne 1 then message,"Could not find coefficients!"
        structure.coeffs[positive,*]=complex(float(structure.coeffs[positive,*]+structure.coeffs[negative,*]),$
                                         imaginary(structure.coeffs[positive,*]-structure.coeffs[negative,*]))/2.
        structure.coeffs[negative,*]=conj(structure.coeffs[positive,*])
      endfor
    endfor
    ; Replace coefficients in the input structure itself
    if structure.type eq "decomp" then coeffs=structure.coeffs $
      else structure.coeffs=transpose(coeffs)
  endif
endif

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history, structure, $
  "Object kept real!"

end

