pro shapelets_recentre, structure,               $
                        EXTEND=extend,           $
                        PRECISION=precision,     $
                        ORDER=order,             $
                        MAX_ITER=max_iter,       $
                        TRANSLATION=translation, $
                        SILENT=silent

;$Id: shapelets_recentre.pro, v2$
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
;      SHAPELETS_SHEAR
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Applies a translation to an object to ensure that its centre of light
;      lies exactly on the origin (the centre of the basis functions).
;
; INPUTS:
;      STRUCTURE   - A shapelet decomp structure or shapecat structure.
;
; OPTIONAL INPUTS:
;      EXTEND      - Increses n_max (and sets new coefs to zero) before starting
;                    translation(s). Gives new coefficients somewhere to go.
;      PRECISION   - Desired accuracy of final centroid [pixels].
;                    DEFAULT: 1/1000th pixel.
;      ORDER       - Expansion order to which translations are performed.
;                    DEFAULT: first order.
;      MAX_ITER    - Maximum number of iterations, to avoid an endless loop.
;                    DEFAULT: 100.
;
; KEYWORD PARAMETERS:
;      SILENT       - Operate silently.
;
; OUTPUTS:
;      STRUCTURE   - The same structure, but with its objects(s) recentred.
;      TRANSLATION - The toal translation applied in the end.
;
; MODIFICATION HISTORY:
;      Apr 05 - Written by Richard Massey
;-

COMPILE_OPT idl2
ON_ERROR,2

; Parse input preferences
if keyword_set(extend) then shapelets_extend_nmax,structure,round(extend)
if not keyword_set(precision) then precision=0.001
if not keyword_set(max_iter) then max_iter=100
warning_level=0.1 ; Warning given if initial offset larger than this # of betas

if structure.type eq "decomp" or structure.type eq "decomp_cartesian" then begin 

  ; Save input in case we need to revert
  structure_in=structure

  ; Determine initial offset
  centroid=shapelets_centroid(structure,/BASIS_FUNCTIONS)
  if total(finite(centroid)) ne 2 then begin
    message,"Could not determine object centroid!",/info,noprint=silent
    return
  endif

  ; Correct offset by iterating
  iter=0
  translation=0.
  mag_shift=!values.f_infinity
  while mag_shift gt precision^2 and iter lt max_iter do begin
    ; Calculate initial offset in units of beta
    shift=-centroid
    ; Check that the iteration is converging
    old_mag_shift=mag_shift
    mag_shift=(centroid[0]^2+centroid[1]^2)
    if mag_shift gt old_mag_shift then begin
      message,"Iteration diverging!",/info,noprint=silent
      return
    endif
    ; Check that the shift is a reasonable task, and reduce it if it looks big
    if mag_shift gt warning_level then begin
      message,"WARNING: Trying to correct for a large centroid offset!",/info,noprint=silent
      shift=shift*((float(iter)/max_iter+0.5)<1)
    endif
    translation=translation+shift
    shapelets_translate,structure,shift
    centroid=(shapelets_centroid(structure)-structure.x)
    old_mag_shift=mag_shift
    iter+=1
  endwhile
 
  ; Update the location of the centre of the basis functions for that object
  structure.x=structure.x-translation

  ; Check that the offset really is now within desired tolerances
  if (centroid[0]^2+centroid[1]^2) gt precision^2 then message,/info,noprint=silent,$
    "Recentering failed to converge in "+strtrim(string(iter),2)+" iterations!"
  
  ;; Recover input decomp if everything's gone horribly wrong
  ;if not finite(shapelets_flux(structure)) then begin
  ;  message,"Recentering failed!",/info,noprint=silent
  ;  structure=structure_in
  ;endif
  
endif else if structure.type eq "shapecat" then begin 

  ; Prepare empty array to report shifts performed
  translation=fltarr(structure.n,2)

  ; Loop over each object in turn
  for i=0L,structure.n-1 do begin
    if (i+1) mod 500 eq 0 then message,"Recentering object #"+strtrim(string(i),2)+$
      "/"+strtrim(string(structure.n),2),/info,noprint=silent
    ; Extract single object from the catalogue
    decomp=shapelets_shapecat2decomp(structure,i)
    ; Shift that object around
    shapelets_recentre, decomp, precision=precision, max_iter=max_iter, $
      translation=indiv_translation, silent=silent
    ; Put it back into the catalogue
    structure.coeffs[i,0:decomp.n_coeffs-1]=decomp.coeffs
    ; Fill in array with which to report the performed shifts
    translation[i,*]=indiv_translation
  endfor
  structure.x=structure.x-translation

endif else message,"Cannot recentre objects in "+structure.type+" structures!"

end

