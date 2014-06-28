function shapelets_concentration,decomp,            $
                                 RATIO=ratio,       $
				 KURTOSIS=kurtosis

;$Id: shapelets_concentration.pro, v2$
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
;      SHAPELETS_CONCENTRATION
;
; PURPOSE:
;      Returns concentration morphology index of a Cartesian decomp structure.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      DECOMP - A Cartesian shapelet decomposition.
;
; KEYWORD PARAMETERS:
;      KURTOSIS - Calculate the object kurtosis.
;      RATIO    - Calculate the ratio of R/beta.
;      ERROR    - Simultaneously calculate the error on the concentration.
;
; OUTPUTS:
;      CONCENTRATION - A floating point variable containing the object's
;                      concentration morphology index.
;                      If /ERROR is set, this becomes a vector
;                      [concentration,concentration_error]
;
; OPTIONAL OUTPUTS:
;      FLUX     - Flux of object (needed to calculate asymmetry anyway).
;      RSQUARED - R^2 size of object (needed to calculate asymmetry anyway).
;      MATRIX   - Cartesian to polar shapelet conversion matrix.
;
; EXAMPLE USE:
;      result = shapelets_concentration(decomp,/error)
;
; TO DO:
;      Add the more complicated option involving integration within
;      various Petrosian radii suggested in Polar Shapelets.
;
; MODIFICATION HISTORY:
;      Jul 05 - Ability to accept polar decomp structures added by RM.
;      Jun 05 - Series that doesn't converge in Polar Shapelets hacked in by RM.
;      Jun 04 - Written by Richard Massey.
;-

COMPILE_OPT idl2

if keyword_set(kurtosis) then begin

  ; Backwardly compatible
  if not tag_exist(decomp,"polar") then decomp=create_struct(decomp,"polar",0B)
  
  ; Remember input variable, so that it can be restored later
  decomp_input=decomp
  
  ; Calculate polar shapelet coefficients
  shapelets_polar_convert,decomp,/polar

  ; Calculate total flux
  m0=where(decomp.m eq 0,n_m0)
  if n_m0 gt 0 then begin
    
    ; Calculate dimensionless equlvalent of flux
    flux=total(float(decomp.coeffs[m0]))

    ; Calculate concentration
    if flux eq 0. then concentration=!values.f_nan else begin

      ; Calculate weights
      weights=fltarr((decomp.n_max+2)/2)
      weights[0]=1.    ; Overall normalisation unimportant
      if decomp.n_max ge 2 then begin
    	weights[1]=weights[0]
    	if decomp.n_max ge 4 then begin
    	  for n=4,decomp.n_max,2 do begin
    	    weights[n/2]=( 10.*weights[n/2-1] + (n-2)*weights[n/2-2] ) /n
    	 endfor
    	endif
      endif
      ; Combine coefficients
      concentration=float(total(decomp.coeffs[m0]*weights))/flux

      ; Rescale to roughly match spread of Chirs Conselice's concentration values
      ;concentration=-50*concentration+5

    endelse
  endif

  ; Restore input variable
  decomp=decomp_input

endif else if keyword_set(ratio) then begin

  ; Calculate R^2 size
  rsquared=shapelets_rsquared(decomp,flux=flux,error=error,cartesian=cartesian,polar=polar)

  ; Calculate concentration
  if decomp.beta eq 0. then concentration=!values.f_nan else $
    concentration=sqrt(rsquared)/decomp.beta

endif else begin

  message,"No method specified!"

endelse

return,concentration

end


