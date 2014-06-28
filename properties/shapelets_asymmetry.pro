function shapelets_asymmetry,decomp,        $
                             FLUX=flux,     $
                             MATRIX=matrix, $
                             ERROR=error
 
;$Id: shapelets_asymmetry.pro, v2$
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
;      SHAPELETS_ASYMMETRY
;
; PURPOSE:
;      Returns asymmetry morphology index of a Cartesian decomp structure.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      DECOMP - A Cartesian shapelet decomposition.
;
; KEYWORD PARAMETERS:
;      ERROR  - Simultaneously calculate the error on the asymmetry.
;
; OUTPUTS:
;      ASYMMETRY - A floating point variable containing the object's asymmetry
;                  morphology index.
;                  If /ERROR is set, this becomes a vector
;                  [asymmetry,asymmetry_error]
;
; OPTIONAL OUTPUTS:
;      FLUX   - Flux of object (needed to calculate asymmetry anyway).
;      MATRIX - Cartesian to polar shapelet conversion matrix.
;
; EXAMPLE USE:
;      result = shapelets_asymmetry(decomp)
;
; MODIFICATION HISTORY:
;      Jun 04 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Backwardly compatible
if not tag_exist(decomp,"polar") then decomp=create_struct(decomp,"polar",0B)

; Remember input variable, so that it can be restored later
decomp_input=decomp

; Convert coefficients to polar shapelet format if not already like that
shapelets_polar_convert,decomp,/POLAR,/SILENT

; Calculate flux
tworootpibeta=2*sqrt(!pi)*decomp.beta
roottwooverpi=sqrt(2.)/!pi
mzero=where(decomp.m eq 0)
modd=where(abs(decomp.m mod 2) eq 1)
flux=tworootpibeta*total(float(decomp.coeffs[mzero]))

; Calculate asymmetry
if flux eq 0 then asymmetry=!values.f_nan else $
  asymmetry=roottwooverpi/flux*total(abs(decomp.coeffs[modd]))

; Calculate error on asymmetry index if required
if keyword_set(error) then begin
  ; Errors will not be quite correct if we do not have the full covariance
  ;   matrix and try to perform the linear tranformation into polar space
  if size(decomp.coeffs_error,/n_dimensions) ne 1 then message,$
    "Warning: error calculation using covariance matrix is not yet written!"
  flux_error = tworootpibeta*sqrt(total((float(decomp.coeffs_error[mzero]))^2))
  if flux eq 0. or asymmetry eq 0. then asymmetry_error=!values.f_infinity else begin
    coeffs_error_sq=(float(decomp.coeffs_error[modd])/float(decomp.coeffs[modd]))^2+$
                    (imaginary(decomp.coeffs_error[modd])/imaginary(decomp.coeffs[modd]))^2
    asymmetry_error=roottwooverpi*sqrt(total(coeffs_error_sq))
    asymmetry_error=abs(asymmetry)*sqrt((asymmetry_error/asymmetry*flux)^2+(flux_error/flux)^2)
  endelse
  flux=[flux,flux_error]
  asymmetry=[asymmetry,asymmetry_error]
endif

; Restore input variable
decomp=decomp_input

return,asymmetry

end


