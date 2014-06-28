function shapelets_rsquared, decomp,        $
                             FLUX=flux,     $
                             MATRIX=matrix, $
                             ERROR=error,   $
                             POLAR=polar,   $
                             CARTESIAN=cartesian

;$Id: shapelets_rsquared.pro, v2$
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
;
;+
; NAME:
;      SHAPELETS_RSQUARED
;
; PURPOSE:
;      Computes the R^2 size measure from a linear summation of shapelet
;      coefficients, read in from a Cartesian decomp structure.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      decomp - A Cartesian shapelet decomp structure, with at least
;               beta and coefficients defined.
;
; KEYWORD PARAMETERS:
;     [/CART] - Perform summation in Cartesian shapelet space. [DEFAULT]
;    [/POLAR] - Perform sunnation in polar shapelet space.
;    [/ERROR] - Simultaneously calculate the error on R^2 (and the flux).
;
; OUTPUTS:
;    rsquared - Shapelet rms size measure, squared.
;      [flux] - Floating point variable optionally returned with the object's 
;               flux (calculation needed to work out the centroid anyway).
;    [matrix] - Cartesian to polar conversion matrix, if /POLAR was set.
;
; EXAMPLE USE:
;      result = shapelets_rsquared(decomp,flux=flux)
;
; NOTES: 
;      Not yet compatible with oversampling, nor full covariance error matrix
;      of coefficients.
;
; MODIFICATION HISTORY:
;      Oct 03 - Calculation of errors added by RM.
;      Feb 03 - Polar/Cartesian switches added by RM. 
;      Apr 02 - Written by R.Massey
;-

COMPILE_OPT idl2

;
; SET DEFAULTS
;
; Backwardly compatible
if not tag_exist(decomp,"polar") then decomp=create_struct(decomp,"polar",0B)
; Decide whether default method should be in Cartesian or polar shapelet space
if ((keyword_set(cartesian)+keyword_set(polar)) mod 2) eq 0 then cartesian=1-decomp.polar
; Remember input variable, so that it can be restored later
decomp_input=decomp

if keyword_set(cartesian) then begin 
  ;
  ; WORK IN CARTESIAN SHAPELET SPACE  
  ;
  shapelets_polar_convert,decomp,/CARTESIAN,/SILENT

  ; Calculate weights
  fweight=fltarr(decomp.n_coeffs)
  rweight=fltarr(decomp.n_coeffs)
  shapelets_make_nvec, decomp.n_max, n_1, n_2, n_coeffs

  for i=0,n_coeffs-1 do begin
    n1=n_1[i] & n2=n_2[i]
    if n1 mod 2 eq 0 and n2 mod 2 eq 0 then begin

      fweight[i]=$
       2.^(.5*(2.-n1-n2)) * $
       sqrt(factorial(n1)*factorial(n2))/factorial(n1/2)/factorial(n2/2)

      rweight[i]=$
       2.^(.5*(4.-n1-n2)) * (1.+n1+n2) * $
       sqrt(factorial(n1)*factorial(n2))/factorial(n1/2)/factorial(n2/2)

    endif
  endfor

  ; Calculate flux and size
  rootpi   = sqrt(!pi)
  flux     = rootpi*decomp.beta*total(fweight*decomp.coeffs)
  rsquared = rootpi*decomp.beta^3/flux*total(rweight*decomp.coeffs)

  ; Calculate errors on the flux and size
  if keyword_set(error) then begin
    flux_error       = rootpi*decomp.beta*sqrt(total((fweight*decomp.coeffs_error)^2))
    rsquared_variance= !pi*decomp.beta^6*total((rweight*decomp.coeffs_error)^2)
    rsquared_error   = abs(rsquared)*$
                       sqrt( rsquared_variance/(rsquared*flux)^2 $
                            +(flux_error/flux)^2 )
    flux             = [flux,flux_error]
    rsquared         = [rsquared,rsquared_error]
  endif

endif else begin
  
  ;
  ; WORK IN POLAR SHAPELET SPACE 
  ;
  ; Convert coefficients to polar shapelet format if not already like that
  shapelets_polar_convert,decomp,/POLAR,/SILENT
  m0=where(decomp.m eq 0)
  rootpi   = sqrt(!pi)
  
  flux     = rootpi*decomp.beta*2*total(float(decomp.coeffs[m0]))
  rsquared = 4*rootpi*decomp.beta^3/flux*total((float(decomp.coeffs*(1+decomp.n)))[m0])

  if keyword_set(error) then begin
    ; Errors will not be quite correct if we do not have the full covariance
    ;   matrix and try to perform the linear transformation into polar space
    if size(decomp.coeffs_error,/n_dimensions) ne 1 then message,$
      "Warning: error calculation using covariance matrix is not yet written!"
    flux_error=rootpi*decomp.beta*2*sqrt(total((float(decomp.coeffs_error)^2)[m0]))
    rsquaredf_variance=16*!pi*decomp.beta^6*total(((float(decomp.coeffs_error)*(1+decomp.n))^2)[m0])
    rsquared_error=rsquared*$
                   sqrt( rsquaredf_variance/(rsquared*flux)^2 $
                        +(flux_error/flux)^2 )

    flux=[flux,flux_error]
    rsquared=[rsquared,rsquared_error]
  endif

endelse

; Restore input variable
decomp=decomp_input

return,rsquared

end
