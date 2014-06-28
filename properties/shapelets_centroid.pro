function shapelets_centroid, decomp,                         $
                             FLUX=flux,                      $
                             ERROR=error,                    $
                             BASIS_FUNCTIONS=basis_functions,$
                             IMAGE_COORDS=image_coords,      $
                             N_MAX=n_max,                    $
                             GAUSSIAN=GAUSSIAN,              $
                             MATRIX=matrix,                  $
                             POLAR=polar,                    $
                             CARTESIAN=cartesian

;$Id: shapelets_centroid.pro, v2$
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
;      SHAPELETS_CENTROID
;
; PURPOSE:
;      Computes the centroid from a linear summation of shapelet coefficients,
;      read in from a Cartesian decomp structure.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      decomp - A Cartesian shapelet decomp structure, with at least
;               beta and coefficients defined.
;
; OPTIONAL INPUTS:
;       N_MAX - Summation performed only to this n_max. Overrides GAUSSIAN
;               keyword. Flux summation still performed to all orders.
;
; KEYWORD PARAMETERS:
;        CART - Perform summation in Cartesian shapelet space. [DEFAULT]
;       POLAR - Perform summation in polar shapelet space.
;       ERROR - Simultaneously calculate the error on the centroid (and flux).
;     BASIS_F - Centroid is returned relative to the centre of the basis fns.
;     IMAGE_C - Centroid is returned in global image coordinates.
;    GAUSSIAN - Gaussian weighted centroid calculated, by setting N_MAX to 2.
;
; OUTPUTS:
;    CENTROID - Shapelet measure of centroid
;
; OPTIONAL OUTPUTS:
;        FLUX - Floating point variable optionally returned with the object's 
;               flux (calculation needed to work out the centroid anyway).
;      MATRIX - Cartesian to polar conversion matrix, if /POLAR was set.
;
; EXAMPLE USE:
;      result = shapelets_centroid(decomp,flux=flux)
;
; NOTES: 
;      Not yet compatible with oversampling, nor full covariance error matrix
;      of coefficients.
;
; MODIFICATION HISTORY:
;      Jan 06 - N_MAX and GAUSSIAN keywords added by RM.
;      Jul 05 - BASIS_FUNCTIONS and IMAGE_COORDS keywords added by RM.
;      Jul 05 - Rendered capable of accepting polar decomp structures by RM.
;      Oct 03 - Calculation of error on flux added by RM.
;      Feb 03 - Written by Richard Massey
;-

COMPILE_OPT idl2

;
; SET DEFAULTS
;
; Backwardly compatible
if not tag_exist(decomp,"polar") then decomp=create_struct(decomp,"polar",0B)
; Decide whether default method should be in Cartesian or polar shapelet space
if ((keyword_set(cartesian)+keyword_set(polar)) mod 2) eq 0 then cartesian=1-decomp.polar
; Decide up to what order the summation should be performed
if n_elements(n_max) eq 0 then if keyword_set(gaussian) then n_max=2 else n_max=decomp.n_max
n_max=0>fix(n_max)<decomp.n_max
; Remember input variable, so that it can be restored later
decomp_input=decomp
rootpi=sqrt(!pi)


if keyword_set(cartesian) then begin 

  ;
  ; WORK IN CARTESIAN SHAPELET SPACE
  ;
  ; Convert coefficients to Cartesian shapelet format if not already like that
  if decomp.polar then shapelets_polar_convert,decomp,/CARTESIAN,/SILENT
  
  ; Calculate weights
  fweight=fltarr(decomp.n_coeffs)
  xweight=fltarr(decomp.n_coeffs)
  yweight=fltarr(decomp.n_coeffs)

  for i=0,decomp.n_coeffs-1 do begin
    n1=decomp.n1[i] & n2=decomp.n2[i]
    ;if n1+n2 le n_max then begin
      if n1 mod 2 eq 0 and n2 mod 2 eq 0 then begin

        fweight[i]=$
         2.^(.5*(2.-n1-n2)) * $
         sqrt(factorial(n1)*factorial(n2))/factorial(n1/2)/factorial(n2/2)

      endif else if n1 mod 2 eq 1 and n2 mod 2 eq 0 and (n1+n2) le n_max then begin

        xweight[i]=sqrt(n1+1.)*2.^(.5*(2.-n1-n2))*$
         sqrt(factorial(n1+1)*factorial(n2))/$
         factorial((n1+1)/2)/factorial(n2/2)

      endif else if n1 mod 2 eq 0 and n2 mod 2 eq 1 and (n1+n2) le n_max then begin
     
        yweight[i]=sqrt(n2+1.)*2.^(.5*(2.-n1-n2))*$
         sqrt(factorial(n1)*factorial(n2+1))/$
         factorial(n1/2)/factorial((n2+1)/2)

      endif
   ; endif
  endfor

  ; Calculate flux and centroid
  flux = rootpi*decomp.beta*total(fweight*decomp.coeffs)
  xc   = rootpi*decomp.beta^2*total(xweight*decomp.coeffs)
  yc   = rootpi*decomp.beta^2*total(yweight*decomp.coeffs)
  if flux ne 0. then centroid=[xc,yc]/flux else centroid=[!values.f_nan,!values.f_nan]
  
  ; Calculate errors on the centroid
  if keyword_set(error) then begin
    flux_error=rootpi*decomp.beta*sqrt(total((fweight*decomp.coeffs_error)^2))
    if flux eq 0. then centroid_error=[!values.f_infinity,!values.f_infinity] else begin
      xc_error = rootpi*decomp.beta^2*sqrt(total((xweight*decomp.coeffs_error)^2))
      yc_error = rootpi*decomp.beta^2*sqrt(total((yweight*decomp.coeffs_error)^2))
      centroid_error=fltarr(2,/nozero)
      if xc eq 0 then xc_error=!values.f_infinity else $
        centroid_error[0] = abs(centroid[0]) * sqrt((xc_error/xc)^2+(flux_error/flux)^2)
      if yc eq 0 then yc_error=!values.f_infinity else $
        centroid_error[1] = abs(centroid[0]) * sqrt((yc_error/yc)^2+(flux_error/flux)^2)
    endelse
    flux     = [flux,flux_error]
    centroid = [[centroid],[centroid_error]]
  endif

endif else begin

  ;
  ; WORK IN POLAR SHAPELET SPACE
  ;
  ; Convert coefficients to polar shapelet format if not already like that
  shapelets_polar_convert,decomp,/POLAR,/SILENT
  m0=where(decomp.m eq 0,n_m0)
  m1=where(decomp.m eq 1 and decomp.n le n_max,n_m1)

  if n_m1 eq 0 then begin
    if arg_present(flux) then flux=shapelets_flux(decomp,cartesian=cartesian,polar=polar,error=error)
    centroid=[0.,0.]
    if keyword_set(error) then centroid = [[centroid],[0.,0.]]
  endif else if n_m0 eq 0 then begin
    if arg_present(flux) then begin
      flux=[0.,0.]
      if keyword_set(error) then flux=[[flux],[0.,0.]]
    endif
    centroid=[!values.f_nan,!values.f_nan]
    if keyword_set(error) then centroid=[[centroid],[!values.f_infinity,!values.f_infinity]]
  endif else begin
  
    ; Calculate flux and centroid
    flux  = rootpi*decomp.beta*2*total(float(decomp.coeffs[m0]))
    weight= sqrt(decomp.n+1)
    xc    = rootpi*decomp.beta^2*sqrt(8)*total((float(decomp.coeffs)*weight)[m1])
    yc    = rootpi*decomp.beta^2*sqrt(8)*total((imaginary(decomp.coeffs)*weight)[m1])
    ;				      ^--- 2 from weight calculation, 2 from m=1 and m=-1, 1/root(2) from inside summation
    if flux ne 0. then centroid=[xc,yc]/flux else centroid=[!values.f_nan,!values.f_nan]

    ; Calculate errors on the centroid
    if keyword_set(error) then begin
      flux_error = rootpi*decomp.beta*2*sqrt(total((float(decomp.coeffs_error[m0]))^2))
      if size(decomp.coeffs_error,/n_dimensions) ne 1 then message,$
   	"Warning: error calculation using covariance matrix is not yet written!"
      if flux eq 0. then cen_error=[!values.f_infinity,!values.infinity] else begin
   	xc_error = rootpi*decomp.beta^2*sqrt(8)*sqrt(total((float(decomp.coeffs)*weight)[m1])^2)
   	yc_error = rootpi*decomp.beta^2*sqrt(8)*sqrt(total((imaginary(decomp.coeffs)*weight)[m1])^2)
   	centroid_error=fltarr(2,/nozero)
   	if xc eq 0 then xc_error=!values.f_infinity else $
   	  centroid_error[0] = abs(centroid[0]) * sqrt((xc_error/xc)^2+(flux_error/flux)^2)
   	if yc eq 0 then yc_error=!values.f_infinity else $
   	  centroid_error[1] = abs(centroid[0]) * sqrt((yc_error/yc)^2+(flux_error/flux)^2)
      endelse
      flux     = [flux,flux_error]
      centroid = [[centroid],[centroid_error]]
    endif
  
  endelse

endelse

; Restore input variable
decomp=decomp_input

; Return a position relative to the bottom left of the postage stamp,
;   rather than the centre of the shapelet basis functions
if not keyword_set(basis_functions) then if n_elements(decomp.x) gt 0 then $
  centroid[*,0]+=decomp.x
if keyword_set(image_coords) then $
  if tag_exist(decomp,"sex") then $
    if n_elements(decomp.sex.xmin) gt 0 then $
      centroid[*,0]+=decomp.sex.xmin

return,centroid

end
