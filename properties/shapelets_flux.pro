function shapelets_flux,structure_input, $
                        MATRIX=matrix,   $
                        ERROR=error,     $
                        POLAR=polar,     $
                        CARTESIAN=cartesian
 
;$Id: shapelets_flux.pro, v2$
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
;      SHAPELETS_FLUX
;
; PURPOSE:
;      Returns total flux of a Cartesian decomp or shapecat structure.
;      This does the same as shapelets_moments.pro but in function form.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      structure: A Cartesian shapelet structureosition.
;
; KEYWORD PARAMETERS:
;     [/CART] - Perform summation in Cartesian shapelet space. [DEFAULT]
;    [/POLAR] - Perform sunnation in polar shapelet space.
;    [/ERROR] - Simultaneously calculate the error on the flux.
;
; OUTPUTS:
;        flux - A floating point variable containing objects' total flux.
;               If /ERROR is set, this becomes a vector [[flux],[flux_error]]
;
; EXAMPLE USE:
;      result = shapelets_flux(decomp)
;
; NOTES: 
;      Not yet compatible with a full covariance error matrix of coefficients.
;
; MODIFICATION HISTORY:
;      Jul 05 - Ability to accept many shapelet structures added by RM.
;      Sep 03 - Calculation of error on flux added by RM.
;      Feb 02 - Written by Richard Massey.
;-

COMPILE_OPT idl2

;
; SET DEFAULTS
;
; Backwardly compatible
if not shapelets_structure_type(structure_input,message=message) then message,message
; Don't upset the input variable
structure=structure_input
; Deal with easy cases first
if structure.type eq "sexcat" then begin
  flux=structure.flux
endif else if structure.type eq "pstamp" then begin
  flux=total(pstamp.image)
endif else begin
  ; Decide whether default method should be in Cartesian or polar shapelet space
  if not keyword_set(cartesian) then cartesian=0B
  if not keyword_set(polar) then polar=0B
  if not ((cartesian+polar) mod 2) then begin
    if structure.polar then begin
      cartesian=0B
      polar=1B
    endif else begin
      cartesian=1B
      polar=0B
    endelse
  endif
  ; Convert coefficients to that type
  shapelets_polar_convert,structure,POLAR=polar,CARTESIAN=cartesian,/SILENT
  ; Get coefficients into a local variable, in the desired form
  if structure.type eq "decomp" then begin
    coeffs=structure.coeffs
    coeffs_error=structure.coeffs_error
    n_max=structure.n_max
    n_objects=1
  endif else if structure.type eq "shapecat" then begin
    coeffs=transpose(structure.coeffs)
    coeffs_error=transpose(structure.coeffs_error)
    n_max=structure.maxn_max
    n_objects=structure.n
  endif else message,"Cannot calculate flux of "+structure.type+" structures!"

  if polar then begin 
    
    ;
    ; WORK IN POLAR SHAPELET SPACE
    ;
    
    ; Sum m=0 coefficients to calculate flux
    shapelets_make_nvec,n_max,n,m,/POLAR
    m0=where(m eq 0,n_m0)
    constant=2*sqrt(!pi)*structure.beta
    if n_m0 gt 0 then flux=constant*total(float(coeffs[m0,*]),1) else flux=fltarr(n_objects)
  
    ; Calculate error on flux
    if keyword_set(error) then begin
      if n_m0 gt 0 then flux_error = constant*sqrt(total((float(coeffs_error[m0,*]))^2,1)) else flux_error=fltarr(n_objects)
      flux=[[flux],[flux_error]]
    endif
  
  endif else begin
  
    ;
    ; WORK IN CARTESIAN SHAPELET SPACE  
    ;
  	
    ; Calculate weights
    shapelets_make_nvec,n_max,n1,n2,n_coeffs
    constant=sqrt(!pi)*structure.beta
    weight=fltarr(n_coeffs)
    eveneven=where(n1 mod 2 eq 0 and n2 mod 2 eq 0,n_eveneven)
    if n_eveneven gt 0 then weight[eveneven]=2.^(.5*(2.-n1[eveneven]-n2[eveneven])) $
  			 *sqrt(factorial(n1[eveneven])*factorial(n2[eveneven]))     $
  			 /factorial(n1[eveneven]/2)/factorial(n2[eveneven]/2)
    weight=weight#replicate(1,n_objects)

    ; Calculate flux
    flux=constant*total(weight*coeffs,1)
    
    ; Calculate error on flux
    if keyword_set(error) then begin
      flux_error=constant*sqrt(total((weight*coeffs_error)^2,1))
      flux=[[flux],[flux_error]]
    endif
  
  endelse
endelse
  
return,flux

end



