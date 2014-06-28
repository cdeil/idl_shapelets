pro shapelets_shear, structure, gammas,   $
                     CARTESIAN=cartesian, $
                     POLAR=polar,         $
                     EXTEND=extend,       $
                     ORDER=order,         $
                     NOHISTORY=nohistory, $
		                 MAINTAIN=maintain

;$Id: shapelets_shear.pro, v2$
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
;      Applies a (pure) shear to one object in a decomp structure, or to all
;      objects in a shapecat catalogue.
;
; INPUTS:
;      STRUCTURE - A shapelet decomp or shapecat structure.
;      GAMMAS    - Shear in [gamma_1,gamma_2] or complex(gamma1,gamma2) format.
;
; OPTIONAL INPUTS:
;      EXTEND    - Increses n_max before starting to shear, and set the new
;                  coefficients to zero. This gives the shear somewhere to go.
;      ORDER     - Include additional orders of shear operator by setting this
;                  to an integer higher than 1. Default is to apply only those
;                  terms of order gamma.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      STRUCTURE - the structure is returned, sheared.
;
; TO DO:
;      Accept sexcats as an input. This would be an interesting exercise.
;      Properly transform the errors on the shapelet coefficients.
;
; MODIFICATION HISTORY:
;      Jul 05 - Sign error fixed in polar implementation by RM.
;      Apr 05 - Shapecats accepted as input by RM.
;      Apr 05 - POLAR keyword added by RM.
;      Apr 05 - Higher order transformation added properly by RM.
;      Feb 05 - IDL2 notation (square brackets, etc) compatability by RM.
;      Aug 04 - Approximate higher order shear calculation added by RM.
;      Apr 02 - Written by Richard Massey
;-

COMPILE_OPT idl2

; Maintain backwards compatibility
if not shapelets_structure_type(structure,message=message) then message,message

; Parse inputs
if not keyword_set(order) then order=1
if not keyword_set(cartesian) then cartesian=0B
if not keyword_set(polar) then polar=0B
if size(gammas,/type) eq 6 then begin
  gamma1=float(gammas[0])
  gamma2=imaginary(gammas[0])
  gammac=gammas[0]
endif else begin
  gamma1=gammas[0]
  if n_elements(gammas) ge 2 then gamma2=gammas[1] else gamma2=0
  gammac=complex(gamma1,gamma2)
endelse
polar_input=structure.polar

; Decide whether default method should be in Cartesian or polar shapelet space
if not ((cartesian+polar) mod 2) then begin
  if structure.polar then begin
    cartesian=0B
    polar=1B
  endif else begin
    cartesian=1B
    polar=0B
  endelse
endif

; Perform shear to higher than first order
if fix(order) gt 1 then begin

  shapelets_polar_convert,structure,CARTESIAN=cartesian,POLAR=polar,/SILENT
  shapelets_exponentiate_operations, "shapelets_shear", $
    structure, gammas, order, CARTESIAN=cartesian, POLAR=polar, EXTEND=extend
  if keyword_set(maintain) then $
    shapelets_polar_convert,structure,CARTESIAN=1-polar_input,POLAR=polar_input,/SILENT

endif else begin

  ; Increase n_max to contain new coefficients that might be created
  if keyword_set(extend) then shapelets_extend_nmax,structure,extend
  
  ; Perform shear to first order
  if cartesian then begin
  
    ; Perform shear in Cartesian shapelet space
  
    ; Obtain initial shapelet coefficients
    shapelets_polar_convert,structure,/CARTESIAN,/SILENT
    if structure.type eq "shapecat" then begin
      old_coeffs=transpose(structure.coeffs)
      shapelets_make_nvec,structure.maxn_max,n1,n2,n_coeffs
    endif else if structure.type eq "decomp" then begin
      old_coeffs=structure.coeffs
      shapelets_make_nvec,structure.n_max,n1,n2,n_coeffs
      n1=structure.n1
      n2=structure.n2
      n_coeffs=structure.n_coeffs
    endif else message,"Structure type not recognised!"
    new_coeffs=old_coeffs

    ; Loop over all of the shapelet coefficients
    for i=0,n_coeffs-1 do begin
  
      ; S_1
      ; Expand leftwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i]+2 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-gamma1/2.*sqrt((n1[i]+1)*(n1[i]+2))*old_coeffs[j,*]
  
      ; Expand downwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i] and n2 eq n2[i]+2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+gamma1/2.*sqrt((n2[i]+1)*(n2[i]+2))*old_coeffs[j,*]
  
      ; Expand rightwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i]-2 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+gamma1/2.*sqrt(n1[i]*(n1[i]-1))*old_coeffs[j,*]
  
      ; Expand upwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i] and n2 eq n2[i]-2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-gamma1/2.*sqrt(n2[i]*(n2[i]-1))*old_coeffs[j,*]
     
      ; S_2
      ; Expand diagonally down-left in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i]+1 and n2 eq n2[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-gamma2*sqrt((n1[i]+1)*(n2[i]+1))*old_coeffs[j,*]
  
      ; Expand diagonally up-right in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i]-1 and n2 eq n2[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+gamma2*sqrt(n1[i]*n2[i])*old_coeffs[j,*]
  
    endfor  
   
  endif else begin
  
    ; Perform shear in polar shapelet space
  
    ; Obtain polar shapelet coefficients
    shapelets_polar_convert,structure,/POLAR,/SILENT
    if structure.type eq "shapecat" then begin
      old_coeffs=transpose(structure.coeffs)
      shapelets_make_nvec,structure.maxn_max,n,m,n_coeffs,/POLAR
    endif else if structure.type eq "decomp" then begin
      old_coeffs=structure.coeffs
      ;shapelets_make_nvec,structure.n_max,n,m,n_coeffs,/POLAR
      n=structure.n
      m=structure.m
      n_coeffs=structure.n_coeffs
    endif else message,"Structure type not recognised!"
    new_coeffs=old_coeffs
  
    ; Loop over all of the shapelet coefficients
    for i=0,n_coeffs-1 do begin

      ; Expand up-right in polar shapelet coefficient space
      j=where(n eq n[i]-2 and m eq m[i]-2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+gammac/4.*sqrt((n[i]+m[i])*(n[i]+m[i]-2))*old_coeffs[j,*]
  
      ; Expand down-right in polar shapelet coefficient space
      j=where(n eq n[i]+2 and m eq m[i]-2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-gammac/4.*sqrt((n[i]-m[i]+2)*(n[i]-m[i]+4))*old_coeffs[j,*]
  
      ; Expand up-left in polar shapelet coefficient space
      j=where(n eq n[i]-2 and m eq m[i]+2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+conj(gammac)/4.*sqrt((n[i]-m[i])*(n[i]-m[i]-2))*old_coeffs[j,*]
  
      ; Expand down-left in polar shapelet coefficient space
      j=where(n eq n[i]+2 and m eq m[i]+2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-conj(gammac)/4.*sqrt((n[i]+m[i]+2)*(n[i]+m[i]+4))*old_coeffs[j,*]

    endfor  
 
  endelse

  ; Reinsert new coefficients into the old structure
  if structure.type eq "shapecat" then begin
    new_coeffs=transpose(new_coeffs)
    if tag_exist(structure,"moments") then if structure.moments then $
      message,"The moments in your shapecat need updating!",/info
  endif
  structure.coeffs=new_coeffs
   
  ; Recover input format
  if keyword_set(maintain) then $
    shapelets_polar_convert,structure,POLAR=polar_input,CARTESIAN=1-polar_input,/SILENT

endelse

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history, structure, $
  "Sheared by {"+strmid(strtrim(gamma1,2),0,5)+","+strmid(strtrim(gamma2,2),0,5)+"}"

end

