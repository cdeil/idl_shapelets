pro shapelets_flex, structure, flexion1, flexion2, $
                    CARTESIAN=cartesian,           $
                    POLAR=polar,                   $
                    EXTEND=extend,                 $
                    ORDER=order,                   $
                    CENTROID=centroid,             $
                    DERIVATIVES=derivatives,       $
                    NOHISTORY=nohistory,           $
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
;      SHAPELETS_FLEXION
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Applies a flexion (the slight bending due to gradients in a shear field)
;      to one object in a decomp structure, or to all objects in a shapecat
;      catalogue.
;
; INPUTS:
;      STRUCTURE - A shapelet decomp or shapecat structure.
;      FLEXION1  - Spatial derivatives of shear in complex number format.
;                  Flexion1=(dgamma1/dx-dgamma2/dy)+i(dgamma1/dx-dgamma2/dy).
;      FLEXION2  - Spatial derivatives of shear in complex number format.
;                  Flexion2=(dgamma1/dx-dgamma2/dy)+i(dgamma1/dx-dgamma2/dy).
;
;[gamma_1,gamma_2] or complex(gamma1,gamma2) format.
;
; OPTIONAL INPUTS:
;      EXTEND    - Increses n_max before starting to shear, and set the new
;                  coefficients to zero. This gives the shear somewhere to go.
;      ORDER     - Include additional orders of shear operator by setting this
;                  to an integer higher than 1. Default is to apply only those
;                  terms of order gamma.
;
; KEYWORD PARAMETERS:
;      CENTROID  - Keep the centroid the same. If this is not set, applying
;                  flexion can shift the centroid of an object. Set this to
;                  provide a counter-shift. It can be equivalently seen as
;                  operating in the post-flexion coordinate system.
;
; OUTPUTS:
;      STRUCTURE - the structure is returned, with flexion.
;
; MODIFICATION HISTORY:
;      Dec 05 - Factor of two fixed, and CENTROID option added, by RM.
;      Nov 05 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Maintain backwards compatibility
if not shapelets_structure_type(structure,message=message) then message,message

; Parse inputs
if not keyword_set(order) then order=1
if not keyword_set(cartesian) then cartesian=0B
if not keyword_set(polar) then polar=0B
if keyword_set(derivatives) then begin
  S11=flexion1[0]
  S12=flexion1[1]
  S21=flexion2[0]
  S22=flexion2[1]
  F=complex(S11+S22,S21-S12)
  G=complex(S11-S22,S21+S12)
endif else begin
  F=flexion1
  G=flexion2
  S11=float(F+G)/2.
  S12=imaginary(G-F)/2.
  S21=imaginary(F+G)/2.
  S22=float(F-G)/2.
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
  shapelets_exponentiate_operations, "shapelets_flex", $
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
    tworoottwo=2*sqrt(2)
    fourroottwo=4*sqrt(2)
    for i=0,n_coeffs-1 do begin
  
      ; S^(2)_11
      j=where(n1 eq n1[i]-3 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S11/tworoottwo*sqrt(n1[i]*(n1[i]-1)*(n1[i]-2))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+3 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S11/tworoottwo*sqrt((n1[i]+1)*(n1[i]+2)*(n1[i]+3))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+1 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S11/tworoottwo*(1-n1[i])*sqrt(n1[i]+1)*old_coeffs[j,*]

      j=where(n1 eq n1[i]-1 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S11/tworoottwo*(2+n1[i])*sqrt(n1[i])*old_coeffs[j,*]
  
      ; S^(2)_12
      j=where(n1 eq n1[i] and n2 eq n2[i]-3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S12/tworoottwo*sqrt(n2[i]*(n2[i]-1)*(n2[i]-2))*old_coeffs[j,*]

      j=where(n1 eq n1[i] and n2 eq n2[i]+3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S12/tworoottwo*sqrt((n2[i]+1)*(n2[i]+2)*(n2[i]+3))*old_coeffs[j,*]

      j=where(n1 eq n1[i] and n2 eq n2[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S12/tworoottwo*(1-n2[i])*sqrt(n2[i]+1)*old_coeffs[j,*]

      j=where(n1 eq n1[i] and n2 eq n2[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S12/tworoottwo*(2+n2[i])*sqrt(n2[i])*old_coeffs[j,*]

      ; S^(2)_21
      j=where(n1 eq n1[i] and n2 eq n2[i]-3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S21/fourroottwo*sqrt(n2[i]*(n2[i]-1)*(n2[i]-2))*old_coeffs[j,*]

      j=where(n1 eq n1[i] and n2 eq n2[i]+3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S21/fourroottwo*sqrt((n2[i]+1)*(n2[i]+2)*(n2[i]+3))*old_coeffs[j,*]

      j=where(n1 eq n1[i]-2 and n2 eq n2[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S21/fourroottwo*3*sqrt(n2[i]*n1[i]*(n1[i]-1))*old_coeffs[j,*]

      j=where(n1 eq n1[i]-2 and n2 eq n2[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S21/fourroottwo*sqrt((n2[i]+1)*n1[i]*(n1[i]-1))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+2 and n2 eq n2[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S21/fourroottwo*sqrt(n2[i]*(n1[i]+1)*(n1[i]+2))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+2 and n2 eq n2[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S21/fourroottwo*3*sqrt((n2[i]+1)*(n1[i]+1)*(n1[i]+2))*old_coeffs[j,*]

      j=where(n1 eq n1[i] and n2 eq n2[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S21/fourroottwo*(2-n2[i]-2*n1[i])*sqrt(n2[i]+1)*old_coeffs[j,*]

      j=where(n1 eq n1[i] and n2 eq n2[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S21/fourroottwo*(5+n2[i]+2*n1[i])*sqrt(n2[i])*old_coeffs[j,*]

      ; S^(2)_22
      j=where(n1 eq n1[i]-3 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S22/fourroottwo*sqrt(n1[i]*(n1[i]-1)*(n1[i]-2))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+3 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S22/fourroottwo*sqrt((n1[i]+1)*(n1[i]+2)*(n1[i]+3))*old_coeffs[j,*]

      j=where(n1 eq n1[i]-1 and n2 eq n2[i]-2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S22/fourroottwo*3*sqrt(n1[i]*n2[i]*(n2[i]-1))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+1 and n2 eq n2[i]-2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S22/fourroottwo*sqrt((n1[i]+1)*n2[i]*(n2[i]-1))*old_coeffs[j,*]

      j=where(n1 eq n1[i]-1 and n2 eq n2[i]+2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S22/fourroottwo*sqrt(n1[i]*(n2[i]+1)*(n2[i]+2))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+1 and n2 eq n2[i]+2,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-S22/fourroottwo*3*sqrt((n1[i]+1)*(n2[i]+1)*(n2[i]+2))*old_coeffs[j,*]

      j=where(n1 eq n1[i]+1 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S22/fourroottwo*(2-n1[i]-2*n2[i])*sqrt(n1[i]+1)*old_coeffs[j,*]

      j=where(n1 eq n1[i]-1 and n2 eq n2[i],n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+S22/fourroottwo*(5+n1[i]+2*n2[i])*sqrt(n1[i])*old_coeffs[j,*]
  
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
    sixteenroottwo=16*sqrt(2)
    for i=0,n_coeffs-1 do begin

      ; FIRST FLEXION (curly F)
      ; Expand upwards in polar shapelet coefficient space, from lower m coefficients
      j=where(n eq n[i]-3 and m eq m[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+F/sixteenroottwo*3*sqrt((n[i]-m[i])*(n[i]+m[i]-2)*(n[i]+m[i]))*old_coeffs[j,*]

      j=where(n eq n[i]-1 and m eq m[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+F/sixteenroottwo*(10+3*n[i]-m[i])*sqrt(n[i]+m[i])*old_coeffs[j,*]
      
      j=where(n eq n[i]+1 and m eq m[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+F/sixteenroottwo*(4-3*n[i]-m[i])*sqrt(n[i]-m[i]+2)*old_coeffs[j,*]
      
      j=where(n eq n[i]+3 and m eq m[i]-1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-F/sixteenroottwo*3*sqrt((n[i]-m[i]+4)*(n[i]-m[i]+2)*(n[i]+m[i]+2))*old_coeffs[j,*]

      ; Expand downwards in polar shapelet coefficient space, from higher m coefficients
      j=where(n eq n[i]-3 and m eq m[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+conj(F)/sixteenroottwo*3*sqrt((n[i]-m[i]-2)*(n[i]-m[i])*(n[i]+m[i]))*old_coeffs[j,*]

      j=where(n eq n[i]-1 and m eq m[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+conj(F)/sixteenroottwo*(10+3*n[i]+m[i])*sqrt(n[i]-m[i])*old_coeffs[j,*]
      
      j=where(n eq n[i]+1 and m eq m[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+conj(F)/sixteenroottwo*(4-3*n[i]+m[i])*sqrt(n[i]+m[i]+2)*old_coeffs[j,*]
      
      j=where(n eq n[i]+3 and m eq m[i]+1,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-conj(F)/sixteenroottwo*3*sqrt((n[i]-m[i]+2)*(n[i]+m[i]+4)*(n[i]+m[i]+2))*old_coeffs[j,*]


      ; SECOND FLEXION (curly G)
      ; Expand upwards in polar shapelet coefficient space, from lower m coefficients
      j=where(n eq n[i]-3 and m eq m[i]-3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+G/sixteenroottwo*sqrt((n[i]+m[i])*(n[i]+m[i]-2)*(n[i]+m[i]-4))*old_coeffs[j,*]

      j=where(n eq n[i]-1 and m eq m[i]-3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+G/sixteenroottwo*sqrt((n[i]+m[i])*(n[i]+m[i]-2)*(n[i]-m[i]+2))*old_coeffs[j,*]
      
      j=where(n eq n[i]+1 and m eq m[i]-3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-G/sixteenroottwo*sqrt((n[i]+m[i])*(n[i]-m[i]+2)*(n[i]-m[i]+4))*old_coeffs[j,*]
      
      j=where(n eq n[i]+3 and m eq m[i]-3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-G/sixteenroottwo*sqrt((n[i]-m[i]+2)*(n[i]-m[i]+4)*(n[i]-m[i]+6))*old_coeffs[j,*]

      ; Expand downwards in polar shapelet coefficient space, from higher m coefficients
      j=where(n eq n[i]-3 and m eq m[i]+3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+conj(G)/sixteenroottwo*sqrt((n[i]-m[i])*(n[i]-m[i]-2)*(n[i]-m[i]-4))*old_coeffs[j,*]

      j=where(n eq n[i]-1 and m eq m[i]+3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]+conj(G)/sixteenroottwo*sqrt((n[i]-m[i])*(n[i]-m[i]-2)*(n[i]+m[i]+2))*old_coeffs[j,*]
      
      j=where(n eq n[i]+1 and m eq m[i]+3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-conj(G)/sixteenroottwo*sqrt((n[i]-m[i])*(n[i]+m[i]+2)*(n[i]+m[i]+4))*old_coeffs[j,*]
      
      j=where(n eq n[i]+3 and m eq m[i]+3,n_j)
      if n_j gt 0 then new_coeffs[i,*]=new_coeffs[i,*]-conj(G)/sixteenroottwo*sqrt((n[i]+m[i]+2)*(n[i]+m[i]+4)*(n[i]+m[i]+6))*old_coeffs[j,*]

    endfor  
 
  endelse

  ; Reinsert new coefficients into the old structure
  if structure.type eq "shapecat" then begin
    new_coeffs=transpose(new_coeffs)
    if tag_exist(structure,"moments") then if structure.moments then $
      message,"The moments in your shapecat need updating!",/info
  endif
  structure.coeffs=new_coeffs

  ; Keep the object centred in the same place
  if keyword_set(centroid) then begin
    beta=structure.beta
    if n_elements(beta) gt 1 then message,"shapelets_translate is going to fail! It needs to be rewritten to cope with different shifts for each object in a shapecat."
    junk=shapelets_quadrupole(structure,ellipticity=ellipticity,rsquared=rsquared)
    shift=-rsquared/4./beta*(6*F+5*conj(F)*ellipticity+G*conj(ellipticity))
    shapelets_translate,structure,[float(shift),imaginary(shift)],POLAR=polar,CARTESIAN=cartesian
 endif
   
  ; Recover input format
  if keyword_set(maintain) then $
    shapelets_polar_convert,structure,POLAR=polar_input,CARTESIAN=1-polar_input,/SILENT

endelse

; Add operation to object's history record
;if not keyword_set(nohistory) then shapelets_update_history, structure, $
;  "Flexion of {"+strmid(strtrim(gamma1,2),0,5)+","+strmid(strtrim(gamma2,2),0,5)+"} applied."

end

