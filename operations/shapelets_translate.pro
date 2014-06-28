pro shapelets_translate, structure, translation,$
                         EXTEND=extend,         $
                         CARTESIAN=cartesian,   $
                         POLAR=polar,           $
                         ORDER=order,           $
                         NOHISTORY=nohistory,   $
                         MAINTAIN=maintain

;$Id: shapelets_translate.pro, v2$
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
;      SHAPELETS_TRANSLATE
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Translates an object (to first order, using ladder operators) some
;      number of pixels.
; 
; INPUTS:
;      STRUCTURE   - A shapelet decomp or shapecat structure.
;      TRANSLATION - [x,y] floating point amount by which to shift object, in
;                    units of pixels.
;
; OPTIONAL INPUTS:
;      EXTEND      - Integer value by which to increase n_max by padding with
;                    zeros before performing shft. If extend=1, default is 
;                    n_max=n_max+4.
;      ORDER       - Include additional orders of translation operator by setting
;                    this to an integer higher than 1. Default is to apply only
;                    terms of first order in the shift.
;
; KEYWORD PARAMETERS:
;      POLAR       - Perform translation in polar shapelet space. 
;                    DEFAULT: Cartesian.
;
; OUTPUTS:
;      STRUCTURE   - Shifted shapelet decomp or shapecat structure.
; 
; MODIFICATION HISTORY:
;      Apr 05 - Input now needs to be specified in units of pixels by RM.
;      Apr 05 - Higher order transformation added by RM.
;      Apr 02 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Maintain backwards compatibility
if not shapelets_structure_type(structure,message=message) then message,message
if not keyword_set(order) then order=1
polar_input=structure.polar

; Perform translation

; Determine structure type
if structure.type eq "shapecat" then begin 

  ; Extend n_max first if requested
  if keyword_set(extend) then shapelets_extend_nmax,structure,extend

  structure.x[*,0]=structure.x[*,0]+translation[0]
  structure.x[*,1]=structure.x[*,1]+translation[1]
  if structure.moments then begin
    structure.centroid[*,0]=structure.centroid[*,0]+translation[0]
    structure.centroid[*,1]=structure.centroid[*,1]+translation[1]
  endif
  ;if structure.sextractor then begin
  ;  structure.sexx[*]=structure.sexx[*]+translation[0]*structure.beta
  ;  structure.sexy[*]=structure.sexy[*]+translation[1]*structure.beta
  ;endif

endif else if structure.type eq "decomp" or structure.type eq "decomp_cartesian" then begin 

  ; Decide which method to use
  if ((keyword_set(cartesian)+keyword_set(polar)) mod 2) eq 0 then begin
    polar=structure.polar
    cartesian=1-polar
  endif
  
  ; Perform translation to higher than first order
  if fix(order) gt 1 then begin
    if keyword_set(cartesian) and keyword_set(polar) then begin
      polar=structure.polar
      cartesian=1-structure.polar
    endif else shapelets_polar_convert,structure,cartesian=cartesian,polar=polar,/SILENT
    shapelets_exponentiate_operations, "shapelets_translate", $
      structure, translation, order, cartesian=cartesian, polar=polar, extend=extend
    if keyword_set(maintain)  and structure.polar ne polar_input then $
      shapelets_polar_convert,structure,cartesian=1-polar_input,polar=polar_input,/SILENT
  endif else begin
  
    ; Extend n_max first if requested
    if keyword_set(extend) then shapelets_extend_nmax,structure,extend
    
    ; Perform translation
    if not keyword_set(polar) then begin
    
      ; Obtain Cartesian shapelet coefficients
      shapelets_polar_convert,structure,/CARTESIAN,/SILENT
    
      ; Perform translation using Cartesian step operators
      shapelets_make_nvec,structure.n_max,n1,n2
      new_coeffs=structure.coeffs
      for i=0,structure.n_coeffs-1 do begin
    
  	; delta_x
  	; Expand rightwards in Cartesian shapelet coefficient space
  	j=where(n1 eq n1[i]+1 and n2 eq n2[i])
  	if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]-translation[0]/structure.beta*sqrt((n1[i]+1.)/2.)*structure.coeffs[j[0]]
    
  	; Expand leftwards in Cartesian shapelet coefficient space
  	j=where(n1 eq n1[i]-1 and n2 eq n2[i])
  	if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]+translation[0]/structure.beta*sqrt(float(n1[i])/2.)*structure.coeffs[j[0]]
    
  	; delta_y
  	; Expand upwards in Cartesian shapelet coefficient space
  	j=where(n1 eq n1[i] and n2 eq n2[i]+1)
  	if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]-translation[1]/structure.beta*sqrt((n2[i]+1.)/2.)*structure.coeffs[j[0]]
    
  	; Expand downwards in Cartesian shapelet coefficient space
  	j=where(n1 eq n1[i] and n2 eq n2[i]-1)
  	if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]+translation[1]/structure.beta*sqrt(float(n2[i])/2.)*structure.coeffs[j[0]]
  
      endfor  
    
      ; Insert new coefficients back into the decomp structure
      structure.coeffs=new_coeffs
     
      ; Recover input format
      if keyword_set(maintain) and polar_input then $
        shapelets_polar_convert,structure,/POLAR,/SILENT
  
    endif else begin
    
      ; Obtain polar shapelet coefficients
      shapelets_polar_convert,structure,/POLAR,/SILENT
      shapelets_make_nvec,structure.n_max,n,m,/POLAR
      new_polar_coeffs=structure.coeffs
      
      ; Perform translation in polar shapelet space
      for i=0,structure.n_coeffs-1 do begin
  
        j=where(m eq m[i]-1 and n eq n[i]-1)
        if j[0] ne -1 then new_polar_coeffs[i] = new_polar_coeffs[i] + $
             ( complex(translation[0],translation[1])/structure.beta * $
               sqrt(float(n[i]+m[i])/8.) * structure.coeffs[j[0]] ) 
    
        j=where(m eq m[i]+1 and n eq n[i]-1)
        if j[0] ne -1 then new_polar_coeffs[i] = new_polar_coeffs[i] + $
             ( complex(translation[0],-1*translation[1])/structure.beta * $
               sqrt(float(n[i]-m[i])/8.) * structure.coeffs[j[0]] ) 
    
        j=where(m eq m[i]-1 and n eq n[i]+1)
        if j[0] ne -1 then new_polar_coeffs[i] = new_polar_coeffs[i] - $
             ( complex(translation[0],translation[1])/structure.beta * $
               sqrt(float(n[i]-m[i]+2.)/8.) * structure.coeffs[j[0]] ) 
  
        j=where(m eq m[i]+1 and n eq n[i]+1)
        if j[0] ne -1 then new_polar_coeffs[i] = new_polar_coeffs[i] - $
             ( complex(translation[0],-1*translation[1])/structure.beta * $
               sqrt(float(n[i]+m[i]+2.)/8.) * structure.coeffs[j[0]] ) 
  
      endfor
    
      ; Insert new coefficients back into the decomp structure
      structure.coeffs=new_polar_coeffs
   
      ; Convert back to Cartesian shapelet coefficients
      if keyword_set(maintain) and not polar_input then $
        shapelets_polar_convert,structure,/CARTESIAN,/SILENT
   
    endelse

  endelse
  
endif else message,"Cannot translate objects in "+structure.type+" structures!"

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,structure,$
  "Translated by ["+strtrim(translation[0],2)+","+strtrim(translation[1],2)+"]."

end

