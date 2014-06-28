pro shapelets_extend_nmax, structure, n_increase, SILENT=silent

;$Id: shapelets_extend_nmax.pro, v2$
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
;      SHAPELETS_EXTEND_NMAX
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Increases the value of n_max in a structure by n_increase, by padding
;      the coefficient arrays with zeros. Can also decrease n_max by truncating
;      the coefficient arrays.
; 
; INPUTS:
;      STRUCTURE  - A shapelet decomp structure or shapecat structure.
;      N_INCREASE - (Positive or negative) integer amount by which to change
;                   n_max.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      SILENT     - Operate silently.
;
; OUTPUTS:
;      STRUCTURE  - The same structure, but with a new n_max.
; 
; MODIFICATION HISTORY:
;      Sep 05 - Bug fixed in creation of m array for polar decomps by RM.
;      Apr 05 - Single routine accepts decomp or shapecat structures by RM.
;      Aug 04 - Rendered compatible with polar shapelet catalogues by RM.
;      Jun 04 - Made robust to changes of tag names in structures by RM.
;      May 04 - Shapecat version written by RM.
;      Apr 02 - Decomp version written by Richard Massey.
;-

COMPILE_OPT idl2

; Parse input requirements
if not shapelets_structure_type(structure,message=message) then message,message
if not keyword_set(n_increase) then begin
  message,'Not changing n_max!',/info,noprint=silent
  return
  n_increase=0
endif

; Determine structure type
if structure.type eq "decomp" then begin 

  ; Determine how many coefficients we want to end up with
  n_max = structure.n_max+n_increase
  shapelets_make_nvec,n_max,n1,n2,n_coeffs

  ; Apply a rather ad hoc rule to increase the number of pixels in any subsequent recomp
  sizeinc = sqrt((n_max+1)/(structure.n_max+1))  
  sizeinc = 0
  
  ; Pad coefficient arrays with zeros and change n_max 
  structure_new={name:structure.name,type:structure.type}
  zero=0. & if tag_exist(structure,"polar") then if structure.polar then zero=complex(0.,0.)
  for i=0,n_tags(structure)-1 do begin
    tagname=strupcase((tag_names(structure))[i])
    case strlowcase((tag_names(structure))[i]) of
      "name":
      "type":
      "coeffs":       begin
                        coeffs=structure.coeffs
                        while n_elements(coeffs) lt n_coeffs do coeffs = [coeffs,zero]
                        if n_elements(coeffs) gt n_coeffs then coeffs = coeffs[0:n_coeffs-1]
                        structure_new=create_struct(structure_new,"coeffs",coeffs)
                      end
      "coeffs_error": begin
                        coeffs_error=structure.coeffs_error
                        while n_elements(coeffs_error) lt n_coeffs do coeffs_error = [coeffs_error,zero]
                        if n_elements(coeffs_error) gt n_coeffs then coeffs_error = coeffs_error[0:n_coeffs-1]
                        structure_new=create_struct(structure_new,"coeffs_error",coeffs_error)
                      end
      "error":        begin
                        error=structure.error
                        while n_elements(error) lt n_coeffs do error = [error,zero]
                         if n_elements(error) gt n_coeffs then error = error[0:n_coeffs-1]
                       structure_new=create_struct(structure_new,"error",error)
                      end
      "n_max": structure_new=create_struct(structure_new,"n_max",n_max)
      "n_coeffs": structure_new=create_struct(structure_new,"n_coeffs",n_coeffs)
      "n1": structure_new=create_struct(structure_new,"n1",n1)
      "n2": structure_new=create_struct(structure_new,"n2",n2)
      "n": structure_new=create_struct(structure_new,"n",n1+n2)
      "m": structure_new=create_struct(structure_new,"m",n1-n2)
      "nl": structure_new=create_struct(structure_new,"nl",n2)
      "nr": structure_new=create_struct(structure_new,"nr",n1)
      "x":structure_new=create_struct(structure_new,"x",structure.x + round(sizeinc)*[1.,1.])
      "n_pixels":structure_new=create_struct(structure_new,"n_pixels",structure.n_pixels + 2*round(sizeinc)*[1,1])
      else: structure_new=create_struct(structure_new,tagname,structure.(i))
    endcase
  endfor
  structure=structure_new

endif else if structure.type eq "shapecat" then begin 

  ; Determine how many coefficients we want to end up with
  n_max = structure.maxn_max+n_increase
  shapelets_make_nvec,n_max,n1,n2,n_coeffs
  n=n1+n2

  if n_max lt structure.maxn_max then begin

    ; Truncate coefficients to reduce n_max
    message,'Truncating coefficients beyond '+$
            'n_max='+strtrim(string(n_max),2),/info,noprint=silent
    keep=where(n le n_max,n_coeffs)
    coeffs=structure.coeffs[*,keep]
    coeffs_error=structure.coeffs_error[*,keep]
  
  endif else begin

    ; Pad with zeros to increase n_max
    message,'Padding coefficients with zeros to increase from '+$
            'n_max='+strtrim(string(structure.maxn_max),2)+$
            ' to n_max='+strtrim(string(n_max),2),/info,noprint=silent
    if structure.polar then begin
      coeffs=[[structure.coeffs],[complexarr(structure.n,n_coeffs-structure.n_coeffs)]]
      coeffs_error=[[structure.coeffs_error],[complexarr(structure.n,n_coeffs-structure.n_coeffs)]]
    endif else begin
      coeffs=[[structure.coeffs],[fltarr(structure.n,n_coeffs-structure.n_coeffs)]]
      coeffs_error=[[structure.coeffs_error],[fltarr(structure.n,n_coeffs-structure.n_coeffs)]]
    endelse

  endelse
  
  ; Adjust shapecat struture to reflect new length of coefficient list
  structure_tagnames=tag_names(structure)
  structure_new={name:structure.name,type:structure.type}
  for i=0,n_tags(structure)-1 do begin
    case strupcase(structure_tagnames[i]) of
      "NAME":
      "TYPE":
      "COEFFS": structure_new=create_struct(structure_new,"coeffs",coeffs)
      "COEFFS_ERROR": structure_new=create_struct(structure_new,"coeffs_error",coeffs_error)
      else: structure_new=create_struct(structure_new,structure_tagnames[i],structure.(i))
    endcase
  endfor
  structure=structure_new
  structure.maxn_max=n_max
  structure.n_coeffs=n_coeffs

endif else message,"Cannot extend n_max in "+structure.type+" structures!"




end
