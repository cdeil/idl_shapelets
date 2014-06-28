pro shapelets_add, structure1, structure2, SILENT=silent, NOHISTORY=nohistory

;$Id: shapelets_add.pro, v2$
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
;      SHAPELETS_ADD
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Superimposes two images (stored as shapelet decompositions),
;      append a shapelet decomposition to a shapelet catalogue, or
;      concatenates two catalogues.
;
; INPUTS:
;      STRUCTURE1 - A decomp, shapecat or sexcat structure.
;      STRUCTURE2 - Another decomp, shapecat or sexcat structure.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      SILENT     - Operate silently.
;
; OUTPUTS:
;      The added or concatenated result is returned as STRUCTURE1.
;
; MODIFICATION HISTORY:
;      Jan 06 - Shapecat concatenation bug fixed by RM.
;      Jul 05 - Polar decomp structures incorporated by RM.
;      Apr 05 - Sexcat and shapecat inputs incorporated by RM.
;      Apr 05 - Converted from a function to a procedure by RM.
;      Apr 02 - Writen by Richard Massey.
;-

ON_ERROR,2
COMPILE_OPT idl2

; Parse input variables
recognised1=shapelets_structure_type(structure1,message=message1,/silent)
recognised2=shapelets_structure_type(structure2,message=message2,/silent)
if not (recognised1 and recognised2) then message,"Input structure(s) not recognised: "+message1+" "+message2

; Concatenate SExtractor catalogues
if structure1.type eq "sexcat" then begin 

  if structure2.type ne "sexcat" then message,"Cannot concatenate a SExtractor catalogue with a "+structure2.type+" structure!"
    
  ; Check that the number of tags match
  n_tags1=n_tags(structure1)
  n_tags2=n_tags(structure2)
  if (n_tags1 ne n_tags2) then begin
    message,"The two sexcats do not contain the same number of tags and are therefore not compatible!",/info
  endif else begin

    ; Check that the tag names match
    names1=tag_names(structure1)
    names2=tag_names(structure2)
    for i=0,n_tags1-1 do if (names1[i] ne names2[i]) then begin
      message,"The two sexcats' tags do not match, so they are therefore not compatible",/info
    endif else begin

      ; Begin a new (combined) catalogue
      new_sexcat={name:structure1.name+" and "+structure2.name, type:"sexcat"} 
      
      ; Compile a concatenated list of variables
      for i=0,n_tags1-1 do begin
        case strupcase(names1[i]) of
          "NAME":
          "TYPE":
          "HISTORY":
          "N": new_sexcat=create_struct(new_sexcat,"n",long(structure1.n)+long(structure2.n))
          "SEEING": new_sexcat=create_struct(new_sexcat,"seeing",(structure1.seeing+structure2.seeing)/2.)
          else: new_sexcat=create_struct(new_sexcat,names1[i],[structure1.(i),structure2.(i)])
        endcase
      endfor

      ; Return new catalogue
      structure1=new_sexcat

    endelse
  endelse

; Concatenate shapelet catalogues
endif else if structure1.type eq "shapecat" and structure2.type eq "shapecat" then begin 
 
  ; Check to see the two catalogues are compatible (makes things much simpler later on)
  if structure1.maxn_max gt structure2.maxn_max then begin
    shapelets_extend_nmax,structure2,structure1.maxn_max-structure2.maxn_max
  endif else if structure1.maxn_max lt structure2.maxn_max then begin
    shapelets_extend_nmax,structure1,structure2.maxn_max-structure1.maxn_max
  endif
  
  ; Check that the structure types (polar or Cartesian) match
  if ((structure1.polar+structure2.polar) mod 2 ne 0) then begin
    shapelets_polar_convert, structure2, polar=structure1.polar, cartesian=1-structure1.polar
  endif
  
  ; Check that the number of tags match
  n_tags1=n_tags(structure1)
  n_tags2=n_tags(structure2)
  if (n_tags1 ne n_tags2) then begin
    message,"The two catalogues do not contain the same number of tags and are therefore not compatible!",/info
  endif else begin

    ; Check that the tag names match
    names1=tag_names(structure1)
    names2=tag_names(structure2)
    for i=0,n_tags1-1 do if (names1[i] ne names2[i]) then begin
      message,"The two catalogues' tags do not match, so they are therefore not compatible",/info
    endif else begin

      ; Begin a new (combined) catalogue
      new_shapecat={name:structure1.name+" and "+structure2.name, type:"shapecat"} 
      
      ; Compile a concatenated list of variables
      for i=0,n_tags1-1 do begin
        case strupcase(names1[i]) of
          "NAME":
          "TYPE":
          "HISTORY": new_shapecat=create_struct(new_shapecat,"history",structure1.history)
          "N": new_shapecat=create_struct(new_shapecat,"n",long(structure1.n)+long(structure2.n))
          "N_COEFFS": new_shapecat=create_struct(new_shapecat,"n_coeffs",structure1.n_coeffs)
          "MAXN_MAX": new_shapecat=create_struct(new_shapecat,"maxn_max",structure1.maxn_max)
          "SEEING": new_shapecat=create_struct(new_shapecat,"seeing",(structure1.seeing+structure2.seeing)/2.)
          "POLAR": new_shapecat=create_struct(new_shapecat,"polar",structure1.polar)
          "MOMENTS": new_shapecat=create_struct(new_shapecat,"moments",structure1.moments)
          "SEXTRACTOR": new_shapecat=create_struct(new_shapecat,"sextractor",structure1.sextractor)
          else: new_shapecat=create_struct(new_shapecat,names1[i],[structure1.(i),structure2.(i)])
        endcase
      endfor

      ; Return new catalogue
      structure1=new_shapecat

    endelse
  endelse
  
  ; Update history tag
  if not keyword_set(nohistory) then shapelets_update_history,structure1,"Concatenated with "+structure2.name

; Superpose shapelet decompositions
endif else if structure1.type eq "decomp" and structure2.type eq "decomp" then begin 

  ; Check that the two objects' scale sizes match
  if structure1.beta ne structure2.beta then message,"Scale sizes of the two objects do not match!"

  ; Check that the two objects' truncation orders match
  if structure1.n_max gt structure2.n_max then begin
    shapelets_extend_nmax,structure2,structure1.n_max-structure2.n_max
  endif else if structure1.n_max lt structure2.n_max then begin
    shapelets_extend_nmax,structure1,structure2.n_max-structure1.n_max
  endif

  ; Check that the structure types (polar or Cartesian) match
  if ((structure1.polar+structure2.polar) mod 2 ne 0) then begin
    shapelets_polar_convert, structure2, polar=structure1.polar, cartesian=1-structure1.polar
  endif

  ; Copy decomp structure
  new_decomp=structure1

  ; Add the other's shapelet coefficients (the model is linear)
  new_decomp.coeffs=structure1.coeffs+structure2.coeffs

  ; Do the same to the (statistically independent) errors on the coefficients
  new_decomp.coeffs_error=new_decomp.coeffs*sqrt((structure1.coeffs_error/structure1.coeffs)^2+(structure2.coeffs_error/structure2.coeffs)^2)
  
  ; Update history tag
  if not keyword_set(nohistory) then shapelets_update_history,shapecat,"Added to "+structure2.name

  ; Return answer
  structure1=new_decomp

; Append a shapelet decomposition to a shapecat
endif else if (structure1.type eq "decomp" and structure2.type eq "shapecat") or $
              (structure1.type eq "shapecat" and structure2.type eq "decomp") then begin 

  ; Get into local variables
  if structure1.type eq "decomp" then begin
    decomp=structure1 & shapecat=structure2
  endif else begin
    decomp=structure2 & shapecat=structure1
  endelse

  ; Check that the two structures' truncation orders match
  if decomp.n_max gt shapecat.maxn_max then begin
    shapelets_extend_nmax,shapecat,decomp.n_max-shapecat.maxn_max
  endif else if decomp.n_max lt shapecat.maxn_max then begin
    shapelets_extend_nmax,decomp,shapecat.maxn_max-decomp.n_max
  endif

  ; Check that the structure types (polar or Cartesian) match
  if ((shapecat.polar+decomp.polar) mod 2 ne 0) then begin
    shapelets_polar_convert, decomp, polar=shapecat.polar, cartesian=1-shapecat.polar
  endif

  ; Begin a new catalogue
  new_shapecat={name:shapecat.name,type:"shapecat"} 
  n_tags=n_tags(shapecat)
  names=tag_names(shapecat)
       
  ; Append decomp structure to the new catalogue
  if shapecat.sextractor then message,"WARNING: SExtractor entries have been stripped from the shapecat!",/info,noprint=silent
  for i=0,n_tags-1 do begin
    case strupcase(names[i]) of
      "NAME":
      "TYPE":
      "HISTORY": new_shapecat=create_struct(new_shapecat,"history",shapecat.history)
      "N": new_shapecat=create_struct(new_shapecat,"n",shapecat.n+1)
      "MAXN_MAX": new_shapecat=create_struct(new_shapecat,"maxn_max",shapecat.maxn_max)
      "N_COEFFS": new_shapecat=create_struct(new_shapecat,"n_coeffs",shapecat.n_coeffs)
      "POLAR": new_shapecat=create_struct(new_shapecat,"polar",shapecat.polar)
      "SEEING": new_shapecat=create_struct(new_shapecat,"seeing",shapecat.seeing)
      "X": new_shapecat=create_struct(new_shapecat,"x",[shapecat.x,transpose(decomp.x)])
      "BETA": new_shapecat=create_struct(new_shapecat,"beta",[shapecat.beta,decomp.beta])
      "N_MAX": new_shapecat=create_struct(new_shapecat,"n_max",[shapecat.n_max,decomp.n_max])
      "COEFFS": new_shapecat=create_struct(new_shapecat,"coeffs",[shapecat.coeffs,transpose(decomp.coeffs)])
      "COEFFS_ERROR": new_shapecat=create_struct(new_shapecat,"coeffs_error",[shapecat.coeffs_error,transpose(decomp.coeffs_error)])
      "FLAG": new_shapecat=create_struct(new_shapecat,"flag",[shapecat.flag,intarr(1,2)])
      "CHISQ": new_shapecat=create_struct(new_shapecat,"chisq",[shapecat.chisq,decomp.chisq[1]])
      ; Strip SExtractor entries from catalogue, as they are no longer current
      "SEXTRACTOR": new_shapecat=create_struct(new_shapecat,"sextractor",0B)
      else: if strmid(names[i],0,3) ne "SEX" and strmid(names[i],0,3) ne "OBJ" $
            then new_shapecat=create_struct(new_shapecat,names[i],shapecat.(i))
    endcase
  endfor
  
  ; Update history tag
  if not keyword_set(nohistory) then shapelets_update_history,shapecat,"Object "+decomp.name+" appended"
  
  ; Update moments in the catalogue
  ;if shapecat.moments then shapelets_shapecat_moments, shapecat
  ;if shapecat.shear_estimates then shapelets_shapecat_shear, shapecat

  ; Set variable for output
  structure1=new_shapecat

endif else message,"A "+structure2.type+" structure cannot be added to a "+structure1.type+" structure!"

end
