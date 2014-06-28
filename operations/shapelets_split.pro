pro shapelets_split, catalogue, selected, rest=rest, NOHISTORY=nohistory

;$Id: shapelets_split.pro, v2$
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
;      SHAPELETS_SPLIT
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Splits a (IDL shapecat or sexcat structure) catalogue.
;
; INPUTS:
;      CATALOGUE  - Shapelet shapecat or sexcat structure. 
;      SELECTED   - Array of numbers, containing the catalogue IDs of those
;                   objects that should be kept in the output catalogue. 
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      CATALOGUE - Shapelet catalogue structure, containing only the selected 
;                  objects from the input catalogue.
;
; OPTIONAL OUTPUTS:
;      REST      - Shapelet catalogue structure, containing the rest of the
;                  objects from the input catalogue.
;
; MODIFICATION HISTORY:
;      Jan 06 - Errors when only one unique object is selected caught by JB.
;      Jul 05 - Separate splitters combined into one general routine by RM .
;      Mar 05 - Shapecat splitter written by RM.
;      Feb 05 - Sexcat splitter adapted by Joel Berge.
;      May 03 - Sexcat splitter written by Richard Massey.
;-

ON_ERROR,2
COMPILE_OPT idl2

; Error catching
if not shapelets_structure_type(catalogue,message=message) then message,message
if not keyword_set(selected) then begin
  if tag_exist(catalogue,"good") then begin
    selected=where(catalogue.good)
  endif else message,"Objects to keep not selected!"
endif else if min(selected) lt 0 or max(selected) gt catalogue.n-1 then begin
  message,"Objects specified out of possible range"
endif

; Remove duplicate entries (these can be forced using shapelets_add later)
selected_sort=selected[sort(selected)]
selected_uniq=selected_sort[uniq(selected_sort)]

; Deal with selections of only a single object
if size(selected_uniq,/n_elements) eq 1 then begin
    selected_uniq_tmp=lonarr(1)
    selected_uniq_tmp[0]=selected_uniq
    selected_uniq=selected_uniq_tmp
endif

; Begin a new (split) catalogue
catalogue_out={name:"Part of "+catalogue.name, type:catalogue.type} 

; Compile a concatenated list of variables
names=tag_names(catalogue)
for i=0,n_tags(catalogue)-1 do begin
  case strupcase(names[i]) of
    "NAME":
    "TYPE":
    "N": catalogue_out=create_struct(catalogue_out,"n",n_elements(selected_uniq))
    else: begin
            case size(catalogue.(i),/N_DIMENSIONS) of
              0: catalogue_out=create_struct(catalogue_out,names[i],catalogue.(i))
              1: catalogue_out=create_struct(catalogue_out,names[i],(catalogue.(i))[selected_uniq])
              2: catalogue_out=create_struct(catalogue_out,names[i],(catalogue.(i))[selected_uniq,*])
              3: catalogue_out=create_struct(catalogue_out,names[i],(catalogue.(i))[selected_uniq,*,*])
              4: catalogue_out=create_struct(catalogue_out,names[i],(catalogue.(i))[selected_uniq,*,*,*])
              else: message,"Very high number of dimensions for one entry. Something's gone wrog!"
            endcase
          end
  endcase
endfor

; Add remark to object's history
if not keyword_set(nohistory) and catalogue.type eq "shapecat" then shapelets_update_history,catalogue_out,$
  "Catalogue split, keeping "+strtrim(catalogue_out.n,2)+" of "+strtrim(catalogue.n,2)+" objects"

; Sift out the remaining objects
if arg_present(rest) then begin
  not_selected=bytarr(catalogue.n)
  for i=0,catalogue.n-1 do begin
    match=where(selected_uniq eq i,n_match)
    not_selected[i]=1B-byte(n_match<1)
  endfor
  shapelets_split,catalogue,where(not_selected)
  rest=catalogue
endif

; Return new catalogue to external scope
catalogue=catalogue_out


end
