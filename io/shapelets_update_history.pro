pro shapelets_update_history, structure, new_history

;$Id: shapelets_update_history.pro, v1$
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
;      SHAPELETS_UPDATE_HISTORY
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Append a string to an object's history record
;
; INPUTS:
;      STRUCTURE   - Shapelets shapecat or decomp structure.
;      NEW_HISTORY - String, describing something that has just been done.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      STRUCTURE is returned with an updated history tag.
;
; TO DO:
;      Make the history tag an array of strings.
;
; MODIFICATION HISTORY:
;      Jul 05 - Written by Richard Massey.
;-

; Parse input.
if not keyword_set(new_history) then new_history=""

; Append a full stop. There's no excuse for bad punctuation!
if strlen(new_history) gt 0 and strmid(new_history,0,/REVERSE_OFFSET) ne "." then new_history=new_history+"."

if not tag_exist(structure,"history") then begin
  ; Create a new history tag.
  structure=create_struct(structure,"history",new_history)
endif else begin
  ; Include this string in the existing history tag.
  new_structure={name:structure.name, type:structure.type} 
  names=tag_names(structure)
  for i=0,n_tags(structure)-1 do begin
    case strupcase(names[i]) of
      "NAME":
      "TYPE":
      "HISTORY": new_structure=create_struct(new_structure,"history",[structure.history,new_history])
      else: new_structure=create_struct(new_structure,names[i],structure.(i))
    endcase
  endfor
  structure=new_structure
  ;structure.history=[structure.history+" "+new_history]
endelse

end
