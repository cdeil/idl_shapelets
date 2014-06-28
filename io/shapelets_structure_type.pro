function shapelets_structure_type, structure,       $
                                   MESSAGE=message, $
                                   SILENT=silent,   $
                                   VERBOSE=verbose

;$Id: shapelets_structure_type.pro, v2$
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
;       SHAPELETS_STRUCTURE_TYPE
;
; PURPOSE:
;       Determine whether a variable is a known shapelets structure, in a
;       very robust manner.
;
; CATEGORY:
;       Shapelets.
;
; INPUTS:
;       STRUCTURE - A shapelet data structure (e.g. decomp, pstamp, image).
;
; OPTIONAL INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;       SILENT    - Operates silently.
;       VERBOSE   - Operates noisily.
;
; OUTPUTS:
;       0 if the input is not a recognised shapelet structure, 1 if it is.
;
; OPTIONAL OUTPUTS:
;       MESSAGE   - A message containing a description of why the structure
;                   type was recognised. 
;
; MODIFICATION HISTORY:
;       Apr 05 - Written by Richard Massey
;-

COMPILE_OPT idl2
ON_ERROR,2

if not keyword_set(verbose) then verbose=0B

; Test anything that might go wrong
if not keyword_set(structure) then begin
  recognised=0B
  message="Variable is not defined!"
endif else begin
  if size(structure,/TYPE) ne 8 then begin
    recognised=0B
    message="Variable is not an IDL structure!"
  endif else begin
    if not tag_exist(structure,"type") then begin
      recognised=0B
      message="Variable is not of a known shapelets structure type!"
    endif else begin
      if size(structure.type,/TYPE) ne 7 then begin
        recognised=0B
        message="Variable is not of a known shapelets structure type!"
      endif else begin
        case strlowcase(structure.type) of
          "image":            begin & recognised=1B & message="image" & end
          "pstamp":           begin & recognised=1B & message="pstamp" & end
          "sexcat":           begin & recognised=1B & message="sexcat" & end
          "decomp":           begin & recognised=1B & message="decomp" & if not tag_exist(structure,"polar") then structure=create_struct(structure,"polar",0B) & end
          "decomp_cartesian": begin & recognised=1B & structure.type="decomp" & message="decomp" & if not tag_exist(structure,"polar") then structure=create_struct(structure,"polar",0B) & end
          "decomp_polar":     begin & recognised=0B & message="Structure type reserved!" & end
          "shapecat":         begin & recognised=1B & message="shapecat" & end
          "focus":            begin & recognised=1B & message="focus" & end
          "chisq_grid":       begin & recognised=1B & message="chisq_grid" & end
          else: begin
                  recognised=0B
                  message="Structure type not recognised!"
                end
        endcase
      endelse
    endelse
  endelse
endelse

; Print message to screen
if verbose or not recognised then message,/info,noprint=silent,$
  (not recognised)?message:"Structure recognised as a "+message

; Report findings
return,recognised

end
