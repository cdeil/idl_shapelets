pro shapelets_ascii2idl_shapecat, filename, shapecat_ascii, shapecat_idl

;$Id: shapelets_ascii2idl_shapecat.pro, v2$
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
; ************************************************************************
; ************************************************************************
;
; NAME:
;       SHAPELETS_ASCII2IDL_SHAPECAT
;
; PURPOSE:
;       Converts a shapecat stored on disc from a shapelets .shape file
;       (in ASCII format) to a shapelets .shapelet file (in IDL SAVE format).
;
; CATEGORY:
;       Shapelet catalogue manipulation.
;
; INPUTS:
;       filename - Filename, without path or extension.
;
; OUTPUTS:
;       A new IDL-format catalogue written to disc.
;       NB: The old ASCII-format catalogue is not overwritten, because they
;       have different file extensions.
;
; OPTIONAL OUTPUTS:
;       shapecat_ascii - Old catalogue structure can be returned.
;       shapecat_idl   - New catalogue structure can be returned.
;
; MODIFICATION HISTORY:
;       May 04 - Written by Richard Massey

COMPILE_OPT idl2, OBSOLETE

; Read in ASCII-format catalogue
shapelets_read_ascii_shapecat, shapecat_ascii, filename;, /moments

;Create new catalogue structure in memory
shapecat_idl={name:shapecat_ascii.name,type:shapecat_ascii.type}

; Insert additional tags that are required by shapelets code versions 2.0 and above
names=tag_names(shapecat_ascii)
for i=0,n_tags(shapecat_ascii)-1 do begin
  case names[i] of
    "NAME":
    "TYPE":
    "ERROR":
    "FLUX":  shapecat_idl=create_struct(shapecat_idl,"MOMENTS",1B,names[i],shapecat_ascii.(i))
    "SEXID": shapecat_idl=create_struct(shapecat_idl,"SEXTRACTOR",1B,names[i],shapecat_ascii.(i))
    else:    shapecat_idl=create_struct(shapecat_idl,names[i],shapecat_ascii.(i))
  endcase
endfor

; Calculate object moments
resolve_routine, "shapelets_read_shapecat"
shapelets_shapecat_moments, shapecat_idl

; Write out new IDL-format catalogue to disk
shapelets_write_shapecat, shapecat_idl, filename

end



