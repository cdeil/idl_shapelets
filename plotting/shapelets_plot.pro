pro shapelets_plot, structure, structure2, CARTESIAN_COEFFICIENTS=cartesian_coefficients,   $
                                           POLAR_COEFFICIENTS=polar_coefficients,           $
                                           STATISTICS=statistics,                           $
                                           _ref_extra = ex

;$Id: shapelets_plot.pro, v2$
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
;       SHAPELETS_PLOT
;
; PURPOSE:
;       Generic plotting routine. Acts as a wrapper for many more specific
;       plotting utilities. First decides what the structure type is, rather
;       as you would expect from object-oriented code, then calls the relevant
;       plotting subroutine.
;
; CATEGORY:
;       Shapelets.
;
; INPUTS:
;       STRUCTURE - A shapelet data structure (e.g. decomp, pstamp, image).
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       Most keywords are passed to subroutines, where most graphics keywords
;       are accepted.
;
; OUTPUTS:
;       Draws a plot or plots to the current output device. Depending upon the
;       type of plots, it may wait for user input, then draw some more.
;
; MODIFICATION HISTORY:
;       Apr 05 - Robust error catching implented by RM.
;       Feb 05 - Written by Richard Massey.
;-

COMPILE_OPT idl2
ON_ERROR,2

recognised=shapelets_structure_type(structure,message=message,/silent)
if not recognised then begin
  if size(structure,/n_dimensions) eq 2 then begin
    ; Deal with simple case of an image array.
    shapelets_plot_image, structure, _extra= ex
    return
  endif else message,message
endif

case strlowcase(structure.type) of
  "pstamp":   begin
                shapelets_plot_pstamp, structure, _extra= ex
                success=1B
              end
  "focus":    begin
                if not keyword_set(structure2) then message,"Usage: shapelets_plot,focus,pstamp"
                shapelets_plot_focus, structure, structure2, _extra= ex
                success=1B
              end
  "image":    begin
                if keyword_set(statistics) then begin
                  shapelets_plot_image_statitics, structure;, _extra= ex
                endif else begin
                  shapelets_plot_image, structure.image, _extra= ex
                  !p.region=0
                  ;atv,image.image
                endelse
                success=1B
              end
  "sexcat":   begin
                shapelets_plot_sexcat, structure, structure2, _extra= ex
                success=1B
              end
  "shapecat": begin
                shapelets_plot_shapecat, structure;, _extra= ex
                success=1B
              end
  "chisq_grid": begin
                  shapelets_plot_chisq_grid, structure;, _extra= ex
                  success=1B
                end
  else: success=0B
endcase
switch strlowcase(structure.type) of
  "decomp_polar": message,"Plotting routines not yet properly coded for this structure type."
  "decomp_cartesian": 
  "decomp": begin
              if keyword_set(cartesian_coefficients) then coef=1B 
              if keyword_set(polar_coefficients) then polar=polar_coefficients
              shapelets_plot_decomp,structure, coef=coef, polar=polar, _extra= ex
              success=1B
            end
endswitch

if not success then message,"Structure is not of shapelets type, or of an older version!"

end
