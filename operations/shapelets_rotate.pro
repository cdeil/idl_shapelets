pro shapelets_rotate, structure, angle,   $
                      CENTRE=centre,      $
                      POLAR=polar,        $
		                  CARTESIAN=cartesian,$
		                  ORDER=order,        $
                      EXTEND=extend,      $
		                  NOHISTORY=nohistory,$
		                  MAINTAIN=maintain

;$Id: shapelets_rotate.pro, v2$
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
;      SHAPELETS_ROTATE
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Rotates an object by converting it's Cartesian shapelet coefficients to
;      polar shapelet coefficients, then applying the rotation matrix - which,
;      in this basis, is a simple multiplication.
; 
; INPUTS:
;      STRUCTURE - Shapelet decomp or shapecat structure.
;      ANGLE     - Anticlockwise, in degrees.
;
; OPTIONAL INPUTS:
;      CENTRE - Centre of rotation [default: shapelet centre]
;               Apply Translation operator before & after?
;      MATRIX - Cartesian to polar shapelet conversion matrix, if it has
;               already been calculated.
;
; KEYWORD PARAMETERS:
;      POLAR     - Perform translation in polar shapelet space. 
;                  Default: Cartesian.
;      MAINTAIN  - Restore the (polar/Cartesian) input format on output.
;      NOHISTORY - Do not record this operation in the object's history tag.
;
; OPTIONAL OUTPUTS:
;      MATRIX - Cartesian to polar shapelet conversion matrix, so that it 
;               doesn't have to be calculated again in the future, for speed.
;
; OUTPUTS:
;      STRUCTURE - Rotated shapelet model(s).
; 
; MODIFICATION HISTORY:
;      Jun 05 - Routine condensed by RM.
;      Apr 05 - Generalised to accept shapecats as input by RM.
;      Feb 02 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Backwardly compatible
if not shapelets_structure_type(structure,message=message) then message,message

; Parse input 
if keyword_set(extend) then shapelets_extend_nmax,structure,round(extend)
if keyword_set(centre) then message,"Code for CENTRE keyword not yet written!"
case message of
  "decomp":   shapelets_make_nvec, structure.n_max, n, m, /POLAR
  "shapecat": begin
                shapelets_make_nvec, structure.maxn_max, n, m, /POLAR
		m=m##replicate(1,structure.n)
	        message,"Need to calculate transformation of objects' global coordinates!",/info
              end
  else: message,"Cannot rotate objects in "+message+" structures!"
endcase
if keyword_set(order) then if order ne 1 then $
  message,"Arbitrary rotations are trivial in polar shapelet space. Specifying the order is not necessary!",/info

; Remember input format for later use
polar_input=structure.polar

; Convert to polar coefficients
shapelets_polar_convert,structure,/POLAR,/SILENT

; Convert polar coefficients to {mag,phase}
phase	   = atan(imaginary(structure.coeffs),double(structure.coeffs))
modulus  = abs(structure.coeffs)
ephase   = atan(imaginary(structure.coeffs_error),double(structure.coeffs_error))
emodulus = abs(structure.coeffs_error)

; Perform the actual rotation, which simply involves a change in the phase
; of shapelet coefficients, modulated by their spin
phase  += m*angle/!radeg
ephase += m*angle/!radeg

; Convert polar coefficients back to {real,imaginary}
structure.coeffs       = modulus*complex(cos(phase),sin(phase))
structure.coeffs_error = emodulus*complex(cos(ephase),sin(ephase))

; Convert back to Cartesian shapelets if necessary
if keyword_set(maintain) then $
  shapelets_polar_convert,structure,polar=polar_input,cartesian=1-polar_input,/SILENT

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,structure,$
  "Rotated by "+strtrim(angle,2)+" degrees"

end
