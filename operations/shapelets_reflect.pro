pro shapelets_reflect, structure,          $
                       XAXIS=xaxis,        $
                       YAXIS=yaxis,        $
                       ANGLE=angle,        $
                       POLAR=polar,        $
		       CARTESIAN=cartesian,$
		       POSITION=position,  $
                       ORDER=order,        $
		       EXTEND=extend,      $
                       NOHISTORY=nohistory,$
		       MAINTAIN=maintain

;$Id: shapelets_reflect.pro, v2$
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
;      SHAPELETS_REFLECT
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Makes a mirror image of (an) object(s) by swapping parity.
; 
; INPUTS:
;      STRUCTURE - Shapelet decomp or shapecat structure.
;
; OPTIONAL INPUTS:
;      ANGLE     - Angle of mirror line [degrees, a/c/w from x axis].
;      POSITION  - Position of mirror line relative to centre of basis
;                  functions [x,y], in units of beta.
;      EXTEND    - Increse n_max by this amount (and set new coefs to zero) 
;                  before flipping. Only useful if POSITION is also set.
;
; KEYWORD PARAMETERS:
;      XAXIS     - Perform flip in the x axis (overrides everything).
;                  This is also the default behaviour.
;      YAXIS     - Perform flip in the y axis (overrides ANGLE).
;      POLAR     - Perform rotation in polar shapelet space. 
;                  Default: Cartesian.
;      MAINTAIN  - Restore the (polar/Cartesian) input format on output.
;      NOHISTORY - Do not record this operation in the object's history tag.
;
; OUTPUTS:
;      STRUCTURE - Flipped shapelet decomp or shapecat structure.
; 
; MODIFICATION HISTORY:
;      Jul 05 - History tags updated by RM.
;      Apr 05 - Generalised to accept shapecats as input by RM.
;      Apr 05 - Ability to use arbitrary mirror lines added by RM.
;      Feb 02 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Maintain backwards compatibility
if not shapelets_structure_type(structure,message=message) then message,message

; Increase n_max if requested (only useful if POSITION is also specified).
if keyword_set(extend) then shapelets_extend_nmax,structure,extend

; Decide whether default method should be in Cartesian or polar shapelet space
if not keyword_set(cartesian) then cartesian=0B
if not keyword_set(polar) then polar=0B
if not ((cartesian+polar) mod 2) then begin
  if structure.polar then begin
    cartesian=0B
    polar=1B
  endif else begin
    cartesian=1B
    polar=0B
  endelse
endif
polar_input=structure.polar

; Pre-rotate image so that refection is possible along any line, using a simple
; flip in the x axis. It would screw up the overall positions of objects in a
; shapecat if we did this as a double-rotation after the flip.
new_history="Reflected in a horizontal line"
if not keyword_set(xaxis) then begin
  if keyword_set(yaxis) then begin
    shapelets_rotate,structure,-90,/NOHISTORY
    new_history="Reflected in a vertical line"
  endif else if keyword_set(angle) then begin
    shapelets_rotate,structure,-angle,/NOHISTORY
    new_history="Reflected in a line at "+strmid(strtrim(angle,2),0,5)+" degrees from x axis"
  endif
endif

; Pre-shift image so that refection is possible along any line, using a simple
; flip in the x axis.
if keyword_set(position) then begin
  ; Simply need to rotate the position via a rotation matrix...
  message,"POSITION keyword won't work properly if ANGLE is set too!",/info
  shapelets_translate,structure,-position,order=order,polar=polar,cartesian=cartesian,/NOHISTORY
  new_history=new_history+" passing through ("+strmid(strtrim(position[0],2),0,5)+","+strmid(strtrim(position[1],2),0,5)+")"
  if structure.type eq "shapecat" then message,$
    "Need to calculate transformation of objects' global coordinates!",/info
endif

; Flip object(s)
if polar then begin
  ; Polar shapelet method.
  shapelets_polar_convert,structure,/POLAR,/SILENT
  structure.coeffs=conj(structure.coeffs)
endif else begin
  ; Cartesian shapelet method.
  shapelets_polar_convert,structure,/CARTESIAN,/SILENT
  if structure.type eq "decomp" then begin 
    structure.coeffs=structure.coeffs*(1-2*(structure.n2 mod 2))
  endif else if structure.type eq "shapecat" then begin
    shapelets_make_nvec,structure.maxn_max,n1,n2,n_coeffs
    oddn2=where(n2 mod 2)
    structure.coeffs[*,oddn2]=-structure.coeffs[*,oddn2]
    if tag_exist(structure,"moments") then if structure.moments then $
      message,"The moments in your shapecat need updating!",/info
  endif else if structure.type eq "sexcat" then begin
    message,"TO DO!"
  endif else message,"Cannot reflect objects in "+structure.type+" structures!"
endelse

; Shift image back to original position.
if keyword_set(position) then $
  shapelets_translate,structure,position,order=order,polar=polar,cartesian=cartesian,/NOHISTORY

; Rotate image back to original orientation.
if not keyword_set(xaxis) then begin
  if keyword_set(yaxis) then begin
    shapelets_rotate,structure,90,/NOHISTORY
  endif else if keyword_set(angle) then begin
    shapelets_rotate,structure,angle,/NOHISTORY
  endif
endif

; Convert back to input format if necessary
if keyword_set(maintain) then $
  shapelets_polar_convert,structure,polar=polar_input,cartesian=1-polar_input,/SILENT

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,structure,new_history

end
