pro shapelets_circularise, structure,           $
                           NOHISTORY=nohistory, $
                           MAINTAIN=maintain

;$Id: shapelets_circularise.pro, v2$
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
;+
; NAME:
;      SHAPELETS_CIRCULARISE
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Circularises objects by setting to zero all of their polar shapelet
;      coefficients where m is nonzero. When averaged around all angles,
;      the negative parts of these basis states cancel out the positive
;      parts.
;
; INPUTS:
;      STRUCTURE - A shapelet decomp or shapecat structure.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      MAINTAIN  - Restore the (polar/Cartesian) input format on output.
;      NOHISTORY - Do not record this operation in the object's history tag.
;
; OUTPUTS:
;      STRUCTURE - Circulaised version of the input object(s).
;
; MODIFICATION HISTORY:
;      Jul 05 - Rendered capable of accepting shapecat structures by RM.
;      Feb 02 - Written by Richard Massey.
;-

; Backwardly compatible
if not shapelets_structure_type(structure,message=message) then message,message

; Convert to polar shapelet coefficients
if structure.polar then cartesian_input=0B else cartesian_input=1B
shapelets_polar_convert,structure,/POLAR,/SILENT

; Set coefficient values away from the diagonal to zero
if structure.type eq "decomp" then begin
  mnonzero=where(structure.m ne 0,n_mnonzero)
  if n_mnonzero gt 0 then structure.coeffs[mnonzero]=complex(0.,0.)
endif else if structure.type eq "shapecat" then begin
  shapelets_make_nvec,structure.maxn_max,n,m,/POLAR
  mnonzero=where(m ne 0,n_mnonzero)
  if n_mnonzero gt 0 then structure.coeffs[*,mnonzero]=complex(0.,0.)
endif else if structure.type eq "sexcat" then begin
  message,"Now, why would you want to do that?"
endif else message,"Cannot circularise objects in "+structure.type+" structures!"

; Convert back to Cartesian shapelets if necessary
if keyword_set(maintain) and cartesian_input then $
  shapelets_polar_convert,structure,/CARTESIAN,/SILENT

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,structure,"Circularised"

end
