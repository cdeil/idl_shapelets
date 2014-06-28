function shapelets_shapecat2decomp, shapecat, id,        $
                                    OVER=over,           $
                                    INTEGRATE=integrate, $
                                    V1=v1,               $
                                    N_PIXELS=n_pixels

;$Id: shapelets_shapecat2decomp.pro, v2$
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
;      SHAPELETS_SHAPECAT2DECOMP
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Extracts one (decomp structure) object from a (shapecat structure)
;      catalogue.
;
; INPUTS:
;      SHAPECAT - Shapelet catalogue structure. Can contain thousands of
;                 objects.
;      ID       - If the catalogue contains more than one object, specify
;                 which one of them is to be converted.
;
; OPTIONAL INPUTS:
;      OVER     - Oversampling amount in decomp structure. Default: 1
;
; KEYWORD PARAMETERS:
;      INTEGRATE- Integrate basis functions within pixels. Default: 1
;      V1       - Assume shapelet coefficients are stored as a_ij/a_00, as
;                 in version 1 of shapelet code. Default is now just a_ij.
;
; OUTPUTS:
;      DECOMP   - Cartesian decomp structure containing all of the information
;                 that it is possible to reasonably extract from a catalogue.
;
; MODIFICATION HISTORY:
;      Jul 05 - External decomp creation routine incorporated by RM.
;      Jun 05 - Polar flag added to decomp structures by RM.
;      Apr 05 - Bug fixed in determination of n_pixels by JB.
;      Jan 05 - N_pixels/centroid consistent with shapecats implemented by RM.
;      Oct 04 - N_PIXELS keyword added by Joel Berge.
;      Aug 04 - Inter-pixel centroid position conserved by RM.
;      Jun 04 - Integers changed to long integers by RM.
;      Apr 04 - Convention on storage of aij/a00 to aij changed by RM.
;      Mar 04 - Tidied up by RM.
;      Apr 02 - Can now take polar catalogues/objects as input.
;      Dec 01 - Written by Richard Massey.
;-

; Initiliase variables
COMPILE_OPT idl2
on_error,2

; Maintain backwards compatibility
if not shapelets_structure_type(shapecat,message=message) then message,message
if not tag_exist(shapecat,"sextractor") then shapecat=create_struct(shapecat,"sextractor",0B)

; Parse inputs
if not keyword_set(id) then id=0
if id ge shapecat.n or id lt 0 then $
  message,"Object #"+strtrim(string(id),2)+" does not exist in this catalogue of only "+strtrim(string(shapecat.n),2)+" objects!"
n_max=shapecat.n_max[id]<shapecat.maxn_max

; Make a new decomp structure
decomp=shapelets_create_decomp(n_max,polar=shapecat.polar)

; Extract shapelet parameters from catalogue
decomp.name=shapecat.name+"_"+strtrim(string(id),2)
decomp.beta=shapecat.beta[id]
decomp.coeffs=reform(shapecat.coeffs[id,0:decomp.n_coeffs-1])
decomp.coeffs_error=reform(shapecat.coeffs_error[id,0:decomp.n_coeffs-1])
if n_elements(over) eq 0 then decomp.over=1
if n_elements(integrate) eq 0 then decomp.integrate=1B

; Rescale coefficients from relative to a_00 overall brightness
; (STORAGE CONVENTION CHANGED APRIL 2004)
if keyword_set(v1) and decomp.n_coeffs gt 1 and not shapecat.polar then coeffs[1:*]=coeffs[1:*]*coeffs[0]

; Estimate a decent number of pixels to be used if the decomp ever gets
; (re)pixellated
if keyword_set(n_pixels) then begin
  if n_elements(n_pixels) ge 2 then decomp.n_pixels=n_pixels[0:1] else decomp.n_pixels=replicate(n_pixels[0],2)
  decomp.x=fix(decomp.n_pixels)/2+(shapecat.x[id,*] mod 1)
endif else begin
  if shapecat.sextractor then begin
    decomp.n_pixels=[shapecat.objx_max[id]-shapecat.objx_min[id]+1,shapecat.objy_max[id]-shapecat.objy_min[id]+1]
    decomp.x=reform(shapecat.x[id,*])-[shapecat.objx_min[id],shapecat.objy_min[id]]
  endif else begin  ; Not specified, so guess something arbitraily
    n_pix_min=30
    n_pix_max=100
    th_max=decomp.beta*sqrt(n_max+0.5)
    n_pixels=n_pix_min>round([1.,1.]*1.5*th_max^1.1)<n_pix_max
    n_pixels=n_pixels + (n_pixels mod 2)         ; Make them even numbers
    decomp.x=n_pixels/2+(shapecat.x[id,*] mod 1) ; Put centre of basis functions in central pixel
    decomp.n_pixels=round(2*decomp.x)            ; Make sure that the central pixel really is in the centre
  endelse
endelse

; Recover sextractor information from shapelet catalogue if available
if shapecat.sextractor then begin
  sex={seeing:shapecat.seeing,  			      $
       x:[shapecat.sexx[id],shapecat.sexy[id]], 	      $
       xmin:[shapecat.objx_min[id],shapecat.objy_min[id]],    $
       xmax:[shapecat.objx_max[id],shapecat.objy_max[id]],    $
       id:shapecat.sexid[id],				      $
       class:shapecat.sexclass[id],			      $
       a:shapecat.sexa[id],				      $
       b:shapecat.sexb[id],				      $
       theta:shapecat.sextheta[id],			      $
       e:sqrt(shapecat.sexe1[id]^2+shapecat.sexe2[id]^2),     $
       e1:shapecat.sexe1[id],				      $
       e2:shapecat.sexe2[id],				      $
       flux:shapecat.sexflux[id],			      $
       mag:shapecat.sexmag[id], 			      $
       fwhm:shapecat.sexfwhm[id],			      $
       area:shapecat.sexarea[id] }
  decomp=create_struct(decomp,"sex",sex)
endif

; Initialise object's history tag
shapelets_update_history,decomp,$
  "Object #"+strtrim(string(id),2)+' from catalogue "'+shapecat.name+'"'

return,decomp

end
