function shapelets_create_decomp, n_max, POLAR=polar

;$Id: shapelets_create_decomp.pro, v1$
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
;      SHAPELETS_CREATE_DECOMP
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Initiliase a band new decomp structure.
;      Works like an IDL __DEFINE procedure, but for an anonymous structure
;      (which we need because they will shrink or expand to accomodate
;      varying amounts of data).
;
; INPUTS:
;      N_MAX  - The desired truncation order of the shapelet expansion.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      POLAR  - Whether or not the expansion should be in polar (set) or
;               Cartesian (unset) shapelet coefficients.
;
; OUTPUTS:
;      DECOMP - Returns the new decomp structure.
;
; MODIFICATION HISTORY:
;      Jul 05 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Parse input
n_max=fix(n_max)
n_coeffs=fix(((n_max+1)*(n_max+2))/2)

; Create shiny new decomp structure
decomp=create_struct("name","","type","decomp")
if keyword_set(polar) then begin
  polar=1B
  coeffs=complexarr(n_coeffs)
  coeffs_error=complexarr(n_coeffs)
  shapelets_make_nvec,n_max,n,m,/POLAR
  shapelets_make_nvec,n_max,nr,nl
  decomp=create_struct(decomp,"polar",polar,"coeffs",coeffs,"coeffs_error",coeffs_error,"n",fix(n),"m",fix(m),"nl",fix(nl),"nr",fix(nr))
endif else begin
  polar=0B
  coeffs=fltarr(n_coeffs)
  coeffs_error=fltarr(n_coeffs)
  shapelets_make_nvec,n_max,n1,n2,n_coeffs
  decomp=create_struct(decomp,"polar",polar,"coeffs",coeffs,"coeffs_error",coeffs_error,"n1",fix(n1),"n2",fix(n2))
endelse

; Add meta-parameters
decomp=create_struct(decomp,"beta",0.,"n_max",n_max,"n_coeffs",n_coeffs,"x",fltarr(2))
; Add pixellisation parameters
decomp=create_struct(decomp,"n_pixels",intarr(2),"over",1.,"integrate",1B,"chisq",fltarr(2),"flag",fix(0),"sky_level",0.,"sky_slope",fltarr(2))
; Add history record
systime=systime()
date=strmid(systime,0,4)+strmid(systime,8,3)+strmid(systime,4,4)+strmid(systime,20,4)
time=strmid(systime,11,5)
decomp=create_struct(decomp,"history","Created at "+time+" on "+date+".")
; Add flags to show whether calculations have been pre-done
;decomp=create_struct(decomp,"sextractor",0B,"moments",0B)

; Return the shiny new decomp structure
return,decomp

end
