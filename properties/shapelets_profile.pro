function shapelets_profile, structure, r, centre=centre

;$Id: shapelets_profile.pro, v2$
;
; Copyright © 2005 Richard Massey and Alexandre Refregier.
;
; This file is a part of the Shapelets analysis code.
; www.ast.cam.ac.uk/~rjm/shapelets/
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
;       SHAPELETS_PROFILE
;
; PURPOSE:
;       Compute the azymuthally averaged profile f(r) of an object from 
;       its shapelet coefficients.
;
; INPUTS:
;       STRUCTURE - Cartesian shapelet decomp structure
;       r         - radii at which to evaluate the profile [pixels]
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       f(r).
;
; MODIFICATION HISTORY:
;       Jul 05 - Completely rewritten and combined by Richard Massey.
;       Mar 02 - Decomp parser written by AR.
;       Jul 99 - Image parser written by Alexandre Refregier.
;-

COMPILE_OPT idl2

if not shapelets_structure_type(structure,message=message,/silent) then begin
  if size(structure,/n_dimensions) eq 2 then begin

    message,"TO DO!"

    ; Deal with simple case of an image array.
    n_pixels=size(structure,/DIMENSIONS)
    n1=si[1] & n2=si[2]
    if not keyword_set(centre) then centre=[n1/2.,n2/2.] 	    
    n_th=1		 ; number of theta bins
    
    ;; define angular bins - linear spacing
    ;th_ran = [0., sqrt(2.)*float(min([n1,n2])/2.)]
    ;dth = (th_ran[1]-th_ran[0])/float(n_th)
    ;th_l = findgen(n_th)*dth+th_ran[0]
    ;th_m = th_l+.5*dth
    ;th_h = th_l+dth

    ;; define angular bins - log spacing
    ;th_ran_log = [1., sqrt(2.)*float(min([n1,n2])/2.)]
    ;dlogth =  (alog10(th_ran_log[1])-alog10(th_ran_log[0]))/float(n_th)
    ;logth_l = findgen(n_th)*dlogth+alog10(th_ran_log[0])
    ;logth_m = logth_l+.5*dlogth 
    
    ; define arrays
    ith = fltarr(n_th)
    ith_log = fltarr(n_th)
    nth = lonarr(n_th)
    nth_log = lonarr(n_th)
    thbin=fltarr(n1,n2)
    logthbin=fltarr(n1,n2)
    ; compute profile
    for i=0, n1-1 do begin
       for j=0, n2-1 do begin
  	  th = sqrt(float(i-x0[0])^2+float(j-x0[1])^2)
  	  thi = fix((th-th_ran[0])/dth)
  	  thbin[i,j]=thi
  	  thi_log = fix((alog10(th)-alog10(th_ran_log[0]))/dlogth)
  	  logthbin[i,j]=thi_log      
  	  if thi lt n_th-1 and thi ge 0 then begin
  	    ith[thi] = ith[thi]+im[i,j]
  	    nth[thi] = nth[thi]+1l
  	  endif   
  	  if thi_log lt n_th-1 and thi_log ge 0 then begin
  	    ith_log[thi_log] = ith_log[thi_log]+im[i,j]
  	    nth_log[thi_log] = nth_log[thi_log]+1L
  	  endif 
  	 endfor
    endfor
    gd=where(nth gt 0)
    ith = ith[gd]/float(nth[gd])
    gd_log=where(nth_log gt 0)
    ith_log = ith_log[gd_log]/float(nth_log[gd_log])
    
    
    ; store in structure
    profile={r:th_m[gd],f:ith[gd],r_log:10.^(logth_m[gd_log]),$
      f_log:ith_log[gd_log]}
    return,profile


  endif else message,message
endif else if structure.type eq "decomp" then begin

  ; Stash away the input variable
  decomp=structure
  
  ; Calculate around a different centre
  if keyword_set(centre) then shapelets_translate,centre-shapelets_centre(decomp),order=4

  ; Convert shapelet coefficients to polar format if necessary
  if tag_exist(decomp,"polar") then if decomp.polar then polar_coeffs=decomp.coeffs
  if n_elements(polar_coeffs) eq 0 then polar_coeffs=shapelets_polar_matrix(decomp.n_max,/c2p)#decomp.coeffs

  ; Extract the m=0 coefficients
  shapelets_make_nvec,decomp.n_max,n,m,/POLAR
  polar_coeffs=float(polar_coeffs[where(m eq 0)])
  
  ; Compute the 1D polar shapelet radial functions
  BasisF=shapelets_chi(decomp.n_max-decomp.n_max mod 2,r,beta=decomp.beta,/ARRAY)
  
  ; Reconstruct the radial profile
  radial_profile=polar_coeffs##BasisF

endif else if structure.type eq "shapecat" then begin

  message,"Calculate mean radial profile"
  ; Convert to polar shapelet coeffs
  ; Dilate small objects?
  ; Average coefficients (weighted by 1/flux)
  ; feed that into this routine as a decomp

endif else message,"Can't calculate the radial profile of "+structure.type+" structures!"

; Tell the world
return,radial_profile

end
