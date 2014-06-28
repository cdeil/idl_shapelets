pro shapelets_idl2ascii_shapecat, filename, shapecat_idl, shapecat_ascii

;$Id: shapelets_idl2ascii_shapecat.pro, v2$
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
;       SHAPELETS_IDL2ASCII_SHAPECAT
;
; PURPOSE:
;       Converts a shapecat stored on disc from a shapelets .shapelet file
;       (in IDL SAVE format) to a shapelets .shape file (in ASCII format).
;
; CATEGORY:
;       Shapelet catalogue manipulation.
;
; INPUTS:
;       filename - Filename, without path or extension.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;       A new IDL-format catalogue written to disc.
;       NB: The old ASCII-format catalogue is not overwritten, because they
;       have different file extensions.
;
; OPTIONAL OUTPUTS:
;       shapecat_idl   - Old catalogue structure can be returned.
;       shapecat_ascii - New catalogue structure can be returned.
;
; MODIFICATION HISTORY:
;       Jun 04 - Written by Richard Massey

COMPILE_OPT idl2, OBSOLETE

; Read in IDL-format catalogue
shapelets_read_shapecat, shapecat_idl, filename, /CARTESIAN

; Check catalogue is suitable at a most basic level
if not tag_exist(shapecat_idl,"coeffs") or not tag_exist(shapecat_idl,"x") then $
  message,'Catalogue '+filename+' does not contain even the most basic shapelet parameters!'

; Invent header information
if tag_exist(shapecat_idl,"n") then n=shapecat_idl.n else n=(size(shapecat_idl.coeffs,/dimensions))[0]
if tag_exist(shapecat_idl,"n_pixels") then n_pixels=shapecat_idl.n_pixels else n_pixels=fix([max(shapecat_idl.x[*,0]),max(shapecat_idl.x[*,1])]+[1,1])
if tag_exist(shapecat_idl,"seeing") then seeing=shapecat_idl.seeing else seeing=0.
if tag_exist(shapecat_idl,"pixel_scale") then seeing=shapecat_idl.pixel_scale else pixel_scale=0.
if tag_exist(shapecat_idl,"n_coeffs") then n_coeffs=shapecat_idl.n_coeffs else n_coeffs=(size(transpose(shapecat_idl.coeffs),/dimensions))[0]
if tag_exist(shapecat_idl,"maxn_max") then maxn_max=shapecat_idl.maxn_max else maxn_max=round(-3+sqrt(9-8*(1-n_coeffs)))/2
shapelets_make_nvec,maxn_max,n1,n2,n_coeffs

; Open output file:
openw,lun,shapelets_paths(3)+filename+'.shape',/get_lun
message,'Writing shapelet coefficients to '+filename+'.shape',$
  /info,noprint=silent

; Write out catalogue header (need to change the 999 if n_max>43):
printf,lun,'# Shapelet coefficients of objects in '+filename+'.fits: '
printf,lun,'#   num objs       out of     n_coeffs     n_coeffs       seeing  n_pixels(x)  n_pixels(y)  pixel scale'
printf,lun,format='("#",I11,3(I13),F13.6,2(I13),F13.6)',$
           n,n,maxn_max,n_coeffs,seeing,n_pixels[0],n_pixels[1],pixel_scale
printf,lun,'#'
printf,lun,'#         x0           y0         beta     SEx fwhm      SEx mag        n_max'+$
          '         Flag       SEx ID  Chi squared    SEx class        SEx A        SEx B'+$
	  '    SEx theta       SEx e1       SEx e2     SEx flux        SEx x        SEx y'+$
	  '     SEx area Object:x_min        x_max        y_min        y_max' 
printf,lun,format='($,"# n1",I8," ",(999(I12, :, " ")))',n1 & printf,lun
printf,lun,format='($,"# n2",I8," ",(999(I12, :, " ")))',n2 & printf,lun

; Loop over every object in catalogue, and determine its values
for i=0L,n-1 do begin

  ; Find shapelet coefficients
  coeffs = reform(shapecat_idl.coeffs[i,*])
  ; Renormalise shapelet coefficients to standard way of writing old ASCII catalogues 
  ;if not keyword_set(v1) then 
  coeffs[1:*]=coeffs[1:*]/coeffs[0]

  if tag_exist(shapecat_idl,"coeffs_error") then begin 
    coeffs_error = reform(shapecat_idl.coeffs_error[i,*])
  endif else if tag_exist(shapecat_idl,"error") then begin 
    coeffs_error = reform(shapecat_idl.error[i,*])
  endif else coeffs_error = fltarr(n_coeffs)
 ;if not keyword_set(v1) then coeffs_error[1:*]=coeffs_error[1:*]/coeffs_error[0]
  if tag_exist(shapecat_idl,"flag") then flag=reform(shapecat_idl.flag[i,*]) else flag=[0,0]
  if tag_exist(shapecat_idl,"beta") then beta=shapecat_idl.beta[i] else beta=0.
  if tag_exist(shapecat_idl,"n_max") then n_max=shapecat_idl.n_max[i] else n_max=maxn_max
  if tag_exist(shapecat_idl,"chisq") then chisq=shapecat_idl.chisq[i] else chisq=0.
  if shapecat_idl.sextractor eq 1B then begin
    sexid=shapecat_idl.sexid[i]
    sexfwhm=shapecat_idl.sexfwhm[i]
    sexflux=shapecat_idl.sexflux[i]
    sexmag=shapecat_idl.sexmag[i]
    sexclass=shapecat_idl.sexclass[i]
    sexa=shapecat_idl.sexa[i]
    sexb=shapecat_idl.sexb[i]
    sextheta=shapecat_idl.sextheta[i]
    sexe1=shapecat_idl.sexe1[i]
    sexe2=shapecat_idl.sexe2[i]
    sexx=shapecat_idl.sexx[i]
    sexy=shapecat_idl.sexy[i]
    sexarea=shapecat_idl.sexarea[i]
    im_ran=[shapecat_idl.objx_max[i],shapecat_idl.objx_min[i]shapecat_idl.objy_min[i],,shapecat_idl.objy_max[i]]
  endif else begin
    sexid=0.
    sexfwhm=0.
    sexflux=0.
    sexmag=0.
    sexclass=0.
    sexa=0.
    sexb=0.
    sextheta=0.
    sexe1=0.
    sexe2=0.
    sexx=0.
    sexy=0.
    sexarea=0
    im_ran=[0,0,0,0]
  endelse
  otherdata    = [reform(shapecat_idl.x[i,*]), $
                 beta,                         $
	         sexfwhm,		       $
	         sexmag,		       $
	         n_max, 		       $
	         10*flag[0]+flag[1],	       $
	         sexid, 		       $
	         chisq, 		       $
                 sexclass,                     $
                 sexa,                         $
                 sexb,                         $
                 sextheta,                     $
                 sexe1,                        $
                 sexe2,                        $
                 sexflux,                      $
                 sexx,                         $
                 sexy,                         $
                 sexarea,                      $
                 im_ran[0],                    $
                 im_ran[1],                    $
                 im_ran[2],                    $
                 im_ran[3] ]
  
  ; Write data to catalogue
  printf,lun,format='($,(5(F12.5, :, " ")),(3(I12, :, " ")),(10(F12.4, :, " ")),(5(I12, :, " ")))',otherdata & printf,lun
  printf,lun,format='($,(999(E12.5, :, " ")))',coeffs & printf,lun
  printf,lun,format='($,(999(E12.5, :, " ")))',coeffs_error & printf,lun
endfor

; Write end flag to bottom of shapelet catalogue
printf,lun,'# End'
close,lun
free_lun,lun

; Read new ASCII catalogue into memory if required
if arg_present(shapecat_ascii) then begin
  if not keyword_set(shapecat_ascii) then begin
    shapelets_read_ascii_shapecat,shapecat_ascii,filename
  endif
endif

end



