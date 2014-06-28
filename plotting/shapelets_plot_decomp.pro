pro shapelets_plot_decomp, decomp, recomp,            $
                           POLAR=polar,               $
                           CARTESIAN=cartesian,       $
			                     COEFFICIENTS=coefficients, $
			                     ERRORS=errors,             $
			                     NRAN=nran,                 $
			                     TOP=top,                   $
			                     NOOVER=noover,             $
			                     CRANGE=crange,             $
			                     CLOG=clog,                 $
			                     CBAR=cbar,                 $
                           ISOTROPIC=isotropic,       $
                           REAL=real,                 $
                           IMAGINARY=imaginary,       $
                           COMPOSITE_RI=composite_ri, $
                           MODULUS=modulus,           $
                           ARGUMENT=argument,         $
                           COMPOSITE_MA=composite_ma, $
                           PHASE=phase,               $
                           TITLE=title,               $
			                     XTITLE=xtitle,             $
			                     YTITLE=ytitle,             $
			                     FRAME=frame,               $
			                     CROSSHAIRS=crosshairs

;$Id: shapelets_plot_decomp.pro, v2$
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
;      SHAPELETS_PLOT_DECOMP
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Repixellate a shapelet model (decomp structure) and display it.
;
; INPUTS:
;      DECOMP - structure produced by shapelets_decomp
;
; OPTIONAL INPUTS:
;      NRAN  - Ignore coefficients with n outsdide of nran
;      TOP   - Keep top largest coefficients
;      Plus many more of the usual plotting parameters.
;      See shapelets_plot_image.pro for details.
;
; KEYWORD PARAMETERS:
;      (DEFAULT) - Plot the reconstructed image.
;      CART      - Plot the Cartesian shapelet coefficients.
;      POLAR     - Plot the polar shapelet coefficients (these are complex -
;                  other keywords specify which aspect of them to plot.
;      MODULUS   - (DEFAULT) Plot the modulus of the polar shapelet coeffs.
;      ARGUMENT  - Plot the phases of the polar shapelet coefficients.
;      REAL      - Plot the real parts of the polar shapelet coefficients.
;      IMAGINARY - Plot the imaginary parts of the polar shapelet coefficients.
;      COEF      - A synonym for CARTESIAN.
;      PHASE     - A synonym for ARGUMENT.
;      ERROR     - 1: Plot S/N of object or coefficients.
;                  2: Plots absolute errors.
;      CROSS     - Overlays crosshairs on the centre of the basis functions.
;
; OUTPUTS:
;      Plot drawn to STDOUT.
;
; OPTIONAL OUTPUTS:
;      RECOMP - reconstructed image array
;
; MODIFICATION HISTORY:
;      Jan 06 - Keywords for different plots of complex coefficients added by RM.
;      Jul 05 - Ability to cope with polar decomp structures incorporated by RM.
;      Mar 02 - Polar shapelet options (using shapelets_plot_decomp_polar.pro) added by RM.
;      Dec 01 - Power spectrum, errors etc. incorporated by Richard Massey.
;      Jul 99 - Written by Alexandre Refregier
;-

COMPILE_OPT idl2, HIDDEN

; Test that the input is a decomp structure
if not shapelets_structure_type(decomp,message=message) then message,message
if decomp.type ne "decomp" then message,"This routine is for plotting decomp structures only!"
polar_input=decomp.polar

;
; Plot errors on coefficients if requested
;
if keyword_set(errors) then begin

  if not keyword_set(title) then title=""
  temp_decomp=decomp
  if errors eq 1 then begin
    temp_decomp.coeffs=abs(temp_decomp.coeffs)/temp_decomp.coeffs_error
    title="S/N on coefficients "+title
  endif else begin
    temp_decomp.coeffs=temp_decomp.coeffs_error
    title="Errors of coefficients "+title
  endelse
  shapelets_plot_decomp, temp_decomp, recomp, nran=nran, coef=coef, top=top,$
    noover=noover,cbar=cbar,clog=clog,crange=crange, title=title,$
     polar=polar, frame=frame

endif else begin

  ;
  ; Set range of coefficients to include
  ;
  if not keyword_set(nran) then nran=[0,decomp.n_max]
  nran=0>nran[sort(nran)]<decomp.n_max

  ;
  ; Make plot
  ;
  if keyword_set(cartesian) or keyword_set(coefficients) then begin
    ;
    ; Plot Cartesian shapelet coefficient array
    ;
    ; Convert to Cartesian shapelet coefficients
    decomp_temp=decomp
    shapelets_polar_convert,decomp_temp,/CARTESIAN,/SILENT

    ; Top and tail coefficients to display only a certain n-range, if required
    if keyword_set(top) then begin
      sorted=reverse(sort(abs(decomp_temp.coeffs)))
      a_top=abs(decomp_temp.coeffs[sorted[top-1]])
    endif else begin
      a_top=0.
    endelse

    ; Construct a suitable image
    image=fltarr(nran[1]+1,nran[1]+1)
    for i=0,decomp_temp.n_coeffs-1 do begin
      ni=decomp_temp.n1[i]+decomp_temp.n2[i]
      if ni ge nran[0] and ni le nran[1] and abs(decomp_temp.coeffs[i]) ge a_top then $
         image[decomp_temp.n1[i],decomp_temp.n2[i]]=decomp_temp.coeffs[i]
    endfor
    xtitle='!6n!i1!n'
    ytitle='!6n!i2!n'
    frame=[-.5,nran[1]+.5,-.5,nran[1]+.5]  
    
  endif else if keyword_set(polar) then begin
    ;
    ; Plot polar shapelet coefficient array
    ;
    ; Convert to polar shapelet coefficients
    decomp_temp=decomp
    shapelets_polar_convert,decomp_temp,/POLAR,/SILENT

    ; Rotate so major axis is horizontal
    if keyword_set(rotate) then begin
      e=shapelets_ellipticity(decomp_temp)
      angle=atan(imaginary(e),float(e))/2.
      shapelets_rotate,decomp_temp,-angle
    endif

    ; Top and tail coefficients to display only a certain n-range, if required
    shapelets_extend_nmax,decomp_temp,nran[1]-decomp.n_max,/SILENT
    low_n=where(decomp_temp.n lt nran[0],n_low_n)
    if n_low_n gt 0 then decomp_temp.coeffs[low_n]=complex(0.,0.)
       
    ; Construct a suitable image
    image=dcomplexarr(decomp_temp.n_max+1,2*decomp_temp.n_max+2)
    shapelets_make_nvec,decomp_temp.n_max,n,m,/POLAR
    image[n,m+decomp_temp.n_max]=decomp_temp.coeffs
    image[n,m+decomp_temp.n_max+1]=decomp_temp.coeffs
    frame=[-.5,decomp_temp.n_max+.5,-decomp_temp.n_max-1,decomp_temp.n_max+1]
    xtitle='!6n'
    ytitle='!6m'

    ; There are various options that we could plot
    if keyword_set(real) then begin ; Plot the real part of complex coefficients
      image=float(image)
    endif else if keyword_set(imaginary) then begin ; Plot the imaginary part of complex coefficients
      image=imaginary(image)
    endif else if keyword_set(phase) or keyword_set(argument) then begin ; Plot the argument of complex coefficients
      image=atan(imaginary(image),float(image))
      crange=[-1,1]*!pi
    endif else if keyword_set(split) then begin ; Plot a mixture of the two
      image[*,decomp_temp.n_max+1:*]=abs(image[*,0:decomp_temp.n_max])
      image[*,0:decomp_temp.n_max]=atan(imaginary(image[*,0:decomp_temp.n_max]),float(recomp[*,0:decomp_temp.n_max]))
      image=float(image)
      if not keyword_set(crange) then crange=[-1,1]*!pi
    endif else begin ; Plot the modulus of complex coefficients
      image=abs(image)
    endelse
    if not keyword_set(crange) then crange=[min(image,max=max_image),max_image]
    
  endif else begin
    ;
    ; Plot reconstructed image
    ;
    shapelets_recomp,decomp,recomp,nran=nran,top=top,noover=noover
    if keyword_set(crosshairs) and not keyword_set(frame) then frame=1B
    image=recomp
    
  endelse

  ; Produce actual plot
  shapelets_plot_image,image,cbar=cbar,csize=csize,clog=clog,crange=crange,$
    frame=frame,title=title,xtitle=xtitle,ytitle=ytitle,isotropic=isotropic
  
  ; Overlay crosshairs
  if keyword_set(crosshairs) then begin
    oplot,[0,decomp.n_pixels[0]],[decomp.x[1],decomp.x[1]],psym=-3
    oplot,[decomp.x[0],decomp.x[0]],[0,decomp.n_pixels[1]],psym=-3
  endif
  
endelse

end

