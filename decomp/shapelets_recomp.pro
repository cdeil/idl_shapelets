pro shapelets_recomp, decomp, recomp, $
                      PSF=psf,        $
                      NRAN=nran,      $
                      TOP=top,        $
                      BOTTOM=bottom,  $
                      COMPLEX=complex,$
                      NOOVER=noover,  $
                      SKY=sky

;$Id: shapelets_recomp.pro, v2$
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
;      SHAPELETS_RECOMP
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Compute the recomposed image corresponding to a set of shapelet
;      coefficients calculated using shapelets_decomp.pro. The input decomp
;      structure also contains a few meta-parameters, including the shapelet
;      scale size beta, and whether or not the basis functions should be
;      integrated within pixels, or merely evaluated at the centre of each
;      pixel.
;
; INPUTS:
;      decomp - structure produced by shapelets_decomp.pro
;
; OPTIONAL INPUTS:
;      psf    - A PSF as a decomp structure type. The model is (re-)convolved
;               with this PSF in shapelet space before it is pixellated. 
;      nran   - Integer array [n_low,n_high]. Ignore coefficients with n 
;               outsdide this range.
;      top    - Integer. Use only this number of (Cartesian) shapelet 
;               coefficients during recomposition, selected as the ones with 
;               the largest absolute values.
;      bottom - Integer. Use only this number of (Cartesian) shapelet 
;               coefficients during recomposition, selected as the ones with 
;               the smallest absolute values.
;
; KEYWORD PARAMETERS:
;      /NOOVER - If set, do not oversample the basis functions.
;      /COMPLEX- If set, return a complex array (complex part is numerical
;                error for a Cartesian shapelet model, but can be populated
;                with polar shapelets).
;      /SKY    - If set, include fit to sky background in reconstructed image.
;
; OUTPUTS:
;      recomp - Reconstructed (2D floating point) image array.
;
; TO DO:
;      Could simplify sky replacement etc. by using shapelets_make_ls_matrix.pro
;
; MODIFICATION HISTORY:
;      Apr 05 - COMPLEX keyword added by RM.
;      Apr 04 - PSF reconvolution added by RM.
;      Nov 01 - Modified by R. Massey to allow numerical integration
;               of basis functions
;      Feb 01 - Modified by AR to allow oversampling of the basis functions
;      Jul 99 - Written by A. Refregier
;-

COMPILE_OPT idl2

; Declarations
if not keyword_set(nran) then nran = [0,decomp.n_max]
shapelets_make_nvec,decomp.n_max,n1,n2,n_coeffs

; Reconvolve with PSF, if necessary
decomp_r=decomp
if keyword_set(psf) then shapelets_convolve,decomp_r,psf

; Eliminate top few coefficients, if requested
if keyword_set(top) then begin
  sorted=reverse(sort(abs(decomp_r.coeffs)))
  coeffs_top=abs(decomp_r.coeffs[sorted[top-1]])
  decomp_r.coeffs[where(abs(decomp_r.coeffs) lt coeffs_top)]=0.
endif

; Eliminate bottom few coefficients, if requested
if keyword_set(bottom) then begin
  sorted=sort(abs(decomp_r.coeffs))
  coeffs_bottom=abs(decomp_r.coeffs[sorted[bottom-1]])
  decomp_r.coeffs[where(abs(decomp_r.coeffs) gt coeffs_bottom)]=0.
endif

; Oversample, if necessary
if keyword_set(noover) then ov=1 else ov=decomp_r.over
ovf=float(ov)

; Create x arrays
shapelets_make_xarr, [decomp_r.n_pixels[0]*ov, decomp_r.n_pixels[1]*ov], x1, x2, x0=decomp_r.x*ovf

; Convert to Cartesian shapelet representation, if necessary
if tag_exist(decomp_r,"polar") then if decomp_r.polar then $
  shapelets_polar_convert,decomp_r,/CARTESIAN

; construct the basis functions
Basis=shapelets_phi([decomp_r.n_max,0], x1/(decomp_r.beta*ovf), x2/(decomp_r.beta*ovf), $
      integrate=decomp_r.integrate, /array) / (decomp_r.beta*ovf)

; Eliminate coeffs outside n_mnax range if requested
if keyword_set(nran) then begin
  n=decomp_r.n1+decomp_r.n2
  outside_n_range=where(n lt nran[0] or n gt nran[1],n_outside_n_range)
  if n_outside_n_range gt 0 then decomp_r.coeffs[outside_n_range]=0.
endif

; Reconstruct the image by summing the (weighted) basis functions
recomp = complexarr(decomp_r.n_pixels*ov)
for i=0,n_coeffs-1 do recomp += decomp_r.coeffs[i]*Basis[*,*,i]

; Resample --- er, does this work?
if ov ne 1 then recomp=rebin(recomp,decomp_r.n_pixels[0],decomp_r.n_pixels[1])

; Add sky background
if keyword_set(sky) then begin                              ; 
  recomp=recomp+decomp_r.skyfit[0]                          ; Add sky background
  if n_elements(decomp_r.skyfit) gt 1 then begin            ; 
    recomp=recomp+decomp_r.skyfit[1]*x1/(decomp_r.beta*ovf) ; Add sky gradient
    recomp=recomp+decomp_r.skyfit[2]*x2/(decomp_r.beta*ovf) ;
  endif                                                     ; 
endif  
  
; Make array complex
if keyword_set(complex) then recomp=complex(recomp) else recomp=float(recomp)

end










