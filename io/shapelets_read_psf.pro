pro shapelets_read_psf,DECOMP,FWHM,                       $
                       N_MAX=n_max,                       $
                       RECOMP=recomp,                     $
                       PIXSIZE=pixsize,                   $
                       SNAP=snap,                         $
                       HST=hst,                           $
                       GEMS=gems,                         $
                       SUBARU=subaru, SINDEX=sindex,      $
                       DIFFSNAP=diffsnap,                 $
                       CHARGE_DIFFUSION=charge_diffusion, $
                       SILENT=silent

;$Id: shapelets_read_psf.pro, v2$
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
;         SHAPELETS_READ_PSF
;
; PURPOSE:
;         Creates a shapelet decomp structure representing a PSF.
;
; CATEGORY:
;         Shapelets.
;
; INPUTS:
;         None are absolutely necessary. Select one of the options below.
;
; OPTIONAL INPUTS:
;         FWHM     - FWHM of the Gaussian.
;         N_MAX    - Truncation parameter of the returned decomp structure.
;         CHARGE_DIFFUSION - Charge diffusion length [microns] for SNAP CCDs.
;         PIXSIZE  - Pixel scale [arcseconds].
;
; KEYWORD PARAMETERS:
;         SILENT   - Operates silently.
;         HST      - Generate a (TinyTim) model of the HST WFPC2 PSF.
;         SNAP     - Generate a model of the predicted SNAP PSF.
;         DIFFSNAP - Generate a model of the SNAP PSF, deconvolved from the
;                    HST PSF. Convolution with this kernel will change a
;                    HST image into a SNAP image.
;
; OUTPUTS:
;         DECOMP   - A decomp structure containing the PSF is returned.
;
; OPTIONAL OUTPUTS:
;         RECOMP   - A pixellated image of the PSF model.
;
; EXAMPLE USE:
;         With no keywords, returns a unit-normalised Gaussian (00 state) with width
;         "fwhm" pixels. eg "shapelets_read_psf,psf,4.0"
;         With either /HST, /PC, /ACS, /SNAP or /GAIA, will read in a raytraced PSF 
;           from that telescope,decompose it into shapelets and then return that. 
;           NB: /HST is the WFPC L and /PC is the (higher resolution) Planetary Camera.
;         A pixellised image of the PSF model is optionally returned in recomp.
;         
; TO DO:
;         Can I model the GAIA PSF?
;
; MODIFICATION HISTORY:
;         Mar 05 - SUBARU PSFs relocated to centres of postage stamps by RM.
;         Jan 05 - SUBARU PSF added by Molly Peeples.
;         Aug 04 - GEMS PSF from Catherine Heymans added by RM.
;         Jul 04 - Variable charge diifusion length for SNAP added by WH.
;         Jul 04 - FWHM/beta bug fixed by Will High.
;         Apr 03 - Variable pixel size from command line added by RM.
;                  This renders HDF,PC,ACS options obsolescent.
;         Mar 03 - ACS option added by RM.
;         Nov 02 - Difference between SNAP and HST added by RM
;                  (convolve HDF galaxies with this to get SNAP galaxies).
;         Oct 02 - HST Planetary camera PSF added by RM.
;                  (same as HST WFPC, but with a different pixel scale).
;         Jun 02 - SNAP PSF added by RM.
;         Apr 02 - HST WFPC (TinyTim) PSF incorporated by RM.
;         Apr 02 - Written by R.Massey
;-

COMPILE_OPT idl2

name=""                                         ; Default telescope name
if not keyword_set(pixsize) then pixsize=0.04   ; Default "/pixel

if keyword_set(hst) then begin
  
  ; Print intention to screen
  message,'Reading in TinyTim model of HST PSF',/info,noprint=silent
  
  ; Set pixel scale (for backwards compatability)
  if keyword_set(hdf) then pixsize=0.04  ; DITHERed HDF has 0.04"/pix
  if keyword_set(pc)  then pixsize=0.046 ; Planetary Camera (in centre) has 0.046"/pix 
  if keyword_set(acs) then pixsize=0.05  ; ACS wide field imager has 0.05"/pix 
  
  ; Read in HST TinyTim PSF image and decompose it into shapelets
  fits_read,shapelets_paths(4,silent=silent)+'TinyTim_PSF.fits',psf; This has 0.01"/pix (10x oversampled WFPC)

  ; Set optimum parameters for decomposition
  x0=[150.50,150.48]
  if not keyword_set(n_max) then n_max=8
  beta=5. ; 3.89 minimises chisq but 5 also captures a bit of the 3rd diffraction wing
  decomp=shapelets_decomp(psf,beta,n_max,recomp=recomp,centre=x0,/silent)
  decomp.beta=decomp.beta/10.   ; Now has 0.1"/pix (regular WFPC)
  decomp.n_pixels=[32,32]        ; Arbitrary image viewing parameters
  decomp.x=[16,16]               ; 

  ; Add in effect of charge diffusion to neighbouring pixels
  shapelets_read_psf,ch_diff,.355*2.3548,/silent
  shapelets_convolve,decomp,ch_diff,gamma=0.6,n_max_g=n_max
  
  ; Set final pixel scale
  decomp.beta=decomp.beta*(0.1/pixsize)
  
  ; Renormalise to unity explicitly so that convolution conserves flux exactly
  shapelets_recomp,decomp,recomp
  weight=total(recomp) 
  decomp.coeffs[*]=decomp.coeffs[*]/weight
  recomp=recomp/weight
  decomp.name='HST TinyTim'

endif else if keyword_set(snap) then begin
  
  ; Print intention to screen
  message,'Reading in predicted model of SNAP PSF',/info,noprint=silent

  ; Read in Jason`s SNAP PSF and decompose it into shapelets
  fits_read,shapelets_paths(4,silent=silent)+'SNAP_psf_800nm_01radians.fits',psf ; This has 0.016 arcsec/pixel
  psf=psf[923:1122,923:1122]
  
  ; Set optimum parameters for decomposition
  x0=[99.95,101.5]
  if not keyword_set(n_max) then n_max=12
  w=5.3 ; 5.8 with nmax of 8, 7.5 with 12 to get 4th ring
  decomp=shapelets_decomp(psf,w,n_max,recomp=recomp,centre=x0,/silent)
  decomp.beta=decomp.beta*(0.016/pixsize) ; }
  decomp.n_pixels=[30,30]                 ; } Image file was oversampled
  decomp.x=[15.5,15.5]                    ; }

  ; Add in effect of charge diffusion to neighbouring pixels
  if not keyword_set(diff_length) then diff_length = 4.0
  shapelets_read_psf,ch_diff,diff_length*2.3548*0.10,/silent
  shapelets_convolve,decomp,ch_diff,gamma=0.6,n_max_g=n_max
 
  ; Renormalise so that convolution conserves flux
  shapelets_recomp,decomp,recomp
  weight=total(recomp)
  decomp.coeffs[*]=decomp.coeffs[*]/weight
  ;recomp=recomp/weight
  decomp.name='SNAP full'

endif else if keyword_set(diffsnap) then begin
  
  ; Print intention to screen
  message,'Calculating difference between HST and SNAP PSFs',/info,noprint=silent

  ; Read in Jason`s SNAP PSF model
  fits_read,shapelets_paths(4,silent=silent)+'SNAP_psf_800nm_01radians.fits',psf ; This has 0.016 arcsec/pixel
  psf=psf[923:1122,923:1122]
  
  ; Read in TinyTim model of HST PSF to deconvolve from the SNAP PSF
  shapelets_read_psf,hst,/hst,n_max=n_max,pixsize=0.016,/silent
 
  ; Set optimum parameters for decomposition
  x0=[99.95,101.5]
  if not keyword_set(n_max) then n_max=12
  w=4.8   ; This is still not really optimised!
  decomp=shapelets_decomp(psf,w,n_max,recomp=recomp,centre=x0,psf=hst,/silent)
  decomp.beta=decomp.beta*(0.016/pixsize) ; }
  decomp.n_pixels=[30,30]                 ; } Image file was oversampled
  decomp.x=[15.5,15.5]                    ; }
  decomp.name='SNAP difference'
 
  ; Renormalise so that convolution conserves flux
  shapelets_recomp,decomp,recomp2
  weight=total(recomp2)
  decomp.coeffs[*]=decomp.coeffs[*]/weight
  recomp2=recomp2/weight

endif else if keyword_set(gems) then begin
  
  ; Print intention to screen
  message,'Reading in F606W ACS PSF from GEMS images',/info,noprint=silent
  
  ; Set pixel scale (for backwards compatability)
  if not keyword_set(pixsize) then pixsize=0.03

  ; Set optimum parameters for decomposition
  if not keyword_set(n_max) then n_max=14
  beta=2.13
  x0=[63.55,63.62]

  ; Read in ACS image of stacked stars and decompose it into shapelets
  fits_read,shapelets_paths(4,silent=silent)+'GEMS_ACS_F606W.fits',psf ; This already has 0.03"/pix
  decomp=shapelets_decomp(psf,beta,n_max,recomp=recomp,centre=x0,psf=hst,/silent)
  
  ; Set final pixel scale
  decomp.beta=decomp.beta*(0.03/pixsize)
  
  ; Renormalise to unity explicitly so that convolution conserves flux exactly
  shapelets_recomp,decomp,recomp
  weight=total(recomp) 
  decomp.coeffs[*]=decomp.coeffs[*]/weight
  recomp=recomp/weight
  decomp.name='ACS F606W'

endif else if keyword_set(subaru) then begin

  if keyword_set(sindex) then begin  ; Subaru PSF index
    ; should put checks on sindex
    c=sindex[0]
    r=sindex[1]
  endif else begin
    print,'PSF selection is out of bounds!, setting c=17,r=3'
    c=17
    r=3
  endelse



  message,'Reading in Subaru PSF catalog',/info,noprint=silent
  shapelets_read_shapecat,all_psf,'PSF/Subaru_PSFs_Molly',silent=silent
  psf_index=all_psf.grid[c,r]
  decomp=shapelets_shapecat2decomp(all_psf,psf_index)
  decomp.name='Subaru Grid #'+strtrim(string(c),2)+','+strtrim(string(r),2)

  ; Place in centre of postage stamp
  decomp.x=decomp.n_pixels/2.
  ;  if keyword_set(pixsize) then shapelets_dilate,decomp,pixsize/0.20,/radius
                               ; pixel size for subaru/PSF is 0.20"

  ; Renormalise so that convolution conserves flux
  shapelets_recomp,decomp,recomp
  weight=total(recomp)
  decomp.coeffs[*]=decomp.coeffs[*]/weight
  decomp.n_pixels=decomp.n_pixels-decomp.n_pixels mod 2+1
  decomp.x=decomp.n_pixels/2
  recomp=recomp/weight


endif else begin
  
  ; Print intention to screen
  message,'Creating a Gaussian PSF',/info,noprint=silent

  ; Create a (unit normalised) Gaussian PSF - set FWHM << beta for a delta function
  n_max=0
  shapelets_make_nvec, n_max, n1, n2, n_coeffs
  coeffs=fltarr(n_coeffs) & coeffs[0]=.5/fwhm*2.3548/sqrt(!pi) ; FWHM in pixels, so pixsize=1 by default
  n_pixels=max([fix(fwhm*5+.5),10])
  decomp = {name:"Gaussian_"+strtrim(string(fwhm),2), $
       	  type:"decomp_cartesian",      	         $
       	  history:strarr(1),	       	         $
       	  x:[n_pixels,n_pixels]/2., 	              $
       	  beta:fwhm/2.3548, 			         $
       	  n_max:0,      		 	  	         $
       	  n_coeffs:n_coeffs,		  	         $
       	  coeffs:coeffs,		 	  	         $
       	  coeffs_error:fltarr(n_coeffs),	         $
       	  n_pixels:[n_pixels,n_pixels],             $
       	  moments:0B,  		 	  	         $
       	  sex:0B, 			 	  	         $
       	  n1:n1,  				  	         $
       	  n2:n2,  				  	         $
       	  over:1,			     	  	         $
       	  integrate:1,	          	  	         $
       	  error:fltarr(n_coeffs), 	  	         $
       	  chisq:[0.,0.], 		                   $
       	  skyfit:0}

  if keyword_set(recomp) then shapelets_recomp,decomp,recomp

endelse

end







