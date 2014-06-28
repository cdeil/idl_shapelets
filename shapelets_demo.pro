pro shapelets_demo, image, sexcat, decomp, psf, shapecat

;$Id: shapelets_demo.pro, v2$
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
;       SHAPELETS_DEMO
;
; PURPOSE:
;       Runs through the main shapelet routines to demonstrate
;       their use and to check that they are installed correctly.
;
; CATEGORY:
;       For demonstration purposes only.
;
; CALLING PROCEDURE:
;       shapelets_demo, image, sexcat, decomp, psf, shapecat
;
; INPUTS:
;       Accepts previously defined values of parameters, to save calculating
;       them anew each time this routine is executed. None are necessary.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       Various plots to screen, a SExtractor catalogue and a shapelet
;       catalogue to disc.
;
; DATA REQUIRED:
;       Requires shapelets_demo data available for download from the shapelets
;       website http://www.astro.caltech.edu/~rjm/shapelets/.
;       This includes a small section of the HDF-N and an oversampled image
;       of the TinyTim model of the HST PSF. These should all be placed in
;       the subdirectories described by shapelets_paths.pro
;
; PROCEDURES USED:
;       Virtually all of the shapelet routines!
;
; MODIFICATION HISTORY:
;       Feb 05 - Floating point underflow warning eliminated by RM.
;       Jan 05 - Names of plotting routines rationalised by RM.
;       Sep 03 - Finalised for public release of shapelet code by RM.
;       Mar 02 - Written by Richard Massey
;-

COMPILE_OPT idl2

; Check that IDL can correctly locate the shapelets_demo image
imagefile=shapelets_paths(3)+"example.fits"
if not file_test(imagefile) then $
  message,"Please download the example data to "+imagefile+$
          ", and correct shapelets_paths.pro"


; Read in shapelets_demo image and pixel noise map (image structure)
if not keyword_set(image) then shapelets_read_image,image,"example"

window,0,title="Example image (part of HDF-N)",$ ; Open graphics window
       xsize=512,ysize=512                       ;
shapelets_plot,alog10(image.image>1e-5)          ; Show image (on log scale)
read,"Type a number to continue...",junk         ; Pause for user input
wdelete                                          ; Close graphics window


; Read in SExtractor catalogue for shapelets_demo image (sexcat structure)
if not keyword_set(sexcat) then begin
  if not file_test(shapelets_paths(3)+"example.sex") then $
    sex,"example",telescope="HDF",filter="gauss_2.0_5x5"
  shapelets_read_sexcat,sexcat,"example"
  shapelets_read_image,image,"example"
endif


window,0,title="Size-magnitude diagram for example image"
plot,sexcat.mag+22.1,sexcat.fwhm*0.04,psym=2,$   ; Plot size-magnitude diagram
     xtitle="Magnitude",ytitle="FWHM (arcsec)"   ;
help,sexcat,/structure                           ; Lists elements of structure
read,"Type a number to continue...",junk         ; Pause for user input
wdelete                                          ; Close window


; Extract small, "postage-stamp" image around one object (pstamp structure)
b_id=(sort(sexcat.mag))[0]                       ; Select brightest in field
pstamp=shapelets_sexcat2pstamp(image,sexcat,b_id); Extract postage stamp
window,0,title="Spiral galaxy",xsize=500,ysize=555
shapelets_plot,pstamp.image,title="ORIGINAL - Pretty spiral!!",$
 /frame,/isotropic,/cbar,crange=crange


; Decompose it into shapelets (decomp structure)
if not keyword_set(decomp) then $
 decomp=shapelets_decomp(pstamp,5.3,20,name="HDF Spiral",$
  centre=[357.7-pstamp.im_ran[0],178.8-pstamp.im_ran[2]]);,psf=psf
window,1,title="Spiral galaxy",xsize=500,ysize=555
shapelets_plot,decomp,title="SHAPELET MODEL - Pretty good reconstruction!!",$
 /frame,/isotropic,/cbar,crange=crange
read,"Type a number to continue...",junk


; Display shapelet coefficient matrix
shapelets_plot,decomp,/coef,/cbar,title="Cartesian shapelet coefficients"
read,"Type a number to continue...",junk
shapelets_plot,decomp,/polar,/cbar,title="Polar shapelet coefficients"
read,"Type a number to continue...",junk


; Perform some image manipulations in shapelet space
decomp_manip=decomp                              ; Don"t alter our good model!
shapelets_rotate,decomp_manip,45                 ; Rotate 45 degrees a/c/w
shapelets_shear,decomp_manip,[0.04,0.03],extend=5; Shear (and increase n_max)
shapelets_dilate,decomp_manip,0.05               ; Enlarge by factor 1.1
;shapelets_plot,decomp_manip                     ; Equivalent to next 2 lines
shapelets_recomp,decomp_manip,recomp_manip
shapelets_plot_image,recomp_manip,/frame,/cbar,$
          title="SHAPELET MODEL - After some manipulation"
read,"Type a number to continue...",junk
wdelete,0
wdelete,1


; More clever shapelet decomposition (optimising n_max and beta automatically)
decomp_optimal=shapelets_focus(pstamp,focus=focus,recomp=recomp_optimal, $
      n_max_range=[2,14],theta_min_geom=0.2,/verbose);,psf=psf

	 
; Read in PSF model for the example image (decomp structure)
if not keyword_set(psf) then shapelets_read_psf,psf,/HST,pixsize=0.04,n_max=12
window,0,title="TinyTim PSF",xsize=500,ysize=555
shapelets_plot,psf
read,"Type a number to continue...",junk
wdelete,0


; Create and read in shapelet coefficient catalogue (shapecat structure)
if not keyword_set(shapecat) then begin
  if not file_test(shapelets_paths(3)+"example.shapecat") then $
    shex,"example",sexcat=sexcat,image=image,n_max=15;,/plot
  shapelets_read_shapecat,shapecat,"example",/moments
  decomp=shapelets_shapecat2decomp(shapecat,(sort(shapecat.mag))[0])
endif


; Match up objects in the SExtractor and shapelet catalogues
match=intarr(shapecat.n)
for i=0L,n_elements(match)-1 do match[i]=where(sexcat.id eq shapecat.sexid[i])


; Display some plots comparing SExtractor and shapelets output
window,0,title="SExtractor vs shapelets",xsize=500,ysize=555
zpoint=22.08
plot,sexcat.mag[match]+zpoint,shapecat.mag+zpoint,psym=1,/isotropic,$
     xran=[20,32],/xstyle,yran=[20,32],/ystyle,$
     title="SExtractor vs shapelets",$
	xtitle="SExtractor magnitude",ytitle="Shapelet magnitude"
read,"Type a number to continue...",junk
plot,[0,0],[0,0],xran=[-0.1,0.1]*2,yran=[-0.1,0.1]*2,/nodata,/isotropic,$
     title="Difference between SExtractor and shapelet centroids",$
	xtitle="x offset [arcsec]",ytitle="y offset [arcsec]"
plots,transpose(shapecat.centroid-sexcat.x[match,*])*0.04,psym=1
read,"Type a number to continue...",junk
plot,sexcat.fwhm[match],shapecat.beta,psym=1,$
     xran=[0,40],/xstyle,yran=[0,8],/ystyle,$
     title="SExtractor vs shapelets",$
	xtitle="SExtractor FWHM",ytitle="!6Optimal shapelet !7b!6"
read,"Type a number to continue...",junk
wdelete,0


; Print a friendly "good bye" message
underflow=check_math(mask=32) ; Ignore floating point UNDERflows
print & print
print,"Thank you for running the shapelets demonstration routine!"
print
print,"We hope that this has provided some useful template code."
print,"For further information, please see the shapelets web page"
print,"at http://www.astro.caltech.edu/~rjm/shapelets/"
print
print,"Richard Massey and Alexandre Refregier."
print

end
