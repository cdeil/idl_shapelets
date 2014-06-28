pro shapelets_wl_pipeline, filename,            $
                           STEP=step,           $
                           TELESCOPE=telescope, $
                           SEEING=seeing,       $
                           PLOTIT=plotit

;$Id: shapelets_wl_pipeline.pro, v1$
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
;      SHAPELETS_WL_PIPELINE
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      As automatically as possible, do everything needed to create a 
;      PSF-corrected shapelet catalogue from an image, to measure weak
;      lensing.
;
; INPUTS:
;      FILENAME - The name of the image file, without path or extension.
;
; OPTIONAL INPUTS:
;      
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      Lots of stuff, all written to disc.
;
; NOTES:
;      Under construction...
;
; MODIFICATION HISTORY:
;      Nov 05 - Tidied up a little by RM.
;      Jun 05 - Written by Richard Massey.
;-


; Decide which file to work on
if not keyword_set(filename) then begin
  filename=""
  read,"Which file? ",filename
endif
if file_test(shapelets_paths(1,/SILENT)+filename+".fits") eq 0 then $
  message,"Image "+filename+".fits not found!"


; Decide what still needs doing for this image
if not keyword_set(step) then begin
  if file_test(shapelets_paths(3,/SILENT)+filename+".sexcat") eq 0 then begin
    step=1
  endif else begin
    if file_test(shapelets_paths(3,/SILENT)+filename+"_deconv_galaxies.shapecat") eq 0 then begin
      step=5
    endif else step=!values.f_infinity
  endelse
endelse


; Do it!
switch step of

  1: begin               ; (1) Run SExtractor
       if not file_test(shapelets_paths(3)+filename+".sexcat") then $
         sex,filename,telescope=telescope,seeing=seeing,/seg;,sexcat=sexcat
       ; Don't hang around; this was probably run in batch mode.
       return
     end


  2: begin               ; (2) Interactively select stars, galaxies etc
       ; Read in SExtractor catalogue
       if not keyword_set(sexcat) then shapelets_read_sexcat,sexcat,filename
       ; Begin interactive session
       shapelets_interactive_select_stars, sexcat, iselect
       ; Terminate program, because this will probably be run up to here for many images
       return
     end


  3: begin               ; (3) Mask image
       shapelets_make_image_mask,filename,ISELECT=iselect,SEXCAT=sexcat,IMAGE=image,/STARS_CIRCLE,/PLOTIT
     end

   
  4: begin               ; (4) Run shex on stars
       ; Read in image and its SExtractor catalogue
       if not keyword_set(sexcat) then shapelets_read_sexcat,sexcat,filename
       if not keyword_set(image) then shapelets_read_image,image,filename,/ESTIMATE_NOISE
       ; Read in star selection
       if not keyword_set(iselect) then restore,shapelets_paths(3,/SILENT)+sexcat.name+".iselect"
       ; Decompose stars into shapelets
       shex,filename,image=image,sexcat=sexcat,index=iselect.stars,$
         output_filename=filename+"_stars1";,verb=2,/plot
     end


  5: begin               ; (5) Run shex on stars, using optimum beta and n_max
       ; Read in catalogues
       if not keyword_set(sexcat) then shapelets_read_sexcat,sexcat,filename
       shapelets_read_shapecat,shapecat,filename+"_stars1"
       ; Calculate the range of scales required by the PSF
       theta_min=fltarr(shapecat.n)
       theta_max=fltarr(shapecat.n)
       for i=0,shapecat.n-1 do begin
         shapelets_theta_minmax, shapecat.beta[i], shapecat.n_max[i], theta
         theta_min[i]=theta[0]
         theta_max[i]=theta[1]
       endfor
       ; Make a pretty plot
       good=where(shapecat.flag[*,0] lt 3,n_good)
       plothist,theta_max[good],bin=0.1,xran=[0.1,10],/xlog
       plothist,theta_min[good],bin=0.01,/over
       plothist,theta_max,/linestyle,bin=0.1,/over,col=150
       plothist,theta_min,/linestyle,bin=0.01,/over,col=150
       mean_theta_min=median(theta_min[good])
       mean_theta_max=median(theta_max[good])
       oplot,mean_theta_min*[1,1],[0,1000]
       oplot,mean_theta_max*[1,1],[0,1000]
       ; Decide what n_max and beta to use
       n_max = round( mean_theta_max / mean_theta_min - 0.5 )
       n_max = n_max + (n_max mod 2)
       beta  = mean([mean_theta_max/sqrt( n_max + 0.5 ),mean_theta_min*sqrt( n_max + 0.5 )])
       ; Update positions of objects in the SExtractor catalogue
       sexcat_newpos=sexcat & sexcat_newpos.x[shapecat.sexid,*]=shapecat.x
       ; Run shex with fixed beta and n_max
       if not keyword_set(image) then shapelets_read_image,image,filename,/ESTIMATE_NOISE
       shex,filename,image=image,sexcat=sexcat_newpos,index=shapecat.sexid,$
            n_min=n_max,n_max=n_max,fixed_beta=beta,$
            output_filename=filename+"_stars2";,verb=2,/plot
     end


  6: begin               ; (6) Interpolate PSF to the positions of the galaxies
       ; Read in catalogue with fixed n_max/beta   
       shapelets_read_shapecat,shapecat,filename+"_stars2",$
         /moments,/cartesian;,n_moments=2
       ; Make some catalogue cuts 
       good=where(sqrt(shapecat.ellipticity[*,0]^2+shapecat.ellipticity[*,1]^2)$
                        lt 0.5,n_good)
       if n_good lt shapecat.n then shapecat=shapelets_split_shapecat(shapecat,good)
       ; Read in original SExtractor catalogue, containing the positions of 
       ; galaxies, at which we would like to know the PSF
       if not keyword_set(sexcat) then shapelets_read_sexcat,sexcat,filename
       ; Interpolate the PSF to the positions of the galaxies
       shapecat_interpolated=shapelets_interpolate_psf(shapecat,$
                             xint=sexcat.x[*,0],yint=sexcat.x[*,1],$
                             plotit=plotit,degree=7,max_p=1,pca=0,model=model)   
       shapelets_write, model, filename+"_stars3"
       shapelets_write, shapecat_interpolated, filename+"_stars4"
       ; Make sure all stars' centres of light are over the origin
       shapecat_interpolated_recentred=shapecat_interpolated
       shapelets_recentre, shapecat_interpolated_recentred, TRANSLATION=translation
       shapecat_interpolated_recentred.x=shapecat_interpolated.x
       shapelets_shapecat_moments,shapecat_interpolated_recentred
       ; Write catalogue to disc
       shapelets_write, shapecat_interpolated_recentred, filename+"_stars5"
     end


  7: begin               ; (7) Run shex on everything, simultaneously deconvolving from the PSF
       ; Read everything in to memory
       shapelets_read_shapecat,psf,filename+"_stars5",/moments
       if not keyword_set(sexcat) then shapelets_read_sexcat,sexcat,filename
       if not keyword_set(image) then shapelets_read_image,image,filename,/ESTIMATE_NOISE
       ; Decompose it into shapelets
       shex,filename,image=image,sexcat=sexcat,psf=psf,/diamond,n_min=4,$
         output_filename=filename+"_galaxies1";,verb=2,/plot
     end


  8: begin               ; (8) Rotate objects so that North is up
       ; Look in image header for the PA of the exposure
       header=headfits(shapelets_paths(3)+filename+".fits")
       cdelt1=fxpar(header,"CDELT1")
       cdelt2=fxpar(header,"CDELT2")
       cd1_1=fxpar(header,"CD1_1")/cdelt1
       cd1_2=fxpar(header,"CD1_2")/cdelt1
       cd2_1=fxpar(header,"CD2_1")/cdelt2
       cd2_2=fxpar(header,"CD2_2")/cdelt2
       cd_matrix=[[cd1_1,cd1_2],[cd2_1,cd2_2]]
       rotation=mean([acos(mean([cd1_1,cd2_2])),asin(mean([cd2_1,-cd1_2]))])*!radeg
       ; Rotate shapelet models of all galaxies by this amount
       if not keyword_set(shapecat_deconv) then shapelets_read_shapecat,shapecat_deconv,filename+"_galaxies1_deconv",/moments
       if rotation ne 0 then begin
         message,"Rotating shapecat by "+strtrim(rotation,2)+" degrees",/info
         shapelets_rotate,shapecat_deconv,rotation
       endif
       shapelets_write,shapecat_deconv,filename+"_galaxies2_deconv"
     end


  9: begin               ; (9) Compute shear estimators
       ; Calculate shape moments of galaxies in the shapelet catalogue
       shapelets_read_shapecat,shapecat_deconv,filename+"_galaxies2_deconv",/moments
       shapelets_write,shapecat_deconv,filename+"_galaxies3_deconv"
     end

  else: begin
          print,"Finished!"
        end

endswitch

end
