;$Id: shapelets_sexcat2pstamp.pro, v2$
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
;       SHAPELETS_SEXCAT2PSTAMP
;
; PURPOSE:
;       Extract a postage stamp of image data around an object detected by
;       SExtractor data. Uses the fact that the inverse variance
;       is zero on saturated stars and around the edge of the HDF.
;
; CATEGORY:
;       Shapelets/ROI extraction.
;
; INPUTS:
;       IMAGE  - Image structure (read in with shapelets_read_image.pro)
;       SEXCAT - SExtractor catalogue of image 
;                (read in with shapelets_read_sexcat.pro)
;       ID     - ID number of the object of interest in the SEx catalogue
;
; OPTIONAL INPUTS:
;       BACK_SIZE - Size of postage stamp used for noise estimation [pixels]
;                   (DEFAULT: 120)
;       N_GROW    - Number of times SExtractor segmentation map is "grown"
;                   to mask objects during noise estimation (DEFAULT: 8)
;       NFWHM     - Size of final postage stamp in units of SExtractor's
;                   major axis (DEFAULT: 4)
;
; KEYWORD PARAMETERS:
;       NOISE_MAP - Use noise map from image structure (rather than local
;                   background estimation).
;       SHOT_NOISE- Set to include photon shot noise in estimated noise maps
;                   (by adding in quadrature to sky background noise).
;       SEG_MAP   - Use locally-determined segmentation map (rather than the
;                   one from SExtractor, stored in the image structure).
;       SQUARE    - Keep the postage stamp square, rather than chopping off
;                   the corners to get a circular region.
;       BORDER    - Increase border around the edges of postage stamp, to
;                   make pretty plots in papers but for useless elsewhere.
;       NEIGHBOUR - treatment of neighbours: 0(DEFAULT): infinite
;                   errors, i.e. unconstrained fit, 1: set neighbour
;                   pixels to the background level and associated errors
;       LAZY      - Don't bother calculating the noise level if the postage
;                   stamp looks bad at an early stage. Setting this can
;                   avoid triggering crashes in SKY.pro and MMM.pro near zero
;                   padding around the edge of an image.
;       PLOTIT    - Plot the postage stamp to the screen.
;       SILENT    - Operate silently.
;
; OUTPUTS:
;       PSTAMP - Postage stamp structure
;
; MEANINGS OF FLAG: (obects with flags greater than 3 should be discarded)
;       0 - OK
;       1 - Nearby object
;       2 - Severe overlapping with nearby object
;       3 - Object is near a saturated pixel
;       4 - Object is near a masked region
;       5 - Object is near the edge of the image
;       6 - Object is itself masked out
;       7 - Object has zero FWHM
;       8 - Too few background pixels around object
;       9 - Object entirely overlapped by neighbours
;      10 - This routine crashed
;
; MODIFICATION HISTORY:
;       Apr 06 - Bug srahing noise estimation without a se map fixed by RM.
;       Jan 06 - Masking of pixels during noise calculation improved by RM.
;       Jan 06 - Photon shot noise incorporated into noise maps by BR.
;       Jan 06 - LAZY keyword avoids computing noise for high flags by JB.
;       Jan 06 - Errors from very negative image pixels caught by JB.
;       Nov 05 - Integers used for image pixel coordinates by Barnaby Rowe.
;       Apr 05 - Routine restructured and made more robust by RM.
;       Apr 05 - Format of tags used in image structures updated by RM.
;       Apr 05 - Bug fixed in image masking by RM.
;       Jan 05 - Default postage stamp geometry changed to square by RM.
;       Jan 05 - Bug fixed for neighbours in circular postage stamp by AR.
;       Jan 05 - NEIGHBOUR keyword added by AR. Spelling corrected by RM!
;       Jul 04 - Default behaviour of /seg_map flipped by RM.
;       Jul 04 - sexid field added to the pstamp structure by AR.
;       May 04 - Flags improved, and only the largest one kept by Joel Berge.
;       Mar 04 - Need for external noise and segmentaion maps removed by AR.
;       Mar 04 - Option to use circular postage stamps by Alexandre Refregier.
;       Dec 02 - Border flag added by RM to make nice plots in papers.
;       Apr 02 - Shapelet structures incorporated by RM.
;       Dec 01 - Written by Richard Massey.
;-

function shapelets_sexcat2pstamp, IMAGE, SEXCAT, ID,                 $
                                  NOISE_MAP=noise_map,               $
                                  SEG_MAP=seg_map,                   $
                                  SATURATION_LEVEL=saturation_level, $
                                  SQUARE=square,                     $
                                  BORDER=border,                     $
                                  BACK_SIZE=back_size,               $
                                  N_GROW=n_grow,                     $
                                  VERY_LOCAL_NOISE_CALC=very_local_noise_calc,$
                                  SHOT_NOISE=shot_noise,             $
                                  LAZY=lazy,                         $
                                  NFWHM=nfwhm,                       $
                                  NEIGHBOUR=neighbour,               $
                                  N_PIXELS=n_pixels,                 $
                                  SILENT=silent,                     $
                                  PLOTIT=plotit

COMPILE_OPT idl2

; Initialise default parameters
if not keyword_set(nfwhm) then nfwhm=5.         ; Size of final postage stamp in units of SExtractor major axis length
if not keyword_set(back_size) then back_size=200; Size of postage stamp used for noise estimation [pixels]
if n_elements(n_grow) eq 0 then n_grow=8        ; Number of times SExtractor segmentation map is "grown" to mask objects during noise estimation
nfwhm = nfwhm + keyword_set(border)             ; Add a border around postage stamp, to look pretty in plots?
nbgpix_min = 10                                 ; Minimum number of pixels required to evaluate background
mask_neigh = 2.75                               ; Size of region to excise from neighbours in units of Sextractor a and b
on_error,2 & pstamp={flag:10}                   ; If routine crashes, return failure flag
flag=0                                          ; Otherwise, start by assuming success!
flag_interpret=["OK", $                         ; Meanings of flags ; 0
                "Nearby object",                                  $ ; 1
                "Severe overlapping with nearby object",          $ ; 2
                "Object is near a saturated pixel",               $ ; 3
                "Object is near a masked region",                 $ ; 4
                "Object is near the edge of the image",           $ ; 5
                "Object is itself masked out",                    $ ; 6
                "Object has zero FWHM",                           $ ; 7
                "Too few background pixels around object",        $ ; 8
                "Object entirely overlapped by neighbours",       $ ; 9 
                "Routine shapelets_sexcat2pstamp crashed"]          ; 10
flag_interpret_mini=["OK","Nearby object","Severe overlap",$
                     "Near saturation","Near mask","Near edge",$
                     "Masked out","FWHM=0?",$
		                 "Little bkgrd","No backgd","Fatal crash"] 


; Parse input
if shapelets_structure_type(image,/silent) and shapelets_structure_type(sexcat,/silent) then begin
  ; I keep accidentally swapping image and sexcat around. Cope with this.
  imagetype=strupcase(image.type)
  sexcattype=strupcase(sexcat.type)
  if imagetype ne "IMAGE" or sexcattype ne "SEXCAT" then begin
    if imagetype eq "SEXCAT" and sexcattype eq "IMAGE" then begin
      junk=temporary(image)
      image=temporary(sexcat)
      sexcat=temporary(junk)
    endif else message,"Inputs need to be a sexcat and an image structure!"
  endif
endif else begin
  ; Cope with the image being input as just an array, rather than a shapelets structure.
  if not shapelets_structure_type(sexcat,/silent) then begin
    if not shapelets_structure_type(image,/silent) then begin
      message,"Structure types not recognised!"
    endif else begin
      if strupcase(image.type) eq "SEXCAT" then begin
        junk=temporary(image)
        image=temporary(sexcat)
        sexcat=temporary(junk)
      endif else begin
        message,"Structure types not recognised!"
      endelse
    endelse
  endif
  if not shapelets_structure_type(image,/silent) then begin
    if size(image,/N_DIMENSIONS) eq 2 and size(image,/TYPE) eq 4 then begin
      image={name:"Unknown",                   $
    	       type:"image",                     $
    	       history:strarr(1),                $
      	     image:image,                      $
      	     n_pixels:size(image,/DIMENSIONS), $
      	     pixel_scale:1.,                   $
    	       header:strarr(1),                 $
    	       mask:0B,                          $
    	       noise:0.,                         $
    	       seg:0B}
      message,"Input image array converted into a shapelets image structure!",/info
    endif else begin
      message,"Structure types not recognised!"
    endelse
  endif
endelse
; Decide which number object to extract a postage stamp around.
case n_elements(id) of
  0: if sexcat.n eq 0 then id=0 else message,"No object specified to extract from the catalogue!" 
  1: 
  else: begin
          message,"More than one object specified to extract from the catalogue! The list has been truncated.",/info
          id=id[0]
        end
endcase


; Objects detected in the noisy regions around the edges of drizzled images,
; where there has only been a single exposure are often just noise.
; SExtractor typically gives these (and cosmic rays) zero FWHM
if sexcat.fwhm[id] eq 0. then begin
  flag=[flag,7]
  message,"This object has zero FWHM according to Sextractor! Is it just noise?",/info,noprint=silent
endif


; Check if the object lies within the image
if sexcat.x[id,0] lt 0 or sexcat.x[id,0] gt image.n_pixels[0] or $ 
   sexcat.x[id,1] lt 0 or sexcat.x[id,1] gt image.n_pixels[1] then $
   message,"Object lies outside the image!"


; Guess at bounds of postage stamp using SExtractor information
r_pstamp=round(nfwhm*sexcat.a[id]+4.) ; Postage stamp radius [pixels]
xmin = 0 > (fix(sexcat.x[id,0])-r_pstamp) < (image.n_pixels[0]-1)
xmax = 0 > (fix(sexcat.x[id,0])+r_pstamp) < (image.n_pixels[0]-1)
ymin = 0 > (fix(sexcat.x[id,1])-r_pstamp) < (image.n_pixels[1]-1)
ymax = 0 > (fix(sexcat.x[id,1])+r_pstamp) < (image.n_pixels[1]-1)


; Check if the object is big enough to be worth bothering with
; (also because the code may crash if it's just a single pixel)
n_bit = [xmax-xmin+1,ymax-ymin+1]
n_in_bit = n_bit[0]*n_bit[1]
if n_in_bit le 1 then begin 
  if n_in_bit le 0 then begin
    message,"Postage stamp contains no pixels!",/info,noprint=silent
  endif else begin
    message,"Postage stamp contains only one pixel. Erratic behaviour possible!",/info,noprint=silent
  endelse
  return,{flag:max([flag,7])}
endif


; Check if the object is near an edge
if xmin eq 0. or ymin eq 0. or xmax eq (image.n_pixels[0]-1) or $
   ymax eq (image.n_pixels[1]-1) then begin
  message,"Object is near the edge of the image!",/info,noprint=silent
  flag=[flag,5]
endif


; Extract region of interest
bit = image.image[xmin:xmax,ymin:ymax]


; Make postage stamp circular by masking off pixels in corners
if keyword_set(square) then begin
  n_corner=0
  geometry="square"
  not_corner=indgen(n_bit)
endif else begin
  shapelets_make_xarr,n_bit,x1_c,x2_c,x0=sexcat.x[id,*]-[xmin,ymin] ; construct x-y map centered on the object
  corner=where(round(sqrt(x1_c^2+x2_c^2)-r_pstamp) gt 0,n_corner)
  not_corner=where(round(sqrt(x1_c^2+x2_c^2)-r_pstamp) le 0,n_not_corner)
  if n_corner gt 0 then bit[corner]=0.
  geometry="circle"
endelse


; Construct segmented pixel map
message,"Constructing local segmentation map",/info,noprint=silent
bits=fltarr(n_bit[0],n_bit[1])
objsize=mask_neigh*sexcat.a
if keyword_set(square) then begin
  neigh=where(sexcat.id ne id and sexcat.x[*,0]+objsize ge xmin and sexcat.x[*,0]-objsize le xmax and $
       sexcat.x[*,1]+objsize ge ymin and sexcat.x[*,1]-objsize le ymax, n_neigh)
endif else begin ; circle (default)
  neigh=where(sexcat.id ne id and $ 
    sqrt((sexcat.x[*,0]-sexcat.x[id,0])^2+(sexcat.x[*,1]-sexcat.x[id,1])^2) lt r_pstamp+objsize,n_neigh)
endelse
if n_neigh lt 1 then objin=[id] else objin=[neigh,id]	   ; all objects in the postage stamp
n_objin=n_elements(objin)     ; number of neighbouring objects including object of interest
shapelets_make_xarr,n_bit,x1,x2,x0=[.5,.5]   ; construct x-y map
; Should there be a response here if the xbugfixed keyword is set?
x1=x1+xmin & x2=x2+ymin
for i=0,n_objin-1 do begin
  a=sexcat.a[objin[i]] & b=sexcat.b[objin[i]]
  theta=sexcat.theta[objin[i]]/!radeg
  dx1=x1-sexcat.x[objin[i],0] & dx2=x2-sexcat.x[objin[i],1]
  r2_ellip=(dx1*cos(theta)+dx2*sin(theta))^2/a^2+(-dx1*sin(theta)+dx2*cos(theta))^2/b^2
  gd=where(r2_ellip lt mask_neigh^2,n_gd)
  if n_gd gt 0 then begin
    if objin[i] eq id then begin       ; the object we're interested in is the last one
      if max(bits[gd]) gt 0 then begin ; are there any overlaps?
        flag=[flag,2]
        message,"Severe overlapping with nearby object",/info,noprint=silent
      endif
    endif
    bits[gd]=objin[i]+1 ; assign pixels to this object. Following Sextractor indexing convention to start from 1 not 0
  endif
endfor
if keyword_set(seg_map) then begin
  message,"Using locally estimated segmentation map",/info,noprint=silent
endif else begin
  if not tag_exist(image,"seg") then begin
    message,"Segmentation map not supplied in image structure",/info,noprint=silent  
  endif else begin
    if n_elements(image.seg) le 1 then begin
      message,"Segmentation map not supplied in image structure",/info,noprint=silent  
    endif else if (not keyword_set(seg_map)) then begin
      message,"Using externally supplied segmentation map",/info,noprint=silent
      bits=image.seg[xmin:xmax,ymin:ymax]
    endif
  endelse
endelse
if n_corner gt 0 then bits[corner]=0.


; Construct image mask
mask=bytarr(n_bit[0],n_bit[1]) ; by default, assume everything is good
if tag_exist(image,"mask") then if n_elements(image.mask) gt 1 then begin
  mask=image.mask[xmin:xmax,ymin:ymax]
  if n_corner gt 0 then mask[corner]=1B
endif else if image.mask[0] eq 1 then mask[*]=1
if total(mask[not_corner]) gt 0 then flag=[flag,4]


; Check that we can see our object, from its pixels on the segmentation map 
obpix = where((bits eq (id+1)),nobpix)
if nobpix eq 0 then begin
  flag=[flag,9]
  message,"Object is entirely covered by neighbours",/info,noprint=silent
endif else if max(mask[obpix]) eq 1 then begin
  flag=[flag,6]
  message,"Object is (at least partly) behind the mask",/info,noprint=silent
endif


; Construct noise map (but only if the postage stamp looks decent so far)
if keyword_set(lazy) and max(flag) ge 5 then begin  
  ; Don't bother if the postage stamp already looks useless
  ; Create a placeholder noise map
  noise=replicate(1.,n_bit[0],n_bit[1])
  back_mean_ext=0
  back_rms_ext=0
  back_mean_local=0
  back_rms_local=0
endif else begin


  ; METHOD #1: Use supplied noise map to find the background level and rms
  back_mean_ext=0
  ; Find background pixels in the postage stamp itself
  bgpix = where(bits eq 0,nbgpix)
  ; Read in the noise (background rms) from the image structure
  if tag_exist(image,"noise") then begin
    if n_elements(image.noise) eq 1 then begin
      ; Constant inverse variance
      noise_ext=replicate(image.noise[0],n_bit[0],n_bit[1])
      if image.noise[0] eq 0 then begin
        message,"External noise map is not provided, or there is no noise in image!",/info,noprint=silent
        back_rms_ext=0 
      endif else back_rms_ext=1./sqrt(abs(image.noise[0]))
      message,"External noise estimate is "+strtrim(string(back_rms_ext),2),/info,noprint=silent
    endif else begin
      ; Spatially varying inverse variance
      noise_ext=image.noise[xmin:xmax,ymin:ymax]
      if nbgpix gt nbgpix_min then begin
        back_rms_ext=mean(1./sqrt(noise_ext[bgpix]))
        message,"External noise map has rms of "+strtrim(string(back_rms_ext),2),/info,noprint=silent
      endif else begin
        back_rms_ext=mean(1./sqrt(noise_ext))
        message,"WARNING! Only "+strtrim(string(nbgpix),2)+$
          " background pixels in the postage stamp!",/info,noprint=silent
      endelse
    endelse
  endif else begin
    message,"External noise map is not provided!",/info,noprint=silent
    back_rms_ext=0
    noise_ext=replicate(1.,n_bit[0],n_bit[1]) 
  endelse


  ; METHOD #2: Use uncontaminated pixels to locally estimate the background level and rms
  if keyword_set(very_local_noise_calc) then begin
    if nbgpix gt nbgpix_min then begin
      ; Estimate background naively using the background pixels in the science image  
      back_mean_local=mean(bit[bgpix])
      back_rms_local=stddev(bit[bgpix])
      message,"Direct local background estimation [mean, rms]:"+$
        strcompress(string(back_mean_local)+string(back_rms_local)),$
        /info,noprint=silent
      ; Estimate background using iterative 3sigma clipping, as really done by SExtractor  
      bitbg=bit[bgpix]
      ave=back_mean_local & rms=back_rms_local
      mode=2.5*median(bitbg)-1.5*ave
      conv=0B		  ; converged?
      niter=0B
      while not conv do begin
        gd=where(abs(bitbg-mode) lt 3.*rms,n_gd)
        if n_gd lt 10 then conv=1B else begin
          ave_new=mean(bitbg[gd])
          rms_new=stddev(bitbg[gd])
          mode_new=2.5*median(bitbg[gd])-1.5*ave_new
          if abs(rms_new-rms)/rms lt 0.01 then conv=1B
          ave=ave_new
          rms=rms_new
          mode=mode_new
          niter=niter+1
        endelse
      endwhile
      if ave eq 0 then back_mean_local=mean(bit[bgpix]) else back_mean_local=ave
      if rms eq 0 then back_rms_local=stddev(bit[bgpix]) else back_rms_local=rms
      ;noise_local=replicate(1./sqrt(back_rms_local),n_bit[0],n_bit[1])
      noise_local=replicate(back_rms_local^(-2),n_bit[0],n_bit[1])
      message,"Local background estimation after "+strtrim(string(niter),2)+$
        " iteration(s) [mean, sigma(<3rms)]:"+$
        strcompress(string(ave)+string(rms)),/info,noprint=silent
    endif else begin
      flag=[flag,8]
      back_mean_local=0.
      back_rms_local=0.
      noise_local=replicate(1.,n_bit[0],n_bit[1])
      message,"Not enough pixels to evaluate local background",/info,noprint=silent
    endelse
  endif else begin

    ; Extract a (much larger) square postage stamp around the object
    left=round(sexcat.x[id,0]-((back_size/2)>r_pstamp))>0
    right=round(sexcat.x[id,0]+((back_size/2)>r_pstamp))<(image.n_pixels[0]-1)
    bottom=round(sexcat.x[id,1]-((back_size/2)>r_pstamp))>0
    top=round(sexcat.x[id,1]+((back_size/2)>r_pstamp))<(image.n_pixels[1]-1)
    noise_pstamp=image.image[left:right,bottom:top]

    ; Detect objects then grow the region of interest to mask nearby pixels
    if tag_exist(image,"seg") then begin
      if n_elements(image.seg) gt 1 then begin  
        objects=abs(image.seg[left:right,bottom:top])<1
        roi=where(objects,n_roi)
        if n_roi gt 0 then begin
          grown_objects=-objects
          for i=-n_grow,n_grow do begin
            jj=round(sqrt(n_grow^2-i^2))
            for j=-jj,jj do begin
              grown_objects+=transpose(shift(transpose(shift(objects,i)),j))
            endfor
          endfor
          objects=grown_objects<1
          roi=where(objects,n_roi)
          if n_roi gt 0 then begin
            highbad=fix(max(noise_pstamp)+2)-1
            noise_pstamp[roi]=highbad+1
          endif
        endif
      endif else n_roi=0
    endif else n_roi=0

    ; Iterate as done by DAOPHOT, but catching any of the ways SKY.pro or 
    ; MMM.pro temd to crash
    zeronoise=where(noise_pstamp eq 0,n_zero)
    if (4*n_zero gt n_elements(noise_pstamp)) and ($
       n_elements(noise_pstamp)-n_roi) lt nbgpix_min then begin
      if (4*n_zero) gt n_elements(noise_pstamp) then begin
        message,"Zero-valued pixels occupy more than a quarter of the noise postage stamp!",/info;,noprint=silent
        flag=[flag,8]
      endif else begin
        message,"Not enough pixels to evaluate local background",/info,noprint=silent
        flag=[flag,8]
      endelse
      back_mean_local=0.
      back_rms_local=0.
      noise_local=replicate(1.,n_bit[0],n_bit[1])
    endif else begin
      sky,noise_pstamp,back_mean_local,back_rms_local,highbad=highbad,/SILENT
      noise_local=replicate(back_rms_local^(-2),n_bit[0],n_bit[1])
    endelse

  endelse

  ; Decide which estimate (external or local) to use as the default noise map
  if keyword_set(noise_map) and back_rms_ext gt 0. then begin
    message,"Using external noise estimation",/info,noprint=silent
    noise=noise_ext
    back_rms=back_rms_ext
    back_mean=back_mean_ext
  endif else begin
    message,"Using local noise estimation",/info,noprint=silent
    noise=noise_local
    back_rms=back_rms_local
    back_mean=back_mean_local
  endelse

  ; Add photon shot noise
  if keyword_set(shot_noise) then begin
     if max(flag) ge 8 then message,"WARNING: Not adding photon shot noise",/info,noprint=silent else begin
       ; Only add shot noise where there are objects in the segmentation map
       ;add_shot=where(bits ne 0,n_shot)
       ; Only add shot noise where signal is significant at (arbitrary) 3 sigma level
       case image.units of
         "cps": image.image=image.image*exposure_time
         "counts":
         else: message,'Image units "'+image.units+'" are not recognised!',/INFO
       endcase
       add_shot=where(bit gt (3*back_rms+back_mean),n_shot)
       if n_shot gt 0 then noise[add_shot]=(1./noise[add_shot]+(bit[add_shot]-back_mean)>0)^(-1)
       message,"Photon shot noise added to "+strtrim(n_shot,2)+$
         "/"+strtrim(not_corner,2)+" pixels",/info,noprint=silent
    endelse
  endif
  
endelse


; Check segmented pixel map for any other objects nearby
othpix=where((bits ne 0) and (bits ne (id+1)) and (bits ne -1) and (mask eq 0),nothpix)  ; pixels in other sources
if nothpix gt 0 then begin
   flag=[flag,1]       
   n_neigh=n_elements(uniq(bits[not_corner],sort(bits[not_corner])))-2
   message,"Found "+strtrim(string(n_neigh),2)+" neighbour(s)",/info,noprint=silent
endif else n_neigh=0


; Mask out neighbouring objects
if nothpix gt 0 then begin
  if keyword_set(neighbour) then begin
    ; Set the value of neighbours' pixels to the background level
    bit[othpix]=back_mean_local
    print,'bml',back_mean_local
    ; Set the noise in neighbours' pixels to background
    if keyword_set(noise_map) and back_rms_ext gt 0. then begin
      noise[othpix]=1./back_rms_ext^2
    endif else if back_rms_local gt 0. then begin
      noise[othpix]=1./back_rms_local^2 
    endif else begin
      noise[othpix]=0.
    endelse
  endif else begin
    ; Set the noise in neighbours' pixels to be infinite (then their values do not matter)
    noise[othpix]=0.
  endelse 
endif
if not keyword_set(square) then if n_corner gt 0 then noise[corner]=0.


; Convert centroid guess into postage-stamp coordinates
x0=[sexcat.x[id,0],sexcat.x[id,1]]-[xmin,ymin]


; Store in an "pstamp" structure. 
; This definition contains lots of redundant information and is a lottle sloppy.
pstamp={name:image.name+"_#"+strtrim(string(id),2), $
        type:"pstamp",                    	        $ ; Structure type, for object-oriented approach.
        image:bit,                        	        $ ; Image data,
        noise:noise,                      	        $ ; Pixel weight map
        seg:bits,                         	        $ ; Image segmentation map
        mask:mask,                        	        $ ; Image mask (0=good pixel, 1=bad pixel)
        n_pixels:n_bit,                   	        $ ; Size of image arrays
        units:image.units,                          $ ; Image units, inherited from image
        pixel_scale:image.pixel_scale,              $ ; Linear pixel size, inherited from image
        photo_zp:image.photo_zp,                    $ ; Photometric zero point, inherited from image
        exposure_time:image.exposure_time,          $ ; Exposure time, inherited from image
        flag:max(flag),                   	        $ ; 0=Fail, 1=OK, 2=Nearby obj, 4=Noise normalised
        flag_interpret:flag_interpret,    	        $ ; Meanings of the flags
        flag_interpret_mini:flag_interpret_mini,    $ ; Meanings of the flags
        im_ran:[xmin,xmax,ymin,ymax],     	        $ ; Coordinates of BL and TR of postage stamp
        geometry:geometry,                	        $ ; Is it a square, a circle, an ellipse, or something else?
        xo:x0[0],                         	        $ ; Centroid guess in postage stamp coordinate
        yo:x0[1],                         	        $ ;  -----"-----
        back_mean_local:back_mean_local,  	        $ ; Local background estimation - mean
        back_rms_local:back_rms_local,    	        $ ;                             - rms
        back_mean_ext:0.,                 	        $ ; background extimation from an  external map (if available) - mean
        back_rms_ext:back_rms_ext,        	        $ ;                                                            - rms
        sexid:id+1,                       	        $ ; SExtractor ID number
        x:sexcat.x[id,0],                 	        $ ; SExtractor centroid guess in image coordinates
        y:sexcat.x[id,1],                 	        $ ;  -----"-----
        xpeak:sexcat.xpeak[id,0],         	        $ ; \
        ypeak:sexcat.xpeak[id,1],         	        $ ;  |
        xmin:sexcat.xmin[id,0],           	        $ ;  |
        ymin:sexcat.xmin[id,1],           	        $ ;  |
        xmax:sexcat.xmax[id,0],           	        $ ;  |
        ymax:sexcat.xmax[id,1],           	        $ ;  | Other SExtractor
        area:sexcat.area[id],             	        $ ;   > parameters, all in
        a:sexcat.a[id],                   	        $ ;  |  image coordinates
        b:sexcat.b[id],                   	        $ ;  |
        theta:sexcat.theta[id],           	        $ ;  |
        class:sexcat.class[id],           	        $ ;  |
        nu:sexcat.flux[id],               	        $ ;  |
        fwhm:sexcat.fwhm[id],             	        $ ;  |
        mag:sexcat.mag[id]}               	          ; /


; Report flag status to screen
message="Pstamp flag: "+strtrim(pstamp.flag,2)
if n_elements(flag) eq 1 then message=message+", ["+strjoin(strtrim(flag,2),",")+"]"
message,message,/info,noprint=silent


; Plot the postage stamp (and SExtractor parameters) to screen
if keyword_set(plotit) then shapelets_plot_pstamp,pstamp,/mask,title="!6Postage stamp #"+strtrim(id,1);,/oframe


; Return postage stamp structure
return,pstamp

end
