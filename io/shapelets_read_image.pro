pro shapelets_read_image, IMAGE, FILENAME,                 $
                          FULL_PATH=full_path,             $
                          FITS_READ=fits_read,             $
                          PIXEL_SCALE=pixel_scale,         $
                          UNITS=units,                     $
                          PHOTO_ZP=photo_zp,               $
                          EXPOSURE_TIME=exposure_time,     $
                          NO_MASK=no_mask,                 $
                          NO_SEGMENTATION=no_segmentation, $
                          NO_NOISE=no_noise,               $
                          NOISE_LEVEL=noise_level,         $
                          ESTIMATE_NOISE=estimate_noise,   $
                          N_GROW=n_grow,                   $
			                    SKY_SUBTRACT=sky_subtract,       $
			                    SILENT=silent

;$Id: shapelets_read_image.pro, v2$
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
;         SHAPELETS_READ_IMAGE
;
; PURPOSE:
;         Read in fits image, plus segmentation and inverse variance pixel
;         map if available (it looks for the same filename, but with
;         _objects and _weight added before the extension).
;
; CATEGORY:
;         Data I/O.
;
; INPUTS:
;         FILENAME    - Name of .fits image file.
;
; OPTIONAL INPUTS:
;         NOISE_LEVEL - An estimate of the rms background noise.
;         N_GROW      - Number of times SExtractor segmentation map is "grown"
;                       to mask objects during noise estimation.
;         PIXEL_SCALE - [arcsec] Size of pixels, in case not found in header.
;         UNITS       - Units of the image, in case not found in header.
;         PHOTO_ZP    - Photometric zeropoint, in case not found in header.
;         EXPOSURE_TIM- Exposure time, in case not found in header.
;                          
; KEYWORD PARAMETERS:
;         FULL_PATH   - If this keyword is set, the routine expects "filename"
;                       to contain the full path and extension of the image file.
;                       Without it, the default is to look in the directory
;                       specified within shapelets_paths.pro.
;         ESTIMATE_NOI- Estimate the global noise level from the entire image.
;         NO_MASK     - Don't even attempt to read in the mask image.
;         NO_SEGMENTAT- Don't even attempt to read in the segmentation map.
;         NO_NOISE    - Don't even attempt to read in the weight image.
;         SILENT      - Operate silently.
;
; OUTPUTS:
;         Returns an image structure, containing image, header, noise map, etc.
;
; EXAMPLE USE:
;         shapelets_read_image,image,"example"
;
; MODIFICATION HISTORY:
;         Jan 06 - Memory efficiency improved by RM. 
;         Jan 06 - Image units/zeropoints/exposure times etc. stored by RM. 
;         Jan 06 - Switch from fits_read to readfits() for better headers by RM. 
;         Jan 06 - Noise/seg/mask images that don't exist now completely 
;                  excluded from image structure, rather than set to 0B by RM.
;         Sep 05 - Seg maps with >2^16 objects accomodated by Joel Berge.
;         Jul 05 - Negative seg maps from 32 bit SExtractor fixed by Matt Ferry.
;         Jul 05 - Global noise estimation improved by RM.
;         May 05 - FULLPATH keyword changed to FULL_PATH by RM.
;         Apr 05 - Format of noise tag in image structure changed by RM.
;         Mar 05 - Mask file added by RM.
;         Jul 04 - Segmentation map file names changed from _objects.fits to
;                  _seg.fits by AR.
;         May 04 - Weight maps changed from _noise.fits to _weight.fits by RM.
;         Mar 04 - Noise map clarified by AR.
;         Apr 02 - Generalised to accept different filenames by Richard Massey.
;         Mar 02 - Adapted from RM's shex.pro by Alexandre Refregier.
;-

COMPILE_OPT idl2

;
; Decide file name to look for
;
if keyword_set(filename) then begin
  file=filename ; Don't disrupt input variables
endif else message,"You must specify an image file name!"
if keyword_set(full_path) then begin
  ; Strip file extension
  if strpos(strlowcase(file),".fit") ne -1 then $
    file=strmid(file,0,strpos(strlowcase(file),'.'))
endif else begin
  ; Add relevant path name
  file=shapelets_paths(3,silent=silent)+file
endelse
if file_test(file+".fits") then begin
  message,"Reading in image data from "+file+".fits",/info,noprint=silent
endif else message,"Cannot find image "+file+".fits!"


;
; Read fits image
;
if keyword_set(fits_read) then begin
  fits_read,file+".fits",image,header
endif else begin
  image=readfits(file+".fits",header)
endelse


;
; Determine key image properties
;
n_pixels=fix(size(image,/dimensions))
if not keyword_set(units) then begin
  units=strlowcase(strtrim(sxpar(header,"BUNIT",count=count),2))
  if count eq 0 then units="unknown"
endif
if not keyword_set(pixel_scale) then begin
  pixel_scale=float(sxpar(header,"PIXSIZE",count=count))
  if count eq 0 then begin ; Various different ways to store the pixel scale in FITS headers
    cd1_1=sxpar(header,"CD1_1",count=count1_1)*3600.
    cd1_2=sxpar(header,"CD1_2",count=count1_2)*3600.
    cd2_1=sxpar(header,"CD2_1",count=count2_1)*3600.
    cd2_2=sxpar(header,"CD2_2",count=count2_2)*3600.
    if count1_1+count1_2+count2_1+count2_2 eq 4 then begin
      pixel_scale=float(sqrt([cd1_1^2+cd1_2^2,cd2_1^2+cd2_2^2]))
      ;if stddev(pixel_scale)/mean(pixel_scale) lt 0.01 then $ ; rectangular pixels?
      pixel_scale=mean(pixel_scale)
    endif else begin
      pixel_scale=1.
    endelse
  endif
endif
if not keyword_set(photo_zp) then begin
  photo_zp=float(sxpar(header,"PHOTZP",count=count))
  if count eq 0 then photo_zp=float(sxpar(header,"PHOTOZP",count=count))
endif
if not keyword_set(exposure_time) then begin
  exposure_time=sxpar(header,"EXPTIME",count=count)
  if count eq 0 or size(exposure_time,/TYPE) eq 7 then exposure_time=1. else begin
    exposure_time=float(exposure_time)
  endelse
endif
mean_image=mean(image)
rms_image=stddev(image)


;
; Subtract constant sky background level
;
if keyword_set(sky_subtract) then begin
  ;sky_background=2.5*median(image.image)-1.5*mean(image.image)
  highbad=max(image)
  sky,image,sky_background,highbad=highbad,/SILENT
  message,"Subtracting a sky background level of "+strtrim(sky_background,2),/info,noprint=silent
  image=temporary(image)-sky_background
endif



;
; Store in image-type IDL structure
;
image={name:filename,           $
       type:"image",            $
       history:strarr(1),       $
       image:image,             $
       header:header,           $
       n_pixels:n_pixels,       $
       mean_image:mean_image,   $
       rms_image:rms_image,     $
       units:units,             $
       pixel_scale:pixel_scale, $
       photo_zp:photo_zp,       $
       exposure_time:exposure_time}


;
; Read in "segmented image" created by SExtractor
; (this flags every pixel with the SEx ID of the object it belongs to -
;  or zero, if it is part of the background.)
;
found_seg_map=0B
if not keyword_set(no_segmentation) then begin
  if file_test(file+"_seg.fits") then begin
    found_seg_map=1B
    message,"Reading in segmented pixel map",/info,noprint=silent
    fits_read,file+"_seg.fits",seg
    ; If 32 bit versions of SExtractor find too many objects, the integer ID
    ; system first jumps from positive to negative, then starts repeating
    ; values. This fix is needed to make unique ID numbers for all objects.
    if min(seg) lt 0  then begin
      nonzero=where(seg ne 0)
      n_objects=n_elements(uniq(seg[nonzero],sort(seg[nonzero])))
      if n_objects eq 65535L then begin ; Suspicious that there are exactly 2e16 objects
        message,"Unwrapping repeated integers in segmentation map",/info,noprint=silent
        seg=long(seg)
        seg_width=(size(seg,/DIMENSIONS))[0]
        negatives=fix(median(where(seg lt 0,n_negative)/seg_width))
        wraparound=where(seg[*,negatives:*] gt 0)
        seg[seg_width*negatives+wraparound]=temporary(seg[seg_width*negatives+wraparound])+65536L
      endif
      message,"Correcting negative indices in segmentation map",/info,noprint=silent
      seg=(seg lt 0)*65536L+seg
      ; Save memory
      if n_objects lt 65535L then seg=uint(seg) else seg=ulong(seg)
    endif
    image=create_struct(temporary(image),"seg",seg)
    delvarx,seg
  endif else message,"No segmentation map available!",/info,noprint=silent
endif


;
; Read in or construct noise map
;
if not keyword_set(no_noise) then begin
  if keyword_set(noise_level) then begin
    message,"Storing input noise level",/info,noprint=silent
    image=create_struct(temporary(image),"noise",1./((noise_level[0]^2)>1e-10))
  endif else if keyword_set(estimate_noise) then begin
    message,"Estimating noise level",/info,noprint=silent
    ;background=2.5*median(image.image)-1.5*mean(image.image)
    ;noise=1./(mean(abs(image.image-background)))^2
    if found_seg_map then begin
      objects=abs(image.seg)<1
      roi=where(objects,n_roi)
      if n_roi gt 0 then begin
        if not keyword_set(n_grow) then n_grow=2
        grown_objects=-objects
        for k=-n_grow,n_grow do grown_objects=+shift(objects,k)
        tr_objects=transpose(objects)
        for k=-n_grow,n_grow do grown_objects+=transpose(shift(tr_objects,k))
        for k=1,(n_grow+1)/2 do begin
     	    for i=-1,1,2 do begin
   	        for j=-1,1,2 do begin
     	        grown_objects+=transpose(shift(transpose(shift(objects,i*k)),j*k))
     	      endfor
          endfor
        endfor
        objects=grown_objects<1
        roi=where(objects,n_roi)
        delvarx,objects
        temp_image=image.image
        if n_roi gt 0 then begin
          highbad=max(temp_image)+1
          temp_image[roi]=highbad
        endif
        delvarx,roi
        sky,temp_image,mean_background,noise,highbad=highbad,/SILENT
        delvarx,temp_image
      endif else sky,image.image,mean_background,noise,/SILENT
    endif else sky,image.image,mean_background,noise,/SILENT
    noise=1./noise^2
    image=create_struct(temporary(image),"noise",noise)
    delvarx,noise
  endif else if file_test(file+"_weight.fits") then begin
    message,"Reading in inverse variance (weight) map",/info,noprint=silent
    fits_read,file+"_weight.fits",noise
    image=create_struct(temporary(image),"noise",noise)
    delvarx,noise
  endif else message,"No weight map available",/info,noprint=silent
endif


;
; Read in or construct image mask
;
if not keyword_set(no_mask) then begin
  if file_test(file+"_mask.fits") then begin
    message,"Reading in image mask",/info,noprint=silent
    fits_read,file+"_mask.fits",mask
      if size(mask,/TYPE) ne 1 then mask=byte(temporary(mask))
    image=create_struct(temporary(image),"mask",mask)
    delvarx,mask
  endif else message,"No mask available",/info,noprint=silent
endif


;
; Report our success story
;
message,"Image "+file+".fits now in memory",/info,noprint=silent

end
