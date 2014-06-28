pro shapelets_make_image_mask, filename,                  $
                               ISELECT=iselect,           $
                               SEXCAT=sexcat,             $
                               IMAGE=image,               $
                               SPIKES_ONLY=spikes_only,   $
                               STARS_SQUARE=stars_square, $
                               STARS_CIRCLE=stars_circle, $
                               STARS_CROSS=stars_cross,   $
                               PLOTIT=plotit

;+
; NAME:
;      SHAPELETS_MAKE_IMAGE_MASK
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Create a binary mask of an image that is zero everywhere, except for
;      ones around saturated or potentially corrupted regions.
;
; INPUTS:
;      FILENAME     - The name of the image file, without its extension.
;
; OPTIONAL INPUTS:
;      ISELECT      - A structure containing information about the saturation
;                     level of an image, and a list of any saturated objects.
;      SEXCAT       - A Source Extractor catalogue of the image, in shapelets
;                     structure format.
;      IMAGE        - The image and segmentation map in a shapelets structure.
;
; KEYWORD PARAMETERS:
;      SPIKES_ONLY  - Masks the stars using only extended vertical rectangles.
;      STARS_CIRCLE - Masks the stars using circles as well [DEFAULT].
;      STARS_SQUARE - Masks the stars using squares as well as the spikes.
;      STARS_CROSS  - Masks the stars using a + pattern as well as the spikes.
;      PLOTIT       - Plot a pretty picture to STDOUT.
;
; OUTPUTS:
;      Writes an image mask to disk, and incorporates it into the input
;      image structure if there was one.
;
; MODIFICATION HISTORY:
;      Aug 05 - Cleaned up and generalised by Richard Massey.
;      May 05 - Written by Will High.
;-

; Read in data
if not keyword_set(sexcat) then shapelets_read_sexcat,sexcat,filename
if not keyword_set(iselect) then begin
  iselect_filename=shapelets_paths(3)+filename+".iselect"
  if not file_test(iselect_filename) then $
    shapelets_interactive_select_stars, sexcat, iselect
  restore,iselect_filename,/verbose
endif
if not keyword_set(image) then shapelets_read_image,image,filename
if n_elements(image.seg) eq 1 and image.seg[0] eq 0 then message,"Segmentaion map not available!"


; Create a blank mask, and index its pixels
mask=bytarr(image.n_pixels[0],image.n_pixels[1])
x=rebin(indgen(image.n_pixels[0],1),image.n_pixels[0],image.n_pixels[1],/sample)
y=rebin(indgen(1,image.n_pixels[1]),image.n_pixels[0],image.n_pixels[1],/sample)

; Loop over each saturated object
for j=0L,iselect.n_saturated-1 do begin
  
  ; Obtain object indices (segmentation map integers are one higher)
  obj=iselect.saturated[j] 
  segobj=iselect.saturated[j]+1
  ;if (j mod 100 eq 0) then $
    print,strtrim(j,2)+'/'+strtrim(iselect.n_saturated,2)

  ; Get relevant quantities into local variables
  xmin=sexcat.xmin[obj,0]>0
  xmax=sexcat.xmax[obj,0]<image.n_pixels[0]-1
  ymin=sexcat.xmin[obj,1]>0
  ymax=sexcat.xmax[obj,1]<image.n_pixels[1]-1

  ; Get the length of the postage stamp's sides, and the centre of the object
  xside=xmax-xmin
  yside=ymax-ymin
  cent=round(sexcat.x[obj,*])

  ; Mask the spike.
  xside=float(xside)
  yside=float(yside)
  if ((abs(yside-xside)/xside gt 0.1) and (yside gt xside)) then begin
    
    ; It's probably a saturated object with a large spike, so mask it

    ; Cut out a postage stamp from the larger image and segmentaion arrays
    pstamp=image.image[xmin:xmax,ymin:ymax]
    seg=image.seg[xmin:xmax,ymin:ymax]
    
    ; Mask out other objects using the segmentation map
    seg[where(seg ne segobj)]=0
    if (where(seg ne 0))[0] eq -1 then begin
    	message,"Object #"+strtrim(segobj,2)+" not in segmentation map!"
    endif else seg[where(seg ne 0)]=1
    pstamp=pstamp*seg
    
    ; Get the index arrays for postage stamps
    centpstamp=cent-[xmin,ymin]
    xpstamp=rebin(indgen(xside+1,1),xside+1,yside+1,/sample)-centpstamp[0]
    ypstamp=rebin(indgen(1,yside+1),xside+1,yside+1,/sample)-centpstamp[1]
    
    ; Find the widest row
    ixx=fltarr(yside+1)
    for momind=0,yside do begin
    	flux=total(pstamp[*,momind])
    	if where(flux eq 0) ne -1 then flux[where(flux eq 0)]=1
    	ixx[momind]=total(pstamp[*,momind]*abs(xpstamp[*,momind]))
    endfor
    
    ; Find the tallest column (not useful right now...)
    iyy=fltarr(xside+1)
    for momind=0,xside do begin
    	flux=total(pstamp[momind,*])
    	if where(flux eq 0) ne -1 then flux[where(flux eq 0)]=1
    	iyy[momind]=total(pstamp[momind,*]*abs(ypstamp[momind,*]))
    endfor
    widest_ind=where(ixx eq max(ixx))
    tallest_ind=where(iyy eq max(iyy))

    ; Check for errors
    if n_elements(widest_ind) ne 1 or $
      n_elements(tallest_ind) ne 1 then begin
    	message,"Multiple peaks in moments"
    endif
    widest_ind=widest_ind[0]
    tallest_ind=tallest_ind[0]
    if widest_ind eq -1 or tallest_ind eq -1 then message,"Moment measurement failed"
    
    ; Plot the pstampage and the moments with centroids
    if keyword_set(plotit) then begin
      pmulti=!p.multi
      !p.multi=[0,1,3]
      pstamp[centpstamp[0],*]=min(pstamp)
      pstamp[*,centpstamp[1]]=min(pstamp)
      pstamp[*,widest_ind]=max(pstamp)
      if centpstamp[1] eq widest_ind then $
    	pstamp[*,widest_ind]=10^mean(alog10(pstamp > 10))
      shapelets_plot_image,alog10(pstamp > 10)
      plot,ixx,charsize=1.5
      oplot,[1,1]*centpstamp[1],[-1,1]*1e10
      oplot,[1,1]*widest_ind,[-1,1]*1e10,line=2
      plot,iyy,charsize=1.5
      oplot,[1,1]*centpstamp[0],[-1,1]*1e10
      wait,1
      !p.multi=pmulti
    endif

    ; Make the spike mask
    cent=[sexcat.x[obj,0],widest_ind+ymin]
    cent=round(cent)
    half_longside=round(yside/2.+$
    			abs(sexcat.x[obj,1]-cent[1]))
    half_longside=round(1.3*half_longside)
    
    ; Make the width of the spike a tenth of the xside
    half_shortside=round(xside/10.)
    bottom=cent[1]-half_longside>0
    top=cent[1]+half_longside<image.n_pixels[1]-1
    left=cent[0]-half_shortside>0
    right=cent[0]+half_shortside<image.n_pixels[0]-1
    mask[left:right,bottom:top]=1B
  endif else if keyword_set(spikes_only) then continue

  ; Mask the bulky center
  if keyword_set(stars_square) then begin
    
    ; Mask the stars using a rectangle
    halfside=round(xside/2.)
    bottom=(cent[1]-halfside)>0
    top=(cent[1]+halfside)<(image.n_pixels[1]-1)
    left=(cent[0]-halfside)>0
    right=(cent[0]+halfside)<(image.n_pixels[0]-1)
    mask[left:right,bottom:top]=1B

  endif else if keyword_set(stars_cross) then begin

    ; Mask the stars using a cross shape
    message,"Not yet written!"

  endif else begin

    ; Mask the stars using a circle
    radius=round(xside/2.)
    diameter=2*radius

    ; Make a little mask
    bottom=cent[1]-radius
    top=cent[1]+radius
    left=cent[0]-radius
    right=cent[0]+radius
    maskpstamp=mask[left>0:right<(image.n_pixels[0]-1),bottom>0:top<(image.n_pixels[1]-1)]

    ; Cut out circle
    xmaskpstamp=rebin(indgen(diameter+1,1),diameter+1,diameter+1,/sample)
    ymaskpstamp=rebin(indgen(1,diameter+1),diameter+1,diameter+1,/sample)
    centmaskpstamp=[1,1]*radius+$
        	 [left<0,bottom<0]
    here=where((xmaskpstamp-centmaskpstamp[0])^2+$
               (ymaskpstamp-centmaskpstamp[1])^2 lt radius^2)
    maskpstamp[here]=1B

    ; Put the little mask in the big mask
    mask[left>0:right<(image.n_pixels[0]-1),bottom>0:top<(image.n_pixels[1]-1)]=maskpstamp

  endelse

endfor

; Store mask in the image structure, and write it to disk
if tag_exist(image,"mask") then begin
  new_image={name:image.name,type:image.type}
  for i=0,n_tags(image)-1 do begin
    tag_name=strupcase((tag_names(image))[i]) 
    if tag_name ne "NAME" and tag_name ne "TYPE" and tag_name ne "MASK" then $
      new_image=create_struct(temporary(new_image),tag_name,image.(i))
  endfor
  image=new_image
endif
image=create_struct(temporary(image),"mask",mask)
shapelets_write,image,filename,/MASK

end
