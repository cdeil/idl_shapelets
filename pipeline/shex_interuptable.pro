;$Id: shex.pro, v2$
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
; ***********************************************************************
; ***********************************************************************
;
;+
; NAME:
;        SHEX
;
; PURPOSE:
;        Having located objects within an image using SExtractor, this routine
;        now decomposes all of the objects into shapelets. It uses the
;        shapelets focus suite of routines to optimise the nmax, beta and
;        centroid parameters. All of the decompositions are written to disc
;        in a shapelets catalgoue. This can be later read to memory using
;        shapelets_read_sexcat.pro
;
; CATEGORY:
;        Shapelets/ROI modelling.
;
; CALLING PROCEDURE:
;        shex,'example',sexcat=sexcat,n_max=15,/plot
;
; INPUTS:
;        FILENAME_INPUT  - file name of input image.
;
; OPTIONAL INPUTS:
;        FILENAME_OUTPUT - File name of output catalogue. Default is the
;                          same as filename_input, but with .shapecat extension.
;        IMAGE           - Can set in advance so it doesn't have to be loaded in.
;        SEXCAT          - Can set in advance so it doesn't have to be loaded in.
;        INDEX    	 - Vector of SExtractor index number to process 
;                 	   (Default: all of the SExtractor catalogue is used). Note
;                 	   that the SExtractor index is equal to sextractor ID-1
;        PSF      	 - Shapelet decompstructure containing PSF to deconvolve from.
;                 	   Currently, this has to be the same everywhere on an image.
;        N_MIN           - Minimum value of n_max=n1+n2 allowed for decomposition.
;        N_MAX           - Highest order of n_max=n1+n2 allowed for decomposition.
;        BETA_TOLERANCE  - Accuracy with which the optimum beta is obtained.
;        FIXED_BETA      - Force this value of beta for all objects. Still iterates
;                          to find the best centroid and n_max.
;	 CHISQ_TARGET    - Ideal value of reduced chi^2 for the residual image.
;        CHISQ_TOLERANCE - Acceptable accuracy for chi^2.
;	 CHISQ_FLATNESS  - Minimum difference in chi^2 between two decompositions
;                          with n_max differing by two to trigger the flatness
;                          constraint in shapelets_focus_nmax
;        THETA_MIN_GEOM  - Minimum scale on which it is possible for the image to
;                          contain data.
;        FULL_FOCUS      - Record every single attempt made at decomposition during
;                          iteration. Makes things a little slower, but plots nicer. 
;        NFWHM    	 - Size of postage stamp in units of a_sextractor
;        SG_CUT   	 - Threshold level for star/galaxy separation, using the 
;                 	   SExtractor stellar clasification from {0,1}. If this is
;                 	   positive, only galaxies with SExtractor values lower
;                 	   than sg_cut are decomposed. If this is negative, only
;                 	   stars with SExtractor values greater than -sg_cut are 
;                 	   decomposed. If not set or zero, everything is decomposed.
;        TRIM_FAILURES   - Discard objects with bad postage stamps or failed
;                          iterations to save space in the catalogue. (Default: keep
;                          everything).
;        TOO_BIG  	 - Maximum radius of postage stamp in pixels (objects that
;                 	   look as if they are bigger than this get skipped). Can set
;                 	   so that IDL doesn't run out of memory.
;        SKY      	 - Allow sky background subtraction by fitting during image
;                 	   decomposition.
;                          1: locally-determined constant, 
;                 	   2: locally-determined plane
;                          3: global constant
;                          4: 3+1 &c.
;        TRUE_IMAGE      - Additional input image for plotting purposes only (e.g. 
;                 	   noise-free simulated image, or image in a second colour)
;        SCAT_IN         - Input scat catalogue for plotting purposes only
;                          (e.g. noise free input catalogue from simulations)
;
; KEYWORD PARAMETERS:
;        SILENT     	 - Operates silently.
;        VERBOSE    	 - Operates noisily (VERBOSE=2 to operate very noisily).
;        PLOTIT     	 - Display pretty pictures to screen.
;        SEG_MAP    	 - Use locally-determined segmentation map, rather than that
;                   	   supplied in the image structure.
;        NOISE_MAP  	 - Use the noise map supplied in the postage stamp, rather
;                   	   than local noise esimation.
;        FULL_FOCUS      - Record every single attempt made at decomposition during
;                          iteration. Makes things a little slower, but plots nicer. 
;        SQUARE     	 - Use a square postage stamp (default is circular).
;        NEIGHBOUR  	 - Treatment of neighbours: 0(default): infinite
;                   	   errors, i.e. unconstrained fit, 1: set neighbour
;                   	   pixels to the background level and associated errors.
;        NON1       	 - Force shapelet coefficients with n=1 to be zero, as 
;                   	   suggested by Kuijken.
;
; OUTPUTS:
;        Shapelet catalogue written to disc, with the same file name and path
;        as the input image, but with the extension .shapecat.
;
; OPTIONAL OUTPUTS:
;        SHAPECAT  	 - Shapelet catalogue data structure stored to memory.
;                  	   View contents by typing "help,shapecat,/structure"
;        IMAGE     	 - Image data structure, if none was entered originally.
;        SEXCAT    	 - SExtractor catalogue data structure, if none was entered.
;
; DATA REQUIRED:
;        Needs an image and a SExtractor catalogue to initially locate the 
;        objects for decomposition.
;
; PROCEDURES USED:
;        This combines most of the decomposition and focus routines in the
;        shapelets code release.
;
; NOTES:
;        The first of the two flags for each object in a shapelet catalogue describes
;        its segmentation into a small postage stamp. The second describes its
;        decomposition into shapelets. Their precise meanings are given below, but see
;        the IDL headers in shapelets_image2pstamp.pro and shapelets_focus.pro
;        for further information.
;
; MODIFICATION HISTORY:
;        Oct 05 - Creation of temporary shapecat by JB
;        Jul 05 - Significant load of bug fixes and clean up by RM.
;        Apr 05 - NON1 keyword added by RM.
;        Mar 05 - image_in input changed to true_image by RM to avoid conflict 
;        Feb 05 - Framework for variable PSFs added by Joel Berge 
;        Jan 05 - Neighbour and nfwhm keywords added by AR & RM
;        Aug 04 - Beta and filename_output keywords added by RM
;        Jul 04 - index, image_in and scat_in keywords added by AR
;        Jul 04 - Too_big paremeter added by RM
;        Jun 04 - Long/integer fixed to allow images with >32768 objects by RM
;                 note: will run into problems with >32768 objects if the
;                       (single precision) SExtractor segmentation map is used
;        May 04 - Sky background subtraction enabled by RM
;        May 04 - Catalogues written to disc as IDL save files by RM
;        Apr 04 - Noise_map and seg_map keywords added by Alexandre Refregier
;        Oct 03 - Bug fix on centroiding by RM
;        Sep 03 - Converted for use on non-UNIX platforms by RM
;        Mar 03 - Significantly tidied up by RM
;        Aug 01 - Written by Richard Massey
;-

pro shex, FILENAME_INPUT,                  $ 
          FILENAME_OUTPUT=filename_output, $
          FULL_PATH=full_path,             $
	  INDEX=index,                     $
          IMAGE=image,                     $
	  SEXCAT=sexcat,                   $
          SHAPECAT=shapecat,               $
	  BETA_TOLERANCE=beta_tolerance,   $
	  FIXED_BETA=fixed_beta,           $
          CHISQ_TARGET=chisq_target,	   $
	  CHISQ_TOLERANCE=chisq_tolerance, $
          CHISQ_FLATNESS=chisq_flatness,   $
	  THETA_MIN_GEOM=theta_min_geom,   $
	  FULL_FOCUS=full_focus,           $
	  N_MIN=n_min,                     $
	  N_MAX=n_max,                     $
	  NOISE_MAP=noise_map,             $
	  SEG_MAP=seg_map,                 $
	  TOO_BIG=too_big,                 $
	  SG_CUT=sg_cut,                   $
          SQUARE=square,                   $
          NEIGHBOUR=neighbour,             $
          NFWHM=nfwhm,                     $
	  TRIM_FAILURES=trim_failures,     $
	  PSF=psf,                         $
	  SKY=sky,                         $
	  NON1=non1,                       $
	  TRUE_IMAGE=true_image,           $
	  SCAT_IN=scat_in,                 $
	  PLOTIT=plotit,                   $
          SILENT=silent,                   $
	  VERBOSE=verbose

COMPILE_OPT idl2
;ON_ERROR,2

; Initialise default parameters
if not keyword_set(verbose) then verbose=0B         ; How much output to screen?
if keyword_set(silent) then begin & plotit=0 & verbose=0 & endif
if not keyword_set(n_min) then n_min=2              ; Allowed range of n_max
if not keyword_set(n_max) then n_max=20             ;
if not keyword_set(sky)   then sky=0                ; Sky subtraction options
if not keyword_set(filename_output) then filename_output=filename_input ; File name of output catalogue
if not keyword_set(too_big) then too_big=100        ; Skip objects bigger than this radius [pixels]
if keyword_set(sg_cut) then begin                   ; Star/galaxy classification cut
  case sign(sg_cut) of                              ; Write initial message to screen
    1: message,"Decomposing galaxies in "+filename_input+".fits",/info,noprint=silent
   -1: message,"Decomposing stars in "+filename_input+".fits",/info,noprint=silent
  endcase
endif else message,"Decomposing objects in "+filename_input+".fits",/info,noprint=silent

; Read in SExtractor catalogue of image
if not keyword_set(sexcat) then $
  shapelets_read_sexcat,sexcat,filename_input,silent=silent,full_path=full_path
sexcat_tagnames=tag_names(sexcat)
n_sexcat_tagnames=n_elements(sexcat_tagnames)
n_sexcat_fields=0L
for i=0,n_sexcat_tagnames-1 do n_sexcat_fields=n_sexcat_fields+size(sexcat.(i),/N_DIMENSIONS)


; Read in image(s)
if not keyword_set(image) then $
  shapelets_read_image,image,filename_input,silent=silent,full_path=full_path,sky_subtract=(sky ge 3)
sky=sky mod 3

; Select and reorder the objects as requested
if keyword_set(index) then begin
  n=long(n_elements(index))
  reorder=long(index)
endif else begin
  n       = long(n_elements(sexcat.id)) ; Number of objects
  ;reorder = reverse(sort(sexcat.flux)) ; Look at pretty objects first...
  ;reorder = sort(sexcat.flux)          ; ...or low S/N ones first
  reorder = indgen(n,/long)             ; ...or in SExtractor order
  ;reorder = reverse(indgen(n,/long))   ; ...or in reverse SExtractor order
endelse  


; Determine seeing, and alter the filename if we are deconvolving from a PSF
if keyword_set(psf) then begin
  if not shapelets_structure_type(psf,message=message) then message,"PSF type not recognised"
  case psf.type of
    "decomp": seeing=sqrt(shapelets_rsquared(psf)) ; Expected for the case of a constant PSF
    "shapecat":  seeing=mean(sqrt(psf.rsquared)) ; Expected for Joel-style pre-done interpolation
    "image": message,"Please decompose the PSF into shapelets first"
    "interpolation": message,"Reserved for futute use"
    "decomp_polar": message,"Reserved for future use"
    else: message,"PSF structure type not recognised"
  endcase
  filename_output=filename_output+"_deconv"
endif else begin
  if ( max(sexcat.class) gt 0.5 ) then $
       seeing=median(sexcat.fwhm[where(sexcat.class gt 0.5)]) $
  else seeing=median(sexcat.fwhm)
endelse


; Prepare empty arrays to contain output catalogue
shapelets_make_nvec,n_max,n1,n2,n_coeffs
; if no temporary shapecat already exists  ;######################
if not file_test(shapelets_paths(2)+filename_output+'_tmp.shapecat') then begin   ;######################
    x=fltarr(1,2) 
    beta=0.
    n_max_cat=0L
    coeffs=fltarr(1,n_coeffs)
    coeffs_error=fltarr(1,n_coeffs)
    flag=bytarr(1,2)
    chisq=0.
    objx_min=0L
    objx_max=0L
    objy_min=0L
    objy_max=0L
    back_mean_local=0.
    back_rms_local=0.
    back_mean_ext=0.
    back_rms_ext=0.
    chisq_rms=0.
    sexid=0L
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
    sexarea=0L
    sexra=0.
    sexdec=0.
endif else begin ;but if a temporary shapecat exists, initialise arrays from temporary ones ;######################
    shapelets_read_shapecat,tmp_shcat,filename_output+'_tmp.shapecat',/moments,/silent ;######################
    x=[[[0.],[0.]],tmp_shcat.x] ;######################
    beta=[0.,tmp_shcat.beta] ;######################
    n_max_cat=[0.,tmp_shcat.n_max] ;######################
    coeffs=[transpose(replicate(0.,n_coeffs)),tmp_shcat.coeffs] ;######################
    coeffs_error=[transpose(replicate(0.,n_coeffs)),tmp_shcat.coeffs_error] ;######################
    flag=[[[0.],[0.]],tmp_shcat.flag] ;######################
    chisq=[0.,tmp_shcat.chisq] ;######################
    objx_min=[0.,tmp_shcat.objx_min] ;######################
    objx_max=[0.,tmp_shcat.objx_max] ;######################
    objy_min=[0.,tmp_shcat.objy_min] ;######################
    objy_max=[0.,tmp_shcat.objy_max] ;######################
    back_mean_local=[0.,tmp_shcat.back_mean_local] ;######################
    back_rms_local=[0.,tmp_shcat.back_rms_local] ;######################
    back_mean_ext=[0.,tmp_shcat.back_mean_ext] ;######################
    back_rms_ext=[0.,tmp_shcat.back_rms_ext] ;######################
    chisq_rms=[0.,tmp_shcat.chisq_rms] ;######################
    sexid=[0,tmp_shcat.sexid] ;######################
    sexfwhm=0.  ;all these SExtractor parameters are put at the end, and not in the temporary shapecat ;######################
    sexflux=0.  ;so they can be initialised to 0 ;######################
    sexmag=0. ;######################
    sexclass=0. ;######################
    sexa=0. ;######################
    sexb=0. ;######################
    sextheta=0. ;######################
    sexe1=0. ;######################
    sexe2=0. ;######################
    sexx=0. ;######################
    sexy=0. ;######################
    sexarea=0L ;######################
    sexra=0. ;######################
    sexdec=0. ;######################
endelse ;######################

; Loop over all objects in SExtractor catalogue, attempting to decompose them
; into shapelets
n_tried  = 0L
n_done   = 0L
n_stored = 0L
for i=0L,n-1L do begin

    if (reorder[i] le max(sexid) and reorder[i] ne 0) then begin ;######################
        message,'Object #'+strtrim(reorder[i],1)+' has already been decomposed in temporary shapecat',/info ;######################
        continue ;######################
    endif ;######################

  ; Reset flags
  skip=0B
  if verbose ge 1 then print

  ; Skip stars/galaxies as required
  if keyword_set(sg_cut) then begin
    if (sg_cut gt 0 and sexcat.class[reorder[i]] gt sg_cut) then begin
      if keyword_set(verbose) then $
   	message,'Skipping star #'+strtrim(reorder[i],1)+', '+$
   	strtrim(i+1,1)+'/'+strtrim(n,1),/info,noprint=silent
print,'SHEX : Skipping star #'+strtrim(reorder[i],1)+', '+$
   	strtrim(i+1,1)+'/'+strtrim(n,1) ;needed for output .txt file when using nohup
      skip=1B
    endif else if (sg_cut lt 0 and sexcat.class[reorder[i]] lt abs(sg_cut)) then begin
      if keyword_set(verbose) then $
   	message,'Skipping galaxy #'+strtrim(reorder[i],1)+', '+$
   	strtrim(i+1,1)+'/'+strtrim(n,1),/info,noprint=silent
print,'SHEX : Skipping galaxy #'+strtrim(reorder[i],1)+', '+$
   	strtrim(i+1,1)+'/'+strtrim(n,1)
        skip=1B
    endif
  endif


  ; Skip really big objects to avoid running out of memory and crashing!
  if (sexcat.a[reorder[i]] gt too_big) then begin
    if keyword_set(verbose) then $
      message,'Skipping bloody great big object #'+strtrim(reorder[i],1)+', '+$
      strtrim(i+1,1)+'/'+strtrim(n,1),/info,noprint=silent
print,'SHEX : Skipping bloody great big object #'+strtrim(reorder[i],1)+', '+$
      strtrim(i+1,1)+'/'+strtrim(n,1)
    skip=1B
  endif
  
  if skip then begin
  
    ; Fake a pstamp structure
    pstamp={im_ran:fltarr(4),$
            back_mean_local:0.,back_rms_local:0.,$
	    back_mean_ext:0.,back_rms_ext:0.,$
            flag:10}

    ; Fake a decomp structure
    decomp=shapelets_create_decomp(n_max)
    decomp=create_struct(decomp,"flag",10)
   endif else begin
    ; Write progress report to screen
    message,'Decomposing object #'+strtrim(reorder[i],1)+', '+$
      strtrim(i+1,1)+'/'+strtrim(n,1)+$
      ' at pixel ('+strtrim(round(sexcat.x[reorder[i],0]),1)+ $
      ','+strtrim(round(sexcat.x[reorder[i],1]),1)+')' $
      ,/info,noprint=silent
print,'SHEX : Decomposing object #'+strtrim(reorder[i],1)+', '+$
      strtrim(i+1,1)+'/'+strtrim(n,1)+$
      ' at pixel ('+strtrim(round(sexcat.x[reorder[i],0]),1)+ $
      ','+strtrim(round(sexcat.x[reorder[i],1]),1)+')'
  
  
    ; Extract postage stamp from larger image
    pstamp=shapelets_sexcat2pstamp(image,sexcat,reorder[i],silent=(1-verbose)>0,$
   				   seg_map=seg_map,noise_map=noise_map,square=square,$
   				   nfwhm=nfwhm,neighbour=1)

  
    if pstamp.flag ge 3 then begin
  
      ; Fake a decomp structure
      decomp=shapelets_create_decomp(n_max)
      focus={flag:10}
      
      ; Skip objects near edge of image, etc.
      message,'Skipping object with postage stamp flag='+(strtrim(pstamp.flag,2)),/info,noprint=silent
print,'SHEX : Skipping object with postage stamp flag='+(strtrim(pstamp.flag,2))
      skip=1B

    endif else begin

      ; Increment counter
      n_tried+=1L

      ; Set up PSF into required form to pass to decomposition routines
      if keyword_set(psf) then begin
   	case psf.type of
   	  "decomp": psf_obj=psf ; Expected for the case of a constant PSF
   	  "shapecat": begin		  ; Expected for pre-done interpolation
   			psf_matched=0B
   			if psf.n ge reorder[i]+1 then begin
   			  ; Find matching PSF
   			  psf_obj=shapelets_shapecat2decomp(psf,reorder[i])
   			  ; Check that it really is in the same place as the galaxy
   			  psf_gal_dist=sqrt((psf.x[reorder[i],0]-pstamp.x)^2+$
   					    (psf.x[reorder[i],1]-pstamp.y)^2)
   			  if psf_gal_dist lt 1 then psf_matched=1B else $
   			    message,"WARNING: Interpolated PSF does not match galaxy catalogue!",/info
   			endif
   			if not psf_matched then begin
   			  ; Find nearest PSF
   			  psf_gal_dists=(psf.x[*,0]-pstamp.x)^2+(psf.x[*,1]-pstamp.y)^2
   			  psf_gal_dist=sqrt(min(psf_gal_dists,closest_psf))
   			  message,"Closest PSF at a distance of "+strtrim(string(psf_gal_dist),2)+" pixels",/info
   			  psf_obj=shapelets_shapecat2decomp(psf,closest_psf)
   			endif
   		      end
   	  else: message,"PSF structure type not recognised"
   	endcase
      endif
      
      
      ; Optimise ("focus") shapelet decomposition parameters n_max, beta and x_c:
      decomp=shapelets_focus(PSTAMP,			      $
   			     FOCUS=focus,		      $
   			     RECOMP=recomp,		      $
   			     psf=psf_obj,		      $
   			     BETA_TOLERANCE=beta_tolerance,   $
   			     N_MAX_RANGE=[n_min,n_max],       $
   			     CHISQ_TARGET=chisq_target,       $
   			     CHISQ_TOLERANCE=chisq_tolerance, $
   			     CHISQ_FLATNESS=chisq_flatness,   $
   			     THETA_MIN_GEOM=theta_min_geom,   $
   			     FULL_FOCUS=full_focus,	      $
   			     MAX_LOOPS=max_loops,	      $
   			     FIXED_BETA=fixed_beta,	      $
   			     VERBOSE=(verbose ge 2),	      $
   			     SILENT=silent,		      $
   			     sky=sky,			      $
   			     non1=non1)
      x0=focus.x+[pstamp.x,pstamp.y]-[pstamp.xo,pstamp.yo]


      ; Check again that all is well 
      if shapelets_flux(decomp) le 0 then begin
   	if shapelets_flux(decomp) eq 0 then message,"Model has zero flux",/info,noprint=silent else $
   	message,"Model has negative flux",/info,noprint=silent
        focus.flag=10
      endif
      if focus.flag ge 10 then skip=1B
      
      ; Plot to screen if required
      if not skip and keyword_set(plotit) then begin
   	shapelets_plot,focus,pstamp,decomp=decomp,recomp=recomp,image_in=true_image,scat_in=scat_in,psf=psf_obj,noise_map=noise_map
   	if i ne n-1 then read,junk
      endif
    
    endelse

  endelse 
  
  ; Append individual coefficient array to storage array
  if skip then begin
    if keyword_set(trim_failures) then continue
    coeffs=[coeffs,fltarr(1,n_coeffs)]
    coeffs_error=[coeffs_error,fltarr(1,n_coeffs)]
    x=[x,fltarr(1,2)]
    beta=[beta,0.]
    n_max_cat=[n_max_cat,0]
    chisq=[chisq,!values.f_infinity]
    chisq_rms=[chisq_rms,0.]
  endif else begin
    ; Pad object with zero coefficients to fit in the catalogue
    shapelets_extend_nmax,decomp,n_max-decomp.n_max
    n_done+=1L
    coeffs=[coeffs,transpose(decomp.coeffs)]
    coeffs_error=[coeffs_error,transpose(decomp.coeffs_error)]
    x=[x,transpose(x0)]
    beta=[beta,decomp.beta]
    n_max_cat=[n_max_cat,focus.n_max]
    chisq=[chisq,decomp.chisq[1]]
    chisq_rms=[chisq_rms,focus.chisq_abs_target-focus.chisq_target]
    if n_done eq 1 then begin
      ; Record strings that explain the various flags from the first successful decomposition
      message,"Recording flag definitions",/info,noprint=silent
      n_flags=n_elements(pstamp.flag_interpret)>n_elements(pstamp.flag_interpret_mini)>$
 	      n_elements(focus.flag_interpret)>n_elements(focus.flag_interpret_mini)
      flag_interpret=strarr(n_flags,2)
      flag_interpret[0:n_elements(pstamp.flag_interpret)-1,0]=pstamp.flag_interpret
      flag_interpret[0:n_elements(focus.flag_interpret)-1,1]=focus.flag_interpret
      flag_interpret_mini=strarr(n_flags,2)
      flag_interpret_mini[0:n_elements(pstamp.flag_interpret_mini)-1,0]=pstamp.flag_interpret_mini
      flag_interpret_mini[0:n_elements(focus.flag_interpret_mini)-1,1]=focus.flag_interpret_mini
      chisq_target=focus.chisq_target
    endif
  endelse
  flag=[flag,[[pstamp.flag],[focus.flag]]]
  objx_min=[objx_min,pstamp.im_ran[0]]
  objx_max=[objx_max,pstamp.im_ran[1]]
  objy_min=[objy_min,pstamp.im_ran[2]]
  objy_max=[objy_max,pstamp.im_ran[3]]
  back_mean_local=[back_mean_local,pstamp.back_mean_local]
  back_rms_local=[back_rms_local,pstamp.back_rms_local]  
  back_mean_ext=[back_mean_ext,pstamp.back_mean_ext]
  back_rms_ext=[back_rms_ext,pstamp.back_rms_ext]
  for j=0,n_sexcat_tagnames-1 do begin
    sexcat_tagname=strupcase(sexcat_tagnames[j])
    case sexcat_tagname of
      "NAME":
      "TYPE":
      "SEEING":
      "X": begin
             sexx=[sexx,sexcat.x[reorder[i],0]]
             sexy=[sexy,sexcat.x[reorder[i],1]]
           end
      "ID": sexid=[sexid,sexcat.id[reorder[i]]]
      "FWHM": sexfwhm=[sexfwhm,sexcat.fwhm[reorder[i]]]
      "FLUX": sexflux=[sexflux,sexcat.flux[reorder[i]]]
      "MAG": sexmag=[sexmag,sexcat.mag[reorder[i]]]
      "CLASS": sexclass=[sexclass,sexcat.class[reorder[i]]]
      "A": sexa=[sexa,sexcat.a[reorder[i]]]
      "B": sexb=[sexb,sexcat.b[reorder[i]]]
      "THETA": sextheta=[sextheta,sexcat.theta[reorder[i]]]
      "E1": sexe1=[sexe1,sexcat.e1[reorder[i]]]
      "E2": sexe2=[sexe2,sexcat.e2[reorder[i]]]
      "AREA": sexarea=[sexarea,sexcat.area[reorder[i]]]
      "RA": sexra=[sexra,sexcat.ra[reorder[i]]]
      "DEC": sexdec=[sexdec,sexcat.dec[reorder[i]]]
      else:
    endcase
  endfor

  ; Increment counter
  message,'Adding object to catalogue',/info,noprint=(1-verbose)>0
  n_stored+=1L
  
;create and save temporary shapecat, growing loop after loop (saved
;every 1000 objects)
  if (i ne 0 and i mod 1000 eq 0) then begin ;######################
      if n_stored eq 0 then message,"WARNING: No objects successfully decomposed into shapelets!",/info else begin ;######################

  ; Create structure ;######################
          tmp_shapecat={name:filename_output,                    $ ;######################
                        type:"shapecat",                         $ ;######################
                        n:n_stored,                              $ ;######################
                        maxn_max:n_max,                          $ ;######################
                        n_coeffs:n_coeffs,                       $ ;######################
                        polar:0B,                                $ ;######################
                        x:x[1:*,*],                              $ ;######################
                        beta:beta[1:*],                          $ ;######################
                        n_max:n_max_cat[1:*],                    $ ;######################
                        coeffs:coeffs[1:*,*],                    $ ;######################
                        coeffs_error:coeffs_error[1:*,*],        $ ;######################
                        flag:flag[1:*,*],                        $ ;######################
                        flag_interpret:flag_interpret,           $ ;######################
                        flag_interpret_mini:flag_interpret_mini, $ ;######################
                        chisq:chisq[1:*],                        $ ;######################
                        objx_min:objx_min[1:*,*],                $ ;######################
                        objx_max:objx_max[1:*,*],                $ ;######################
                        objy_min:objy_min[1:*,*],                $ ;######################
                        objy_max:objy_max[1:*,*],                $ ;######################
                        back_mean_local:back_mean_local[1:*],    $ ;######################
                        back_rms_local:back_rms_local[1:*],      $ ;######################
                        back_mean_ext:back_mean_ext[1:*],        $ ;######################
                        back_rms_ext:back_rms_ext[1:*],          $ ;######################
                        chisq_target:chisq_target,               $ ;######################
                        chisq_rms:chisq_rms[1:*],                $ ;######################
                        pixel_scale:image.pixel_scale,           $ ;######################
                        photo_zp:image.photo_zp,                 $ ;######################
                        seeing:seeing,                           $ ;######################
                        sexid:sexid[1:*,*]} ;######################
      
; Write to disc ;######################
          shapelets_write, tmp_shapecat, filename_output+'_tmp', full_path=full_path ;######################

      endelse ;no need to write SExtractor information before the whole process has run ;######################
; => need sexid !!! ;######################

      if not keyword_set(silent) then message,'Temporary shapecat written to disk',/continue ;######################
  endif ;######################

endfor


; Write catalogue to disc
if not keyword_set(silent) then print
if n_stored eq 0 then message,"WARNING: No objects successfully decomposed into shapelets!",/info else begin

  ; Create structure
  shapecat={name:filename_output,                    $
            type:"shapecat",                         $
            n:n_stored,                              $
            maxn_max:n_max,                          $
            n_coeffs:n_coeffs,                       $
            polar:0B,                                $
            x:x[1:*,*],                              $
            beta:beta[1:*],                          $
            n_max:n_max_cat[1:*],                    $
            coeffs:coeffs[1:*,*],                    $
            coeffs_error:coeffs_error[1:*,*],        $
            flag:flag[1:*,*],                        $
            flag_interpret:flag_interpret,           $
            flag_interpret_mini:flag_interpret_mini, $
            chisq:chisq[1:*],                        $
            back_mean_local:back_mean_local[1:*],    $
            back_rms_local:back_rms_local[1:*],      $
            back_mean_ext:back_mean_ext[1:*],        $
            back_rms_ext:back_rms_ext[1:*],          $
            chisq_target:chisq_target,               $
            chisq_rms:chisq_rms[1:*],                $
	    pixel_scale:image.pixel_scale,           $
	    photo_zp:image.photo_zp,                 $
	    seeing:seeing}

  ; Append information from SExtractor catalogue
  shapecat=create_struct(temporary(shapecat),"sextractor",1B)
  if n_elements(sexid)    gt 1 then shapecat=create_struct(temporary(shapecat),"sexid",sexid[1:*])
  if n_elements(sexx)     gt 1 then shapecat=create_struct(temporary(shapecat),"sexx",sexx[1:*])
  if n_elements(sexy)     gt 1 then shapecat=create_struct(temporary(shapecat),"sexy",sexy[1:*])
  if n_elements(sexra)    gt 1 then shapecat=create_struct(temporary(shapecat),"sexra",sexra[1:*])
  if n_elements(sexdec)   gt 1 then shapecat=create_struct(temporary(shapecat),"sexdec",sexdec[1:*])
  if n_elements(sexfwhm)  gt 1 then shapecat=create_struct(temporary(shapecat),"sexfwhm",sexfwhm[1:*])
  if n_elements(sexflux)  gt 1 then shapecat=create_struct(temporary(shapecat),"sexflux",sexflux[1:*])
  if n_elements(sexmag)   gt 1 then shapecat=create_struct(temporary(shapecat),"sexmag",sexmag[1:*])
  if n_elements(sexclass) gt 1 then shapecat=create_struct(temporary(shapecat),"sexclass",sexclass[1:*])
  if n_elements(sexa)     gt 1 then shapecat=create_struct(temporary(shapecat),"sexa",sexa[1:*])
  if n_elements(sexb)     gt 1 then shapecat=create_struct(temporary(shapecat),"sexb",sexb[1:*])
  if n_elements(sextheta) gt 1 then shapecat=create_struct(temporary(shapecat),"sextheta",sextheta[1:*])
  if n_elements(sexarea)  gt 1 then shapecat=create_struct(temporary(shapecat),"sexarea",sexarea[1:*])
  if n_elements(sexe1)    gt 1 then shapecat=create_struct(temporary(shapecat),"sexe1",sexe1[1:*])
  if n_elements(sexe2)    gt 1 then shapecat=create_struct(temporary(shapecat),"sexe2",sexe2[1:*])
  if n_elements(objx_min) gt 1 then shapecat=create_struct(temporary(shapecat),"objx_min",fix(objx_min[1:*]))
  if n_elements(objx_max) gt 1 then shapecat=create_struct(temporary(shapecat),"objx_max",fix(objx_max[1:*]))
  if n_elements(objy_min) gt 1 then shapecat=create_struct(temporary(shapecat),"objy_min",fix(objy_min[1:*]))
  if n_elements(objy_max) gt 1 then shapecat=create_struct(temporary(shapecat),"objy_max",fix(objy_max[1:*]))
  shapecat=create_struct(temporary(shapecat),"morphology",0B,"shear_estimates",0B)

  ; Write to disc
  shapelets_write, shapecat, filename_output, full_path=full_path
  spawn,'rm -f '+shapelets_paths(2)+filename_output+'_tmp.shapecat' ;######################

  ; Display final message proclaiming success!
  if not keyword_set(silent) then begin
    message,'Finished '+filename_output+'.shapecat',/info
    print,"Success rate: "+strtrim(string(n),2)+" objects found by SExtractor"
    print,"              "+strtrim(string(n_stored),2)+" recorded in shapelet catalogue "+$
          "("+strmid(strtrim(string(100.*n_stored/n),2),0,5)+"%)"
    print,"              "+strtrim(string(n_tried),2)+" decompositions attempted "+$
          "("+strmid(strtrim(string(100.*n_tried/n),2),0,5)+"%)"
    print,"              "+strtrim(string(n_done),2)+" good shapelet models "+$
          "("+strmid(strtrim(string(100.*n_done/n_tried),2),0,5)+"%)"
    print & print
    print,'Deleted temporary shapecat : '+filename_output+'_tmp.shapecat' ;######################
    print & print ;######################
  endif

endelse

end
