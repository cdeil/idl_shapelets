;$Id: shapelets_plot_focus.pro, v2$
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
; ***********************************************************************
; ***********************************************************************
;
;+
; NAME:
;      SHAPELETS_IMAGE_PROFILE
;
; CATEGORY:
;      Subroutine of shapelets_plot_focus.pro.
;
; PURPOSE:
;      Computes the radial profile of an image array.
;
; INPUTS:
;      IM - 2D image array.
;
; OPTIONAL INPUTS:
;      X0 - Object"s centre (DEFAULT: centre of the array)
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      PROFILE
;
; MODIFICATION HISTORY:
;      Apr 05 - Modified by Richard Massey.
;      Mar 02 - Modified by AR to be more general 
;      Jul 99 - Written by A. Refregier
;-

function shapelets_image_profile, im, x0=x0

COMPILE_OPT idl2, HIDDEN

; declarations
si=size(im)
n1=si[1] & n2=si[2]
if not keyword_set(x0) then begin    ; set object center
   x0=[n1/2.,n2/2.]             
endif
n_th=15              ; number of theta bins

; define angular bins - linear spacing
th_ran = [0., sqrt(2.)*float(min([n1,n2])/2.)]
dth = (th_ran[1]-th_ran[0])/float(n_th)
th_l = findgen(n_th)*dth+th_ran[0]
th_m = th_l+.5*dth
th_h = th_l+dth

; define angular bins - log spacing
th_ran_log = [1., sqrt(2.)*float(min([n1,n2])/2.)]
dlogth =  (alog10(th_ran_log[1])-alog10(th_ran_log[0]))/float(n_th)
;th_l = findgen(n_th)*dth+th_ran_log[0]
;th_m = th_l+.5*dth
;th_h = th_l+dth
logth_l = findgen(n_th)*dlogth+alog10(th_ran_log[0])
logth_m = logth_l+.5*dlogth 

; define arrays
ith = fltarr(n_th)
ith_log = fltarr(n_th)
nth = lonarr(n_th)
nth_log = lonarr(n_th)
thbin=fltarr(n1,n2)
logthbin=fltarr(n1,n2)
; compute profile
for i=0, n1-1 do begin
   for j=0, n2-1 do begin
      th = sqrt(float(i-x0[0])^2+float(j-x0[1])^2)
      thi = fix((th-th_ran[0])/dth)
      thbin[i,j]=thi
      thi_log = fix((alog10(th)-alog10(th_ran_log[0]))/dlogth)
      logthbin[i,j]=thi_log      
      if thi lt n_th-1 and thi ge 0 then begin
        ith[thi] = ith[thi]+im[i,j]
        nth[thi] = nth[thi]+1l
      endif   
      if thi_log lt n_th-1 and thi_log ge 0 then begin
        ith_log[thi_log] = ith_log[thi_log]+im[i,j]
        nth_log[thi_log] = nth_log[thi_log]+1L
      endif 
     endfor
endfor
gd=where(nth gt 0)
ith = ith[gd]/float(nth[gd])
gd_log=where(nth_log gt 0)
ith_log = ith_log[gd_log]/float(nth_log[gd_log])


; store in structure
profile={r:th_m[gd],f:ith[gd],r_log:10.^(logth_m[gd_log]),$
  f_log:ith_log[gd_log]}
return,profile

end



; ***********************************************************************
; ***********************************************************************
;
;+
; NAME:
;       SHAPELETS_PLOT_FOCUS
; 
; PURPOSE:
;       Plot the route taken through possible n_max and beta values to
;       obtain a final decomposition.
;
; CATEGORY:
;       Shapelets/plotting.
;
; INPUTS:
;       FOCUS    - Focus history structure.
;       PSTAMP   - Postage stamp structure.
;
; OPTIONAL INPUTS:
;       DECOMP   - Shapelet decomp structure of final decomposition.
;       IMAGE_IN - Second input image to plot (eg. noise free simulated image). 
;       SCAT_IN  - Input sextractor catalog to plot input ellipse as well.
;       PSF      - PSF for convolution, as a decomp structure type.
;
; KEYWORDS:
;       None.
;
; OUTPUTS:
;       Several plots are drawn to STDOUT.
;
; MODIFICATION HISTORY:
;       Aug 05 - Upper geometic constraints changed by RM to match focus.
;       Apr 05 - Updated by RM.
;       Jan 05 - idl2 compatability ([]s, etc) implemented by RM.
;       Jan 05 - Color scale ranges normalised, and PSF keyword added by AR.
;       Nov 04 - Window positions generalised to my screen by Richard Massey.
;       Jul 04 - Modified by AR to allow more flexibility in indexing.
;       Aug 03 - Written by Alexandre Refregier.
;-

pro shapelets_plot_focus, FOCUS, PSTAMP,     $
                          DECOMP=decomp,     $
                          RECOMP=recomp,     $
                          PSF=psf,           $
                          IMAGE_IN=image_in, $
                          SCAT_IN=scat_in,   $
                          NOISE_MAP=noise_map

COMPILE_OPT idl2

; Parse required inputs
if n_params() ne 2 then message,"Usage: shapelets_plot_focus,focus,pstamp"
if not shapelets_structure_type(focus,message=message) then message,message
if focus.type ne "focus" then message,"Not a shapelets focus structure!"
if not shapelets_structure_type(pstamp,message=message) then message,message
if pstamp.type ne "pstamp" then message,"Not a shapelets pstamp structure!"


; Display parameters
device,get_screen_size=dsize  ;Ask IDL for the display's resolution...
if keyword_set(psf) then n_windows=4 else n_windows=3
wborder=[5,20]                ; Number of pixels used by window manager as border around [L/R/B,T] of windows
csize=0.09+0.03*n_windows     ; Fractional size of colour bar
charsize=1.3                  ; Size of label text, realtive to IDL default
wsize=[(dsize[0])/n_windows-2*wborder[0],$                   ; Window width 
       (dsize[0])/n_windows*(1+csize)-wborder[0]-wborder[1]] ; Window height
c_ran=[min(pstamp.image),max(pstamp.image)]                  ; Colour scale range


;
; Plot the original image
;
window,0,title="Original object",$
  xsize=wsize[0],ysize=wsize[1],$
  xpos=0,        ypos=dsize[1]-wsize[1]
shapelets_plot_pstamp,pstamp,title="!6Original object",$
  /mask,cran=c_ran,csize=csize,linestyle=1


;
; Plot the reconstructed image
;
window,1,title="Shapelets reconstruction",$
  xsize=wsize[0],            ysize=wsize[1],$
  xpos=wsize[0]+2*wborder[0],ypos=dsize[1]-wsize[1]
if keyword_set(decomp) then begin
  if not shapelets_structure_type(decomp,message=message) then message,message
  if decomp.type ne "decomp" then message,"Not a shapelets decomp structure!"
endif else decomp=shapelets_decomp(pstamp,focus.beta,focus.n_max,recomp=recomp,centre=focus.x,psf=psf)
if not keyword_set(recomp) then shapelets_recomp,decomp,recomp
shapelets_plot_image,recomp,/fr,/col,title="!6Reconstructed object",csize=csize;,cran=c_ran
plt_ellipse,pstamp.xo,pstamp.yo,$                 ; Sextractor ellipse
    pstamp.a,pstamp.b,(pstamp.theta+90.)/!radeg,linestyle=1
if keyword_set(focus) then begin ; plot focus history if available
  ;oplot,[focus.x_guess[0]],[focus.x_guess[1]],psym=1     ; guess
  oplot,focus.x_history[0,*],focus.x_history[1,*],psym=-3 ; whole history
  oplot,[focus.x[0]],[focus.x[1]],psym=4,symsize=sqrt(2.) ; final 
endif
; Plot unweighted ellipticity
q=shapelets_quadrupole(decomp,flux=flux)
if n_elements(flux) eq 0 then flux=0
a=sqrt((q[0,0]+q[1,1]+sqrt((q[0,0]-q[1,1])^2+4.*q[0,1]^2))/flux>1e-10)  ; principal axes [pixels]
b=sqrt((q[0,0]+q[1,1]-sqrt((q[0,0]-q[1,1])^2+4.*q[0,1]^2))/flux>1e-10)
alpha=.5*atan(2.*q[0,1],q[0,0]-q[1,1])   ; position angle, counter-clockwise from x-axis
plt_ellipse,decomp.x[0],decomp.x[1],a,b,alpha-!pi/2.
centroid=shapelets_centroid(decomp)
; Plot centre of light
oplot,[centroid[0]],[centroid[1]],psym=7,symsize=3,color=0
; Plot minimum and maximum scales
theta_min=focus.beta/sqrt(focus.n_max+1)
theta_max=focus.beta*(focus.n_max+1)^0.59
plt_ellipse,focus.x[0],focus.x[1],theta_min,theta_min,0,linestyle=2
plt_ellipse,focus.x[0],focus.x[1],theta_max,theta_max,0,linestyle=2


;
; Plot histogram of pixel values
;
window,4,title="Histogram",$
    xsize=wsize[0],ysize=wsize[1]*2/3,$ 
    xpos=0,        ypos=dsize[1]-wsize[1]*5/3-wborder[0]-wborder[1]
!p.region=0
gd=where(pstamp.seg ne -1)
im_ran=c_ran
binsize=(im_ran[1]-im_ran[0])/50
plothist,pstamp.image[gd],x_values,y_values,$
   xrange=c_ran,/xstyle,bin=binsize,linestyle=0,$
   xtitle="!6Pixel value !8I!6",ytitle="!6Number of pixels",title="!6Postage stamp histogram"
yrange=[0,10*max(y_values)]
plothist,pstamp.image[where(pstamp.seg eq 0)],linestyle=2,/over,bin=binsize
if keyword_set(noise_map) then begin
  ext_linestyle=2
  local_linestyle=1
  xyouts,.9,.85,/norm,charsize=charsize,align=1,$
        "!3- - ext back ave= "+string(pstamp.back_mean_ext,format="(e9.2)")+"!c"+$
          "- - ext back rms= "+string(pstamp.back_rms_ext,format="(e9.2)")+"!c"+$
        ". . . loc back ave= "+string(pstamp.back_mean_local,format="(e9.2)")+"!c"+$
        ". . . loc back rms= "+string(pstamp.back_rms_local,format="(e9.2)")
  ;xyouts,.9,.85,/norm,charsize=charsize,align=0,"!3. . .!c. . .!c- -!c- -"
endif else begin
  ext_linestyle=1
  local_linestyle=2
  xyouts,.9,.85,/norm,charsize=charsize,align=1,$
        "!3- - loc back ave= "+string(pstamp.back_mean_local,format="(e9.2)")+"!c"+$
          "- - loc back rms= "+string(pstamp.back_rms_local,format="(e9.2)")+"!c"+$
        ". . . ext back rms= "+string(pstamp.back_rms_ext,format="(e9.2)")
  ;xyouts,.9,.85,/norm,charsize=charsize,align=0,"!3- -!c- -!c. . ."
endelse
oplot,replicate(pstamp.back_mean_local,2),yrange,linestyle=local_linestyle,psym=-3
oplot,replicate(pstamp.back_mean_local-pstamp.back_rms_local,2),yrange,linestyle=local_linestyle,psym=-3
oplot,replicate(pstamp.back_mean_local+pstamp.back_rms_local,2),yrange,linestyle=local_linestyle,psym=-3
oplot,replicate(pstamp.back_mean_ext,2),yrange,linestyle=ext_linestyle,psym=-3
oplot,replicate(pstamp.back_mean_ext-pstamp.back_rms_ext,2),yrange,linestyle=ext_linestyle,psym=-3
oplot,replicate(pstamp.back_mean_ext+pstamp.back_rms_ext,2),yrange,linestyle=ext_linestyle,psym=-3


;
; Plot focus history if available
;
if keyword_set(focus) then begin
  case n_windows of
    3: xsize=wsize[0]
    4: xsize=2*wsize[0]+2*wborder[0]
  endcase
  window,3,title="Focus history",$
    xsize=xsize,                   ysize=wsize[1]*2/3,$
    xpos=2*(wsize[0]+2*wborder[0]),ypos=dsize[1]-wsize[1]*5/3-wborder[0]-wborder[1]
  ;th_max=focus.beta*sqrt(focus.n_max+1.)
  ;th_min=focus.beta/sqrt(focus.n_max+1.)
  th_maxmax=abs(min([focus.x[0],$  ; imposed by image boundaries
       pstamp.im_ran[1]-pstamp.im_ran[0]-focus.x[0],$ 
       focus.x[1],pstamp.im_ran[3]-pstamp.im_ran[2]-focus.x[1]]))
  beta_ran=[0.,th_maxmax]
  n_max_ran=[0,max(focus.n_max_history)+2]
  !p.region=0
  plot,[0,0],/nodata,yran=beta_ran,xran=n_max_ran,$
    xtitle="!8n!6!dmax!n",ytitle="!6!7b!6",title="!6Focus history"
  if focus.n_iterations gt 1 then begin
    oplot,focus.n_max_history,focus.beta_history,psym=-3
    oplot,focus.n_max_history,focus.beta_history,psym=1
  endif
  usersym, cos(2*!pi*findgen(21)/20), sin(2*!pi*findgen(21)/20)
  oplot,[focus.n_max_guess<(n_max_ran[1])],[focus.beta_guess],psym=8,symsize=3
  oplot,[focus.n_max],[focus.beta],psym=7,symsize=3
  flux=shapelets_flux(decomp,/error)
  snr=flux[0]/flux[1]
  xyouts,.9,.85,/normal,charsize=charsize,align=1,$
       "!7b!3="+string(focus.beta,format="(f6.2)")+", !8n!3!dmax!n="+strtrim(string(focus.n_max),2)+$
       "!C!7v!3!e2!n="+strtrim(string(focus.chisq),2)+", snr="+strtrim(string(snr),2)+$
       "!Cpstamp: "+pstamp.flag_interpret_mini[pstamp.flag]+$;" ("+strtrim(pstamp.flag,2)+")!c"+$
       "!Cfocus: "+focus.flag_interpret_mini[focus.flag];+" ("+strtrim(focus.flag,2)+")"
  ; Plot geometrical constraints
  nn_max=xgen(0.,n_max_ran[1]+4)
  oplot,nn_max,focus.theta_min_geom*sqrt(nn_max+1.),linestyle=3,psym=-3
  oplot,nn_max,th_maxmax/(nn_max+1.)^0.59,linestyle=3,psym=-3
  ; Plot singular matrix constraint
  if keyword_set(psf) then begin
    N_alpha=indgen(100)+1
    N_beta=psf.n_max+1.
    N_gamma=N_alpha-0.5
    singular_alpha=psf.beta*sqrt((  N_alpha*(N_beta^2-N_gamma^2)/$
  				   (N_beta*(N_gamma^2-N_alpha^2))    )>0)
    oplot,indgen(100),singular_alpha,linestyle=3,psym=-3
  endif
endif


;
; Plot reconvolved image if deconvolution was applied
;
if keyword_set(psf) then begin
  decomp_r=decomp
  shapelets_convolve,decomp_r,psf
  decomp_r.x=decomp.x
  decomp_r.n_pixels=decomp.n_pixels
  shapelets_recomp,decomp_r,recomp_psf;,nran=[0,0]
  window,7,title="Reconvolved reconstruction",$
    xsize=wsize[0],                ysize=wsize[1],$
    xpos=3*(wsize[0]+2*wborder[0]),ypos=dsize[1]-wsize[1]
  shapelets_plot_image,recomp_psf,/fr,/col,title="!6Reconvolved Reconstruction",csize=csize;,cran=c_ran
  plt_ellipse,pstamp.xo,pstamp.yo,$                 ; Sextractor ellipse
      pstamp.a,pstamp.b,(pstamp.theta+90.)/!radeg,linestyle=2
endif



;
; Plot object profile
;
prof=shapelets_image_profile(pstamp.image,x0=[pstamp.XO,pstamp.YO])
window,6,title="Profile",$
    xsize=wsize[0],            ysize=wsize[1]*2/3,$
    xpos=wsize[0]+2*wborder[0],ypos=dsize[1]-wsize[1]*5/3-wborder[0]-wborder[1]
!p.region=0
r=xgen(0.,max(prof.r))
fsh_r=shapelets_profile(decomp,r)
profile_ylog=1B
if profile_ylog then begin
  plot,prof.r,prof.f,yran=[max(fsh_r)*1e-4,max(fsh_r)],/ystyle,psym=6,$
    xtitle="!8r!6 [pixels]",ytitle="!8I!6 (!8r!6)",title="!6Profile",/ylog
endif else begin
  plot,prof.r,prof.f,yran=[c_ran[0],max(fsh_r)],/ystyle,psym=6,$
    xtitle="!8r!6 [pixels]",ytitle="!8I!6 (!8r!6)",title="!6Profile"
endelse
oplot,[0,1000],replicate(pstamp.back_mean_local,2),linestyle=local_linestyle,psym=-3
oplot,[0,1000],replicate(pstamp.back_mean_local-pstamp.back_rms_local,2),linestyle=local_linestyle,psym=-3
oplot,[0,1000],replicate(pstamp.back_mean_local+pstamp.back_rms_local,2),linestyle=local_linestyle,psym=-3
oplot,[0,1000],replicate(pstamp.back_mean_ext,2),linestyle=ext_linestyle,psym=-3
oplot,[0,1000],replicate(pstamp.back_mean_ext-pstamp.back_rms_ext,2),linestyle=ext_linestyle,psym=-3
oplot,[0,1000],replicate(pstamp.back_mean_ext+pstamp.back_rms_ext,2),linestyle=ext_linestyle,psym=-3
oplot,r,fsh_r,linestyle=0,psym=-3
if keyword_set(image_in) then begin   ; plot input porfile if available
  bit_ref=image_in.image[pstamp.im_ran[0]:pstamp.im_ran[1],$
                       pstamp.im_ran[2]:pstamp.im_ran[3]] 
  mom=shapelets_image_moments(bit_ref)
  prof_in=shapelets_image_profile(bit_ref,x0=mom.xc)
  oplot,prof_in.r,prof_in.f,psym=1
endif
if keyword_set(psf) then begin
   oplot,r,shapelets_profile(decomp_r,r),linestyle=3,psym=-3
   oplot,r,shapelets_profile(psf,r)*$
         (shapelets_profile(decomp_r,0))[0]/(shapelets_profile(psf,0))[0],linestyle=4,psym=-3
endif



;
; Plot residual map
;
if keyword_set(psf) then residual=pstamp.image-recomp_psf $
    else residual=pstamp.image-recomp
window,2,title="Residuals",$
  xsize=wsize[0],                ysize=wsize[1],$
  xpos=2*(wsize[0]+2*wborder[0]),ypos=dsize[1]-wsize[1]
shapelets_plot_image,residual,/fr,/col,title="!6Residuals",csize=csize



;
; Plot input image if available
;
if keyword_set(image_in) then begin
  
  ; plot input image
  window,5,title="Input image",$
    xsize=wsize[0],                ysize=wsize[1],$
    xpos=2*(wsize[0]+2*wborder[0]),ypos=dsize[1]-2*wsize[1]-wborder[0]-wborder[1]
  bit_ref=image_in.image[pstamp.im_ran[0]:pstamp.im_ran[1],$
                       pstamp.im_ran[2]:pstamp.im_ran[3]]
  shapelets_plot_image,bit_ref,cran=c_ran,/frame,/col,csize=csize,$
    title="!6Input image"
  plt_ellipse,decomp.x[0],decomp.x[1],dmom.a,dmom.b,dmom.alpha-!pi/2.
  
  ; compute image moments after removing neighboring sources
  othpix= where((pstamp.seg ne 0) and (pstamp.seg ne pstamp.sexid),nothpix)  
  if nothpix gt 0 then  bit_ref[othpix]=0.
  ;cmp_mom,bit_ref,mom                  ; Alexandre's old version
  mom=shapelets_image_moments(bit_ref)  ; Richard's standardised routine name
  plt_ellipse,mom.xc[0],mom.xc[1],mom.a,mom.b,mom.alpha-!pi/2.,linestyle=1
  oplot,[mom.xc[0]],[mom.xc[1]],psym=1
  
  ; plot input parameter ellipse if input scat catalog is provided
  if keyword_set(scat_in) then begin
    centroid=shapelets_centroid(decomp)+[pstamp.im_ran[0],pstamp.im_ran[2]]
    dist_in=sqrt((scat_in.x[*,0]-centroid[0])^2+$
               (scat_in.x[*,1]-centroid[1])^2)
    dist=min(dist_in,match)
    plt_ellipse,scat_in.x[match,0]-pstamp.im_ran[0],$
              scat_in.x[match,1]-pstamp.im_ran[2],$
              scat_in.a[match]*2.,scat_in.b[match]*2.,$
              scat_in.theta[match]/!radeg-!pi/2.,linestyle=3
  endif
endif 

end
