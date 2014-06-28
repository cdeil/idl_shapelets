pro plt_objs,obj,focus,decomp,$
   chi2_grid=chi2_grid,centroid=centroid,image_in=image_in

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
;+
; July 2003 - Modified by AR to include plot of residuals
; May 2002   - Minor style adjustments by R. Massey
; April 2002 - Written by A. Refregier
; PURPOSE: plot the results of the focussing of shapelets on an object.
; The history of the search of optimal shapelet parameters is also shown.
; INPUT: pstamp: pstamp structure from shapelets_image2pstamp.pro
;        focus: focus structure from shapelets_focus.pro    
;        decomp: decomposition structure from shapelets_focus.pro  
; OPTIONAL INPUT: chi2_grid: plot the Chi^2(beta,nmax) grid array if provided.
;                   if this keyword is set to a scalar the array is
;                   calculated with ranges around the optimal values
;                 centroid: plot centroid search history
;                 image_in: input image to plot as well (eg. noise
;                      free simulated image)
; OUTPUT: plots
;-

COMPILE_OPT OBSOLETE

; declarations
th_max=focus.beta*sqrt(focus.nmax+1.)
th_min=focus.beta/sqrt(focus.nmax+1.)
th_maxmax=abs(min([focus.x0(0),pstamp.im_ran(1)-pstamp.im_ran(0)-focus.x0(0),$ 
            focus.x0(1),pstamp.im_ran(3)-pstamp.im_ran(2)-focus.x0(1)]))
                  ; imposed by image boundaries
c_ran=[min(pstamp.image),max(pstamp.image)]     ; colar bar range
wsize=[525,550]    ; window x and y sizes
dsize=[1000,1200]   ; display size
flag_pstamp=['OK','','nearby object','severe overlap','near edge','FWHM=0','no backgd',$
             'object masked out']
flag_focus=['OK','!7h!6!imin!n-!7h!6!imax!n bounce','max n!imax!n reached','converged by flatness',$
            'failed to converge','centroid off edge','fatal error']


; plot original object
window,0,title="Original object",$
  xsize=wsize(0),ysize=wsize(1),xpos=0,ypos=dsize(1)-wsize(1)
shapelets_plot_pstamp,pstamp,cran=c_ran,title='Original object',/mask
;oplot,focus.x0_h(0,*),focus.x0_h(1,*)
;oplot,[focus.x0,focus.x0],[0,3000],linestyle=2
;oplot,[0,3000],[focus.x0,focus.x0],linestyle=2

; plot reconstructed object
r=1.
shapelets_decomp,pstamp.image,focus.beta,focus.nmax,decomp,recomp,$
    x0=focus.x0,/ls,/over,noise=pstamp.noise,/integrate
window,1,title="Shapelets reconstruction",$
  xsize=wsize(0),ysize=wsize(1),xpos=wsize(0)+8,ypos=dsize(1)-wsize(1)
shapelets_plot_image,recomp,/fr,/col,cran=c_ran,title='Reconstructed object'
oplot,[focus.x0(0)],[focus.x0(1)],psym=6              ; final 
oplot,focus.x0_h(0,*),focus.x0_h(1,*)                 ; history
oplot,[focus.x0_guess(0)],[focus.x0_guess(1)],psym=1  ; guess
plt_ellipse,pstamp.xo,pstamp.yo,$                 ; Sextractor ellipse
    pstamp.a,pstamp.b,(pstamp.theta+90.)/!radeg,lines=2
;plt_ellipse,focus.x0(0),focus.x0(1),focus.beta,focus.beta,0.
;plt_ellipse,focus.x0(0),focus.x0(1),th_min,th_min,0.,lines=2
;plt_ellipse,focus.x0(0),focus.x0(1),th_max,th_max,0.,lines=2
shapelets_moments,decomp,dmom
plt_ellipse,decomp.x(0),decomp.x(1),dmom.a,dmom.b,dmom.alpha-!pi/2.
oplot,[dmom.xc(0)],[dmom.xc(1)],psym=1

; plot residuals
window,2,title='Residuals',$
  xsize=wsize(0),ysize=wsize(1),xpos=(wsize(0)+8)*2,ypos=dsize(1)-wsize(1)
shapelets_plot_image,pstamp.image-recomp,/fr,/col,title='Residuals'

; plot search history on beta vs nmax plane
window,3,title="Focus history",$
  xsize=wsize(0),ysize=300,xpos=wsize(0)+8,ypos=dsize(1)-wsize(1)-325
!p.region=0
if keyword_set(chi2_grid) then begin
  beta_ran=focus.beta*[.5,1.5] & nmax_ran=[2,(focus.nmax+2)<15]
  s=size(chi2_grid)
  if s(0) eq 0 then begin
    print,'chi^2 grid: nmax_min, nmax_max:',nmax_ran(0),nmax_ran(1)
    print,'chi^2 grid: beta_min, beta_max:',beta_ran(0),beta_ran(1)
    shapelets_make_chi2_grid,obj,chi2_grid,beta_ran=beta_ran,n_ran=nmax_ran,x0=focus.x0
  endif
  shapelets_plot_chisq_grid,chi2_grid;,title='Search history',xtit='n!dmax!n',ytit='beta'
  contour,chi2_grid.chi2,chi2_grid.nmax,chi2_grid.beta,$
    level=focus.chi20+focus.tol_chi2*[-2.,-1.,0.,1.,2.],$
    /over,c_lines=[1,2,0,2,1]
endif else begin
;beta_ran=[focus.beta*.5,focus.beta*1.5]
beta_ran=[0.,th_maxmax]
nmax_ran=[0,max(focus.nmax_h)+2]
plot,focus.nmax_h,focus.beta_h,/nodata,yran=beta_ran,xran=nmax_ran,$
    xtitle='!6n!imax!n',ytitle='!7b!6'
endelse
oplot,focus.nmax_h,focus.beta_h
oplot,focus.nmax_h,focus.beta_h,psym=1
oplot,[focus.nmax_guess],[focus.beta_guess],psym=7
oplot,[focus.nmax],[focus.beta],psym=6
flux=shapelets_flux(decomp,/error)
snr=flux(0)/flux(1)
xyouts,.4,.85,/normal,charsize=1.7,$
       '!7b!6='+string(focus.beta,format='(f6.2)')+', '+$
       'n!imax!n='+string(focus.nmax,format='(i2)')+'!c'+$
       '!7v!6!e2!n='+string(focus.chi2,format='(f7.3)')+', '+$
       'snr='+string(snr,format='(f6.1)')+'!c'+$
       'pstamp: '+string(flag_pstamp(obj.flag))+'!c'+$
       'focus: '+string(flag_focus(focus.flag))
;string(obj.flag)+string(focus.flag)

; plot geometrical constraints
nnmax=xgen(0.,nmax_ran(1)+4)
oplot,nnmax,focus.th_min*sqrt(nnmax+1.),lines=3
oplot,nnmax,th_maxmax/sqrt(nnmax+1.),lines=3

; plot history of search of centroid if requested
if keyword_set(centroid) then begin
;  oplot,focus.x0_h(0,*),focus.x0_h(1,*)
;  oplot,[focus.x0(0)],[focus.x0(1)],psym=6 
  window,4,xsize=500,ysize=500,xpos=130,ypos=50,title='Centroid search'
  !p.region=0
  plot,[reform(focus.x0_h(0,*)),focus.x0(0),focus.x0_guess(0)],$
    [reform(focus.x0_h(1,*)),focus.x0(1),focus.x0_guess(1)],$
    /nodata,/ynozero,xtitle='x',ytitle='y',title='Centroid search'
  oplot,focus.x0_h(0,*),focus.x0_h(1,*)
  oplot,[focus.x0_guess(0)],[focus.x0_guess(1)],psym=1
  oplot,[focus.x0(0)],[focus.x0(1)],psym=6
endif

; plot input image if available
if keyword_set(image_in) then begin
 ; plot input image
  window,5,title='Input image',$
    xsize=wsize(0),ysize=wsize(1),xpos=(wsize(0)+8)*2,$
    ypos=dsize(1)-2*wsize(1)-25
  bit_ref=image_in.image(obj.im_ran(0):obj.im_ran(1),$
                       obj.im_ran(2):obj.im_ran(3))
  shapelets_plot_image,bit_ref,/frame,/col,cran=c_ran
;  ; compute moments after removing neighboring sources
;  othpix= where((obj.seg ne 0) and (obj.seg ne shcat.sexid(io)+1),nothpix)  
;  if nothpix gt 0 then  bit_ref(othpix)=0.
;  cmp_mom,bit_ref,mom
;  plt_ellipse,mom.xc(0),mom.xc(1),mom.a,mom.b,mom.alpha-!pi/2.
;  oplot,[mom.xc(0)],[mom.xc(1)],psym=1
;  ; plot input parameter ellipse if input scat catalog is provided
;  if keyword_set(scat_in) then begin
;    dist_in=sqrt((scat_in.x(*,0)-shcat.centroid(io,0))^2+$
;               (scat_in.x(*,1)-shcat.centroid(io,1))^2)
;    dist=min(dist_in,match)
;    plt_ellipse,scat_in.x(match,0)-obj.im_ran(0),$
;              scat_in.x(match,1)-obj.im_ran(2),$
;              scat_in.a(match)*2.,scat_in.b(match)*2.,$
;              scat_in.theta(match)/!radeg-!pi/2.,lines=2
;;    print,'match: x,y:',scat_in.x(match,*)
;;    print,'match: a,b,theta:',scat_in.a(match),$
;;         scat_in.b(match),scat_in.theta(match)
  endif


end


