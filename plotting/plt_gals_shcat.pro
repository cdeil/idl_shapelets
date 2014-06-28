pro plt_gals_shcat,shcat,class_c=class_c,mag_ran=mag_ran

; August 998 - Written by A. Refregier
;
; PURPOSE: plot various figures relative to the galaxies in a shapelet
; catalog (shcat).
; INPUT: shcat: shapelet catalog read by shapelet_read_shcat.pro
; KEYWORDS: class_c: class cut for star/gal separation (default:.8)
;           mag_ran: magnitude range
; OUTPUT: plots

COMPILE_OPT OBSOLETE

; declarations
if not keyword_set(class_c) then class_c=.8    ; class cut for galaxies
print,'Class cut:',class_c
e_max=2.                            ; max |e| to include in gal catalog
omega_rad=(4096.*.04/3600./!radeg)^2  ; temporary
print,'Warning: default image are is used (arcmin^2):',omega_rad*!radeg^2*3600.
if not keyword_set(mag_ran) then mag_ran=[min(shcat.mag),max(shcat.mag)]
n_mag = 20                  ; number of magnitude bins

; compute ellipticities
e1=shcat.ellipticity(*,0)
e2=shcat.ellipticity(*,1)
ell=sqrt(e1^2+e2^2)
alpha=atan(e2,e1)/2.             ; correct?

; select galaxies
ga=where(shcat.sexclass le class_c and ell le e_max,n_ga) 

; plot mag-size distribution for galaxies
plot,shcat.mag(ga),sqrt(shcat.rsquared(ga)),psym=3,/ytype,xran=mag_ran, $
  xtitle='mag',ytitle='r (pixels)'

; plot ellipticities
read,test
!p.multi = 0
plot,e1(ga),e2(ga),psym=3,$
  xran=[-5,5],yran=[-5,5],$
  xtitle='e1',ytitle='e2'
oplot, [-10., 10.], [0., 0.], lines=1
oplot, [0., 0.], [-10., 10.], lines=1
e1_m = amean(e1(ga)) &  e2_m=amean(e2(ga))
e1_sig = asigma(e1(ga)) & e2_sig = asigma(e2(ga))
e1_err = e1_sig/sqrt(n_ga) &  e2_err= e2_sig/sqrt(n_ga) 
ploterr_xy, e1_m, e2_m, e1_err, e2_err
oplot, [e1_m], [e2_m], psym=6
print,'e1,e2:'
print,'  mean:',e1_m, e2_m
print,'  error:',e1_err, e2_err
print,'  sigma:',e1_sig, e2_sig

; plot ellipticities vs magnitude
read,test
!p.multi=[0,1,2]
;   all gals
plot,shcat.mag(ga),e1(ga),psym=3,yrange=[-2.,2.],xran=mag_ran, $
  xtitle='mag',ytitle='e1'
plot,shcat.mag(ga),e2(ga),psym=3,yrange=[-2.,2.],xran=mag_ran, $
  xtitle='mag',ytitle='e2'
!p.multi=0

; compute statistics as a function of magnitude
mag=bin(mag_ran(0), mag_ran(1),nb=n_mag)
dmag=mag.bsize
dngdm=fltarr(n_mag)          ; differential counts [gals mag^-1 srad^-1]
ng=fltarr(n_mag)             ; cummulative counts [gals srad^-1]
sige1=fltarr(n_mag)          ; ellipticity rms, differential
sige2=fltarr(n_mag)
sige1c=fltarr(n_mag)         ; ellipticity rms, cummulative
sige2c=fltarr(n_mag)

for i=0,n_mag-1 do begin
  ; differential
  gd=where(shcat.sexclass le class_c and ell le e_max and $
               shcat.mag ge mag.l(i) and shcat.mag lt mag.h(i),ngi)
  dngdm(i)=ngi/dmag/omega_rad
  if ngi gt 0 then begin
    sige1(i)=asigma(e1(gd))
    sige2(i)=asigma(e2(gd))
  endif

  ; cummulative
  gd=where(shcat.sexclass le class_c and ell le e_max and $
               shcat.mag gt mag_ran(0) and shcat.mag lt mag.m(i),ngi)
  ng(i)=ngi/omega_rad
  if ngi gt 0 then begin
    sige1c(i)=asigma(e1(gd))
    sige2c(i)=asigma(e2(gd))
  endif
endfor

; plot number counts
read,test
!p.multi=[0,1,2]
plot,mag.m,dngdm/2./!radeg^2,psym=6,/ytype,$
  xtitle='mag',ytitle='dn/dm (gals/deg^2/(.5 mag)'
;oplot,[19.5,21.5,24.5,27.5,29.5],10.^[2.8,3.6,4.6,5.2,6.0],$
;  psym=1   ; Pozzetti et al. 

plot,mag.m,ng/!radeg^2,psym=6,/ytype,$
  xtitle='mag',ytitle='N(>m) (gals/deg^2)'

; plot ellipticity dispersion vs mag
read,test
plot,mag.m,sige1,psym=5,$
  xtitle='mag',ytitle='sigma(ei)',title='differential'
oplot,mag.m,sige2,psym=7
legend,['e1','e2'],psym=[5,7],box=0
plot,mag.m,sige1c,psym=5,$
  xtitle='m!ii!n',ytitle='sigma(ei)',title='cummulative'
oplot,mag.m,sige2c,psym=7
!p.multi=0

; plot ellipticity accross the field
read,test
xy_ran = [min(shcat.centroid), max(shcat.centroid)]
plot,shcat.centroid(ga,0),shcat.centroid(ga,1),$
  psym=3,xran=xy_ran, yran=xy_ran, /xstyle, /ystyle, $
  xtitle='x (pixels)',ytitle='y (pixels)'
scale=30.
plt_evec,shcat.centroid(ga,0),shcat.centroid(ga,1),$
  ell(ga),alpha(ga)-!pi/2.,xs=scale,ys=scale

end

