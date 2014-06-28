pro plt_gals,cat,class_c=class_c,mag_ran=mi_ran

; August 998 - Written by A. Refregier
;
; PURPOSE: plot various figures relative to the galaxies in a catalog.
; INPUT: cat: catalog read by shapelet_read_scat.pro
; KEYWORDS: class_c: class cut for star/gal separation (default:.8)
;           mag_ran: magnitude range
; OUTPUT: plots

; declarations
if not keyword_set(class_c) then class_c=.8    ; class cut for galaxies
print,'Class cut:',class_c
;a_psf=1.57                  ; psf rms radius, assumed circular [pixels]
a_psf = 0.
;omega_rad=cat.omega/3600./!radeg^2
;omega_rad=(4096.*.04/3600./!radeg)^2  ; temporary
omega_rad=1./!radeg^2          ; default
print,'Warning: default image area is used (deg^2):',$
   omega_rad*!radeg^2
;th_x=8./60./!radeg              ; field linear sizes [rad]
;th_y=16./60./!radeg  
;th_x=1.78/60./!radeg
;th_y=3.27/60./!radeg
;th_pix=cat.th_pix/3600./!radeg        ; pixel size [rad]
;n_x=2000.                   ; number of pixels accross the field
;n_y =4000.
;n_x=713.
;n_y=1309.
mag_ran = [min(cat.mag),max(cat.mag)] ; magnitude range to consider
if keyword_set(mi_ran) then mag_ran=mi_ran
n_mag = 20                  ; number of magnitude bins
d=sqrt(cat.a^2+cat.b^2)

; select galaxies
ga=where(cat.class le class_c)    ; all galaxies
gar=where(cat.class le class_c and cat.a gt a_psf and cat.b gt a_psf, n_gar)
                                  ; resolved galaxies (a,b>a_psf)

; plot mag-size distribution for galaxies
;!p.multi=[0,1,2]
plot,cat.mag(ga),d(ga),psym=3,/ytype,xran=mag_ran, $
  xtitle='m!ii!n',ytitle='d (pixels)'
oplot,[-100,100],[1.,1.]*a_psf,lines=1

; plot maj-minor relation for galaxies
;read,test
;plot,cat.a(ga),cat.b(ga),psym=3,xrange=[0.,25.],yrange=[0.,25.],$
;  xtitle='a (pixels)',ytitle='b (pixels)'
;oplot,[1.,1.]*a_psf,[0.,100.],lines=1
;oplot,[0.,100.],[1.,1.]*a_psf,lines=1
;oplot,[0.,100.],[0.,100.],lines=2

; compute deconvolved ellipticities
alpha = cat.theta/!radeg    ; convert PA to radians
e=(cat.a^2-cat.b^2)/(cat.a^2+cat.b^2-2.*a_psf^2)
e1=e*cos(2.*alpha)
e2=e*sin(2.*alpha)

; plot ellipticities
read,test
!p.multi = 0
plot,e1(gar),e2(gar),psym=3,$
  xtitle='e1',ytitle='e2'
oplot, [-5., 5.], [0., 0.], lines=1
oplot, [0., 0.], [-5., 5.], lines=1
e1_m = mean(e1(gar)) &  e2_m=mean(e2(gar))
e1_sig = stdev(e1(gar)) & e2_sig = stdev(e2(gar))
e1_err = e1_sig/sqrt(n_gar) &  e2_err= e2_sig/sqrt(n_gar) 
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
plot,cat.mag(ga),e1(ga),psym=3,yrange=[-2.,2.],xran=mag_ran, $
  xtitle='m!ii!n',ytitle='e1 (deconvolved)',title='all galaxies'
plot,cat.mag(ga),e2(ga),psym=3,yrange=[-2.,2.],xran=mag_ran, $
  xtitle='m!ii!n',ytitle='e2 (deconvolved)'
;; resolved gals
;read,test
;plot,cat.mi(gar),e1(gar),psym=3,yrange=[-2.,2.],xran=mag_ran, $
;  xtitle='m!ii!n',ytitle='e1 (deconvolved)',$
;  title='resolved galaxies (a,b>a_psf)'
;plot,cat.mi(gar),e2(gar),psym=3,yrange=[-2.,2.],xran=mag_ran, $
;  xtitle='m!ii!n',ytitle='e2 (deconvolved)'

!p.multi=0

; compute statistics as a function of magnitude
mag=bin(mag_ran(0), mag_ran(1),nb=n_mag)
dmag=mag.bsize
dngdm=fltarr(n_mag)          ; differential counts [gals mag^-1 srad^-1]
ng=fltarr(n_mag)             ; cummulative counts [gals srad^-1]
dngrdm=fltarr(n_mag)         ; same but for resolved galaxies (ad,bd>0)
ngr=fltarr(n_mag)
sige1=fltarr(n_mag)          ; ellipticity rms, differential
sige2=fltarr(n_mag)
sige1c=fltarr(n_mag)         ; ellipticity rms, cummulative
sige2c=fltarr(n_mag)


for i=0,n_mag-1 do begin

  ; all galaxies
    ; differential
  gd=where(cat.class le class_c and $
               cat.mag ge mag.l(i) and cat.mag lt mag.h(i),ngi)
  dngdm(i)=ngi/dmag/omega_rad
    ; cummulative
  gd=where(cat.class le class_c and $
               cat.mag gt mag_ran(0) and cat.mag lt mag.m(i),ngi)
  ng(i)=ngi/omega_rad

  ; resolved galaxies
    ; differential
  gd=where(cat.class le class_c and cat.a gt a_psf and cat.b gt a_psf and $
               cat.mag ge mag.l(i) and cat.mag lt mag.h(i),ngi)
  dngrdm(i)=ngi/dmag/omega_rad
  if ngi gt 1 then begin
    sige1(i)=stdev(e1(gd))
    sige2(i)=stdev(e2(gd))
  endif
    ; cummulative
  gd=where(cat.class le class_c and cat.a gt a_psf and cat.b gt a_psf and $
              cat.mag gt mag_ran(0) and cat.mag lt mag.m(i),ngi)
  ngr(i)=ngi/omega_rad
  if ngi gt 1 then begin
    sige1c(i)=stdev(e1(gd))
    sige2c(i)=stdev(e2(gd))
  endif
  
endfor

; plot number counts
read,test
!p.multi=[0,1,2]
plot,mag.m,dngdm/2./!radeg^2,psym=6,/ytype,yrange=[1.,1e6],$
  xtitle='m!ii!n',ytitle='dn/dm (gals/deg^2/(.5 mag)'
oplot,mag.m,dngrdm/2./!radeg^2,psym=5
;oplot,[19.5,21.5,24.5,27.5,29.5],10.^[2.8,3.6,4.6,5.2,6.0],$
;  psym=1   ; Pozzetti et al. 
legend,['all galaxies','resolved galaxies'],$
  psym=[6,5],box=0

plot,mag.m,ng/!radeg^2,psym=6,/ytype,yrange=[1.,1e6],$
  xtitle='m!ii!n',ytitle='N(>m) (gals/deg^2)'
oplot,mag.m,ngr/!radeg^2,psym=5

; plot ellipticity dispersion vs mag
read,test
plot,mag.m,sige1,psym=5,$
  xtitle='m!ii!n',ytitle='sigma(ei)',title='differential'
oplot,mag.m,sige2,psym=7
legend,['e1','e2'],psym=[5,7],box=0
plot,mag.m,sige1c,psym=5,$
  xtitle='m!ii!n',ytitle='sigma(ei)',title='cummulative'
oplot,mag.m,sige2c,psym=7

; plot ellipticity accross the field
read,test
!p.multi=0
xy_ran = [min(cat.x), max(cat.x)]
plot,cat.x(gar,0),cat.x(gar,1),psym=3,$
  xran=xy_ran, yran=xy_ran, /xstyle, /ystyle, $
  xtitle='x (pixels)',ytitle='y (pixels)'
scale=30.
plt_evec,cat.x(gar,0),cat.x(gar,1),sqrt(cat.e1(gar)^2+cat.e2(gar)^2),$
  alpha(gar)-!pi/2.,xs=scale,ys=scale

; plot ellipse as function of position on the chip
;read, test
;plot,cat.x(gar),cat.y(gar),psym=3,xran=xy_ran, yran=xy_ran, /xstyle,$
;  /ystyle, xtitle='x (pixels)',ytitle='y (pixels)', /nodata
;plt_ellipse,cat.x(gar),cat.y(gar), cat.a(gar), cat.b(gar), $
;  cat.alpha(gar)-!pi/2., xs=10., ys=10.


!p.multi=0

end

