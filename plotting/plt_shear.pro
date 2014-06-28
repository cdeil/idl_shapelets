pro plt_shear,shcat,class_c=class_c,mag_ran=mag_ran,e_max=e_max,$
  flag_max=flag_max,snr_min=snr_min

; August 998 - Written by A. Refregier
;
; PURPOSE: plot shear estimator statistics
; INPUT: shcat: shapelet catalog read by shapelet_read_shcat.pro
; KEYWORDS: class_c: class cut for star/gal separation (default:.8)
;           mag_ran: magnitude range
;           e_max: maximum value of the ellipticity to keep
;           flag_max: maximum flag values (2 values)
; OUTPUT: plots

COMPILE_OPT OBSOLETE

; declarations
if not keyword_set(class_c) then class_c=.8    ; class cut for galaxies
print,'Class cut:',class_c
if not keyword_set(e_max) then e_max=2. 
print,'e_max:',e_max                           ; ellipticity cut for galaxies
if not keyword_set(mag_ran) then mag_ran=[min(shcat.mag),max(shcat.mag)] 
if not keyword_set(flag_max) then flag_max=[2,1]
print,'max flag1, flag2:',flag_max(0),flag_max(1)
snr=shcat.flux/shcat.flux_error  ; signal-to-noise ratio
if not keyword_set(snr_min) then snr_min=10.
tek_color
col=[!p.color,2,3]

; compute ellipticities
e1=shcat.ellipticity(*,0)
e2=shcat.ellipticity(*,1)
ell=sqrt(e1^2+e2^2)
alpha=atan(e2,e1)/2.             ; correct?

; select galaxies
;ga=where(shcat.sexclass le class_c and ell le e_max,n_ga) 
ga=where(shcat.sexclass le class_c and ell le e_max $
          and shcat.mag ge mag_ran(0) and shcat.mag lt mag_ran(1) $
          and snr gt snr_min $
          and shcat.flag(*,0) le flag_max(0) and shcat.flag(*,1) le flag_max(1),n_ga) 
print,'No of galaxies: total, kept:',n_elements(shcat.mag),n_ga

; plot mag-size distribution for galaxies
plot,shcat.mag,sqrt(shcat.rsquared),psym=3,/ytype, $
  xtitle='mag',ytitle='r (pixels)'
oplot,shcat.mag(ga),sqrt(shcat.rsquared(ga)),psym=6

; plot ellipticities
read,test
!p.multi = 0
plot,e1(ga),e2(ga),psym=3,$
  xran=[-5,5],yran=[-5,5],$
  xtitle='e1',ytitle='e2'
;oplot,e1(gad),e2(gad),psym=6
oplot, [-10., 10.], [0., 0.], lines=1
oplot, [0., 0.], [-10., 10.], lines=1
e1_m = amean(e1(ga)) &  e2_m=amean(e2(ga))
e1_sig = asigma(e1(ga)) & e2_sig = asigma(e2(ga))
e_sig=sqrt(e1_sig^2+e2_sig^2)
e1_err = e1_sig/sqrt(n_ga) &  e2_err= e2_sig/sqrt(n_ga) 
ploterr_xy, e1_m, e2_m, e1_err, e2_err
oplot, [e1_m], [e2_m], psym=6
print,'e1,e2:'
print,'  mean:',e1_m, e2_m
print,'  error:',e1_err, e2_err
print,'  sigma:',e1_sig, e2_sig

; compute and print shear
p_gamma=2.-e_sig^2
print,'P_gamma:',p_gamma
g1=e1/p_gamma & g2=e2/p_gamma
g1_m=amean(g1(ga)) & g2_m=amean(g2(ga))
g1_sig=asigma(g1(ga)) & g2_sig=asigma(g2(ga))
g1_err=g1_sig/sqrt(n_ga) &  g2_err= g2_sig/sqrt(n_ga)
print,'gamma1, gamma2:'
print,'  mean:',g1_m,g2_m
print,'  error:',g1_err,g2_err
print,'  sigma:',g1_sig,g2_sig

; plot shear distribution
read,test
!p.multi=[0,1,2]
plothist,g1(ga),bin=.02,xran=[-1.,1.],$
  xtitle='gamma1',ytitle='No galaxies'
oplot,[0.,0.],[0,1000],lines=1
oplot,[1.,1.]*g1_m,[0,1000],col=col(2)
plothist,g2(ga),bin=.02,xran=[-1.,1.],$
  xtitle='gamma1',ytitle='No galaxies'
oplot,[0.,0.],[0,1000],lines=1
oplot,[1.,1.]*g2_m,[0,1000],col=col(2)
!p.multi=0

; plot shear vs mag
read,test
!p.multi=[0,1,2]
plot,shcat.mag(ga),g1(ga),psym=3,$
  xtitle='mag',ytitle='gamma1'
oplot,[-100,100],[0,0],lines=0
oplot,[-100,100],[1,1]*g1_m,lines=0,col=col(2)
oplot,[-100,100],[1,1]*g1_m+g1_err,lines=2,col=col(2)
oplot,[-100,100],[1,1]*g1_m-g1_err,lines=2,col=col(2)
slice_fit,shcat.mag(ga),g1(ga),10,5,/plotit,/bars,color=col(1)
plot,shcat.mag(ga),g2(ga),psym=3,$
  xtitle='mag',ytitle='gamma2'
oplot,[-100,100],[0,0],lines=0
oplot,[-100,100],[1,1]*g2_m,lines=0,col=col(2)
oplot,[-100,100],[1,1]*g2_m+g2_err,lines=2,col=col(2)
oplot,[-100,100],[1,1]*g2_m-g2_err,lines=2,col=col(2)
slice_fit,shcat.mag(ga),g2(ga),10,5,/plotit,/bars,color=col(1)
!p.multi=0



; plot unnormalised ellipticities
;u1=shcat.quadrupole(*,0,0)-shcat.quadrupole(*,1,1)
;u2=2.*shcat.quadrupole(*,0,1)
;fr2=shcat.quadrupole(*,0,0)+shcat.quadrupole(*,1,1)
;read,test
;!p.multi=[0,1,2]
;plot,fr2(ga),u1(ga),psym=3,/xtype,$
;  xtitle='FR^2',ytitle='u1'
;oplot,[-50,1000],[0,0]
;plot,fr2(ga),u2(ga),psym=3,/xtype,$
;  xtitle='FR^2',ytitle='u1'
;oplot,[-50,1000],[0,0]
;!p.multi=0

end

