pro plt_shear_scat,scat,class_c=class_c,mag_ran=mi_ran

; August 998 - Written by A. Refregier
;
; PURPOSE: plot shear statistics derived from sextractor catalog. This
; is particularly useful if the sextractor catalog is is scat_in input catalog
; of a simulation
; INPUT: scat: catalog read by shapelet_read_scat.pro (esp. scat_in:
;              input catalog)
; KEYWORDS: class_c: class cut for star/gal separation (default:.8)
;           mag_ran: magnitude range
; OUTPUT: plots

COMPILE_OPT OBSOLETE

; declarations
if not keyword_set(class_c) then class_c=.8    ; class cut for galaxies
print,'Class cut:',class_c
;a_psf=1.57                  ; psf rms radius, assumed circular [pixels]
mag_ran = [min(scat.mag),max(scat.mag)] ; magnitude range to consider
if keyword_set(mi_ran) then mag_ran=mi_ran
e_max=1.                     ; maximum ellipticity to include

; compute deconvolved ellipticities
alpha = scat.theta/!radeg    ; convert PA to radians
ell=scat.e
e1=scat.e1
e2=scat.e2

; select galaxies
ga=where(scat.class le class_c and ell le e_max and $
      scat.mag ge mag_ran(0) and scat.mag lt mag_ran(1),n_ga)    ; all galaxies

; plot mag-size distribution for galaxies
;!p.multi=[0,1,2]
plot,scat.mag(ga),scat.fwhm(ga),psym=3,/ytype,xran=mag_ran, $
  xtitle='m!ii!n',ytitle='d (pixels)'

; plot ellipticities
read,test
!p.multi = 0
plot,e1(ga),e2(ga),psym=3,$
  xtitle='e1',ytitle='e2'
oplot, [-5., 5.], [0., 0.], lines=1
oplot, [0., 0.], [-5., 5.], lines=1
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

; plot shear vs mag
read,test
!p.multi=[0,1,2]
plot,scat.mag(ga),g1(ga),psym=3,$
  xtitle='mag',ytitle='gamma1'
oplot,[-100,100],[0,0],lines=2
plot,scat.mag(ga),g2(ga),psym=3,$
  xtitle='mag',ytitle='gamma2'
oplot,[-100,100],[0,0],lines=2
!p.multi=0

end

