pro plt_scat_mom,scat_mom

; April 04 - Written by A. Refregier
;
; PURPOSE: plot statistics for the scat_mom catalog. In particular,
; compare derived moments from the image to the input catalog
; INPUT: scat_mom: from mk_scat_mom.pro
; OUTPUT: plot

COMPILE_OPT OBSOLETE

; declarations
gd=where(scat_mom.flag eq 0)     ; select sources without problems in postage stamp
gd2=where(scat_mom.flag le 2)    ; include sources with unproblematic neighbors 
orig=!p.color
col=[orig,4,2]  
tek_color 

; plot flag statistics
plothist,scat_mom.flag,bin=.5,$
  xtitle='flag',ytitle='No of objects'

; plot fluxes
read,test
plot,scat_mom.flux_orig(gd),scat_mom.flux(gd),psym=6,$
  /xtype,/ytype,/nodata,$
  xtitle='flux original',ytitle='flux moments'
oplot,scat_mom.flux_orig,scat_mom.flux,psym=3,col=col(0)
oplot,scat_mom.flux_orig(gd2),scat_mom.flux(gd2),psym=1,col=col(1)
oplot,scat_mom.flux_orig(gd),scat_mom.flux(gd),psym=6,col=col(2)
oplot,[1e-4,1e8],[1e-4,1e8]

; plot radii
read,test
r2_orig=sqrt((scat_mom.a_orig^2+scat_mom.b_orig^2))
r2=sqrt(scat_mom.a^2+scat_mom.b^2)
plot,r2_orig(gd),r2(gd),$
  /xtype,/ytype,/nodata,$
  xtitle='R original (pixels)',ytitle='R moments (pixels)'
oplot,r2_orig,r2,psym=3
oplot,r2_orig(gd2),r2(gd2),psym=1,col=col(1)
oplot,r2_orig(gd),r2(gd),psym=6,col=col(2)
oplot,[1e-4,1e8],[1e-4,1e8]

; plot ellipticities
read,test
!p.multi=[0,1,2]
plot,scat_mom.e1_orig(gd),scat_mom.e1(gd),/nodata,$
  xran=[-1,1],yran=[-1,1],$
  xtitle='e1 original',ytitle='e1 moments'
oplot,scat_mom.e1_orig,scat_mom.e1,psym=3
oplot,scat_mom.e1_orig(gd2),scat_mom.e1(gd2),psym=1,col=col(1)
oplot,scat_mom.e1_orig(gd),scat_mom.e1(gd),psym=6,col=col(2)
oplot,[-10,10],[-10,10]

plot,scat_mom.e2_orig(gd),scat_mom.e2(gd),/nodata,$
  xran=[-1,1],yran=[-1,1],$
  xtitle='e2 original',ytitle='e2 moments'
oplot,scat_mom.e2_orig,scat_mom.e2,psym=3
oplot,scat_mom.e2_orig(gd2),scat_mom.e2(gd2),psym=1,col=col(1)
oplot,scat_mom.e2_orig(gd),scat_mom.e2(gd),psym=6,col=col(2)
oplot,[-10,10],[-10,10]
!p.multi=0

end
