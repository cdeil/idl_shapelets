pro shapelets_plot_shapecat,shcat

; April 04 - Written by A. Refregier
;
; PURPOSE: plot statistics relative to shapelet catalog
; INPUT: shcat: shapelet catalog from shapelets_read_shcat.pro
; OUTPUT plots

; declarations
n=shcat.n
gd=where(shcat.flag(*,0) eq 0 and shcat.flag(*,1) le 1$
         and finite(shcat.chisq),n_gd)
print,'No objects: total, not flagged:',n,n_gd

; plot size-magnitude distribution
!p.region=0
!p.multi=[0,1,2,0,0]
plot,shcat.mag,sqrt(shcat.rsquared),/ytype,psym=3,$
  xtitle='mag',ytitle='R [pixels]'
oplot,shcat.mag(gd),sqrt(shcat.rsquared(gd)),psym=6

; plot distribution of Chi^2
gd2=where(finite(shcat.chisq))
plothist,shcat.chisq(gd2),xran=[0,3.],bin=.01,$
  xtitle='chi^2',ytitle='No. of objects'
plothist,shcat.chisq(gd),bin=.01,/over,lines=2
oplot,[1,1],[0,10000],lines=2
read,junk

; plot distribution of flags
wait=0B
if n_elements(uniq(shcat.flag[*,0])) gt 1 then begin
  !p.multi=[0,1,2,0,0]
  plothist,shcat.flag[*,0],bin=.5,$
    xtitle='Flag - postage stamp',ytitle='no. of objects',$
    title='0: OK, 2: neighbor, 3: severe overlap, 4: edge, 5: FWHM=0, '+$
          '6: background error, 7: masked out'
  wait=1B
endif else !p.multi=[0,1,1,0,0]
if n_elements(uniq(shcat.flag[*,1])) gt 1 then begin
  plothist,shcat.flag[*,1],bin=.5,$
    xtitle='Flag - focus',ytitle='no. of objects',$
    title='0: OK, 1: bounced, 2: nmax_max, 3: flatness, 4: no converge, '+$
          '5: centroid error, 6: fatal'
  wait=1B
endif
if wait then read,junk

; plot distribution of nmax and beta
!p.multi=[0,1,2,0,0]
plothist,shcat.n_max,$
  xtitle='n!dmax!n',ytitle='No. of objects'
plothist,shcat.n_max[gd],/OVER,LINESTLYE=2,/FILL
plothist,shcat.beta,bin=.1,$
  xtitle='beta (pixels)',ytitle='No of objects'
plothist,shcat.beta[gd],bin=.1,/OVER,LINESTLYE=2,/FILL
read,junk
!p.multi=0

; plot theta_min vs theta_max
;read,test
;th_min=shcat.beta*sqrt(shcat.n_max+1.)
;th_max=shcat.beta/sqrt(shcat.n_max+1.)
;plot,shcat.mag,sqrt(shcat.rsquared),/ytype,psym=3,$
;  xtitle='mag',ytitle='theta (pixels)',xran=[0,10]
;for i=0,n-1,10 do oplot,shcat.mag(i)*[1.,1.],[th_min(i),th_max(i)]
;plot,th_min,th_max,psym=3,$
;  xtitle='theta_min (pixels)',ytitle='theta_max (pixels)'

; plot centroid offset
dx=shcat.centroid(*,0)-shcat.sexx(*)
dy=shcat.centroid(*,1)-shcat.sexy(*)
plot,dx,dy,xran=[-5,5],yran=[-5,5],$
  xtitle='x-x_sex (pixels)',ytitle='y-y_sex (pixels)',psym=3
oplot,[-100,100],[0,0],lines=2
oplot,[0,0],[-100,100],lines=2
oplot,dx(gd),dy(gd),psym=6

dr=sqrt(dx^2+dy^2)
;gd=where(dr lt 5.,n_gd)
print,'delta x, y (pixels): mean:',amean(dx(gd)),amean(dy(gd))
print,'delta x, y (pixels): error:',asigma(dx(gd))/sqrt(n_gd),$
       asigma(dy)/sqrt(n_gd)
print,'delta x, y (pixels): sigma:',asigma(dx(gd)),asigma(dy(gd))

; plot background statistics
bin_back=.005
read,test
!p.multi=[0,1,2]
plothist,shcat.back_mean_local,bin=bin_back,$
  xtitle='background mean',ytitle='No of objects'
plothist,shcat.back_mean_local(gd),/over,lines=2,bin=bin_back
;plothist,shcat.back_mean_ext,/over,line=2
;legend,['local','external'],lines=[0,2],box=0
plothist,shcat.back_rms_local,bin=bin_back,$
  xtitle='background rms',ytitle='No of objects'
plothist,shcat.back_rms_local(gd),/over,lines=2,bin=bin_back
;plothist,shcat.back_rms_ext,/over,line=2
!p.multi=0

end
