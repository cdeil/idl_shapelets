pro plt_stars,cat,simple=simple,class_c=class_c,mi_ran=mi_ran,telescope=telescope


;     Oct 05 - size of pixels automatically read from SExtractor
;              configuration file added by JB
;     Nov 04 - plot of mag-class and of flux_radius added by Joel Berge
; August 998 - Written by A. Refregier
;
; PURPOSE: plot various figures relative to the stars in a catalog.
; This is a measure of the PSF.
; INPUT: cat: catalog read by rd_cat.pro
; OPTIONAL INPUT: simple: simplify by only plotting stellar ellipticity
;      class_c: cut for star/galaxy classifier (default=0.8)
;      mi_ran: magnitude range for selected stars (default=[-15,-11])
;      telescope - telescope used for image, default='SUBARU'
; OUTPUT: plot
  
; declarations
if not keyword_set(telescope) then telescope='SUBARU'
if not keyword_set(class_c) then class_c = 0.8  ; class cut for stars
if not keyword_set(mi_ran) then mi_ran = [-15., -10.] 
                              ; magnitude range for useful stars
th_pix=extract_pixsc(telescope)                    ; pixel size [arcsec]
npix_smooth=400               ; smoothing scale for PSF smoothing [pixels]
print,th_pix
radius_d = sqrt((cat.a^2+cat.b^2)/2.)     ; radius
x_ran=[min(cat.x(*,0)),max(cat.x(*,0))]
y_ran=[min(cat.x(*,1)),max(cat.x(*,1))]
mag_ran=[min(cat.mag)-1,max(cat.mag)+1]

; select useful (unsaturated bright) stars
st=where(cat.class gt class_c and cat.mag gt mi_ran(0) and $
   cat.mag lt mi_ran(1),ns)
print,'No. of stars selected: ',ns


;plot mag-flux_radius relation
plot,cat.mag,cat.flux_radius*1.5*th_pix,xtitle='mag',ytitle='1.5*flux_radius (arcsec)',psym=3, $
  title='mag-1.5*flux_radius',yrange=[0,2],xrange=mag_ran;[13,40]
oplot,cat.mag(st),cat.flux_radius(st)*1.5*th_pix,psym=5
read,test


; plot mag-size relation
if not keyword_set(simple) then begin
plot, cat.mag, cat.fwhm*th_pix, psym=3, yrange=[.1,1000],xrange=mag_ran,$;[13,28],$
   xtitle='mag', ytitle='FWHM (arcsec)',/ytype,title='mag-fwhm relation'
oplot,cat.mag(st), cat.fwhm(st)*th_pix,psym=5
xyouts, /normal, .2, .9, 'No. of stars:'+string(ns)
xyouts, /normal, .2, .86, 'Class cut:'+string(class_c)
xyouts, /normal, .2, .82, 'mag range:'+string(mi_ran(0), mi_ran(1))


;plot mag-class relation
read,test
plot, cat.mag, cat.class, psym=3, yrange=[.1,1000],xrange=mag_ran,$;[13,28],$
   xtitle='mag', ytitle='class',/ytype,title='mag-class relation'
oplot,cat.mag(st), cat.class(st),psym=5
xyouts, /normal, .2, .9, 'No. of stars:'+string(ns)
xyouts, /normal, .2, .86, 'Class cut:'+string(class_c)
xyouts, /normal, .2, .82, 'mag range:'+string(mi_ran(0), mi_ran(1))


; plot size distribution
read, test
plothist,cat.fwhm(st)*th_pix,bin=.1,$
  xtitle='FWHM (arcsec)',ytitle='No. of stars',title='fwhm distribution'
fwhm_mean = mean(cat.fwhm(st))
xyouts, /normal, .2, .9, 'mean FWHM (arcsec)'+$
    string(fwhm_mean*th_pix)
print,'mean FWHM (pixels, arcsec):', fwhm_mean,fwhm_mean*th_pix


;plot 1.5*flux_radius distribution
read, test
plothist,1.5*cat.flux_radius(st)*th_pix,bin=.02,$
  xtitle='1.5*FLUX_RADIUS (arcsec)',ytitle='No. of stars',title='1.5*flux_radius distribution'
fr_mean = mean(1.5*cat.flux_radius(st))
xyouts, /normal, .2, .9, 'mean 1.5*FLUX_RADIUS (arcsec)'+$
    string(fr_mean*th_pix)
print,'mean 1.5*FLUX_RADIUS (pixels, arcsec):', fr_mean,fr_mean*th_pix


; plot ellipticities
read,test
endif
plot,cat.e1(st),cat.e2(st),xrange=[-1., 1.], yrange=[-1., 1.], psym=3,$
  xtitle='e1',ytitle='e2'
oplot, [-5., 5.], [0., 0.], lines=1
oplot, [0., 0.], [-5., 5.], lines=1
e1_m = mean(cat.e1(st)) &  e2_m=mean(cat.e2(st))
e1_sig = stdev(cat.e1(st)) & e2_sig = stdev(cat.e2(st))
e1_err = e1_sig/sqrt(ns) &  e2_err= e2_sig/sqrt(ns) 
xyouts, /normal, .2, .9, 'e1,e2: mean:'+string(e1_m, e2_m)
xyouts, /normal, .2, .86, 'e1,e2: error:'+ string(e1_err, e2_err)
xyouts, /normal, .2, .82, 'e1,e2: rms:'+ string(e1_sig, e2_sig)
;ploterr_xy, e1_m, e2_m, e1_err, e2_err
print,'e1,e2:'
print,'  mean:',e1_m, e2_m
print,'  error:',e1_err, e2_err
print,'  sigma:',e1_sig, e2_sig

; plot ellipticities as a function of position on the chip
read,test
xy_ran = [min(cat.x), max(cat.x)]
plot,cat.x(st,0),cat.x(st,1),psym=3,xran=xy_ran, yran=xy_ran, /xstyle, /ystyle, $
  xtitle='x (pixels)',ytitle='y (pixels)'
scale=300.
plt_evec,cat.x(st,0),cat.x(st,1),sqrt(cat.e1(st)^2+cat.e2(st)^2),$
  cat.theta(st)-!pi/2.,xs=scale,ys=scale

; compute and plot smoothed PSF ellipticity map
n_xs=fix((x_ran(1)-x_ran(0))/npix_smooth)
n_ys=fix((y_ran(1)-y_ran(0))/npix_smooth)
xs=(indgen(n_xs)+.5)*npix_smooth+x_ran(0)
ys=(indgen(n_ys)+.5)*npix_smooth+y_ran(0)
e1_s=fltarr(n_xs,n_ys) & e2_s=fltarr(n_xs,n_ys)
n_gd=fltarr(n_xs,n_ys)
for i=0l,ns-1 do begin
  is=fix((cat.x(st(i),0)-x_ran(0))/npix_smooth)
  js=fix((cat.x(st(i),1)-y_ran(0))/npix_smooth)
  if is ge 0 and is lt n_xs and js ge 0 and js lt n_ys then begin
    n_gd(is,js)=n_gd(is,js)+1l
    e1_s(is,js)=e1_s(is,js)+cat.e1(st(i))
    e2_s(is,js)=e2_s(is,js)+cat.e2(st(i))
  endif
endfor
e1_s=e1_s/float(n_gd)
e2_s=e2_s/float(n_gd)
read,test
plot,[0],[0],psym=3,xran=xy_ran, yran=xy_ran, /xstyle, /ystyle, $
  xtitle='x (pixels)',ytitle='y (pixels)',/nodata
scale=300.
e_scale=npix_smooth/.5
for i=0,n_xs-1 do begin
  for j=0,n_ys-1 do begin
    if n_gd(i,j) gt 0 then begin
      alpha=atan(e2_s(i,j),e1_s(i,j))/2.
      e=sqrt(e1_s(i,j)^2+e2_s(i,j)^2)
      oplot,xs(i)+e*cos(alpha)*e_scale*[-1.,1.],$
           ys(j)+e*sin(alpha)*e_scale*[-1.,1.]
    endif
  endfor
endfor
gds=where(n_gd gt 0,n_gds)
e1s_m = mean(e1_s(gds)) &  e2s_m=mean(e2_s(gds))
e1s_sig = stdev(e1_s(gds)) & e2s_sig = stdev(e2_s(gds))
e1s_err = e1s_sig/sqrt(n_gds) &  e2s_err= e2_sig/sqrt(n_gds) 
print,'smoothed e1,e2:'
print,'  mean:',e1s_m, e2s_m
print,'  error:',e1s_err, e2s_err
print,'  sigma:',e1s_sig, e2s_sig

; plot PSF ellipse as function of position on the chip
if not keyword_set(simple) then begin
read, test
plot,cat.x(st,0),cat.x(st,1),psym=3,xran=xy_ran, yran=xy_ran, /xstyle, /ystyle, $
  xtitle='x (pixels)',ytitle='y (pixels)', /nodata
plt_ellipse,cat.x(st,0),cat.x(st,1), cat.a(st), cat.b(st), cat.theta(st)-!pi/2., $
  xs=30., ys=30.



; plot psf size vs x and y
read, test
!p.multi = [0, 1, 2]
d=sqrt((cat.a*cat.a+cat.b*cat.b)/2)
plot, cat.x(st,0), d(st), psym=3, yrange=[0.,4.],$
  xtitle='x (pixels)',ytitle='d (pixels)'
;slice_fit, cat.x(st,0), d(st),  10,  5, x_sl, d_sl
;oplot, x_sl, d_sl
plot, cat.x(st,1), d(st), psym=3,  yrange=[0.,4.],$
  xtitle='y (pixels)',ytitle='d (pixels)'
;slice_fit, cat.y(st), d(st),  10,  5, y_sl, d_sl
;oplot, y_sl, d_sl
!p.multi = 0


; plot psf 1.5*flux_radius vs x and y
read, test
!p.multi = [0, 1, 2]
plot, cat.x(st,0), 1.5*cat.flux_radius(st)*th_pix, psym=3, yrange=[0.,2.],$
  xtitle='x (pixels)',ytitle='1.5*flux_radius (arcsec)'
;slice_fit, cat.x(st,0), d(st),  10,  5, x_sl, d_sl
;oplot, x_sl, d_sl
plot, cat.x(st,1), 1.5*cat.flux_radius(st)*th_pix, psym=3,  yrange=[0.,2.],$
  xtitle='y (pixels)',ytitle='1.5*flux_radius (arcsec)'
;slice_fit, cat.y(st), d(st),  10,  5, y_sl, d_sl
;oplot, y_sl, d_sl
!p.multi = 0


endif

end
