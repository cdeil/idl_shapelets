pro plt_sgsep,scat,class_c=class_c,mag_ran=mag_ran,telescope=telescope


;    Oct 05 - size of pixels automatically read from SExtractor
;              configuration file added by JB
;    Nov 04 - plots of flux_radius added by Joel Berge
; august 98 - written by A. Refregier
;
; PURPOSE: plot the various figures relative to the star-galaxy
; separation for a catalog of objetcs
; INPUT: scat: catalog read by shapelets_read_scat.pro
; KEYWORDS: mag_ran: magnitude range (default:min,max of magnitudes)
;           class_c: classifier cut off for star/gal separation (default:0.8)
;      telescope - telescope used for image, default='SUBARU'
; OUTPUT: plots

; declarations
if not keyword_set(telescope) then telescope='SUBARU'
n=n_elements(scat.mag)
dc=.02
;mag_ran = [-20., -6.]    ; magnitude range to plot
if not keyword_set(mag_ran) then mag_ran=[min(scat.mag),max(scat.mag)]
if not keyword_set(class_c) then class_c=0.8 ; class cut for separation

th_pix=extract_pixsc(telescope)  ;size of a pixel in arcsec


; compute derived quantities
alpha = scat.theta/!radeg    ; convert PA to radians
e = (scat.a^2-scat.b^2)/(scat.a^2+scat.b^2)  ; ellipticity
e1=e*cos(2.*alpha)
e2=e*sin(2.*alpha)
d = sqrt((scat.a^2+scat.b^2)/2.)     ; radius


; plot class vs mag
plot,scat.mag,scat.class,psym=3,xrange=mag_ran, /xstyle,$
  xtitle='mag',ytitle='class'

; plot mag-size relation with different colors for different
; classification index
read,test
plt_colbar,[1.+dc,0.],title='class',csize=.2
plot,scat.mag,d,psym=3,/ytype,/xstyle,xrange=mag_ran, /noerase,/nodata,$
  xtitle='mag',ytitle='d (pixels)'

for i=0l,n-1 do oplot,[scat.mag(i)],[d(i)],psym=1,symsize=.5,$
  color=(1.+dc-scat.class(i))/(1.+dc)*255.

; plot mag-size relation for stars, and galaxies separately
st=where(scat.class gt class_c, n_st)
ga=where(scat.class le class_c, n_ga)
read,test
!p.region=0
!p.multi=[0,1,2]

d_ran = [min(d), max(d)]
print, 'Class cut: ', class_c
print, 'No. of galaxies, stars: ', n_ga, n_st
plot,scat.mag(ga),d(ga),psym=3,/ytype,xrange=mag_ran,yrange=d_ran,/xstyle,$
  xtitle='mag',ytitle='d (pixels)',title='galaxies - class < '+$
  string(class_c, format='(f4.2)')
xyouts, /normal, .2, .9, 'No. of galaxies:'+string(n_ga)

plot,scat.mag(st),d(st),psym=3,/ytype,xrange=mag_ran,yrange=d_ran,/xstyle,$
  xtitle='mag',ytitle='d (pixels)',title='stars - class > '+$
  string(class_c, format='(f4.2)')
xyouts, /normal, .2, .4, 'No. of stars:'+string(n_st)

!p.multi=0


;plot mag-1.5*flux_radius relation for stars and galaxies with different colors
read,test
plt_colbar,[1.+dc,0.],title='class',csize=.2
plot,scat.mag,1.5*scat.flux_radius*th_pix,psym=3,/ytype,/xstyle,xrange=mag_ran, /noerase,/nodata,$
  xtitle='mag',ytitle='1.5*flux_radius (arcsec)'

for i=0l,n-1 do oplot,[scat.mag(i)],[1.5*scat.flux_radius(i)*th_pix],psym=1,symsize=.5,$
  color=(1.+dc-scat.class(i))/(1.+dc)*255.


; plot mag-1.5*flux_radius relation for stars, and galaxies separately
st=where(scat.class gt class_c, n_st)
ga=where(scat.class le class_c, n_ga)
read,test
!p.region=0
!p.multi=[0,1,2]

fr_ran = [min(1.5*scat.flux_radius*th_pix), max(1.5*scat.flux_radius*th_pix)]
print, 'Class cut: ', class_c
print, 'No. of galaxies, stars: ', n_ga, n_st
plot,scat.mag(ga),1.5*scat.flux_radius*th_pix(ga),psym=3,/ytype,xrange=mag_ran,yrange=fr_ran,/xstyle,$
  xtitle='mag',ytitle='1.5*flux_radius (arcsec)',title='galaxies - class < '+$
  string(class_c, format='(f4.2)')
xyouts, /normal, .2, .9, 'No. of galaxies:'+string(n_ga)

plot,scat.mag(st),1.5*scat.flux_radius*th_pix(st),psym=3,/ytype,xrange=mag_ran,yrange=fr_ran,/xstyle,$
  xtitle='mag',ytitle='1.5*flux_radius (arcsec)',title='stars - class > '+$
  string(class_c, format='(f4.2)')
xyouts, /normal, .2, .4, 'No. of stars:'+string(n_st)

!p.multi=0


; plot positions
read,test
plot,scat.x(ga,0),scat.x(ga,1),psym=3,$
  xtitle='x (pixels)',ytitle='y (arcmin)'
oplot,scat.x(st,0),scat.x(st,1),psym=1,sym=.7

end


