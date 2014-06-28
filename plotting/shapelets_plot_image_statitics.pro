pro shapelets_plot_image_statitics, image


;+
; NAME:
;      SHAPELETS_PLOT_IMAGE_STATITICS
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Plot statistics about an image.
;
; INPUTS:
;      IMAGE - Shapelets image structure.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      Plots to STDOUT.
;
; MODIFICATION HISTORY:
;      Mar 04 - Written by Alexandre Refregier.
;-

; Initial declarations
im=image.image
tol_rms=.05          ; relative tolerance for rms determination

; Compute raw statistics
print,'min, max:',min(im),max(im)
ave=mean(im)
rms=stddev(im)
;mode=2.5*median(im)-1.5*ave
;print,'mean, rms, mode:',ave,rms,mode
print,'mean, sigma:',ave,rms

; Compute approximate mean and rms of background by iterating
; using an alogrythm similar to that used by SEXtractor
conv=0                ; converged?
niter=0
while not conv do begin
  gd=where(abs(im-ave) lt 3.*rms)
  ave_new=mean(im(gd))
  rms_new=stddev(im(gd))
;  mode_new=2.5*median(im(gd))-1.5*ave_new
  if abs(rms_new-rms)/rms lt tol_rms then conv=1
  ave=ave_new
  rms=rms_new
; mode=mode_new
  niter=niter+1
;print,'mean, rms, mode (<3rms):',ave,rms,mode
;print,'mean, sigma, (<3rms):',ave,rms
endwhile
print,'mean, sigma (<3rms), n_iter:',ave,rms,niter

; Compare to noise map
if keyword_set(image.noise) then begin
  noise_rms=mean(1./sqrt(image.noise(gd)))
  print,'noise map: sigma(<3sigma):',noise_rms
endif

; Plot pixel histogram
plothist,image.image<20.*rms,bin=rms/20.,xran=[-3.,10.]*rms+ave,$
  xtitle='pixel value',ytitle='No. of pixels'
oplot,ave*[1.,1.],[0,1e9]
oplot,ave+rms*[1,1],[0,1e9],lines=2
oplot,ave-rms*[1,1],[0,1e9],lines=2
;oplot,ave+3.*rms*[1,1],[0,1e9],lines=1
;oplot,ave-3.*rms*[1,1],[0,1e9],lines=1
;legend,['mean','rms','3rms'],lines=[0,2,1],box=0


end
