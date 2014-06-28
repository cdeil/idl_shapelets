pro shapelets_plot_basis, n_max, $
                          polar=polar, $
                          ps=ps,$
                          color=color,$
                          isotropic=isotropic,$
                          ksb=ksb,$
                          title=title,$
                          real=real,imaginary=imaginary

; Oct 02 - Written by Richard Massey
; Jun 04 - Real/imaginary options added by RM
;
; Makes a plot of all the polar shapelet basis functions out to n_max
;
; SWITCHES:
; /ps   prints to (postscript) file instead of screen
; /colo with /ps makes a colour postscript file
; /ksb  plots only those basis functions which are equivalently used in KSB
; /real plots only the real parts of the basis functions
; /imag plots only the imaginary parts of the basis functions


;
; Parse inputs
;
if keyword_set(real) or keyword_set(imaginary) then polar=1B
if not keyword_set(n_max) then if keyword_set(ksb) then n_max=4 else n_max=6

;
; Set plotting options
;
n_pixels=[40,40]     ; Resolution within each basis function's box
beta=3.6             ; beta [pixels]
n_pixels=[45,45]     ; Resolution within each basis function's box
beta=4.1             ; beta [pixels]
integrate=1          ; Analytically integrate basis functions within pixels?
bits_per_pixel=5     ; Number of greyscale=2^bits
xmargin=!x.margin    ; Remember input margins, so they can be restored later
ymargin=!y.margin    ;
!x.margin=[6,1]      ; Set less conservative margins
!y.margin=[3,1]      ;
col=1B               ; Draw a colour bar
if keyword_set(imaginary) then begin
  ytitle="!8m!6"
  xtitle="!8n!6"
  col=0B
  frame=[-.5,n_max+.5,-n_max-1,n_max+1]
endif else if keyword_set(polar) then begin
  if not keyword_set(title) then title="!6Polar shapelet basis functions"
  ytitle="!8m!6"
  frame=[-.5,n_max+.5,-n_max-1,n_max+1]
endif else begin
  if not keyword_set(title) then title="!6Cartesian shapelet basis functions"
  xtitle="!8n!6!d1!n"
  ytitle="!8n!6!d2!n"
  frame=[-.5,n_max+.5,-.5,n_max+.5]
endelse

;
; Set up basis
;
shapelets_make_nvec, n_max, n1, n2, n_coeffs
nr=n1
nl=n2
n=nr+nl
m=nr-nl
shapelets_make_xarr, n_pixels, x1, x2, x0=x0
x1=x1/beta
x2=x2/beta

;
; Generate images of the shapelet basis functions
;
decomp = {name:"blah",                   $
          type:"decomp",                 $
          history:strarr(1),             $
          x:n_pixels/2.,                 $
          beta:beta,                     $
          n_max:n_max,                   $
          n_coeffs:n_coeffs,             $
          coeffs:fltarr(n_coeffs),       $
          coeffs_error:fltarr(n_coeffs), $
          n_pixels:n_pixels,             $
          sextractor:0B,                 $
          moments:0B,                    $
          n1:n1,                         $
          n2:n2,                         $
          over:1,                        $
          polar:0B,                      $
          integrate:1B,                  $
          chisq:fltarr(2),               $
          skyfit:0}
; Generate images of the Cartesian shapelet basis functions
if not keyword_set(polar) then begin
  basis=fltarr(n_pixels[0],n_pixels[1],n_coeffs)
  for i=0,n_coeffs-1 do begin
    decomp.coeffs[*]=0.
    decomp.coeffs[where(decomp.n1 eq n1[i] and decomp.n2 eq n2[i])]=1
    shapelets_recomp,decomp,recomp
    basis[*,*,i]=recomp
  endfor
; Generate images of the polar shapelet basis functions
endif else begin
  basis=complexarr(n_pixels[0],n_pixels[1],n_coeffs)
  for i=0,n_coeffs-1 do begin
    decomp.coeffs[*]=0.
    shapelets_polar_convert,decomp,/SILENT
    if keyword_set(imaginary) or (keyword_set(polar) and not keyword_set(real) and m[i] lt 0) then begin
      decomp.coeffs[where(decomp.n eq n[i] and decomp.m eq m[i])]=complex(0,-1)
    endif else begin
      decomp.coeffs[where(decomp.n eq n[i] and decomp.m eq m[i])]=complex(1,0)
    endelse
    shapelets_polar_convert,decomp,/SILENT
    shapelets_recomp,decomp,recomp;,/COMPLEX
    basis[*,*,i]=recomp
        
  ;  decomp.coeffs[*]=0.
  ;  shapelets_polar_convert,decomp,/silent
  ;  decomp.coeffs[where(decomp.n eq n[i] and decomp.m eq  m[i])]=complex(0,0.5)
  ;  decomp.coeffs[where(decomp.n eq n[i] and decomp.m eq -m[i])]=complex(0,-0.5)
  ;  shapelets_polar_convert,decomp,/silent
  ;  shapelets_recomp,decomp,recomp
  ;  basis[*,*,i]=basis[*,*,i]+complex(0,1)*recomp
  endfor
endelse

;
; Rescale image from average within pixels to peak value
;
cscale=1./sqrt(!pi)/max(float(basis[*,*,0]))
basis=basis*cscale


;
; Remove all basis functions except those used equivalently in KSB
;
if keyword_set(ksb) then begin
  basis[*,*,1:2]=0
  basis[*,*,6:11]=0
  basis[*,*,13:*]=0
  if not keyword_set(title) then title="!6Shapelet basis functions used by KSB"
endif


;
; Assemble overall image
;
; Slightly wasteful mathod perhaps of pixellating entire plot, rather than just 
; that triangle which will contain data. Still, at least it's then not
; transparent.
xpix=n_pixels[0]*(n_max+1)
ypix=n_pixels[1]*(n_max+1)
image=fltarr(xpix,ypix)
for i=0,n_coeffs-1 do begin

  ; Decide where the basis function should go
  if keyword_set(polar) then begin
    left=n[i]*n_pixels[0]                  & right=left+n_pixels[0]-1
    bottom=(ypix/2)+(m[i]-1)*n_pixels[1]/2 & top=bottom+n_pixels[1]-1
  endif else begin
    left=(n1(i))*n_pixels[0]     & right=left+n_pixels[0]-1
    bottom=(n2(i))*n_pixels[1]   & top=bottom+n_pixels[1]-1
  endelse
  
  ; And place it there
;  if keyword_set(imaginary) or (keyword_set(polar) and not keyword_set(real) and m[i] lt 0) then begin
;    image[left:right,bottom:top]= $
;    image[left:right,bottom:top]+imaginary(basis[*,*,i])
;  endif else begin
    image[left:right,bottom:top]+=float(basis[*,*,i])
;  endelse

endfor

; Open postscript file if requested
csize=0.1


if keyword_set(ps) then begin
  csize=0.1
  set_plot,'ps'
  if keyword_set(imaginary) then begin
    device,filename='polar_basis_imaginary.eps',yoffset=5.,ysize=20.*(1-csize),/encap,bits=bits_per_pixel
  endif else if keyword_set(real) then begin
    device,filename='polar_basis_real.eps',yoffset=5.,ysize=20.,/encap,bits=bits_per_pixel
  endif else if keyword_set(polar) then begin
    device,filename='polar_basis.eps',yoffset=5.,ysize=20.,/encap,bits=bits_per_pixel
  endif else begin
    device,filename='Cartesian_basis.eps',yoffset=5.,ysize=20.,/encap,bits=bits_per_pixel
  endelse
  opub
endif
;loadct,1,/silent


; Set colour scale
crange=[min(image),-1*min(image)]
crange=max(abs(image))*[-1,1]
;crange=[-0.12,0.12]
crange=[-1,1]/sqrt(!pi)



;
; Plot basis functions
;
shapelets_plot_image,image,col=1-keyword_set(imaginary),$
  cs=csize,crange=crange,frame=frame,isotropic=isotropic, $
  xtitle=xtitle,ytitle=ytitle,title=title



; Plot boxes around functions
if keyword_set(polar) then begin
  for i=0,n_max do begin
    oplot,[i-.5,i-.5],[-1*i-1,i+1],psym=-3
    for j=0,i+1 do begin
      oplot,[i-.5,i+.5],[-1*i-1+j*2,-1*i-1+j*2],psym=-3
    endfor
  endfor
  oplot,[n_max+.5,n_max+.5],[-1*n_max-1,n_max+1],psym=-3
endif else begin
  for i=0,n_max-1 do begin
    oplot,[i+.5,i+.5],[-.5,n_max-i+.5],psym=-3
  endfor
  for i=0,n_max-1 do begin
    oplot,[-.5,n_max-i+.5],[i+.5,i+.5],psym=-3
  endfor
endelse

; Print labels
if keyword_set(real) then begin
  xyouts,-0.3,n_max+0.1,"!6Real components",charsize=2.25
endif else if keyword_set(imaginary) then begin
  xyouts,-0.3,n_max+0.1,"!6Imaginary components",charsize=2.25
endif


; Close ps file if appropriate
if keyword_set(ps) then cps


; Restore margins
!x.margin=xmargin
!y.margin=ymargin






end
