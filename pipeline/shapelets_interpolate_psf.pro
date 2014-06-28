function shapelets_interpolate_psf, shapecat,         $
                                    XINTERPOLATE=x,   $
                                    YINTERPOLATE=y,   $
                                    XRANGE=xrange,    $
                                    YRANGE=yrange,    $
                                    SUBSET=subset,    $
                                    PCA=pca,          $
                                    BASIS=basis,      $
                                    MAX_P=max_p,      $
                                    DEGREE=degree,    $
                                    PLOTIT=plotit,    $
				    MODEL=model,      $
                                    N_NEIGHBOURS=n_neighbours
COMPILE_OPT idl2

; Set display options
escale=500
ekey=20
if not keyword_set(xrange) then xrange=[min(shapecat.x[*,0]),max(shapecat.x[*,0])]
if not keyword_set(yrange) then yrange=[min(shapecat.x[*,1]),max(shapecat.x[*,1])]
if keyword_set(x) and keyword_set(y) then n_galaxies=n_elements(x)<n_elements(y)
!P.region=0
!P.multi=0
if n_elements(degree) eq 0 then degree=3
test_object=121<(shapecat.n-1)
n_max=shapecat.n_max[uniq(shapecat.n_max)]
beta=shapecat.beta[uniq(shapecat.beta)]
if n_elements(n_max)>n_elements(beta) gt 1 then message,"Shapecat must have a constant n_max and beta!"
if shapecat.polar then message,"Not yet sure how to handle polar shapelet catalogues!"


; Renormalise each star to unit flux (and simultaneously transpose the matrices)
coeffs=dblarr(shapecat.n_coeffs,shapecat.n)
coeffs_error=dblarr(shapecat.n_coeffs,shapecat.n)
for i=0,shapecat.n-1 do begin
  coeffs[*,i]=reform(shapecat.coeffs[i,*],shapecat.n_coeffs,1)/shapecat.flux[i]
  ; ignore error on the flux
  coeffs_error[*,i]=reform(shapecat.coeffs_error[i,*],shapecat.n_coeffs,1)/shapecat.flux[i]
endfor



; Subtract off average star shape
coeffs_average=total(coeffs,2)/float(shapecat.n)
for i=0,shapecat.n_coeffs-1 do coeffs[i,*]=coeffs[i,*]-coeffs_average[i]
decomp=shapelets_shapecat2decomp(shapecat,0)
decomp.coeffs=coeffs_average
print,"Mean ellipticity is ",shapelets_ellipticity(decomp)



; Make some pretty plots for the delectation of the user
if keyword_set(plotit) then begin
  ; Plot ellipticities
  plot,[0,0],/nodata,xstyle=2,ystyle=2,xrange=xrange,yrange=yrange,$
       title="!6Uncorrected stellar ellipticities",$
       xtitle="!6x [pixels]",ytitle="!6y [pixels]",/iso
  for i=0,shapecat.n-1 do plt_evec,shapecat.x[i,0],shapecat.x[i,1],/e1e2,$
     shapecat.ellipticity[i,0],shapecat.ellipticity[i,1],xscale=escale,yscale=escale
endif



; Find Principal Components of the variation about this average
nvariables=shapecat.n_coeffs
if keyword_set(pca) then begin
  new_coeffs=pcomp(coeffs,                    $
                   coefficients=coefficients, $
                   eigenvalues=eigenvalues,   $
                   variances=percentages,     $
                   nvariables=nvariables)
  ;new_coeffs=transpose(coefficients)#coeffs
  if nvariables lt shapecat.n_coeffs then new_coeffs=[new_coeffs,fltarr(shapecat.n_coeffs-nvariables,shapecat.n)]
  eigenvectors=coefficients/rebin(eigenvalues,shapecat.n_coeffs,shapecat.n_coeffs)
  eigenvectors[where(finite(eigenvectors)) eq 0]=0
  n_significant=n_elements(where(eigenvalues gt 0))
  if pca ne 1 then begin
    cumulative_percentages=fltarr(shapecat.n_coeffs,/nozero)
    for i=0,shapecat.n_coeffs-1 do cumulative_percentages[i]=total(percentages[0:i])
    n_significant=n_elements(where(cumulative_percentages le pca))
  endif
  new_coeffs_error=transpose(coefficients)#coeffs_error
endif else begin
  new_coeffs=coeffs
  new_coeffs_error=coeffs_error
  n_significant=shapecat.n_coeffs
endelse



; Make some pretty plots for the delectation of the user
if keyword_set(plotit) then begin

  ; Draw average star
  read,junk
  decomp=shapelets_shapecat2decomp(shapecat,0)
  decomp.coeffs=coeffs_average
  decomp.n_pixels=50
  decomp.x=decomp.n_pixels/2.
  decomp.beta=decomp.beta*3
  !P.region=0
  !P.multi=[0,3,4]
  shapelets_recomp,decomp,recomp
  shapelets_plot,recomp,cran=cran
  xyouts,0,0,/data,"!6Average star"

  ; Draw principal components
  for i=0,10<(size(new_coeffs))[1] do begin
    decomp.coeffs=reform(coefficients[i,*])
    shapelets_recomp,decomp,recomp
    recomp=recomp*sign(recomp[round(decomp.x[0]),round(decomp.x[1])])
    shapelets_plot,recomp,cran=cran
    xyouts,0,0,/data,"!6PC"+strtrim(string(i),2)
    print,"First "+strtrim(string(i+1),2)+" PCs together contain "+$
          strtrim(string(total(percentages[0:i])*100),2)+"% of the variance"
  endfor
  read,junk
  
  ; Plot the coefficients (in terms of PCs) as a function of position
  usersym, cos(2*!pi*findgen(21)/20), sin(2*!pi*findgen(21)/20),/fill
  for pc=0,10 do begin
    ; Plot using coloured or resized symbols
    !P.multi=0
    plot,[0,0],/nodata,/xstyle,/ystyle,title="!6PC #"+strtrim(string(pc),2),$
         xrange=[min(shapecat.x[*,0]),max(shapecat.x[*,0])],$
         yrange=[min(shapecat.x[*,1]),max(shapecat.x[*,1])],$
         xtitle="!6x [pixels]",ytitle="!6y [pixels]",/iso
    for i=0,shapecat.n-1 do begin
      ;oplot,shapecat.x[i,0]*[1,1],shapecat.x[i,1]*[1,1],$
      ;   symsize=new_coeffs[3,i]/mean(new_coeffs[3,*]),psym=8
      oplot,shapecat.x[i,0]*[1,1],shapecat.x[i,1]*[1,1],symsize=2,psym=8,$
         col=255.*(new_coeffs[pc,i]-min(new_coeffs[pc,*]))/(max(new_coeffs[pc,*])-min(new_coeffs[pc,*]))
    endfor
    
    ; Plot separately as a function of x and as a function of y
    !P.multi=[0,1,2]
    plot,shapecat.x[*,0],new_coeffs[pc,*],psym=1,$
         title="!6PC #"+strtrim(string(pc),2),$
         xtitle="!6x [pixels]",ytitle="!6PC #"+strtrim(string(pc),2)
    plot,shapecat.x[*,1],new_coeffs[pc,*],psym=1,$
         xtitle="!6y [pixels]",ytitle="!6PC #"+strtrim(string(pc),2)
    ;read,junk
  endfor
endif



; Check PCA worked by reconstructing the coefficient array
if keyword_set(pca) then begin

  coeffs_pca=dblarr(shapecat.n_coeffs,shapecat.n)
  for i=0,shapecat.n_coeffs-1 do for j=0,shapecat.n-1 do begin
    coeffs_pca[i,j]=total(coefficients[i,0:n_significant-1]*new_coeffs[0:n_significant-1,j]/eigenvalues[0:n_significant-1])
  endfor


  if keyword_set(plotit) then begin
    print,'PCA reconstruction error:',total((coeffs_pca-coeffs)^2)
    plot,[0,0],/nodata,xrange=.1*[-1,1],yrange=.1*[-1,1],/iso,/xstyle,/ystyle,$
         title="!6PCA reconstruction",$
         xtitle="!6True shapelet coefficient",ytitle="!6Shapelet coefficient after PCA"
    oplot,[0,0],[-100,100]
    oplot,[-100,100],[0,0]
    oplot,[-100,100],[-100,100],/linestyle
    for i=0,shapecat.n_coeffs-1 do oplot,coeffs[i,*],coeffs_pca[i,*],psym=3
    read,junk
  endif

endif


; Make a set of basis functions that we can fit the data to
ls_matrix=shapelets_interpolate_psf_lsmatrix (degree,                                  $
                                              reform(shapecat.x[*,0]),                 $
                                              reform(shapecat.x[*,1]),                 $
                                              XRANGE=xrange,                           $
                                              YRANGE=yrange,                           $
                                              MAX_P=max_p,                             $
                                              PLOTIT=0B,                               $
                                              N_BASIS=n_basis,                         $
                                              BASIS_FUNCTION_VALUE=basis_function,     $
                                              BASIS_FUNCTION_NAME=basis_function_name, $
                                              BASIS_FUNCTION_ORDER=basis_function_order)

; Fit data to basis functions using least-squares inversion, one PC at a time
coeffs_fit=dblarr(shapecat.n_coeffs,shapecat.n)
basis_function_coeffs=fltarr(shapecat.n_coeffs,n_basis)
for i=0,n_significant-1 do begin

  print,i,n_significant
  
;  covariance=fltarr(shapecat.n,shapecat.n)
;  for j=0L,shapecat.n-1 do covariance[j,j]=1;new_coeffs_error[i,j]^(2)
;  covariance=invert(covariance)

;  bf=basis_function
;  ;for j=0L,shapecat.n_coeffs-1 do basis_function[*,j]=basis_function[*,j]/(new_coeffs_error[j,*]^2) ; / since V^-1 / since INVERSE variance 
;  ls_matrix_temp = invert(bf # covariance # transpose(bf)) # bf # covariance 
;
;  window,0
;  plt_image,ls_matrix_temp,/fr,/col
;  window,2
;  ls_matrix_temp = ls_matrix;#covariance
;  plt_image,ls_matrix_temp,/fr,/col




  
;  ls_matrix_temp = ls_matrix#covariance
  ls_matrix_temp = ls_matrix

  basis_function_coeffs[i,*]=float(ls_matrix_temp # reform(new_coeffs[i,*], shapecat.n, 1)) ;Coefficients
  coeffs_fit[i,*]=reform(reform(basis_function_coeffs[i,*],n_basis) # basis_function, 1, shapecat.n)  ;Return the fit
endfor

; Interpolate PSF to the positions of galaxies
if keyword_set(x) and keyword_set(y) then begin

 
  if keyword_set(subset) then begin
    print,"only interpolating to a few positions!"
    good=lindgen(n_galaxies/20)*20 ; Do only one in twenty, for speed
  endif else good=lindgen(n_galaxies)
  x=x[good] & y=y[good]
  n_galaxies=n_elements(good)
  
  ls_matrix=shapelets_interpolate_psf_lsmatrix (degree, x, y, XRANGE=xrange, YRANGE=yrange, MAX_P=max_p, BASIS_FUNCTION_VALUE=basis_function)
  coeffs_int=dblarr(shapecat.n_coeffs,n_galaxies)
  for i=0,n_significant-1 do begin
    coeffs_int[i,*]=reform(reform(basis_function_coeffs[i,*],n_basis) # basis_function, 1, n_galaxies)  ;Return the fit
  endfor
  
  ;
  ;message,'Must now interpolate to the positions of galaxies'

endif


; Revert PCA and recover true shapelet coefficients
if keyword_set(pca) then begin
  array_reconstruct=dblarr(shapecat.n_coeffs,shapecat.n)
  ;n_significant=n_elements(where(eigenvalues gt 0))
  for i=0,shapecat.n_coeffs-1 do begin
    for j=0,shapecat.n-1 do begin
      array_reconstruct[i,j]=total(coefficients[i,0:n_significant-1]*coeffs_fit[0:n_significant-1,j]/eigenvalues[0:n_significant-1])
    endfor
  endfor
  coeffs_fit=array_reconstruct
  
  ; Check on screen
  print,"  average coeffs     true coeffs       after PCA"+$
        "          fitted     true coeffs"
  for i=0,shapecat.n_coeffs-1 do begin
    print,coeffs_average[i],coeffs[i,test_object],coeffs_pca[i,test_object],$
          coeffs_fit[i,test_object],coeffs[i,test_object]
  endfor

  if keyword_set(x) and keyword_set(y) then begin
    help,coeffs_int
    array_reconstruct=dblarr(shapecat.n_coeffs,n_galaxies)
    for i=0,shapecat.n_coeffs-1 do begin
      for j=0,n_galaxies-1 do begin
        array_reconstruct[i,j]=total(coefficients[i,0:n_significant-1]*coeffs_int[0:n_significant-1,j]/eigenvalues[0:n_significant-1])
      endfor
    endfor
    coeffs_int=array_reconstruct
    help,coeffs_int
  endif

endif
print,'Fitting reconstruction error:',total((coeffs_fit-coeffs)^2)




; Add back in average shapelet coefficients
model=shapecat
for i=0,model.n-1 do model.coeffs[i,*]=(coeffs_fit[*,i]+coeffs_average);*model.flux[i]
shapelets_shapecat_moments, model;, n_max=2



if keyword_set(x) and keyword_set(y) then begin
  shapecat_int={name:"Interpolated PSF in "+shapecat.name, $
                type:"shapecat",$
                n:n_galaxies,$
                maxn_max:shapecat.maxn_max,$
                n_coeffs:shapecat.n_coeffs,$
                polar:shapecat.polar,$
                x:[[x],[y]],$
                beta:replicate(beta,n_galaxies),$
                n_max:replicate(n_max,n_galaxies),$
                coeffs:fltarr(n_galaxies,shapecat.n_coeffs),$
                coeffs_error:fltarr(n_galaxies,shapecat.n_coeffs),$
                flag:lonarr(n_galaxies,2),$
                chisq:fltarr(n_galaxies)+1,$
                back_mean_local:fltarr(n_galaxies),$
                back_rms_local:fltarr(n_galaxies),$
                back_mean_ext:fltarr(n_galaxies),$
                back_rms_ext:fltarr(n_galaxies),$
                seeing:shapecat.seeing,$
                sextractor:0B,$
                morphology:0B,$
                shear_estimates:0B,$
                moments:0B,$
                sexid:good}
  shapecat_int.coeffs=transpose(coeffs_int)+(replicate(1,n_galaxies)#coeffs_average)
  for i=0L,shapecat_int.n-1 do shapecat_int.coeffs[i,*]=(coeffs_int[*,i]+coeffs_average);*model.flux[i]
  ; Renormalise each star to unit flux
  shapelets_shapecat_moments, shapecat_int;, n_max=2
  shapecat_int.flux_error=shapecat_int.flux_error/shapecat_int.flux
  shapecat_int.coeffs=shapecat_int.coeffs/(shapecat_int.flux#replicate(1,shapecat_int.n_coeffs))
  shapecat_int.coeffs_error=shapecat_int.coeffs_error/(shapecat_int.flux#replicate(1,shapecat_int.n_coeffs))
  shapecat_int.flux[*]=1.
  help,shapecat_int,/st
endif


; Make some pretty (?) plots
c=0
plot,[0,0],/nodata,xrange=.02*[-1,1]+coeffs_average[c],$
     yrange=.02*[-1,1]+coeffs_average[c],/iso,/xstyle,/ystyle,$
     title="!6Reconstruction after fit to PCs",$
     xtitle="!6True shapelet coefficient",ytitle="!6Fitted shapelet coefficient"
oplot,[coeffs_average[c],coeffs_average[c]],[-100,100],linestyle=2
oplot,[-100,100],[coeffs_average[c],coeffs_average[c]],linestyle=2
oplot,[0,0],[-100,100]
oplot,[-100,100],[0,0]
oplot,[-100,100],[-100,100],/linestyle
oplot,coeffs[c,*]+coeffs_average[c],coeffs_fit[c,*]+coeffs_average[c],psym=3
oplot,coeffs[c,test_object]*[1,1]+coeffs_average[c],coeffs_fit[c,test_object]*[1,1]+coeffs_average[c],psym=1


; Plot ellipticities
if keyword_set(plotit) then begin
  if !d.name eq "X" or !d.name eq "WIN" then window,2,title="Truth"
  plot,[0,0],/nodata,xstyle=2,ystyle=2,xrange=xrange,yrange=yrange,$
       title="!6True stellar ellipticities",$
       xtitle="!6x [pixels]",ytitle="!6y [pixels]",/iso
  for i=0,shapecat.n-1 do plt_evec,shapecat.x[i,0],shapecat.x[i,1],/e1e2,$
     shapecat.ellipticity[i,0],shapecat.ellipticity[i,1],xscale=escale,yscale=escale
  plt_evec,mean(xrange),mean(yrange),ekey/100.,0.0,/e1e2,thick=2,xscale=escale,yscale=escale
  plt_ellipse,mean(xrange),mean(yrange),ekey/100.,ekey/100.,0,thick=2,xscale=escale,yscale=escale
  plothist,shapecat.rsquared,bin=0.05,xran=[2,5],yran=[0,90]

  if !d.name eq "X" or !d.name eq "WIN" then window,0,title="Model"
  plot,[0,0],/nodata,xstyle=2,ystyle=2,xrange=xrange,yrange=yrange,$
       title="!6Model stellar ellipticities",$
       xtitle="!6x [pixels]",ytitle="!6y [pixels]",/iso
  for i=0,model.n-1 do plt_evec,model.x[i,0],model.x[i,1],/e1e2,$
     model.ellipticity[i,0],model.ellipticity[i,1],xscale=escale,yscale=escale
  plt_evec,mean(xrange),mean(yrange),ekey/100.,0.0,/e1e2,thick=2,xscale=escale,yscale=escale
  plt_ellipse,mean(xrange),mean(yrange),ekey/100.,ekey/100.,0,thick=2,xscale=escale,yscale=escale
  plothist,model.rsquared,bin=0.05,xran=[2,5],yran=[0,90]

  if !d.name eq "X" or !d.name eq "WIN" then window,1,title="Residuals"
  plot,[0,0],/nodata,xstyle=2,ystyle=2,xrange=xrange,yrange=yrange,$
       title="!6Residual stellar ellipticities",$
       xtitle="!6x [pixels]",ytitle="!6y [pixels]",/iso
  for i=0,model.n-1 do plt_evec,model.x[i,0],model.x[i,1],/e1e2,$
     shapecat.ellipticity[i,0]-model.ellipticity[i,0],$
     shapecat.ellipticity[i,1]-model.ellipticity[i,1],xscale=escale,yscale=escale
  plt_evec,mean(xrange),mean(yrange),ekey/100.,0.0,/e1e2,thick=2,xscale=escale,yscale=escale
  plt_ellipse,mean(xrange),mean(yrange),ekey/100.,ekey/100.,0,thick=2,xscale=escale,yscale=escale

  if keyword_set(x) and keyword_set(y) then begin
    if !d.name eq "X" or !d.name eq "WIN" then window,0,title="Interpolated"
    plot,[0,0],/nodata,xstyle=2,ystyle=2,xrange=xrange,yrange=yrange,$
  	 title="!6Interpolated stellar ellipticities",$
  	 xtitle="!6x [pixels]",ytitle="!6y [pixels]",/iso
    for i=0L,shapecat_int.n-1 do plt_evec,shapecat_int.x[i,0],shapecat_int.x[i,1],/e1e2,$
       shapecat_int.ellipticity[i,0],shapecat_int.ellipticity[i,1],xscale=escale,yscale=escale
    plt_evec,mean(xrange),mean(yrange),ekey/100.,0.0,/e1e2,thick=2,xscale=escale,yscale=escale
    plt_ellipse,mean(xrange),mean(yrange),ekey/100.,ekey/100.,0,thick=2,xscale=escale,yscale=escale
  endif
endif

return,shapecat_int

end



