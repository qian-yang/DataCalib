pro update_mcal_all_aktest_real_sdss_sigma_newmask_kdot_one_new,uparam,objs, maskornot, offset = offset, objs_out = objs_out, uparam_out = uparam_out
aterm = uparam.aterm
kterm = uparam.kterm
kdot = uparam.kdoterm
aindx = objs.aindx
kindx = objs.kindx
airmass = objs.airmass
a = aterm[aindx]
k = kterm[kindx]
objs.a = a
objs.k = k
dt = objs.dt;objs.mjd - floor(objs.mjd)
objs.mresi = objs.sdss_mag - (objs.instmag + a - (k+kdot*dt) *airmass);objs.mcal
zeropoint = fltarr(16)
for i = 0, 15 do begin
  ind = where((objs.cal eq 1) and (objs.ampind eq i+1))
  djs_iterstat, objs[ind].mresi, mean=mean, median = median, sigrej=5.0, mask = mask, sigma = sigma
  print, "!!!! Compare mean, median", mean, median, sigma
  zeropoint[i] = mean
  index = where(objs.ampind eq i+1)
  objs[index].zeropoint = mean
endfor
  uparam.zeropoint = zeropoint
offset = mean;median(objs[ind].mresi_forflat)
print, "!!! get offset value = ", offset
objs.instmag = objs.instmag; + offset
objs.mcal = objs.zeropoint + objs.instmag + a - (k+kdot*dt) *airmass
objs.mresi = objs.sdss_mag - objs.mcal

file = 'sdssRefStars_DR13_g.fits'
dd = mrdfits(file, 1)
objidarr = dd.objid
nstar = n_elements(dd)
print, 'Updating refmag...', nstar
if maskornot then begin
ns = 0
for i = 0, nstar - 1 do begin
    i1 = dd[i].i1
    i2 = dd[i].i2
    n = dd[i].number
    id = lindgen(n) + i1
    objj = objs[id]
    idd = where(objj.cal eq 1, nthis)
    ; if (nthis gt 0) then begin
    if (nthis gt 5) then begin
        mags = objs[id[idd]].instmag ;magadu
        mcals = objs[id[idd]].mcal
        wei = objs[id[idd]].instmagivar ;1.0/(objs[i1:i2].erradu)^2
        djs_iterstat, mcals, mean=mean, sigrej=5.0, mask = mask, sigma = sigma, maxiter = 2
        ind = where(mask eq 0, nd)
        if ((sigma gt 0.1) or (nd gt 5))  then begin
          tparr = objs[id].flagq mod 17
          tpind = where(tparr gt 0, tpn)
          if (tpn gt 0) then begin
            print, '!!!!Warning: mask whole star (variable)', objidarr[i], mean, sigma, nd ;, mags[ind]
            objs[id[tpind]].flagq *= 17 ;47 ;whole star
            objs[id[tpind]].cal = 0
          endif
        endif else begin
        ; don't mask whole object
            if ((nd le 5) and (nd gt 0)) then begin
              tparr = objs[id[idd[ind]]].flagq mod 13
              tpind = where(tparr gt 0, tpn)
              if (tpn gt 0) then begin
                print, '!Warning: mask outlier', objidarr[i], mean, sigma, nd, mags[ind]
                objs[id[idd[ind[tpind]]]].flagq *= 13 ;mask outlier
                objs[id[idd[ind[tpind]]]].cal = 0
                ; objs[id[idd[ind[tpind]]]].stat = 0 ;?
                ; objs[id[idd[ind]]].remain = 2 ; not decide
              endif
            endif
            ; endif else begin
                ind2 = where(mask eq 1, nd)
                if (nd gt 5) then begin
                    mags = mcals[ind2]
                    wei = wei[ind2]
                    weiall = total(wei)
                    mea = total(mcals*wei/weiall)
                    ; objs[id].refmag = mea
                    if ((abs(mea - objs[id[0]].sdss_mag) gt 0.5)) then begin
                      tparr = objs[id].flagq mod 59
                      tpind = where(tparr gt 0, tpn)
                      if (tpn gt 0) then begin
                      print, "!!!Mask whole star object (sdss offset)", objidarr[i], mean, sigma, 'offset', abs(mea - objs[id[0]].sdss_mag)
                      objs[id[tpind]].flagq *= 59 ; 313
                      objs[id[tpind]].cal = 0
                    endif
                    endif else begin
                        print, "object for calibration!", objidarr[i], mea, objs[id[0]].sdss_mag, sigma
                        ns += 1
                    endelse
                endif else begin
                    if (nd gt 0) then begin
                      tparr = objs[id].flagq mod 7
                      tpind = where(tparr gt 0, tpn)
                      if (tpind gt 0) then begin
                        print, '! Warning: do not bother with object <= 5 times unmasked dections'
                        objs[id[tpind]].flagq *= 7
                        objs[id[tpind]].cal = 0
                        objs[id[tpind]].stat = 0
                      endif
                    endif
                endelse
            ; endelse
        endelse
    endif else begin
        if (nthis gt 0) then begin
          tparr = objs[id].flagq mod 7
          tpind = where(tparr gt 0, tpn)
          if (tpn gt 0) then begin
            objs[id[tpind]].flagq *= 7
            objs[id[tpind]].cal = 0
            objs[id[tpind]].stat = 0
            ; objs[id].remain = 1
          endif
        endif
    endelse

endfor
print, "Number of objects remain for calibration: ", ns
endif
objs_out = objs
uparam_out = uparam
end
