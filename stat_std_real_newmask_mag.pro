pro stat_std_real_newmask_mag
  ; file = 'mag_refmagout_new_175_amp.fits'
file = 'mag_refmagout_new_175_amp_aper1_apercor.fits'
objs = mrdfits(file, 1)
nall = long(n_elements(objs))
file = 'sdssRefStars_DR13_g.fits'
data = mrdfits(file, 1)
nstar = n_elements(data)
objidarr = data.objid
struct_add_field, data, 'std_instmag', fltarr(nstar)
struct_add_field, data, 'std_magadu', fltarr(nstar)
; struct_add_field, data, 'std_mcal', fltarr(nstar)
struct_add_field, data, 'std_mag', fltarr(nstar)
struct_add_field, data, 'sdssmag', fltarr(nstar)
struct_add_field, data, 'refmag', fltarr(nstar)
struct_add_field, data, 'stat', bytarr(nstar)
for i = 0, nstar - 1 do begin
        i1 = data[i].i1
        i2 = data[i].i2
        obj = objs[i1:i2]
        print, i
        index = where((obj.photometric eq 1) and (obj.remain1 eq 1) and ((obj.flagq mod 11) gt 0), nthis)
        help, index
        ; index = where((obj.stat eq 1) and (obj.APERMAG lt 90.0), nthis)
        ; index = where((obj.remain eq 1) and (obj.APERMAG lt 90.0), nthis)
        if (nthis eq 0) then begin
            data[i].stat = 0
        endif else begin
            data[i].stat = 1
            objj = obj[index]
            data[i].sdssmag = objj[0].sdss_mag
            djs_iterstat, objj.mag, mean=mean, median = median, sigrej=5.0, mask = mask, sigma = sigma
            data[i].refmag = mean
            data[i].std_instmag = stddev(objj.instmag1)
            data[i].std_magadu = stddev(objj.APERMAG1)
            ; data[i].std_mcal = stddev(objj.mcal)
            data[i].std_mag = stddev(objj.mag)
        endelse
        index = where((obj.photometric eq 1) and (obj.remain1 eq 1) and ((obj.flagq mod 11) eq 0), nthis)
        help, index
        if (nthis gt 0) then begin
            data[i].stat = 2
            objj = obj[index]
            djs_iterstat, objj.mag, mean=mean, median = median, sigrej=5.0, mask = mask, sigma = sigma
            data[i].refmag = mean
            data[i].sdssmag = objj[0].sdss_mag
            data[i].std_instmag = stddev(objj.instmag1)
            data[i].std_magadu = stddev(objj.APERMAG1)
            ; data[i].std_mcal = stddev(objj.mcal)
            data[i].std_mag = stddev(objj.mag)
        endif
endfor
; file = 'stat_mag_refmagout_new_175_amp.fits'
file = 'stat_mag_refmagout_new_175_amp_aper1_apercor_refmag.fits'
spawn, 'rm -rf ' + file
mwrfits, data, file
end
