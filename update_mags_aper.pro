pro update_mags_aper
; file = 'pars_framesout_new_175_amp.fits'
file = 'pars_framesout_new_175_amp.fits'
frame = mrdfits(file, 1)
nframe = n_elements(frame)
file = 'bokrmphot_g.fits'
data = mrdfits(file, 1)
file = 'refmagout_new_175_amp.fits'
objs = mrdfits(file, 1)
struct_add_field, objs, 'mag', 0.0
struct_add_field, objs, 'photometric', 0
struct_add_field, objs, 'instmag1', 0.0
struct_add_field, objs, 'erradu1', 0.0
struct_add_field, objs, 'apermag1', 0.0
struct_add_field, objs, 'apermagerr1', 0.0
struct_add_field, objs, 'flagq1', 1
struct_add_field, objs, 'remain1', 1
count = data.COUNTS
counterr = data.countsErr
aper = 1
count = count[aper, *]
counterr = counterr[aper, *]
APERMAG = data.APERMAG
APERMAGERR = data.APERMAGERR
objs.apermag1 = transpose(APERMAG[aper, *])
objs.apermagerr1 = transpose(APERMAGERR[aper, *])
struct_add_field, objs, 'count', transpose(count)
struct_add_field, objs, 'count_cor', 0.0
struct_add_field, objs, 'counterr', transpose(counterr)
struct_add_field, objs, 'counterr_cor', 0.0

; mask
flags = data.flags
index = where(flags[aper, *] gt 0)
objs[index].flagq1 *= 2
objs[index].remain1 = 0
index = where(count le 0)
objs[index].flagq1 *= 3
objs[index].remain1 = 0
; mask qso
index = where((objs.flagq mod 11) eq 0)
objs[index].flagq1 *= 11
objs[index].remain1 = 0

; index = where(count[aper, *] gt 0)
; result[index].mag_apercor = -2.5 * alog10(transpose(count[aper, index])/dat[index].exptime)
; result[index].erradu = abs(2.5 / transpose(count[aper, index]) / alog(10.d) * counterr[aper, index])
; file = 'uparamout_new_175_amp.fits'
file = 'uparamout_new_175_amp.fits'
uparam = mrdfits(file, 1)
zeropoint = uparam.zeropoint
file = 'frame_aper_ratio.fits'
ratioarr = mrdfits(file, 1)
ratio = ratioarr.ratio_cal
for i = 0, nframe - 1 do begin
  id = frame[i].frameindex
  print, i, id, ratio[i]
  index = where(objs.FRAMEID eq id)
  if (frame[i].photometric eq 1) then begin
    objs[index].photometric = 1
  endif else if (frame[i].selfmask eq 0) then begin
    objs[index].photometric = 2
  endif
  objs[index].count_cor = objs[index].count * ratio[i]
  objs[index].counterr_cor = objs[index].counterr * ratio[i]
  obj = objs[index]
  ind = where(obj.count gt 0)
  objs[index[ind]].instmag1 = -2.5 * alog10(objs[index[ind]].count_cor/objs[index[ind]].exptime)
  objs[index[ind]].erradu1 = abs(2.5 / objs[index[ind]].count_cor) / alog(10.d) * objs[index[ind]].counterr_cor
  if (i eq 0) then begin
    print, objs[i].instmag, objs[i].instmag1
  endif
  offset = [frame[i].zp1, frame[i].zp2, frame[i].zp3, frame[i].zp4]
  for ccd = 0, 3 do begin
    for amp = 0, 3 do begin
      ind = where(obj.ampind eq ccd*4+amp+1, nt)
      if (nt gt 0) then begin
        objs[index[ind]].mag = objs[index[ind]].instmag1 + zeropoint[ccd*4+amp] + offset[ccd]
      endif
    endfor
  endfor
endfor
; filew = 'mag_refmagout_new_175_amp.fits'
filew = 'mag_refmagout_new_175_amp_aper1_apercor.fits'
spawn, 'rm -rf ' + filew
mwrfits, objs, filew
end
