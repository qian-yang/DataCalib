pro get_qso_curves_cor
; pars
file = 'pars_framesout_new_175_amp_scale_cor.fits'
frame = mrdfits(file, 1)
nframe = n_elements(frame)
; zeropoint
file = 'uparamout_new_175_amp.fits'
uparam = mrdfits(file, 1)
zeropoint = uparam.zeropoint

file = 'frame_aper_ratio.fits'
ratioarr = mrdfits(file, 1)
ratio = ratioarr.ratio_cal

file = 'bokrmphot_rmqso.fits'
objs = mrdfits(file, 1)

; objs = get_data_aper(objs_out, 1)

objs = mrdfits(file, 1)
get_exptime, objs, exptime_out = exptime_out
struct_add_field, objs, 'exptime', exptime_out
get_mjd, objs, mjd_out = mjd_out
struct_add_field, objs, 'mjd', mjd_out
get_airmass, objs, airmass_out = airmass_out
struct_add_field, objs, 'airmass', airmass_out
get_dt, objs, dt_out = dt_out
struct_add_field, objs, 'dt', dt_out
get_amp, objs, amp_out = amp_out
struct_add_field, objs, 'ampind', amp_out
objs_out = objs
objs = get_data_aper(objs_out, 1)
struct_add_field, objs, 'count_cor', 0.0
struct_add_field, objs, 'instmag_cor', 0.0
struct_add_field, objs, 'err_cor', 0.d
struct_add_field, objs, 'mag', 0.0
struct_add_field, objs, 'photometric', 0

for i = 0, nframe - 1 do begin
  id = frame[i].frameindex
  index = where(objs.FRAMEID eq id)
  if (frame[i].photometric eq 1) then begin
    objs[index].photometric = 1
  endif else if (frame[i].selfmask eq 0) then begin
    objs[index].photometric = 2
  endif
  objs[index].count_cor = objs[index].count1 * ratio[i]
  obj = objs[index]
  ind = where(obj.count1 gt 0)
  objs[index[ind]].instmag_cor = -2.5 * alog10(objs[index[ind]].count_cor/objs[index[ind]].exptime)
  if (i eq 0) then begin
    print, objs[i].instmag, objs[i].instmag_cor
  endif
  offset = [frame[i].zp1, frame[i].zp2, frame[i].zp3, frame[i].zp4]
  for ccd = 0, 3 do begin
    for amp = 0, 3 do begin
      ind = where(obj.ampind eq ccd*4+amp+1, nt)
      if (nt gt 0) then begin
        objs[index[ind]].mag = objs[index[ind]].instmag_cor + zeropoint[ccd*4+amp] + offset[ccd]
      endif
    endfor
  endfor
  objs[index].err_cor = frame[i].sigma * objs[index].erradu
endfor
filew = 'lightcurves_bokrm_gi.fits'
spawn, 'rm -rf ' + filew
mwrfits, objs, filew
end
