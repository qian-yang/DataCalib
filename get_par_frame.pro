pro get_par_frame
file = 'framesout_new_175_amp.fits'
; file = 'framesout_new_175_amp_aper1.fits'
frame = mrdfits(file, 1)
n_night = max(frame.nightindex) + 1
ki = 0
file = 'uparamout_new_175_amp.fits'
; file = 'uparamout_new_175_amp_aper1.fits'
uparam = mrdfits(file, 1)
aterm = uparam.aterm
kterm = uparam.kterm
kdot = uparam.kdoterm
zeropoint = uparam.zeropoint
struct_add_field, frame, 'kindx', -1
struct_add_field, frame, 'airmass', 0.0
struct_add_field, frame, 'mjd', 0.0
struct_add_field, frame, 'exptime', 0.0
struct_add_field, frame, 'fileindex', 0
struct_add_field, frame, 'zp1', 0.0
struct_add_field, frame, 'zp2', 0.0
struct_add_field, frame, 'zp3', 0.0
struct_add_field, frame, 'zp4', 0.0
struct_add_field, frame, 'good', 1
file = 'frames_g.fits'
info = mrdfits(file, 1)
nframe = n_elements(frame)
file = 'refmagout_new_175_amp.fits'
; file = 'refmagout_new_175_amp_aper1.fits'
objs = mrdfits(file, 1)
for i = 0, n_night - 1 do begin
  index = where((frame.nightindex eq i) and (frame.photometric eq 1), nf)
  if (nf gt 0) then begin
    frame[index].kindx = ki
    ki += 1
  endif
endfor

for i = 0, nframe - 1 do begin
  index = where(info.frameIndex eq frame[i].frameindex)
  frame[i].fileindex = index
  print, index, frame[i].frameindex
  frame[i].airmass = info[index].airmass
  frame[i].mjd = info[index].mjd
  dt = frame[i].dt;frame[i].mjd - floor(frame[i].mjd)
  frame[i].exptime = info[index].exptime
  if (info[index].good eq 'T') then begin
    if (frame[i].photometric eq 1) then begin
      ; help, info[index].good
      ki = frame[i].kindx
      ai = ki*4
      ; frame[i].zeropoint += kterm[ki] * airmass + kdot * dt + [fltarr(4) + aterm[ai], fltarr(4) + aterm[ai+1], fltarr(4) + aterm[ai+2], fltarr(4) + aterm[ai+3]]
      frame[i].zp1 = - (kterm[ki] + kdot * dt) * frame[i].airmass + aterm[ai]
      frame[i].zp2 = - (kterm[ki] + kdot * dt) * frame[i].airmass + aterm[ai+1]
      frame[i].zp3 = - (kterm[ki] + kdot * dt) * frame[i].airmass + aterm[ai+2]
      frame[i].zp4 = - (kterm[ki] + kdot * dt) * frame[i].airmass + aterm[ai+3]
    endif else begin
      print, frame[i].frameindex
      ; print, info[index].good
      help, info[index].good
      ind = where((objs.remain eq 1) and (objs.FRAMEID eq frame[i].frameindex) and (objs.sdss_mag gt 17.5) and (objs.sdss_mag lt 20) and (objs.erradu lt 0.1) and ((objs.flagq mod 11) gt 0))
      ; print, frame[i].FRAMEID
      ; help, ind
      obj = objs[ind]
      offset = fltarr(4)
      for ccd = 0, 3 do begin
        res = [0.0]
        for amp = 0, 3 do begin
          id = where(obj.ampind eq ccd*4+amp+1, nt)
          if (nt gt 0) then begin
            res = [res, obj[id].sdss_mag - obj[id].instmag - zeropoint[ccd*4+amp]]
          endif
        endfor
        djs_iterstat, res, mean=mean, median = median, sigrej=5.0, mask = mask, sigma = sigma
        offset[ccd] = mean
      endfor
      ; frame[i].zeropoint += [fltarr(4)+ offset[0], fltarr(4) + offset[1],  fltarr(4) + offset[2], fltarr(4) + offset[3]]
      frame[i].zp1 = offset[0]
      frame[i].zp2 = offset[1]
      frame[i].zp3 = offset[2]
      frame[i].zp4 = offset[3]
    endelse
  endif else begin
    frame[i].good = 0
  endelse
endfor
filew = 'pars_framesout_new_175_amp.fits'
; filew = 'pars_framesout_new_175_amp_aper1.fits'
spawn, 'rm -rf ' + filew
mwrfits, frame, filew
end
