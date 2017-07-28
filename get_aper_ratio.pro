pro get_aper_ratio
file = 'bokrmphot_g.fits'
data = mrdfits(file, 1)
nall = long(n_elements(data))
file = 'other_info_g_dt.fits'
dat = mrdfits(file, 1)
sdss_mag = dat.SDSS_MAG
flagq = lonarr(nall)+1
; mask
flags = data.flags
index = where((flags[1, *] gt 0) or (flags[6, *] gt 0))
flagq[index] *= 2
; mask
; file = 'qso_test_cross.fits'
; dd = mrdfits(file, 1)
; nq = n_elements(dd)
; for i = 0, nq - 1 do begin
;     print, i
;     index = where(data.objid eq dd[i].star_id, nid)
;     if (nid gt 0) then begin
;         flagq[index] *=11
;     endif
; endfor
;
frameidarr = data.frameindex
count = data.COUNTS
counterr = data.countsErr
count1 = count[1, *]
count6 = count[6, *]
count1e = counterr[1, *]
count6e = counterr[6, *]
sn1 = count1/count1e
sn6 = count6/count6e
;
file = 'frames_g.fits'
framedata = mrdfits(file, 1)
; aperture_correction, objs, framedata, objs_out = objs_out
nf = n_elements(framedata)
result = replicate({frameid:0l, ratio:0.0, ratio_highsn:0.0, ratio_cal:0.0}, nf)
result.frameid = framedata.frameindex
; struct_add_field, framedata, 'ratio', fltarr(nf)
; struct_add_field, framedata, 'ratio17', fltarr(nf)
; struct_add_field, framedata, 'ratio18', fltarr(nf)
; struct_add_field, framedata, 'ratio19', fltarr(nf)
; struct_add_field, framedata, 'ratio20', fltarr(nf)
; struct_add_field, framedata, 'ratio21', fltarr(nf)
; struct_add_field, framedata, 'ratio22', fltarr(nf)
magarr = findgen(6) + 17
for i = 0, nf - 1 do begin
  frameid = framedata[i].frameindex
  index = where((frameidarr eq frameid) and (flagq eq 1))
  ; sdss = sdss_mag[index]
  ratio = count6[index]/count1[index]
  djs_iterstat, ratio, mean=mean, sigrej=5.0, maxiter = 2, median = median
  result[i].ratio = median
  index = where((frameidarr eq frameid) and (flagq eq 1) and (sn1 gt 10) and (sn6 gt 10))
  ratio = count6[index]/count1[index]
  djs_iterstat, ratio, mean=mean, sigrej=5.0, maxiter = 2, median = median
  result[i].ratio_highsn = median
  index = where((frameidarr eq frameid) and (flagq eq 1) and (sn1 gt 10) and (sn6 gt 10) and (sdss_mag gt 18) and (sdss_mag lt 20))
  ratio = count6[index]/count1[index]
  djs_iterstat, ratio, mean=mean, sigrej=5.0, maxiter = 2, median = median
  result[i].ratio_cal = median
  ; ratioarr = fltarr(6)
  ; for j = 0, 5 do begin
  ;   ind = where((sdss gt magarr[j]) and (sdss le magarr[j]), nd)
  ;   if (nd gt 0) then begin
  ;     rat = ratio[ind]
  ;     djs_iterstat, ratio, mean=mean, sigrej=5.0, maxiter = 2, median = median
  ;     ratioarr[j] = median
  ;   endif
  ; endfor
  ; result[i].ratioarr = ratioarr
  print, result[i]
endfor
filew = 'frame_aper_ratio.fits'
spawn, 'rm -rf ' + filew
mwrfits, result, filew
end
