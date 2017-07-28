pro get_otherinfo
file = 'bokrmphot_g.fits'
objs = mrdfits(file, 1)
nall = long(n_elements(objs))
result = replicate({catalogid:0l, exptime:0.0, ampind:0, $
                    sdss_mag:0.0, sdss_magerr:0.0, dt:0.0}, nall)
result.catalogid = objs.catalogid

get_exptime, objs, exptime_out = exptime_out
result.exptime = exptime_out

get_dt, objs, dt_out = dt_out
result.dt = dt_out

get_amp, objs, amp_out = amp_out
result.ampind = amp_out

objid = objs.objid
objs = 0
file = 'sdssRefStars_DR13.fits'
data = mrdfits(file, 1)
struct_add_field, data, 'i1', 0l
struct_add_field, data, 'i2', 0l
struct_add_field, data, 'number', 0l
nstar = long(n_elements(data))
for i = 0l, nstar - 1l do begin
  index = where(objid eq i, n)
  if (n gt 0) then begin
    print, i
    data[i].number = n
    result[index].sdss_mag = data[i].g
    result[index].sdss_magerr = data[i].err_g
    if ((max(index) - min(index) + 1) eq n) then begin
      data[i].i1 = min(index)
      data[i].i2 = max(index)
    endif else begin
      print, 'Warning: number! ', i
    endelse
  endif
endfor
ind = where(data.number gt 0)
data = data[ind]
mwrfits, data, 'sdssRefStars_DR13_g.fits'
mwrfits, result, 'other_info_g_dt.fits'
end
