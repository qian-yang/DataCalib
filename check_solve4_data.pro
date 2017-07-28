pro check_solve4_data
; data = get_data(6)
; objs = mask_new2(data)
; index = where((objs.flagq mod 9) eq 0, n)
; print, n
; ; print, objs[0]
; print, objs[1117]
; mwrfits, objs, 'data_aper6_175_amp.fits'
file = 'data_aper6_175_amp.fits'
objs = mrdfits(file, 1)
nall = long(n_elements(objs))
struct_add_field, objs, 'OIND', lindgen(nall)
struct_add_field, objs, 'mcal', fltarr(nall)
struct_add_field, objs, 'a', fltarr(nall)
struct_add_field, objs, 'k', fltarr(nall)
struct_add_field, objs, 'kdot', fltarr(nall)
struct_add_field, objs, 'mresi', fltarr(nall)
struct_add_field, objs, 'aindx', lonarr(nall)
struct_add_field, objs, 'kindx', lonarr(nall)
;
sort_index = sort(objs.frameid)
objs = objs[sort_index]
struct_add_field, objs, 'framerank', lindgen(nall)
sort_index = sort(objs.oind)
objs = objs[sort_index]
;
file = 'frames_g.fits'
framedata = mrdfits(file, 1)
; aperture_correction, objs, framedata, objs_out = objs_out
nf = n_elements(framedata)
struct_add_field, framedata, 'NON_NUMBER', lonarr(nf)
struct_add_field, framedata, 'NON_RATIO', fltarr(nf)
struct_add_field, framedata, 'cal_mask', lonarr(nf)
struct_add_field, framedata, 'selfmask', lonarr(nf)
struct_add_field, framedata, 'wholenight', lonarr(nf)
struct_add_field, framedata, 'photometric', lonarr(nf)+1

maxiter = 20
maskornot = 0
breakornot = 0
masknight = 0
init_uparam_one_first, objs, uparam_out = uparam_out
uparam = uparam_out
for iiter = 0L, maxiter-1 do begin
   get_ak_indx_all_newmask_nightly, objs, framedata, objs_out = objs_out
   objs = objs_out
   init_uparam_one, objs, uparam, uparam_out = uparam_out
   uparam = uparam_out

   pfit_update_akmags_all_fast_simple_cubic_noflat_dimension_real_newmask_one_new, objs, uparam, uparam_out=uparam_out, objs_out = objs_out
   uparam=uparam_out
   objs = objs_out

 if (masknight ge 3) then maskornot = 1
   update_mcal_all_aktest_real_sdss_sigma_newmask_kdot_one_new,uparam,objs,maskornot, objs_out = objs_out, uparam_out = uparam_out
   objs = objs_out
   uparam = uparam_out

   clean_frame_newmask_nightly_new, objs, framedata, framedata_out = framedata_out, objs_out = objs_out, nc = nc, masknight = masknight
   framedata = framedata_out
   objs = objs_out

   breakornot += 1
   if (breakornot eq 100) then break
   if ((nc eq 0) and (iiter gt 1)) then breakornot = 99

endfor
   mwrfits, objs, 'refmagout_new_175_amp.fits'
   mwrfits, framedata, 'framesout_new_175_amp.fits'
   mwrfits, uparam, 'uparamout_new_175_amp.fits'
  ;  reminder, context = 'newmask'
end
