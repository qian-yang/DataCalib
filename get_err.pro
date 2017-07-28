pro get_err
file = 'pars_framesout_new_175_amp.fits'
frame = mrdfits(file, 1)
file = 'mag_refmagout_new_175_amp_aper1_apercor_refmag.fits'
objs = mrdfits(file, 1)
ind = where(objs.remain eq 1)
objs = objs[ind]
nframe = n_elements(frame)
struct_add_field, frame, 'mean', 0.0
struct_add_field, frame, 'sigma', -1.0
struct_add_field, frame, 'scale', -1.0
for i = 0, nframe - 1 do begin
  frameid = frame[i].frameindex
  index = where((objs.frameid eq frameid) and (objs.remain1 eq 1) and (objs.STATSTAR eq 1))
  ratio = (objs[index].mag - objs[index].refmag)/objs[index].erradu1
  frame[i].mean = mean(ratio)
  frame[i].sigma = stddev(ratio)
  print, i, frameid, frame[i].mean, frame[i].sigma
  set_plot, 'ps'
  !p.multi=[0,1,2,0,2]
  device, filename = './scales/scale_' + STRCOMPRESS(frameid, /REMOVE_ALL) + '.ps'
  plot, histx(ratio/frame[i].sigma, 0.1), histy(ratio/frame[i].sigma, 0.1), xrange = [-5, 5], xstyle = 1, xtitle = '(mcal - refmag)/erradu * scale_factor'
  plot, objs[index].refmag, abs(objs[index].mag - objs[index].refmag), xrange = [16, 24], yrange = [0, 1.0], xstyle = 1, ystyle = 1, xtitle = textoidl('refmag'), ytitle = textoidl('\sigma(mag)'), psym = 4, symsize = 0.2
  oplot, objs[index].refmag, objs[index].erradu1, psym = 4, color = fsc_color('red'), symsize = 0.2
  device, /close
endfor
filew = 'pars_framesout_new_175_amp_scale_cor.fits'
spawn, 'rm -rf ' + filew
mwrfits, frame, filew
end
