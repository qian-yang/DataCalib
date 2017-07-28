Data Calibration
====
Photometric data calibration for the SDSS-RM project
----
Qian Yang <br>
Feb 11, 2017
----
Reference <br>
Bok 90Prime image data reduction by Ian McGreer (https://github.com/imcgreer/idmrm). The core ubercal code are learnt from the  IDL Ubercal code (by Finkbeiner, Padmanabhan & Schlegel) and python uberpy by Ian (https://github.com/imcgreer/uberpy).
----
Current version <br>
Bok data Ubercal calibration <br>
----
TODO <br>
Bok data i band data error <br>
CFHT data <br>
Spectrophotometry <br>
Mayall data <br>
SDSS data <br>
DRW parameters <br>
Structure function parameters <br>
Short time PSD <br>

----
Main files <br>
get_tref.pro <br>
==> sdssrm-bok_tref.fits (use utDate) <br>

get_otherinfo.pro <br>
==>other_info_g_dt.fits <br>
sdssRefStars_DR13_g.fits <br>

get_nightindex_g.pro <br>
==> sdssrm-bok_tref_nightg.fits <br>

get_aper_data.pro <br>
==> data_aper6_175_amp.fits <br>

get_gframes.pro <br>
==> frames_g.fits <br>
Uniq good frame, 1606 in g band <br>

check_solve4_data.pro (Main) <br>

get_par_frame.pro <br>
==> pars_framesout_new_175_amp.fits <br>

get_aper_ratio.pro <br>
==>  frame_aper_ratio.fits <br>

update_mags_aper.pro <br>
==> mag_refmagout_new_175_amp_aper1_apercor.fits <br>

stat_std_real_newmask_mag.pro <br>
==> stat_mag_refmagout_new_175_amp_aper1_apercor_refmag.fits <br>

get_err.pro <br>
==> pars_framesout_new_175_amp_scale_cor.fits <br>

get_qso_curves_cor.pro <br>
==> lightcurves_bokrm_g_amp_mag_cor.fits <br>

-----
RM mask record <br>
[2, 3, 5, 7, 11, 13, 17, 19, 23, 37, 59, 317, 331] <br>
flag = 0 ==> mask = 2, throw =1 <br>
mag>90 ==> mask=3, throw = 1 <br>
err>0.2 ==> mask=5 <br>
objid, ns <=5 ==> mask = 7 <br>
objid, ns <=0 ==>  mask = 11 (no) <br>
djs_iterstat, 3 sigma outlier ==> mask = 13 (no) <br>
mask, sigma>0.1 ==> mask == 17 (no?) <br>
mag_mean > 20 ==> mask = 19, throw = 0 <br>
frameid, ns<10 ==> mask = 23, throw = 1 <br>
ratio (where > 0.1 or <-0.1) > 0.1==> mask = 29 <br>
remain sigma > 0.3 ==> mask = 31 <br>
objid, ns<=5 ==> mask=37 <br>
objid, ns<=0 ==> mask = 41 <br>
43 no <br>
djs_iterstat, 5sigma, outlier ==> bmask = 47/311 <br>
not mask, nd<=1, mask = 53 <br>
abs(mea - sdss_mag)>0.5 or sigma>0.1 ==> mask = 59/313 <br>
sdss_mag >18 and mask=1 ==> mask = 211 <br>
k<0 ==> mask = 223 <br>
clean_frame, offset>0.1 ==> mask = 317, throw = 1 <br>
-----
Contract <br>
Any question please feel free to email Qian (qianyang.astro@gmail.com).
