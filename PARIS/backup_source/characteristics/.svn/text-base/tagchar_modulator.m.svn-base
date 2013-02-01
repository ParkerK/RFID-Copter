% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: tag - modulation (reflection coefficient)
%    .) full simulator characteristic (SETTINGS.SIMCHAR = true)
%    .) light characteristic not suitable for simulator (SETTINGS.SIMCHAR = false)
%       ... does not include an extrapolated reflection coefficient, impedances, and vectors/matrices
%           concerning the input power Pin (only Pic, Pav)
%       ... takes less time to calculate (by far)
%
% WARNING: Settings are not checked for sanity!
%
%
%
% ***** Copyright / License / Authors *****
% Copyright 2007, 2008, 2009, 2010, 2011 Daniel Arnitz
%   Signal Processing and Speech Communication Laboratory, Graz University of Technology, Austria
%   NXP Semiconductors Austria GmbH Styria, Gratkorn, Austria
% Copyright 2012 Daniel Arnitz
%   Reynolds Lab, Department of Electrical and Computer Engineering, Duke University, USA
%
% This file is part of the PARIS Simulation Framework.
%
% The PARIS Simulation Framework is free software: you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
%
% The PARIS Simulation Framework is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with the PARIS Simulation 
% Framework. If not, see <http://www.gnu.org/licenses/>.
%
%
% ***** Behavior *****
% tagchar_modulator() 
%    Generates a modulator characteristics using the settings defined within this function.
% SETTINGS.SIMCHAR = true:
%    tagchar_modulator(settings) 
%       Generates a modulator characteristics suitable for simulator usage 
%       according to SETTINGS (results are written to file).
% % SETTINGS.SIMCHAR = false:
%    lightchar = tagchar_modulator(settings)
%       Generates a light modulator characteristics according to SETTINGS, results are returned in 
%       LIGHTCHAR. Note that the results are not suitable for simulator usage.
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
%
%
% ***** Todo *****
% - retry to call lsqnonlin in case of a license error
% ? oversampling for calculation / downsampling before saving the results
% ? prevent large steps for pav (better for plots)
% ? switch to complete IIR model (characteristic = coefficients)
%
% *******************************************************************************************************

function lightchar = tagchar_modulator(settings)
version = 'beta 3.0'; % ALSO SET VERSION BELOW !


% *******************************************************************************************************
% version system (uncomment if this is to be used as script)

% just return version number
if nargin == 0
   lightchar = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% initialization when called from command line (i.e., not in batch mode)
if nargin == 0
   clear; close all; clc; pause(0.1)
   version = 'beta 3.0'; % ALSO SET VERSION BELOW !

   % paths
   path = fullfile(fileparts(mfilename('fullpath')), '..');
   globalinit('silent');
   cd(fileparts(mfilename('fullpath')));

   % general settings
   %     basics
   settings.simchar  =  false; % suitability for simulator (true) or for theory only (false)
   settings.writeres =   true; % also write result (characteristic) file or just workspace?
   settings.suffix   = '_test'; % filename suffix (can be empty), start with '_' if not empty
   settings.fs       =    6e9; % Hz sampling frequency
   %     system
   settings.renew_kticket = false; % renew Kerberos ticket?
   %     folders (subfolders)
   settings.folders.log =       'logs'; % log (diary)
   settings.folders.wsp = 'workspaces'; % workspaceg (for plots, etc)
   settings.folders.res =     'tagmod'; % results .mat-file ("tag modulator characteristic")
   
   % data points per dimension
   %     base vectors (esp. warping errors will increase dramatically for small values)
   settings.nf  = 192; % (old: 192, 384)
   settings.nm  =   2; % (old:  10,  12)
   settings.np  = 256; % (old: 128, 192)
   %     for % "assembly/application tolerance" influence plot (settings.nat x settings.nat)
   settings.nat = 250;
      
   % tag characteristics
   %     filenames
   settings.charfiles.meas_tag_zic    = 'meas_tag_zic';      % chip input impedance
   settings.charfiles.meas_tag_zmod   = 'meas_tag_zmod-act'; % modulation impedance (modulated state)
   settings.charfiles.meas_tag_za     = 'meas_tag_za-v6';    % antenna impedance
   settings.charfiles.meas_tag_picmin = 'meas_tag_picmin';   % minimum operational chip input power
   %     modulation impedance
   settings.zmod_max = complex(4e3, -15e3); % "unmodulated" state; modulated is taken from measurements
   settings.rmod_min = 2;% Ohm minimum real part of modulation impedance ("prevent active Zmod")
   %     special treatment of "passive" zmod characteristics (modulation not possible below min Pic)
   %     ... Zic,mod approx Zic,unmod in this case, thus Zmod measurement close to and below min Pic subject
   %         to heavy distortions => modulation impedance below picopmin replaced by zmod_max
   %     ... threshold='fcn-of-f' will introduce distortions (non-differentiable reflection coefficient)
   %         due to limited resolution in combination with interpolation
   %     ... Use picop_tol with care; the char. has the highest resolution close to picmin, thus a high
   %         picop_tol shifts the "step" of zmod out of the highres area
   settings.psv_zmod.mode      = 'act'; % {'psv': keep passive, 'act': make active}
   settings.psv_zmod.threshold = 'max'; % {'max': use max(Picmin), 'fcn-of-f': use Picmin(f)}
   settings.psv_zmod.picop_tol =     0; % dB tolerance area above picopmin (set to -Inf not to replace anything)
   settings.psv_zmod.flen      =     1; % length of smoothing filter to avoid discontinuities
   
   % power characteristics
   %     limits for chip input power (extrapolation if necessary)
   settings.picmax   = 1e-1; % W
   settings.picmin   = 2e-8; % W
   %     limits for available power
   settings.pavmax   =   10; % W (set high enough to get full char. even for mod)
   settings.pavmin   =    0; % W (zero for automatic)
   %     assumed "average" path-loss factor (linear indexing for log-dist model with this PLF)
   settings.plf =  6;
   %     settings for iterations to determine optimal pav bounds
   settings.npr_inop    = 0.95; % ratio of settings.np withing max(pav) through Pmin
   settings.pavit_maxit = 1000; % maximum number of iterations
   settings.pavit_bias  = 1e-1; % bias for next pav bounds 0<bias<=1 (smaller bias: more exact but slower)
   settings.npma        = round(settings.np/4); % length of MA filter (smoothing of power vector resolution)

   % frequency characteristics
   %     maximum/minimum frequency
   settings.fmin = 1e3; % Hz
   settings.fmax = settings.fs/2; % Hz
   %     maximum/minimum operation frequency (regulations) => better resolution here
   %     (will be adapted to match available frequency grid later)
   settings.foptol = 20e6; % better resolution for operation area plus/minus settings.foptol Hz
   settings.fopmax = 960e6 + settings.foptol; % Hz
   settings.fopmin = 860e6 - settings.foptol; % Hz
   %     settings.fopfactor times more data points within area of operation
   settings.fopfactor = 2.5; % ~2-3 is ok
   %     length of MA filter (smoothing of frequency vector resolution)
   settings.nfma = round(settings.nf/5);
   %     frequency resolution outside characteristic
   settings.fres1 = 15; % bins/GHz for zone1 (fmin ... fchmin)
   settings.fres4 =  3; % bins/GHz for zone4 (fmax ... fchmax)
   %     polynomial models for inter/extrapolation of frequency axis
   settings.npoly_zic  = 2; % (unmodulated) chip impedance
   settings.npoly_zmod = 4; % modulation impedance
   
   % assembly tolerance (parallel RC)
   settings.cat = 450e-15; % F
   %     Q-matching (Rat-Calculation)
   settings.fatm          = 915e6; % Hz frequency
   settings.rat_x0        =   1e3; % Ohm initial value for numerical search for Q-matched Rat
   settings.rat_min       =   1e0; % Ohm minimum parallel resistance
   settings.rat_max       =   1e6; % Ohm maximum parallel resistance
   settings.maxerr_qmatch =     1; % percent maximum error in matched Q (check optimization)
   settings.rat_shift     =     0; % percent tolerance (shift) for Rat
   %     min/max for assembly influence plot
   settings.at_plot.rat_min =      1e1; % Ohm
   settings.at_plot.rat_max =      1e4; % Ohm
   settings.at_plot.cat_min =    1e-15; % F
   settings.at_plot.cat_max = 1200e-15; % F
   
   % detuning (antenna impedance)
   %     detuning state
   settings.detuning.res_en = 0.5; % factor for resonance boost -1..1 (-1: near water, 0: none, 1: near-metal)
   settings.detuning.fshift = 100e6; % Hz shift in frequency (towards smaller frequencies)
   %     constants for resonance boost (TEST CHANGES CAREFULLY!)
   settings.detuning.weight_lin_y =   10; % nonlinearity: weight of linear part, R/X-warping
   settings.detuning.weight_lin_f = 1.25; % nonlinearity: weight of linear part, f-warping
   settings.detuning.weight_exp_y =    2; % nonlinearity: weight of exponential part, R/X-warping
   settings.detuning.weight_exp_f =    4; % nonlinearity: weight of exponential part, f-warping
   settings.detuning.weight_exp_w =    3; % nonlinearity: weight for combining trend and resonance
   settings.detuning.fshift_attf  = 5e-9; % attenuation for resonance per Hz frequency shift
   settings.detuning.nfilt        =    3; % samples filter length of MA filter (keep small!)

   % maximum reflection coefficient (rho = 1 is not possible)
   settings.rhomax     = 0.999;
   settings.rhomax_tol =     1; % percent tolerance for IIR fitting

   % random inverse checks rho(f, zmod, pav) -> Zic
   settings.checks = 1e4; % number of checks
   settings.maxerr = 0.01; % maximum relative mismatch for real/imaginary part
   settings.maxavg = settings.maxerr/10; % maximum relative mismatch mean value (real/imaginary part)
   %     demilitarized zone around NaN (border region) => no checks here (outliers)
   settings.dmz_f = 6; % fch
   settings.dmz_p = 6; % pav
   settings.dmz_m = 1; % zmod

   % IIR bandstop filter parameters to try
   %     pronounced V-shape valley
   settings.iir{1}.iirorder =      2; % filter order
   settings.iir{1}.bshift   =    0.2; % 0...1/2; shift of IIR borders from char. border to center (0: none, 1/2: complete)
   settings.iir{1}.cshift   =    0.6; % 0...1;   shift of IIR center towards minimum (0: none, 1: complete)
   settings.iir{1}.dfmin    =   50e6; % Hz minimum initial bandwidth of IIR bandstop
   %     V-U
   settings.iir{2}.iirorder =      3;
   settings.iir{2}.bshift   =   0.15;
   settings.iir{2}.cshift   =    0.5;
   settings.iir{2}.dfmin    =  100e6;
   %     U-shape valley
   settings.iir{3}.iirorder =      3;
   settings.iir{3}.bshift   =   0.25;
   settings.iir{3}.cshift   =    0.6;
   settings.iir{3}.dfmin    =  150e6;
   %     pronounced U-shape valley
   settings.iir{4}.iirorder =      3;
   settings.iir{4}.bshift   =    0.4;
   settings.iir{4}.cshift   =    0.7;
   settings.iir{4}.dfmin    =  150e6;

   % IIR bandstop optimization (only magnitude is optimized)
   %     minimum number of characteristic points to fit IIR (per side)
   settings.minpts = 3; % < 2 is not allowed (slope detection)
   %     options for lsqnonlin
   settings.optimset = optimset('Display','off', 'LargeScale','on',...
      'TolX',1e-9, 'TolFun',5e-9, 'MaxIter',1e3, 'MaxFunEval',1e5);
   %     initial point x0 [zeros, poles, gain]
   settings.optimstart  = [0.98; 0.98; settings.rhomax];
   %     bounds for x (upper and lower bound) [zeros, poles, gain]
   settings.optimset.ub = [1-100*eps; 1-100*eps; 100]; % stable and minimum-phase
   settings.optimset.lb = [0; 0; 0];
   %     maximum number of optimizer retries to get |rho|<settings.rhomax
   settings.maxoptretries = 2;
   %     cost function weights (see function optim_costfcn for explanations)
   settings.cf_w = [1, 3, 2, 0.1, 2, 3, 3, 50, 1]';
   %     multiplier for weight alteration (will increase IIR residual error but improve fitting quality)
   %     (e.g. in case of negative slope => more weight on |rho(f->)|=rhomax and less on slope)
   settings.cf_walt = 3.5; % has to be > 1
   %     multiplier for |rho|>rhomax punishment (multiplied by this factor for each retry)
   settings.cf_wmult = 3; % has to be > 1 ! % 10
   %     median filter to create smooth borders after IIR extrapolation
   settings.bordermed_len = 5; % samples

   % IIR bandstop residual error warning if maximum error above
   settings.reserr_warn     = -30; % dB (bad fitting)
   settings.reserr_critwarn = -10; % dB (very bad fitting)

end


% *******************************************************************************************************
% logging, multithreading, etc

% start logging
if settings.simchar
   settings.diaryfilename = fullfile(settings.folders.log, sprintf('%s%s.log', mfilename, settings.suffix));
   system(sprintf('rm -f %s', settings.diaryfilename));
   diary(settings.diaryfilename);
end

% version, etc.
disp('*******************************************************************************************************');
disp(sprintf('        Created by %s.m, version %s, %s', mfilename, version, datestr(now, 'local'))); %#ok<*DSPS>
disp('*******************************************************************************************************');

% place runtime warning (debug, under construction, ...) here
% disp(sprintf('\n\n ************* WARNING: UNDER CONSTRUCTION *************\n\n'));
 
% renew Kerberos ticket (if requested)
if settings.renew_kticket
   x1 = system('kinit -R');
   x2 = system('aklog');
   if x1==0 && x2==0
      disp(sprintf('\nSuccessfully renewed Kerberos ticket and AFS authentication token.'));
   end
end
clear('x1', 'x2')


% *******************************************************************************************************
% checks (incomplete; added for more tricky settings and on demand)

% minimum consecutive points to fit IIR (per side)
if settings.minpts < 2 % needs at least 2 points per border to determine the slope
   error('Setup error: settings.minpts < 2 is not allowed.');
end
% minimum/maximum bias
if settings.pavit_bias <= 0 || settings.pavit_bias >= 1
   error('Setup error: 0 < settings.pavit_bias < 1 not met.');
end


% *******************************************************************************************************
% antenna impedance (= frequency base)
disp(sprintf('\nPreparing characteristic: antenna impedance vs frequency (using as frequency base)...'));

% load measurement data: antenna impedance vs frequency (typical)
%     load
meas_tag_za = loadmat(settings.charfiles.meas_tag_za, '   ...');
%     extract data
f  = meas_tag_za.fa;
ra = meas_tag_za.ra;
xa = meas_tag_za.xa;

% check if f is (strictly) monotonically increasing
if any(diff(f) <= 0)
   error('Provided vector f (antenna impedance vs frequency) is not monotonically increasing!');
end

% setup for new (non-homogenous) frequency vector
%     maximum/minimum frequency of provided characteristic
settings.fchmin = meas_tag_za.fa(1);
settings.fchmax = meas_tag_za.fa(end);
%     find fopmin and fopmax (operating area); 'extrap' to avoid NaN when out of range
ind_fopmin = interp1(meas_tag_za.fa, [1:1:length(meas_tag_za.fa)], settings.fopmin, 'nearest', 'extrap');
ind_fopmax = interp1(meas_tag_za.fa, [1:1:length(meas_tag_za.fa)], settings.fopmax, 'nearest', 'extrap');
%     correct settings.fopmin, settings.fopmax
settings.fopmin = meas_tag_za.fa(ind_fopmin);
settings.fopmax = meas_tag_za.fa(ind_fopmax);
%     calculate number of samples outside/inside area of operation etc
%         nf_out1      nf_out2       nf_in       nf_out3      nf_out4
%     fmin <---> fchmin <---> fopmin <---> fopmax <---> fchmax <---> fmax
settings.nf_out1 = ceil(settings.fres1 * (settings.fchmin - settings.fmin)/1e9);
settings.nf_out4 = floor(settings.fres4 * (settings.fmax - settings.fchmin)/1e9);
settings.nf_in   = ceil((settings.nf - settings.nf_out1 - settings.nf_out4) * (settings.fopmax - settings.fopmin) / (settings.fchmax - settings.fchmin) * settings.fopfactor);
settings.nf_out2 = round( (settings.nf - settings.nf_out1 - settings.nf_out4 - settings.nf_in) * (settings.fopmin - settings.fchmin) / (settings.fchmax - settings.fopmax + settings.fopmin - settings.fchmin) );
settings.nf_out3 = settings.nf - settings.nf_out1 - settings.nf_out4 - settings.nf_out2 - settings.nf_in;
%     samples inside characteristic
settings.nfch = settings.nf_out2 + settings.nf_in + settings.nf_out3; 
%     safety checks
if settings.nf_out1 < 1
   disp(sprintf('   Error: Less than one point between fmin and fchmin (nout1).'));
end
if settings.nf_out2 < 1
   disp(sprintf('   Error: Less than one point between fchmin and fopmin (nout2).'));
end
if settings.nf_in < 1
   disp(sprintf('   Error: Less than one point between fopmin and fopmax (nin).'));
end
if settings.nf_out3 < 1
   disp(sprintf('   Error: Less than one point between fopmax and fchmax (nout3).'));
end
if settings.nf_out4 < 1
   disp(sprintf('   Error: Less than one point between fchmax and fmax (nout4).'));
end
if settings.nf_in < 1 || settings.nf_out1 < 1 || settings.nf_out2 < 1 || settings.nf_out3 < 1 || settings.nf_out4 < 1
   error('Unable to split frequency according to settings. See messages above.');
end

% create frequency vector
%     for characteristic only
fch_a = linspace(meas_tag_za.fa(1), meas_tag_za.fa(ind_fopmin), settings.nf_out2+1);
fch_b = linspace(meas_tag_za.fa(ind_fopmin), meas_tag_za.fa(ind_fopmax), settings.nf_in); % operation area
fch_c = linspace(meas_tag_za.fa(ind_fopmax), meas_tag_za.fa(end), settings.nf_out3+1);
fch   = [fch_a(1:end-1), fch_b, fch_c(2:end)]'; % assemble (no jumps at borders!)    
%     smoothing (only for characteristic to make sure no extrapolation has to be done for ra and xa)
if 3*settings.nfma > length(fch)
   disp(sprintf('   Warning: Saturating fch smoothing filter length settings.nfma (old: %i, new: %i).',...
      settings.nfma, floor(length(fch)/3)));
   settings.nfma = floor(length(fch)/3);
end
fch = filtfilt(ones(settings.nfma, 1)/settings.nfma, 1, fch);
%     full fmin...fmax
f = [linspace(settings.fmin, meas_tag_za.fa(1), settings.nf_out1+1),...                % not in characteristic
     fch(2:end-1)',...                                                      % in characteristic
     logspace(log10(meas_tag_za.fa(end)), log10(settings.fmax), settings.nf_out4+1)]'; % not in characteristic    

% interpolate ('extrap' just prevents NaN is min/max of f is just outside min/max of meas_tag_za.fa)
ra = interp1(meas_tag_za.fa, ra, fch, 'spline', 'extrap');
xa = interp1(meas_tag_za.fa, xa, fch, 'spline', 'extrap');

% cleanup
clear('ind_fopmin', 'ind_fopmax', 'fch_a', 'fch_b', 'fch_c');


% *******************************************************************************************************
% minimum operational power
disp(sprintf('\nPreparing characteristic: minimum operational chip input power Pic,min...'));

% load measurement data: minimum operational power vs. frequency
meas_tag_picmin = loadmat(settings.charfiles.meas_tag_picmin, '   ...');

% extrapolate/interpolate impedance to frequency vector
%    linear extrap; log domain to avoid negative power levels caused by extrapolation
picopmin = 10.^interp1(meas_tag_picmin.f, log10(meas_tag_picmin.p)+2.5/10, fch, 'linear', 'extrap'); % W
%    spline interp within available data ('extrap' only to cover borders)
ind_f1 = interp1(fch, 1:settings.nfch, meas_tag_picmin.f(1), 'nearest', 'extrap');
ind_f2 = interp1(fch, 1:settings.nfch, meas_tag_picmin.f(end), 'nearest', 'extrap');
picopmin(ind_f1:ind_f2) = 10.^interp1(meas_tag_picmin.f, log10(meas_tag_picmin.p), fch(ind_f1:ind_f2), 'spline');

% figure; hold on;
% plot(meas_tag_picmin.f/1e6, 10*log10(meas_tag_picmin.p)+30, 'bo');
% plot(fch/1e6, 10*log10(picopmin)+30, 'b-');
% hold off; grid on; ylim(xyzlimits(10*log10(meas_tag_picmin.p),10*log10(picopmin))+30);
% setlabels('INTERPOLATION OF MINIMUM OPERATIONAL CHIP POWER LEVELS', 'fch [MHz]', 'P_{ic,min} [dBm]');


% *******************************************************************************************************
% chip input impedance vs input power and frequency (= power base)
disp(sprintf('\nPreparing characteristic: chip input impedance Zic vs chip input power (using as power base) and frequency...'));

% load measurement data: chip input impedance vs frequency and power
meas_tag_zic = loadmat(settings.charfiles.meas_tag_zic, '   ...');

% extrapolation warnings
if min(meas_tag_zic.fic) - min(fch) > 0 || max(meas_tag_zic.fic) - max(fch) > 0
   disp(sprintf('   Warning: Extrapolating frequency axis (old: %.1f to %.1f MHz, new: %.1f to %.1f MHz).',...
      min(meas_tag_zic.fic)*1e-6, max(meas_tag_zic.fic)*1e-6, min(fch)*1e-6, max(fch)*1e-6 ));
end
if min(meas_tag_zic.pic) > settings.picmin || max(meas_tag_zic.pic) < settings.picmax
   disp(sprintf('   Warning: Extrapolating power axis (old: %.1f to %.1f dBm, new: %.1f to %.1f dBm).',...
      10*log10(min(meas_tag_zic.pic))+30, 10*log10(max(meas_tag_zic.pic))+30,...
      10*log10(settings.picmin)+30, 10*log10(settings.picmax)+30 ));
end

% prepare data matrices
%     frequency axis extrapolation
ric_tmp = nan(settings.nfch, length(meas_tag_zic.pic)); % [fch, pic]; note that we've transposed the original matrices
xic_tmp = nan(settings.nfch, length(meas_tag_zic.pic));
%     complete: frequency and power axis extrapolation
ric = nan(settings.nfch, settings.np); % [fch, pic]; note that we've transposed the original matrices
xic = nan(settings.nfch, settings.np);

% frequency axis: inter-/extrapolation using a polynomial model (works only for slight extrapolation)
for i = 1 : length(meas_tag_zic.pic)
   %     standard interpolation
   ric_tmp(:,i) = interp1(meas_tag_zic.fic, meas_tag_zic.ric(i,:)', fch, 'spline', NaN);
   xic_tmp(:,i) = interp1(meas_tag_zic.fic, meas_tag_zic.xic(i,:)', fch, 'spline', NaN);
   %     use polynomial model to extrapolate
   [poly_ric, S_ric, mu_ric] = polyfit(meas_tag_zic.fic, meas_tag_zic.ric(i,:)', settings.npoly_zic);
   [poly_xic, S_xic, mu_xic] = polyfit(meas_tag_zic.fic, meas_tag_zic.xic(i,:)', settings.npoly_zic);
   ric_tmp(isnan(ric_tmp(:,i)),i) = polyval(poly_ric, fch(isnan(ric_tmp(:,i))), [], mu_ric);
   xic_tmp(isnan(xic_tmp(:,i)),i) = polyval(poly_xic, fch(isnan(xic_tmp(:,i))), [], mu_xic);
end

% power axis: extra-/interpolate to new power vector (linear in log scale)
pic = logspace( log10(settings.picmin), log10(settings.picmax), settings.np)';
for i = 1 : settings.nfch
   ric(i,:) = interp1(log10(meas_tag_zic.pic), ric_tmp(i,:), log10(pic), 'linear', 'extrap')';
   xic(i,:) = interp1(log10(meas_tag_zic.pic), xic_tmp(i,:), log10(pic), 'linear', 'extrap')';
end

% check if an extrapolation produced an active circuit
if any(any( ric < 0 ))
      disp(sprintf('   WARNING: Extrapolation created %.2f%% R_ic<0. Settings these values to R_ic=0.',...
         100*sum(sum(ric<0)) / numel(ric) ));
      ric(ric<0) = 0;
end

% chip input capacity
cic = -1 ./ (2*pi*repmat(fch, 1, settings.np) .* xic);

% cleanup
clear('ric_tmp','poly_ric','S_ric','mu_ric', 'xic_tmp','poly_xic','S_xic','mu_xic', ...
   'ind_fic1','ind_fic2', 'fch_tmp');

% figure; hold on;
% surface(meas_tag_zic.fic/1e6, meas_tag_zic.pic, meas_tag_zic.ric, 'FaceColor', 'none');
% surface(fch/1e6, pic, ric', 'EdgeColor', 'none');
% hold off; grid on; view(-25, 30); axis tight; set(gca, 'yScale','log');
% setlabels('CHIP IMPEDANCE INTER/EXTRAPOLATION: Ric', 'f [MHz]', 'Pic [W]', 'Ric [\Omega]');
% setlegend({'original', 'interp2'}, 'NorthEast');
% 
% figure; hold on;
% surface(meas_tag_zic.fic/1e6, meas_tag_zic.pic, meas_tag_zic.xic, 'FaceColor', 'none');
% surface(fch/1e6, pic, xic', 'EdgeColor', 'none');
% hold off; grid on; view(-25, 30); axis tight; set(gca, 'yScale','log'); 
% setlabels('CHIP IMPEDANCE INTER/EXTRAPOLATION: Xic', 'f [MHz]', 'Pic [W]', 'Xic [\Omega]');
% setlegend({'original', 'interp2'}, 'NorthEast');
% 
% figure; hold on;
% ind_fic1 =  ceil(interp1(fch, [1:1:settings.nfch], min(meas_tag_zic.fic), 'linear', 'extrap'));
% ind_fic2 = floor(interp1(fch, [1:1:settings.nfch], max(meas_tag_zic.fic), 'linear', 'extrap'));
% surface(meas_tag_zic.fic/1e6, meas_tag_zic.pic, meas_tag_zic.ric, 'FaceColor', 'none');
% surface(fch(ind_fic1:ind_fic2)/1e6, pic, ric(ind_fic1:ind_fic2,:)', 'EdgeColor', 'none');
% hold off; grid on; view(-25, 30); axis tight; set(gca, 'yScale','log');
% setlabels('CHIP IMPEDANCE FREQUENCY INTERPOLATION DETAIL: Ric', 'f [MHz]', 'Pic [W]', 'Ric [\Omega]');
% setlegend({'original', 'interp2'}, 'NorthEast');
% 
% figure; hold on;
% surface(meas_tag_zic.fic/1e6, meas_tag_zic.pic, meas_tag_zic.xic, 'FaceColor', 'none');
% surface(fch(ind_fic1:ind_fic2)/1e6, pic, xic(ind_fic1:ind_fic2,:)', 'EdgeColor', 'none');
% hold off; grid on; view(-25, 30); axis tight; set(gca, 'yScale','log');
% setlabels('CHIP IMPEDANCE FREQUENCY INTERPOLATION DETAIL: Xic', 'f [MHz]', 'Pic [W]', 'Xic [\Omega]');
% setlegend({'original', 'interp2'}, 'NorthEast');
% 
% figure; hold on;
% ind_pic1 = max(1,            ceil(interp1(pic, [1:1:settings.np], min(meas_tag_zic.pic), 'linear', 'extrap')));
% ind_pic2 = min(settings.np, floor(interp1(pic, [1:1:settings.np], max(meas_tag_zic.pic), 'linear', 'extrap')));
% surface(fch/1e6, pic(ind_pic1:ind_pic2), cic(:,ind_pic1:ind_pic2)', 'EdgeColor', 'none');
% hold off; grid on; view(-25, 30); axis tight; set(gca, 'yScale','log'); 
% setlabels('CHIP IMPEDANCE INTER/EXTRAPOLATION: Cic (EXCLUDING POWER EXTRAP)', 'f [MHz]', 'Pic [W]', 'Cic [F]');


% *******************************************************************************************************
% modulation impedance
% ... large parts copied from Zic part above
disp(sprintf('\nPreparing characteristic: modulation impedance (plus switching characteristic)...'));

% load measurement data: modulation impedance vs power
meas_tag_zmod = loadmat(settings.charfiles.meas_tag_zmod, '   ...');

% extrapolation warnings
if min(meas_tag_zmod.pic) > settings.picmin || max(meas_tag_zmod.pic) < settings.picmax
   disp(sprintf('   Warning: Extrapolating power axis (old: %.1f to %.1f dBm, new: %.1f to %.1f dBm).',...
      10*log10(min(meas_tag_zmod.pic))+30, 10*log10(max(meas_tag_zmod.pic))+30,...
      10*log10(settings.picmin)+30, 10*log10(settings.picmax)+30 ));
end
if min(meas_tag_zmod.fic) - min(fch) > 1 || max(meas_tag_zmod.fic) - max(fch) > 1 % extrap more than 1 Hz
   disp(sprintf('   Warning: Extrapolating frequency axis (old: %.1f to %.1f MHz, new: %.1f to %.1f MHz).',...
      min(meas_tag_zmod.fic)*1e-6, max(meas_tag_zmod.fic)*1e-6, min(fch)*1e-6, max(fch)*1e-6 ));
end

% prepare data matrices
%     frequency axis extrapolation
rmod_tmp = nan(settings.nfch, length(meas_tag_zmod.pic)); % [fch, pic]; note that we've transposed the original matrices
xmod_tmp = nan(settings.nfch, length(meas_tag_zmod.pic));
%     complete: frequency and power axis extrapolation
rmod = nan(settings.nfch, settings.np); % [fch, pic]; note that we've transposed the original matrices
xmod = nan(settings.nfch, settings.np);

% frequency axis: inter-/extrapolation using a polynomial model (works only for slight extrapolation)
for i = 1 : length(meas_tag_zmod.pic)
   %     standard interpolation
   rmod_tmp(:,i) = interp1(meas_tag_zmod.fic, meas_tag_zmod.rmod(i,:)', fch, 'spline', NaN);
   xmod_tmp(:,i) = interp1(meas_tag_zmod.fic, meas_tag_zmod.xmod(i,:)', fch, 'spline', NaN);
   %     use polynomial model to extrapolate
   [poly_rmod, S_rmod, mu_rmod] = polyfit(meas_tag_zmod.fic, meas_tag_zmod.rmod(i,:)', min(size(meas_tag_zmod.rmod,2)-1, settings.npoly_zmod));
   [poly_xmod, S_xmod, mu_xmod] = polyfit(meas_tag_zmod.fic, meas_tag_zmod.xmod(i,:)', min(size(meas_tag_zmod.xmod,2)-1, settings.npoly_zmod));
   rmod_tmp(isnan(rmod_tmp(:,i)),i) = polyval(poly_rmod, fch(isnan(rmod_tmp(:,i))), [], mu_rmod);
   xmod_tmp(isnan(xmod_tmp(:,i)),i) = polyval(poly_xmod, fch(isnan(xmod_tmp(:,i))), [], mu_xmod);
end

% power axis: interpolate to new power vector
% ... linear in double log for rmod for safe exptrapolation (=> active circuits)
% ... linear in log for xmod (negative => extrapolation in loglog scale would be more complex to)
for i = 1 : settings.nfch
   rmod(i,:) = interp1(log(meas_tag_zmod.pic), rmod_tmp(i,:), log(pic), 'linear', 'extrap')';
   xmod(i,:) = interp1(log(meas_tag_zmod.pic), xmod_tmp(i,:), log(pic), 'linear', 'extrap')';
end

% check for active circuits
if any(any( rmod < settings.rmod_min ))
      disp(sprintf('   WARNING: %.2f%% R_mod<rmod_min. Setting these values to minimum.',...
         100*sum(sum(rmod<settings.rmod_min)) / numel(rmod) ));
      rmod(rmod<settings.rmod_min) = settings.rmod_min;
end

% add trajectory during modulation => zmod(fch, mod, pic)
%    ... assume that zmod -> settings.zmod_max for unmodulated state
%    ... exponential characteristic for unmod -> mod -> unmod if linearly indexed
zmod_tran = (logspace(0, 6, settings.nm) - 1) / (1e6-1); 
for i = 1 : settings.nm
   zmod(:, i, :) = complex(rmod, xmod) * (1-zmod_tran(i)) + settings.zmod_max * zmod_tran(i);
end

% if this is a "passive" (mod approx unmod below picopmin) modulation char and we want to have special treatment
% (almost identical mod/unmod => measured modulation impedance consists of noise only)
% => replace zmod by zmod_max below picopmin (plus some tolerance)
if ~meas_tag_zmod.act && ~isinf(settings.psv_zmod.picop_tol)
   disp('   This is a "passive" Zmod characteristic (modulation below Picmin not possible).');
   % mode: keep passive or make active?
   switch(lower(settings.psv_zmod.mode))
      case 'psv'
         fprintf('      ... keeping "passive", but settings Zmod=zmod_max below ');
      case 'act'
         fprintf('      ... making "active" (i.e., modulation below Picmin possible) by holding Zmod @ ');
      otherwise
         error('Unsupported settings.psv_zmod.mode=''%s''.', settings.psv_zmod.mode);
   end
   % special treatment below picmin (more realistic) or max(picmin) (less distortions due to limited resolution)
   switch(lower(settings.psv_zmod.threshold))
      case 'max'
         fprintf('max(Picmin).\n');
         picopmin_zmod = ones(size(picopmin)) * max(picopmin) * 10^(settings.psv_zmod.picop_tol/10);
      case 'fcn-of-f'
         fprintf('Picmin(f). Note that this will create discontinuities due to limited resolution.\n');
         picopmin_zmod = picopmin;
      otherwise
         error('Unsupported settings.psv_zmod.threshold=''%s''.', settings.psv_zmod.threshold);
   end
   % shift picopmin_zmod
   picopmin_zmod = picopmin_zmod * 10^(settings.psv_zmod.picop_tol/10);
   % replace parts below picopmin_zmod
   for i = 1 : settings.nfch
      for j = 1 : settings.nm
         switch(lower(settings.psv_zmod.mode))
            case 'psv' % replace by zmod_max => keep passive
               zmod(i, j, pic < picopmin_zmod(i)) = settings.zmod_max;
            case 'act' % replace by last zmod => make active
               ind_last = ceil( interp1(pic, 1:settings.np, picopmin_zmod(i), 'linear') );
               zmod(i, j, pic < picopmin_zmod(i)) = zmod(i, j, ind_last);
         end
         if settings.psv_zmod.flen > 1 % blunt edges
            zmod(i, j, :) = filtfilt(ones(settings.psv_zmod.flen,1)/settings.psv_zmod.flen, 1, zmod(i, j, :));
         end
      end
   end
end

% cleanup
clear('rmod_tmp','poly_rmod','S_rmod','mu_rmod', 'xmod_tmp','poly_xmod','S_xmod','mu_xmod','zmod_tran',...
   'picopmin_zmod', 'ind_last');

% figure; hold on;
% surface(meas_tag_zmod.fic/1e6, 10*log10(meas_tag_zmod.pic)+30, meas_tag_zmod.rmod, 'FaceColor','none','EdgeColor','k');
% surface(fch/1e6, 10*log10(pic)+30, squeeze(real(zmod(:,  1,:)))', 'FaceColor','interp','EdgeColor','none');
% hold off; grid on; set(gca, 'zScale', 'log');
% xlim(xyzlimits(fch)/1e6); ylim(xyzlimits(10*log10(pic)+30)); setlegend({'meas', 'model'}, 'NorthWest');
% setlabels('INTERPOLATION OF MODULATION IMPEDANCE (MOD. STATE): Re', 'f [MHz]', 'Pic [dBm]', 're(Zmod) [\Omega]');
% 
% figure; hold on;
% surface(meas_tag_zmod.fic/1e6, 10*log10(meas_tag_zmod.pic)+30, abs(meas_tag_zmod.xmod), 'FaceColor','none','EdgeColor','k');
% surface(fch/1e6, 10*log10(pic)+30, abs(squeeze(imag(zmod(:,  1,:))))', 'FaceColor','interp','EdgeColor','none');
% hold off; grid on; set(gca, 'zScale', 'log');
% xlim(xyzlimits(fch)/1e6); ylim(xyzlimits(10*log10(pic)+30)); setlegend({'meas', 'model'}, 'NorthWest');
% setlabels('INTERPOLATION OF MODULATION IMPEDANCE (MOD. STATE): Im', 'f [MHz]', 'Pic [dBm]', '|im(Zmod)| [\Omega]');
% 
% figure; hold on;
% surface(meas_tag_zmod.fic/1e6, 10*log10(meas_tag_zmod.pic)+30, sqrt(meas_tag_zmod.rmod.^2+meas_tag_zmod.xmod.^2), 'FaceColor','none','EdgeColor','k');
% surface(fch/1e6, 10*log10(pic)+30, squeeze(abs(zmod(:,  1,:)))', 'FaceColor','interp','EdgeColor','none');
% hold off; grid on; set(gca, 'zScale', 'log');
% xlim(xyzlimits(fch)/1e6); ylim(xyzlimits(10*log10(pic)+30)); setlegend({'meas', 'model'}, 'NorthWest');
% setlabels('INTERPOLATION OF MODULATION IMPEDANCE (MOD. STATE): Abs', 'f [MHz]', 'Pic [dBm]', 'abs(Zmod) [\Omega]');
% 
% figure; hold on;
% surface(meas_tag_zmod.fic/1e6, 10*log10(meas_tag_zmod.pic)+30, 180/pi*atan(meas_tag_zmod.xmod./meas_tag_zmod.rmod), 'FaceColor','none','EdgeColor','k');
% surface(fch/1e6, 10*log10(pic)+30, 180/pi*squeeze(angle(zmod(:,  1,:)))', 'FaceColor','interp','EdgeColor','none');
% hold off; grid on;
% xlim(xyzlimits(fch)/1e6); ylim(xyzlimits(10*log10(pic)+30)); setlegend({'meas', 'model'}, 'NorthWest');
% setlabels('INTERPOLATION OF MODULATION IMPEDANCE (MOD. STATE): Angle', 'f [MHz]', 'Pic [dBm]', 'angle(Zmod) [deg]');


% *******************************************************************************************************
% assembly (tolerances): matching at minimum power at settings.fatm
disp(sprintf('\nAssembly: Calculating Rat for given Cat by matching Q @ %.0f MHz and Pmin(%.0fMHz)...',...
   settings.fatm*1e-6, settings.fatm*1e-6));

% indices in frequency/power vectors
ind_fatm = interp1(fch, [1:1:settings.nfch], settings.fatm, 'nearest', 'extrap');
ind_patm = interp1(pic, [1:1:length(pic)], picopmin(ind_fatm), 'nearest', 'extrap');

% assembly tolerances 
%     check if tolerances for Rat is integer (percents)
if round(settings.rat_shift) ~= settings.rat_shift
   disp(sprintf('   WARNING: settings.rat_shift is not integer. Rounding from %g to %g percent.', settings.rat_shift, round(settings.rat_shift)));
   settings.rat_shift = round(settings.rat_shift);
end
%     numerically solve the Q-matching
temp.za  = complex(ra(ind_fatm), xa(ind_fatm)); % antenna
temp.zic = complex(ric(ind_fatm, ind_patm), -1/(2*pi*settings.fatm*cic(ind_fatm, ind_patm))); % chip
temp.zl = 1 ./ ( 1./temp.zic + complex(0, 2*pi*settings.fatm*settings.cat) ); % load: chip plus Cat
temp.qa = abs(imag(temp.za))./real(temp.za); % quality factor of antenna
temp.qmatch_cost = @(x) ( abs(imag( (temp.zl*x)./(temp.zl+x) )) ./ real( (temp.zl*x)./(temp.zl+x) ) - temp.qa).^2;
settings.rat = fminsearch(temp.qmatch_cost, settings.rat_x0, struct('Display', 'off'));
if  settings.rat < settings.rat_min
   disp(sprintf('   Warning: Optimization produced an unrealistically low resistance. Truncating to Rat = %.1g Ohm.', settings.rat_min));
   settings.rat = settings.rat_min;
elseif settings.rat > settings.rat_max
   disp(sprintf('   Warning: Optimization produced an unrealistically high resistance. Truncating to Rat = %.1g MOhm.', settings.rat_max/1e6));
   settings.rat = settings.rat_max;
end
%     check this result
temp.zl = temp.zl * settings.rat / (temp.zl + settings.rat);
temp.ql = abs(imag(temp.zl))./real(temp.zl); % quality factor of chip and assembly impedance
if abs(temp.ql - temp.qa) ./ temp.qa * 100 > settings.maxerr_qmatch
   disp(sprintf('   WARNING: Q-matching has failed: Qa = %f, Ql = %f', temp.qa, temp.ql));
end
%     shift from optimal value (detuning)
settings.rat = settings.rat * (1 + settings.rat_shift/100);
%     output
if settings.rat_shift ~= 0
   disp(sprintf('   Shifting Rat by %g percent (from optimum).', settings.rat_shift));
end
disp(sprintf('   Final assembly: %.1g Ohm / %.1g fF parallel', settings.rat, settings.cat*1e15));

% cleanup
clear('temp');


% settings.rat = fminsearch(temp.qmatch_cost, 100, settings.detuning.optimset);

% *******************************************************************************************************
% influence of assembly on rho at Pmin and standard operating frequency fpmeas
if ~settings.simchar
   disp(sprintf('\nSkipping calculation of best assembly impedance.'));
else
   disp(sprintf('\nComputing influence of assembly on rho @ %.0f MHz and Pmin(%.0fMHz)...',...
      settings.fatm*1e-6, settings.fatm*1e-6));
   
   % prepare matrices / vectors
   %     assembly impedances
   rat_at = linspace(settings.at_plot.rat_min, settings.at_plot.rat_max, settings.nat)';
   cat_at = linspace(settings.at_plot.cat_min, settings.at_plot.cat_max, settings.nat)';
   %     resulting reflection coefficient
   rho_at = zeros(2, settings.nat, settings.nat); % unmodulated, modulated
   
   % get antenna and chip impedance at Pmin and fpmeas (approximately)
   za_at  = complex(ra(ind_fatm), xa(ind_fatm));
   zic_at = complex(ric(ind_fatm, ind_patm), -1/(2*pi*settings.fatm*cic(ind_fatm, ind_patm)));
   
   % get modulation impedance at this frequency/power
   [x, y] = meshgrid(pic, fch);
   zmod_min = interp2(x,y, squeeze(zmod(:,  1, :)), pic(ind_patm), settings.fatm, 'linear');
   zmod_max = interp2(x,y, squeeze(zmod(:,end, :)), pic(ind_patm), settings.fatm, 'linear');
   
   % calculate reflection coefficient
   for i = 1 : settings.nat % for Rat
      for j = 1 : settings.nat % for Cat
         % assembly impedance
         zat = 1 / (1/rat_at(i) + complex(0, 2*pi*settings.fatm*cat_at(j)));
         % input impedance (unmodulated)
         zin0 = 1 / (1/zic_at + 1/zat + 1/zmod_max);
         zin1 = 1 / (1/zic_at + 1/zat + 1/zmod_min);
         % reflection coefficient
         rho_at(1,i,j) = (za_at - conj(zin0)) / (za_at + zin0);
         rho_at(2,i,j) = (za_at - conj(zin1)) / (za_at + zin1);
      end
   end
   
   % best assembly (unmodulated)
   %   find minimum (Rat, Cat)
   [min_rho_at, ind_bestrat] = min(abs(squeeze(rho_at(1,:,:))), [], 1);
   [min_rho_at, ind_bestcat] = min(min_rho_at);
   ind_bestrat = ind_bestrat(ind_bestcat);
   %	output
   disp(sprintf('   Best Assembly, Unmod. (from plot): %10.1f Ohm / %6.1f fF parallel', rat_at(ind_bestrat), cat_at(ind_bestcat)*1e15));
end

% cleanup
clear('i', 'j', 'ind_fatm', 'ind_patm', 'za_at', 'zic_at', 'x', 'y', 'zat', 'zmod_min', 'zmod_max',...
   'zin0', 'zin1', 'min_rho_at', 'ind_bestrat', 'ind_bestcat', 'zat_best', 'zat_best_par');



% *******************************************************************************************************
% detuning (empirical)
disp(sprintf('\nDetuning: Warping antenna impedance around its resonance...'));

% extract linear trends and differences to this trend
%     R
det.ra   = ra; % for debugging
det.ra_l = interp1(fch([1,end]), ra([1,end]), fch, 'linear', 'extrap'); % linear trend
det.ra_r = ra - det.ra_l; % difference to this trend
det.ra_d = det.ra_r; % detuned impedance (difference)
%     X
det.xa   = xa; % for debugging
det.xa_l = interp1(fch([1,end]), xa([1,end]), fch, 'linear', 'extrap'); % linear trend
det.xa_r = xa - det.xa_l; % difference to this trend
det.xa_d = det.xa_r; % detuned impedance (difference)

% position of antenna resonance [could be made part of the characterisitc...]
[dummy, det.ind_res] = max(ra);

% enhance resonance
if settings.detuning.res_en ~= 0
   % define nonlinearities (centered at resonance)
   det.nl_x = (fch - fch(det.ind_res))/max(fch - fch(det.ind_res));
   det.nl = @(wl, we, d) (wl-d)/(2*wl)*det.nl_x + (wl+d)/(2*wl) * sign(det.nl_x) .* abs(det.nl_x.^(we*abs(d)+1)); % (weight-linear, weight-exp, detuning)
   det.nl_f = det.nl(settings.detuning.weight_lin_f, settings.detuning.weight_exp_f, settings.detuning.res_en); % for f
   det.nl_y = det.nl(settings.detuning.weight_lin_y, settings.detuning.weight_exp_y, settings.detuning.res_en); % for R, X
   % warp frequency vector (concentrate around resonance point)
   if settings.detuning.res_en > 0
      det.nl_fa = interp1(det.nl_f, fch, det.nl_x, 'nearest', 'extrap'); % this requires a high resolution!
   else
      det.nl_fa = interp1(det.nl_x, fch, det.nl_f, 'nearest', 'extrap'); % this requires a high resolution!
   end
   %     apply
   det.ra_d  = interp1(fch, det.ra_r, det.nl_fa, 'linear', 'extrap');
   det.xa_d  = interp1(fch, det.xa_r, det.nl_fa, 'linear', 'extrap');
   %     re-align extrapolated parts
   det.ra_d = det.ra_d - [linspace(det.ra_d(1), 0, det.ind_res), linspace(0, det.ra_d(end), length(fch)-det.ind_res)]';
   det.xa_d = det.xa_d - [linspace(det.xa_d(1), 0, det.ind_res), linspace(0, det.xa_d(end), length(fch)-det.ind_res)]';
   
   % warp impedance
   %     prepare exponent
   if settings.detuning.res_en > 0
      det.exp_rx = (settings.detuning.res_en+1)^settings.detuning.weight_exp_w;
   else
      det.exp_rx = (settings.detuning.res_en+1)^(1/settings.detuning.weight_exp_w);
   end
   det.ra_d = interp1(det.nl_x * max(abs(det.ra_r)), det.nl_y * max(abs(det.ra_r)), det.ra_d, 'linear','extrap') * det.exp_rx;
   det.xa_d = interp1(det.nl_x * max(abs(det.xa_r)), det.nl_y * max(abs(det.xa_r)), det.xa_d, 'linear','extrap') * det.exp_rx;
   % filter to avoid sharp transients
   det.ra_d = filtfilt(ones(1,settings.detuning.nfilt)/settings.detuning.nfilt, 1, det.ra_d);
   det.xa_d = filtfilt(ones(1,settings.detuning.nfilt)/settings.detuning.nfilt, 1, det.xa_d);
end

% shift resonance frequency
if settings.detuning.fshift > 0
   det.ra_d = interp1(fch, det.ra_d, fch + settings.detuning.fshift, 'linear','extrap');
   det.xa_d = interp1(fch, det.xa_d, fch + settings.detuning.fshift, 'linear','extrap');
end

% strong resonances unlikely for high frequency shifts (massive detuning) => dampening factor
det.att = max(0, 1 - settings.detuning.fshift_attf * settings.detuning.fshift);

% recombine with linear trend (apply)
ra = det.ra_l + det.ra_d * det.att;
xa = det.xa_l + det.xa_d * det.att;
%     antenna impedance over frequency
za = complex(ra, xa);

% output
disp(sprintf('   Resonance boosted by %.0f %%, shifted by %.0f MHz, and attenuated by %.2f dB.',...
   settings.detuning.res_en*100, settings.detuning.fshift/1e6, -20*log10(det.att)));

% figure;
% subplot(2,1,1); hold on;
% plot(fch/1e6, det.ra, 'b-');
% plot(fch/1e6, ra, 'r-');
% hold off; grid on; setlegend({'original', 'detuned'}, 'NorthEast');
% setlabels('DETUNING OF ANTENNA IMPEDANCE', 'f [MHz]', 'R_a [\Omega]');
% subplot(2,1,2); hold on;
% plot(fch/1e6, det.xa, 'b-');
% plot(fch/1e6, xa, 'r-');
% hold off; grid on; setlegend({'original', 'detuned'}, 'SouthEast');
% setlabels('DETUNING OF ANTENNA IMPEDANCE', 'f [MHz]', 'X_a [\Omega]');

% cleanup
clear('dummy', 'det');


% *******************************************************************************************************
% input reflection coefficient
disp(sprintf('\nComputing input reflection coefficient rho(fch, zmod, pic)...'));

% prepare matrices
%     input power
m_pin = zeros(settings.nfch, settings.nm, settings.np); % (fch, zmod, pic)
%     available power
m_pav = zeros(settings.nfch, settings.nm, settings.np); % (fch, zmod, pic)
%     reflection coefficient
rho   = zeros(settings.nfch, settings.nm, settings.np); % (fch, zmod, pic)

% calculate rho(fch, zmod, pin)
for i = 1 : settings.nfch % for: frequency (inside characteristic)
   % input impedance Zin
   %     chip input impedance (unmodulated)
   zic = complex(ric(i,:), -1./(2*pi*fch(i)*cic(i,:)));
   %     application/assembly tolerances (assembly + tuning + detuning)
   zat = 1 ./ ( 1/settings.rat + complex(0, 2*pi*fch(i)*settings.cat) );   
   %     parallel => input impedance (unmodulated)
   zin = zic .* zat ./ (zic + zat);
    
   for j = 1 : settings.nm % for: zmod
      % effective input power Pin out of Pic using efficiency factor
      zat_zmod = zat * squeeze(zmod(i, j, :)).' ./ ( zat + squeeze(zmod(i, j, :)).' ); % zmod uses power... 
      m_pin(i, j, :) = ( (abs(zat_zmod).^2 ./ real(zat_zmod)) + (abs(zic).^2 ./ ric(i,:)) ) ./...
         ( abs(zat_zmod).^2 ./ real(zat_zmod) ) .* pic';
      
      % input impedance (modulated; zin_mod = zin || zmod)
      zin_mod = zin .* squeeze(zmod(i, j, :)).' ./ ( zin + squeeze(zmod(i, j, :)).' );
          
      % reflection coefficient
      rho(i, j, :) = ( za(i) - conj(zin_mod) ) ./ ( za(i) + zin_mod );
      
      % available power Pav(fch, zmod, pic)
      m_pav(i, j, :) = ( squeeze(m_pin(i, j, :)) ./ ( 1 - squeeze(abs(rho(i, j, :))).^2 ) )';
   end
end

% check if m_pav > m_pin > pic (last index is pic for all matrices)
check.pav_pin = 0;
check.pin_pic = 0;
for i = 1 : settings.nfch % for: frequency (inside characteristic)
   for j = 1 : settings.nm % for: zmod
      if any( squeeze(m_pav(i, j, :)) <= squeeze(m_pin(i, j, :)) )
         check.pav_pin = check.pav_pin + 1;
      end
      if any( squeeze(m_pin(i, j, :)) <= pic )
         check.pin_pic = check.pin_pic + 1;
      end
   end
end
%     output
if check.pav_pin > 0
   disp(sprintf('   WARNING: found %.0f combinations fch/zmod with Pav <= Pin', check.pav_pin));
end
if check.pin_pic > 0
   disp(sprintf('   WARNING: found %.0f combinations fch/zmod with Pin <= Pic', check.pin_pic));
end

% cleanup
clear('i', 'zic', 'zat', 'zat_zmod', 'zin', 'zin_mod', 'check');

% figure(1); clf(1); hold on;
% % figure('Visible', 'on'); hold on; 
% surf(fch/1e9, pic, abs(squeeze(rho(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
% surf(fch/1e9, pic, abs(squeeze(rho(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
% hold off; grid on; axis tight; set(gca, 'YScale','log'); view([20, 30]); 
% setlabels('MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED / UNMODULATED)',...
%    'f [GHz]', 'P_{av} [W]', '|\rho|'); setcolorbar();
% 
% figure(2); clf(2); hold on;
% % figure('Visible', 'on'); hold on; 
% surf(fch/1e9, pic, 180/pi*angle(squeeze(rho(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
% surf(fch/1e9, pic, 180/pi*angle(squeeze(rho(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
% hold off; grid on; axis tight; set(gca, 'YScale','log'); view([20, 30]); 
% setlabels('PHASE OF REFLECTION COEFFICIENT (MODULATED / UNMODULATED)',...
%    'f [GHz]', 'P_{av} [W]', 'arg(\rho) [deg]'); setcolorbar();
% 
% % figure('Visible', 'on'); hold on;
% % surf(fch/1e9, pic, abs(squeeze(m_pin(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
% % surf(fch/1e9, pic, abs(squeeze(m_pin(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
% % hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log', 'ZScale','log');
% % setlabels('P_{in} TO P_{ic} LOOKUP TABLE  (MODULATED / UNMODULATED)',...
% %    'f [GHz]', 'P_{ic} [W]', 'P_{in} [W]'); setcolorbar();
% 
% pause(0.1); 


% *******************************************************************************************************
% iteratively generate pav vector and warp rho(fch, zmod, pic) -> rho_pav(fch, zmod, pav)
disp(sprintf('\nGenerating pav vector and warping rho(fch, zmod, pic) -> rho_pav(fch, zmod, pav)...'));

% prepare matrix for reflection coefficient (characteristic only)
rho_pav = zeros(settings.nfch, settings.nm, settings.np); % (fch, zmod, pav)

% logarithmic m_pav to be used for interpolation
m_pav_log = log(m_pav);

% reset counter: intermittent NaNs
intermittend_nans = 0;

% create pav vector with minimal possible unused parts ()
for it_pav = 1 : settings.pavit_maxit
   
%    it_pav
   
   % prepare vector for available power
   % ... linear indexing for log-distance model and a path-loss factor of settings.plf, 
   %     better resolution for operating area pav_max...pav_opmin
   % ... worst-case assumption for pav_opmin (min. picopmin(f), no reflection)
   pav_min   = max(min(m_pav(:)), settings.pavmin)^(-1/settings.plf);
   pav_max   = min(max(m_pav(:)), settings.pavmax)^(-1/settings.plf);
   pav_opmin = max(pav_min^(-settings.plf), min(picopmin))^(-1/settings.plf);
   %     distribute available number of points (out=min:opmin, in=opmin:max)
   settings.np_in  = ceil( settings.np * settings.npr_inop );
   settings.np_out = settings.np - settings.np_in;
   %     check created lengths
   if settings.np_in < 1 || settings.np_out < 0
      error('Unable to split available power according to settings.');
   end
   %     create linear pav vectors
   pav_in  = linspace(pav_max, pav_opmin, settings.np_in)';
   pav_out = linspace(pav_opmin, pav_min, settings.np_out+1)';
   %     assemble, filter and  warp to get linear resolution for log-distance model
   %     ... assembling has to be done mirrored because of ^(-settings.plf)
   pav = filtfilt(ones(settings.npma, 1)/settings.npma, 1, [pav_in; pav_out(2:end)]) .^ (-settings.plf);
   pav = pav(end:-1:1);
      
   % distort reflection coefficient from pic to pav
   % ... power vectors will be closer to log than to linear scale => interp in log domain
   % ... note that this interpolation is not robust against discontinuities in m_pav_log
   % ... interp1c (interpolation abs/angle instead of re/im creates more problems here than it solves:
   %     .) the calculation works with re/im (R, X), so does the reverse check below
   %     .) the grid for rho_pav will be dense enough => only minimal changes
   for i = 1 : settings.nfch % for: frequency (inside characteristic)
      for j = 1 : settings.nm % for: zmod
         rho_pav(i, j, :) = interp1(squeeze(m_pav_log(i, j, :)), squeeze(rho(i, j, :)), log(pav), 'linear', complex(NaN,NaN));
         rho_pav(i, j, abs(rho_pav(i, j, :)) > settings.rhomax) = settings.rhomax; % no active reflection (will introduce discontinuities though)
      end
   end

   % reset bounds if: only NaNs in the frequency dimension for at least one zmod
   ind = findzeros([1; squeeze(all(any(isnan(rho_pav), 2), 1)); 1] - 0.1); ind = ind - 1; % -0.1 => bias towards non-NaN   
   if length(ind) > 2 % handle intermittend NaNs: extract core part
      intermittend_nans = intermittend_nans + 1;
      [dummy, ind_maxdiff] = max(diff(ind));
      ind = ind(ind_maxdiff:ind_maxdiff+1);
   end
   if ind(1) > 1 || ind(2) < settings.np
      settings.pavmin = pav(ind(1));
      settings.pavmax = pav(ind(2));
      continue;
   end
   
   % reset bounds if: not enough support to fit IIR for any zmod
   %     copied check from below
   for i = 1 : settings.np % for: available power
      ind = findzeros([1; squeeze(any(isnan(rho_pav(:,:,i)), 2)); 1] - 0.1); ind = ind - 1; % -0.1 => bias towards non-NaN
      if length(ind) > 4; ind = ind([1,2,length(ind)-1,length(ind)]); end % outermost borders
      if length(ind) == 3; error('One-sided NaN distribution. This shouldn''t happen.'); end % outermost borders
      fit_prob(i) = (length(ind) == 4 && (min(diff(ind)) <   settings.minpts) || ...
                     length(ind) == 2 && (min(diff(ind)) < 2*settings.minpts)); %#ok<*AGROW> % [IIR_CHK] do not change this comment
   end
   %     try to adapt bounds if necessary (can only adapt min/max values)
   if any(fit_prob)
      ind = findzeros([1, fit_prob, 1] - 0.1); ind = ind - 1; % -0.1 => bias towards non-NaN
      if length(ind) > 2 % handle intermittend NaNs: extract core part
         intermittend_nans = intermittend_nans + 1;
         [dummy, ind_maxdiff] = max(diff(ind));
         ind = ind(ind_maxdiff:ind_maxdiff+1);
      end
      if ind(1) > 1 % reset lower bound (in log domain)
         if it_pav < settings.pavit_maxit
            settings.pavmin = 10 ^ (sum( log10(pav(ind(1)-1:ind(1)  )) .*...
               [1-settings.pavit_bias;   settings.pavit_bias] )); % move gradually
         else
            settings.pavmin = pav(ind(1)); % this is the last run => shift to next known value w/o fitting problems
         end
      end
      if ind(2) < settings.np % reset upper bound
         if it_pav < settings.pavit_maxit
            settings.pavmax = 10 ^ (sum( log10(pav(ind(2)  :ind(2)+1)) .*...
               [  settings.pavit_bias; 1-settings.pavit_bias] )); % move gradually
         else
            settings.pavmax = pav(ind(2)); % this is the last run => shift to next known value w/o fitting problems
         end
      end
      continue;
   end
   
   % continue with next step; allow for at least one iteration to get proper resolution for pav
   break;
end

% how many loops did we need?
disp(sprintf('    Needed %i iterations to determine bounds for pav.', it_pav));
if intermittend_nans > 0
   disp(sprintf('    Warning: Skipped %i occurences of intermittend NaNs in pav vector.', intermittend_nans));
end

% cleanup
clear('m_pav_log', 'i', 'j', 'pav_min', 'pav_max', 'pav_opmin', 'ind', 'fit_prob', 'it_pav', 'intermittend_nans');

% figure('Visible', 'on'); hold on;
% surf(fch/1e9, pav, abs(squeeze(rho_pav(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
% surf(fch/1e9, pav, abs(squeeze(rho_pav(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
% hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log');
% setlabels('MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED/UNMODULATED)',...
%    'f [GHz]', 'P_{av} [W]', '|\rho|'); setcolorbar();
% 
% figure('Visible', 'on'); hold on;
% surf(fch/1e9, pav, abs(squeeze(rho_pav(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
% hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log');
% setlabels('MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED)',...
%    'f [GHz]', 'P_{av} [W]', '|\rho|'); setcolorbar();
% 
% figure('Visible', 'on'); hold on;
% surf(fch/1e9, pav, 180/pi*angle(squeeze(rho_pav(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
% hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log');
% setlabels('MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED)',...
%    'f [GHz]', 'P_{av} [W]', 'angle(\rho) [deg]'); setcolorbar();
% 
% return


% ****************************************************************************************************
% extrapolate reflection coefficient using optimized IIR bandstop model
% (phase is not optimized because this is only used to filter noise outside known characteristic)
if ~settings.simchar
   disp(sprintf('\nSkipping extrapolation of reflection coefficient. Results are not usable for simulator.'));
else
   
   disp(sprintf('\nExtrapolating reflection coefficient rhoext_pav(f, zmod, pav) using optimized IIR bandstop models...'));
   
   % prepare matrices
   %     extrapolated reflection coefficient (will be extrapolated below)
   rhoext_pav = nan(settings.nf, settings.nm, settings.np); % (f, zmod, pav)
   rhoext_pav(1+settings.nf_out1:end-settings.nf_out4, :, :) = rho_pav;
   %     residual cost fcn
   reserr = inf(settings.nm, settings.np);
   %     fitting time [s]
   t_fitiir = NaN(settings.nm, settings.np);
   %     best setup
   bestiirset = NaN(settings.nm, settings.np);
   
   % fit IIR and extrapolate
   for i = 1 : settings.nm % for: zmod
      % init
      %     last optimized line (for smoothness)
      lastopt_abs = NaN;
      
      % start optimization at high power values (flatter frequency response)
      %    => better fitting => better results for smoothing
      for j = settings.np : -1 : 1 % for: pav
         % init
         %     start timer
         tic;
         %     initialize result/messages struct
         iiropt.result.z    = [];
         iiropt.result.p    = [];
         iiropt.result.k    = [];
         iiropt.result.x    = [];
         iiropt.result.disp = struct();
         
         % select curve rho(fch)
         rho_ch = squeeze(rho_pav(:, i, j));
         [minrho_ch, ind_minrho_ch] = min(abs(rho_ch));
         
         % we need at least a few points to optimize => skip if not present
         %     obtain borders of NaNs; select outermost positions
         ind = findzeros([1; isnan(rho_ch); 1] - 0.1); ind = ind - 1; % -0.1 => bias towards non-NaN
         if length(ind) > 4; ind = ind([1,2,length(ind)-1,length(ind)]); end
         %     skip optimization if at least one border is too narrow to fit the IIR
         % ... MAKE SURE TO COPY THIS CHECK TO PAV BORDER ITERATION; MARKED BY THE COMMENT "[IIR_CHK]"
         if length(ind) == 4 && (min(diff(ind)) <   settings.minpts) || ... % gap in between => "separate sides"
               length(ind) == 2 && (min(diff(ind)) < 2*settings.minpts) % no gap in between => "linked sides"
            disp(sprintf('    (:,%3i,%3i) No enough characteristic support to fit IIR. Skipping and setting rhoext_pav to NaN.', i, j));
            rhoext_pav(:, i, j) = complex(nan(settings.nf, 1, 1), nan(settings.nf, 1, 1));
            continue;
         end
         
         % optimize for different IIR setups, take best run
         for iirset_ind = 1 : length(settings.iir)
            % initialize loop messages
            %     main
            iiropt.no_support    = false;
            iiropt.reducing_gain = false;
            iiropt.dvio          = NaN;
            %     temp (overwritten at each run)
            iiropt.temp.no_support    = false;
            iiropt.temp.reducing_gain = false;
            iiropt.temp.dvio          = NaN;
            
            % create IIR bandpass filter
            %     set indices of NaN borders
            ind_fmin = ind(1);
            ind_fmax = ind(end);
            %     set minimum/maximum cutoff frequency of IIR plus check
            %     (use log scale of filter to create flatter slope at high frequency)
            fcmin = fch(ind_fmin+floor((ind_fmax-ind_fmin)*settings.iir{iirset_ind}.bshift));
            fcmax = fch(ind_fmax-floor((ind_fmax-ind_fmin)*settings.iir{iirset_ind}.bshift));
            %     check for fcmax < fcmin
            if fcmax < fcmin
               disp(sprintf('    (:,%3i,%3i) WARNING: fcmax < fcmin. Check settings.iir{%i}.bshift for sanity!',...
                  i, j, iirset_ind));
               dummy = fcmax; fcmax = fcmin; fcmin = dummy;
            end
            %     make sure fcmax-fcmin is not smaller than settings.iir{iirset_ind}.dfmin
            if fcmax - fcmin < settings.iir{iirset_ind}.dfmin
               freqdiff = settings.iir{iirset_ind}.dfmin-fcmax+fcmin;
               fcmax = fcmax + freqdiff/2;
               fcmin = fcmin - freqdiff/2;
            end
            %     shift stopband center frequency towards minimum of |rho_ch|
            fshift = fch(ind_minrho_ch) - (fcmax+fcmin)/2;
            fcmin = max(fcmin + fshift*settings.iir{iirset_ind}.cshift, f(1));
            fcmax = min(fcmax + fshift*settings.iir{iirset_ind}.cshift, f(end));
            %     calculate zpk-model of filter
            [iir.z, iir.p, iir.k] = butter(settings.iir{iirset_ind}.iirorder, [fcmin, fcmax]*2/(settings.fs), 'stop');
            
            %          [b,a] = zp2tf(iir.z, iir.p, iir.k);
            %          figure; hold on;
            %          plot(f*1e-9, abs(freqz(b, a, f, settings.fs)), 'r');
            %          plot(f*1e-9, abs(squeeze(rhoext_pav(:, i, j))), 'g-');
            %          plot(fch*1e-9, abs(rho_ch), 'b.-');
            %          hold off; grid on; xlim([0,2]); ylim([0,1]);
            
            % setup cost fcn for optimization
            %     get indices of minimum/maximum frequency not NaN
            settings.optimset.ind_fmin = ind_fmin;
            settings.optimset.ind_fmax = ind_fmax;
            %     desired values (inside characteristic and |rho|=settings.rhomax for start/end points)
            settings.optimset.cf_f = [f(1); fch(settings.optimset.ind_fmin:settings.optimset.ind_fmin+1);...
               fch(ind_minrho_ch); fch(settings.optimset.ind_fmax-1:settings.optimset.ind_fmax); f(end)];
            settings.optimset.cf_desired = [settings.rhomax;...
               abs(rho_ch(settings.optimset.ind_fmin:settings.optimset.ind_fmin+1)); abs(rho_ch(ind_minrho_ch));...
               abs(rho_ch(settings.optimset.ind_fmax-1:settings.optimset.ind_fmax)); settings.rhomax];
            %     set everything to standard (x0: [zero(1); zero(2);... pole(1); pole(2); ...; gain])
            settings.optimset.x0 = [ones(settings.iir{iirset_ind}.iirorder, 1)*settings.optimstart(1);...
               ones(settings.iir{iirset_ind}.iirorder, 1)*settings.optimstart(2); settings.optimstart(3)];
            settings.optimset.cf_w = settings.cf_w;
            optimres.oldresnorm = Inf;
            %     alterations of weights
            if any(isnan( abs(rho_ch(settings.optimset.ind_fmin:settings.optimset.ind_fmax)) )) % gaps in support
               settings.optimset.cf_w(1) = settings.optimset.cf_w(1) * settings.cf_walt; % make |rho(f->0)| more important
               settings.optimset.cf_w(7) = settings.optimset.cf_w(7) * settings.cf_walt; % make |rho(f->end)| more important
            end
            if abs(rho_ch(settings.optimset.ind_fmin)) < abs(rho_ch(settings.optimset.ind_fmin+1)) % negative slope (low frequ)
               settings.optimset.cf_w(1) = settings.optimset.cf_w(1) * settings.cf_walt; % make |rho(f->0)| more important
               settings.optimset.cf_w(3) = settings.optimset.cf_w(3) / settings.cf_walt; % and slope @ low frequ less
            end
            if abs(rho_ch(settings.optimset.ind_fmax)) < abs(rho_ch(settings.optimset.ind_fmax-1)) % negative slope (high frequ)
               settings.optimset.cf_w(7) = settings.optimset.cf_w(7) * settings.cf_walt; % make |rho(f->end)| more important
               settings.optimset.cf_w(5) = settings.optimset.cf_w(5) / settings.cf_walt; % and slope @ high frequ less
            end
            
            % check if desired values are non-NaN (if NaN: not enough characteristic support)
            if any( isnan(settings.optimset.cf_desired) )
               iiropt.temp.no_support = true;
               rhoext_pav(:, i, j) = complex(nan(settings.nf, 1, 1), nan(settings.nf, 1, 1));
               continue;
            else
               iiropt.temp.no_support = false;
            end
            
            % try to optimize pole/zero gain of filter to fit course of |rho_ch| at borders
            for optretries = 0 : settings.maxoptretries
               % run optimizer
               optimres = optimize_iir(settings, iirset_ind, iir, optimres, f, lastopt_abs);
               % calculate frequency response
               iir.h = freqz(optimres.b, optimres.a, f, settings.fs);
               %  check if |iir.h| > settings.rhomax by more than 1% or > 1 (outside characteristic support)
               if any( abs(iir.h([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end])) > min(1, settings.rhomax*(1+settings.rhomax_tol/100)) )
                  % increase punishment for |rho|<=settings.rhomax
                  settings.optimset.cf_w(end) = settings.optimset.cf_w(end) * settings.cf_wmult;
                  % and restart from where we left
                  settings.optimset.x0 = optimres.x;
                  % but with decreased gain
                  settings.optimset.x0(end) = settings.optimset.x0(end) /...
                     max(abs(iir.h([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end]))) * settings.rhomax;
               else
                  break;
               end
            end
            
            % check if |iir.h| > settings.rhomax by more than 1% or > 1 (outside characteristic support)
            if any( abs(iir.h([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end]) ) > min(1, settings.rhomax*(1+settings.rhomax_tol/100)) )
               % ... for message (at end of loop)
               iiropt.temp.reducing_gain = true;
               iiropt.temp.dvio = -min(1, settings.rhomax*(1+settings.rhomax_tol/100)) +...
                  max(abs(iir.h([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end])));
               % decrease gain multiplicator and gain factor
               optimres.x(end) = optimres.x(end) / max(abs(iir.h([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end]))) * settings.rhomax;
               optimres.k = optimres.k / max(abs(iir.h([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end]))) * settings.rhomax;
               % recalculate iir.h
               [optimres.b, optimres.a] = zp2tf(optimres.z, optimres.p, optimres.k);
               iir.h = freqz(optimres.b, optimres.a, f, settings.fs);
               
               %             disp(sprintf('    (:,%3i,%3i) shifted by %g, cost fcn before recalc: %g', i, j, iiropt.temp.dvio, sqrt(optimres.resnorm)));
               
            else
               iiropt.temp.reducing_gain = false;
            end
            
            % residual error: recalculate cost function, use std. weights for comparability
            % ... cost fcn cannot be accessed directly: use optimizer with 0 cost fcn calls
            settings.optimset.cf_w = settings.cf_w;
            settings.optimset.x0 = optimres.x;
            maxfunevals_old = settings.optimset.MaxFunEvals;
            settings.optimset.MaxFunEvals = 0;
            optimres = optimize_iir(settings, iirset_ind, iir, optimres, f, lastopt_abs);
            settings.optimset.MaxFunEvals = maxfunevals_old;
            
            % check if the residual error improved compared to last run (always true for first run)
            if reserr(i,j) > sqrt(optimres.resnorm) % improvement: store results
               
               %             disp(sprintf('    (:,%3i,%3i) run %i is better (old=%g, new=%g)', i, j, iirset_ind, reserr(i,j), sqrt(optimres.resnorm)));
               
               % residual error and index of best setting
               reserr(i,j) = sqrt(optimres.resnorm);
               bestiirset(i,j) = iirset_ind;
               % optimized values: zeros/poles
               iiropt.result.z = optimres.z;
               iiropt.result.p = optimres.p;
               % optimized values: gain
               iiropt.result.k = optimres.k;
               % optimized values: multipliers for zeros/poles [zeros(1:end/2); poles(1:end/2)]
               iiropt.result.x = optimres.x;
               % copy messages (we took this setting)
               iiropt.no_support    = iiropt.temp.no_support;
               iiropt.reducing_gain = iiropt.temp.reducing_gain;
               iiropt.dvio          = iiropt.temp.dvio;
               
            else % no improvement: just continue (or exit loop if this was the last setting)
               
               %             disp(sprintf('    (:,%3i,%3i) run %i is worse (old=%g, new=%g)', i, j, iirset_ind, reserr(i,j), sqrt(optimres.resnorm)));
               
               continue;
            end
            
            %          [b,a] = zp2tf(iir.z, iir.p, iir.k);
            %          figure; hold on;
            %          plot(f*1e-9, abs(freqz(b, a, f, settings.fs)), 'r.-');
            %          plot(fch*1e-9, abs(rho_ch), 'b.-');
            %          hold off; grid on; xlim([0,2]); ylim([0,1]);
            %          continue
            
         end
         
         % display messages
         if iiropt.no_support
            disp(sprintf(...
               '    (:,%3i,%3i) No enough support to fit IIR (NaNs in desired). Skipping and setting rhoext_pav to NaN.', i, j));
         end
         if iiropt.reducing_gain
            disp(sprintf(...
               '    (:,%3i,%3i) |rho|>min(1,rho_max+tol) by %g after %i retries. Reducing gain, keeping poles/zeros.',...
               i, j, iiropt.dvio, settings.maxoptretries));
         end
         
         % recalculate frequency response
         [iiropt.result.b, iiropt.result.a] = zp2tf(iiropt.result.z, iiropt.result.p, iiropt.result.k);
         iiropt.result.h = freqz(iiropt.result.b, iiropt.result.a, f, settings.fs);
         
         % add IIR frequency responses to reflection coefficient
         ind_f_replace = [1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:settings.nf];
         rhoext_pav(ind_f_replace, i, j) = iiropt.result.h(ind_f_replace);
         
         % phase: try to remove "jumps" at borders of interpolation (shift by average phase difference at borders)
         %        cut @ low frequency
         dphase = angle(rhoext_pav(settings.nf_out1+settings.optimset.ind_fmin, i, j)) - angle(rhoext_pav(settings.nf_out1+settings.optimset.ind_fmin-1, i, j));
         rhoext_pav(1:settings.nf_out1+settings.optimset.ind_fmin-1, i, j) =...
            rhoext_pav(1:settings.nf_out1+settings.optimset.ind_fmin-1, i, j) .* exp(complex(0, 1)*dphase);
         %        cut @ high frequency
         dphase = angle(rhoext_pav(settings.nf_out1+settings.optimset.ind_fmax, i, j)) - angle(rhoext_pav(settings.nf_out1+settings.optimset.ind_fmax+1, i, j));
         rhoext_pav(settings.nf_out1+settings.optimset.ind_fmax+1:end, i, j) =...
            rhoext_pav(settings.nf_out1+settings.optimset.ind_fmax+1:end, i, j) .* exp(complex(0, 1)*dphase);
         
         % try to create smooth border, but modify only extrapolated parts
         for k = 1 : settings.bordermed_len
            rhoext_pav(:, i, j) = medfilt1(rhoext_pav(:, i, j), settings.bordermed_len);
            % restore parts that should not be changed
            ind_f_replace = 1 : settings.nf;
            ind_f_replace([1:settings.nf_out1+settings.optimset.ind_fmin, settings.nf_out1+settings.optimset.ind_fmax+1:settings.nf]) = [];
            rhoext_pav(ind_f_replace, i, j) = squeeze(rho_pav(ind_f_replace-settings.nf_out1, i, j));
         end
         
         % store magnitude of this line for the next run
         lastopt_abs = abs(rhoext_pav(:, i, j));
         
         % stop timer, store time
         t_fitiir(i, j) = toc;
         
         %       figure; hold on;
         %       plot(fch*1e-9, abs(rho_ch), 'b.-');
         %       plot(f*1e-9, abs(squeeze(rhoext_pav(:, i, j))), 'ro-');
         %       hold off; grid on; axis tight; xlim([0,4]);
         %
         %       figure; hold on;
         %       plot(fch*1e-9, angle(rho_ch), 'b.-');
         %       plot(f*1e-9, angle(squeeze(rhoext_pav(:, i, j))), 'ro-');
         %       hold off; grid on; axis tight; xlim([0,2]);
         %
         %       return
         
      end
   end
   
   % replace Inf in reserr by NaN
   reserr(isinf(reserr)) = NaN;
   
   % warning if residual error is too high
   %     very bad
   if 20*log10(max(max(reserr))) > settings.reserr_critwarn
      disp(sprintf('   WARNING: Residual error of IIR fitting exceeds critial limits (|%g| > %g [dB]). Fitting might be very bad!',...
         20*log10(max(max(reserr))), settings.reserr_critwarn));
      %     not too bad
   elseif 20*log10(max(max(reserr))) > settings.reserr_warn
      disp(sprintf('   Warning: Residual error of IIR fitting exceeds limits (%g > %g [dB]).',...
         20*log10(max(max(reserr))), settings.reserr_warn));
   end
   
   % cleanup
   clear('i', 'j', 'k', 'dummy', 'rho_ch', 'minrho_ch', 'ind_minrho_ch', 'ind', 'ind_fmin', 'ind_fmax', 'fcmin',...
      'fcmax', 'freqdiff', 'drho_chmin', 'drho_chmax', 'fshift', 'lastopt_abs', 'iir', 'iiropt', 'iirset_ind',...
      'optretries', 'optimres', 'maxfunevals_old', 'dphase', 'ind_f_replace');
   
   % figure('Visible', 'on'); hold on;
   % surf(f/1e9, pav, abs(squeeze(rhoext_pav(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
   % surf(f/1e9, pav, abs(squeeze(rhoext_pav(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
   % hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log');
   % xlim([0, 2*max(fch)]/1e9);
   % setlabels(MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED/UNMODULATED)',...
   %    'f [GHz]', 'P_{av} [W]', '|\rho|'); setcolorbar();
end

% *******************************************************************************************************
% reverse check rho_pav(fch, zmod, pav) -> Zic
disp(sprintf('\nReverse check: rho_pav(fch, zmod, pav) -> Zic...'));

% preparations
err     = zeros(settings.checks, 2); % [abs, angle in deg]
err_ind = zeros(settings.checks, 3); % [fch, zmod, pav]

% settings.cecks random setups
for i = 1 : settings.checks
   
   % create random index [fch, zmod, pav] inside characteristic excluding borders ("demilitarized zone DMZ")
   while(1)
      % random index
      ind_fch = 1 + settings.dmz_f + round(rand * (settings.nfch - 2*settings.dmz_f - 1));
      ind_m   = 1 + settings.dmz_m + round(rand * (settings.nm   - 2*settings.dmz_m - 1));
      ind_pav = 1 + settings.dmz_p + round(rand * (settings.np   - 2*settings.dmz_p - 1));
      % dump current set if any rho_pav in the vicinity (dmz) is NaN
      %     calculate DMZ
      dmz_fch = [max(1, ind_fch - settings.dmz_f) : min(ind_fch + settings.dmz_f, settings.nfch)];
      dmz_m   = [max(1, ind_m   - settings.dmz_m) : min(ind_m   + settings.dmz_m, settings.nm  )];
      dmz_pav = [max(1, ind_pav - settings.dmz_p) : min(ind_pav + settings.dmz_p, settings.np  )];
      %     if no NaNs: keep this set of indices
      if ~any(any(any(isnan( rho_pav(dmz_fch, dmz_m, dmz_pav) ))))
         break;
      end
   end
   %     store indices for later analysis
   err_ind(i, :) = [ind_fch, ind_m, ind_pav];
      
   % first way: calculate Pin using the reflection coefficient
   pin_hat = pav(ind_pav) * ( 1 - abs(rho_pav(ind_fch, ind_m, ind_pav))^2 );
   %     ... its index is also to index for zic => can be used to calculate zic directly
   %         (also checks m_pin)
   ric_hat  = interp1(squeeze(m_pin(ind_fch, ind_m, :)), ric(ind_fch, :), pin_hat, 'spline', 'extrap');
   cic_hat  = interp1(squeeze(m_pin(ind_fch, ind_m, :)), cic(ind_fch, :), pin_hat, 'spline', 'extrap');
   zic_hat1 = complex(ric_hat, -1./(2*pi*fch(ind_fch)*cic_hat));  
   
   % second way: reverse reflection coefficient calculation
   %     get real and imaginary part of reflection coefficient first
   rho_r = real(rho_pav(ind_fch, ind_m, ind_pav));
   rho_i = imag(rho_pav(ind_fch, ind_m, ind_pav));
   %     reverse rho = (za-zin*)/(za+zin) formula
   zin_hat = [1+rho_r, -rho_i; rho_i, rho_r-1] \ [1-rho_r, rho_i; -rho_i, 1-rho_r] * ...
      [real(za(ind_fch)); imag(za(ind_fch))]; % [rin; xin]
   zin_hat = complex(zin_hat(1), zin_hat(2));
   %     get modulation impedance from method #1 (unfortunately we need the power level)
   zmod_hat = interp1(squeeze(m_pin(ind_fch, ind_m, :)), squeeze(zmod(ind_fch, ind_m, :)), pin_hat, 'spline', 'extrap');
   %     calculate assembly impedance at this frequency
   zat_hat = 1 / (1/settings.rat + complex(0, 2*pi*fch(ind_fch)*settings.cat));
   %     zin, zat and zmod are known => calculate zic
   zic_hat2 = 1 / (1/zin_hat - 1/zat_hat - 1/zmod_hat);
   
   % those two values have to match => calculate relative mismatch witch zic_hat2 as reference
   % (relative error is ok because real/imag will not be close to zero)
   err(i, 1) = ( real(zic_hat1) - real(zic_hat2)) / real(zic_hat2);
   err(i, 2) = ( imag(zic_hat1) - imag(zic_hat2)) / imag(zic_hat2); 
   
end

% warning message if error ...
%     ... exceeds limits
if max(abs(err(:,1))) > settings.maxerr
   disp(sprintf('   Warning: Maximum relative error exceeded limits (real part, %g > %g).',...
      max(abs(err(:,1))), settings.maxerr));
end
if max(abs(err(:,2))) > settings.maxerr
   disp(sprintf('   Warning: Maximum relative error exceeded limits (imaginary part, %g > %g).',...
      max(abs(err(:,2))), settings.maxerr));
end
%     ... is not zero-mean
if mean(abs(err(:,1))) > settings.maxavg
   disp(sprintf('   Warning: Relative error is not zero-mean (real part, |%g| > %g).',...
      mean(abs(err(:,1))), settings.maxavg));
end
if mean(abs(err(:,2))) > settings.maxavg
   disp(sprintf('   Warning: Relative error is not zero-mean (imaginary part, |%g| > %g).',...
      mean(abs(err(:,2))), settings.maxavg));
end

% cleanup
clear('i', 'ind_fch', 'ind_m', 'ind_pav', 'dmz_fch', 'dmz_m', 'dmz_pav',...
   'pin_hat', 'ric_hat', 'cic_hat', 'zmod_hat', 'zat_hat', 'zic_hat1', 'zic_hat2', 'rho_r', 'rho_i', 'zin_hat');
%     also rat and cat are not used any more (besides: they can be found in settings anyway)
clear('rat', 'cat');


% *******************************************************************************************************
% fill leftover gaps by fitting the last line that has no gaps (again: only magnitude!)
if settings.simchar
   disp(sprintf('\nFilling leftover gaps (NaNs) in extrapolated characteristic...'));
   
   for i = 1 : settings.nm % for: zmod
      for j = 1 : settings.np % for: pav
         
         % find indices of last non-NaN values [left, right]
         ind_borders = findzeros(isnan(squeeze(rhoext_pav(:, i, j))) - 0.1); % -0.1 => bias towards non-NaN
         %     found no indices => no/all NaN
         if isempty(ind_borders)
            continue;
         end
         %     create binary vector for simple indexing
         ind_nan = false(settings.nf, 1);
         for k = 1 : length(ind_borders) / 2
            ind_nan(ind_borders(2*k-1) : ind_borders(2*k)) = true;
         end
         
         % fill gap using interpolation of magnitude and phase
         %     interpolate
         replacement = interp1c(f(~ind_nan), squeeze(rhoext_pav(~ind_nan, i, j)), f, 'spline');
         %    if > rhomax: switch from spline to cubic interpolation
         if any( abs(replacement(ind_nan)) > settings.rhomax )
            replacement = interp1c(f(~ind_nan), squeeze(rhoext_pav(~ind_nan, i, j)), f, 'cubic');
         end
         %    if > rhomax: switch from cubic to linear interpolation
         if any( abs(replacement(ind_nan)) > settings.rhomax )
            replacement = interp1c(f(~ind_nan), squeeze(rhoext_pav(~ind_nan, i, j)), f, 'linear');
         end
         %     copy (only necessary parts)
         rhoext_pav(ind_nan, i, j) = replacement(ind_nan);
         
      end
   end
   
   
   % provide some information concerning the results
   %     calculate
   info.nans       = sum(sum(sum(isnan(rhoext_pav))));
   info.nofitiir   = sum(sum(nansum(rhoext_pav, 1)==0)); % ==0 unlikely if not all NaNs
   info.unused_pav = sum(nansum(nansum(rhoext_pav, 1), 2)==0); % ==0 unlikely if not all NaNs
   %     output
   disp('   Final extrapolated characteristic contains:');
   disp(sprintf('      %6i (%6.2f %%) NaNs', info.nans, info.nans/numel(rhoext_pav)*100));
   disp(sprintf('      %6i (%6.2f %%) unused rows (not enough support to fit IIR)',...
      info.nofitiir,   info.nofitiir/(size(rhoext_pav,2)*size(rhoext_pav,3))*100));
   disp(sprintf('      %6i (%6.2f %%) unused power values (skipped IIR fitting for all zmod)',...
      info.unused_pav, info.unused_pav/settings.np*100));
   
   % cleanup
   clear('i','j','k', 'ind_borders','ind_nan','replacement');
   
   % figure('Visible', 'on'); hold on;
   % surf(f/1e9, pav, abs(squeeze(rhoext_pav(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
   % surf(f/1e9, pav, abs(squeeze(rhoext_pav(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
   % hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log');
   % xlim([0.5*min(fch), 1.25*max(fch)]/1e9);
   % setlabels('MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED/UNMODULATED)',...
   %    'f [GHz]', 'P_{av} [W]', '|\rho|'); setcolorbar();
   %
   % figure('Visible', 'on'); hold on;
   % surf(f/1e9, pav, angle(squeeze(rhoext_pav(:, 1, :)))',   'EdgeColor','black', 'FaceColor','interp');
   % surf(f/1e9, pav, angle(squeeze(rhoext_pav(:, end, :)))', 'EdgeColor','none',  'FaceColor','interp');
   % hold off; grid on; axis tight; view([20, 30]); set(gca, 'YScale','log');
   % xlim([0.5*min(fch), 1.25*max(fch)]/1e9);
   % setlabels('ANGLE OF REFLECTION COEFFICIENT (MODULATED/UNMODULATED)',...
   %    'f [GHz]', 'P_{av} [W]', '|\rho|'); setcolorbar();
end

% *******************************************************************************************************
% save / return results and reactivate warnings
if settings.simchar
   disp(sprintf('\nSaving results...'));
end

% one last cleanup
clear('ans');

% save characteristic
%     prepare information header
matfilename = sprintf('tagchar_modulator%s', settings.suffix);
if settings.simchar
   characteristic = sprintf('%s\n%s\n%s\n%s',...
      'reflection coefficient rho_pav(fch,zmod,pav)',...
      'extrapolated reflection coefficient rhoext_pav(f,zmod,pav)',...
      'input power vs. chip state pin(fch,zmod,pic)',...
      'chip input impedance zic(fch,pic) and capacity cic(fch,pic)',...
      'minimum operational power picopmin(fch)'); 
else
   characteristic = sprintf('%s\n%s\n%s\n%s',...
      'reflection coefficient rho_pav(fch,zmod,pav)',...
      'available power vs. chip state m_pav(fch,zmod,pic)',...
      'minimum operational power picopmin(fch)'); 
end
createdon    = datestr(now, 0); %#ok<NASGU>
createdby    = sprintf('%s.m, rev.: %s', mfilename, version); 
usedmatfiles = sprintf('%s (%s)\n%s (%s)\n%s (%s)\n%s (%s)',...
   meas_tag_za.matfilename, meas_tag_za.createdby,...
   meas_tag_zic.matfilename, meas_tag_zic.createdby,...
   meas_tag_zmod.matfilename, meas_tag_zmod.createdby,...
   meas_tag_picmin.matfilename, meas_tag_picmin.createdby); 
assembly = sprintf('Used assembly impedance: %.0f fF parallel to %.0f Ohm.', settings.cat*1e15, settings.rat);
detuning = sprintf('Detuning: Resonance enhanced by %.0f %% and shifted by %.0f MHz.', settings.detuning.res_en*100, settings.detuning.fshift/1e6);
%     other data
zic = complex(ric, xic); %#ok<NASGU>
pin = m_pin; %#ok<NASGU>
%     save
if settings.writeres
   if settings.simchar
      % save results to file, do not return anything
      lightchar = struct();
      simchar   = true; %#ok<NASGU> ... mark file as suitable for simulator 
      ads_mapping = load(settings.charfiles.ads_mapping); % also attach the assembly/detuning states
      save(fullfile(settings.folders.res, matfilename), 'matfilename','characteristic','createdby','usedmatfiles','assembly','detuning',...
         'f','fch','zmod','pic','pav', 'picopmin', 'rho_pav','rhoext_pav',...
         'pin','zic','cic', 'settings', 'simchar', 'ads_mapping');
   else
      % do not save results, but return in struct
      lightchar = struct('matfilename',matfilename, 'characteristic',characteristic, 'createdby',createdby,...
         'usedmatfiles',usedmatfiles, 'assembly', assembly, 'detuning', detuning,...
         'fch',fch, 'pic',pic, 'pav',pav, 'picopmin',picopmin, 'rho_pav',rho_pav,...
         'm_pav',m_pav, 'settings',settings);
   end
end

% save entire workspace
% 1) for debugging purposes and plots (it takes hours to run this function)
% 2) Matlab may crash when no x-server is available (detached or -nodisplay) and saveas is called
%    (at least the figures are created with very low quality)
%    => create figures in a separate function
if settings.simchar
   save(fullfile(settings.folders.wsp, sprintf('tagchar_mod%s_workspace', settings.suffix)));
end

% stop recording
diary('off');

end



% *******************************************************************************************************
% *******************************************************************************************************
% optimize IIR bandstop by optimizing gain factors for poles and zeros (=> minimize optimizer_costfcn)
function optimres = optimize_iir(settings, iirset_ind, iir, optimres, f, lastopt_abs)

% complete upper/lower bound for optimization [zero(1); zero(2);... pole(1); pole(2); ...; gain]
ub = [ones(settings.iir{iirset_ind}.iirorder, 1)*settings.optimset.ub(1);...
      ones(settings.iir{iirset_ind}.iirorder, 1)*settings.optimset.ub(2); settings.optimset.ub(3)];
lb = [ones(settings.iir{iirset_ind}.iirorder, 1)*settings.optimset.lb(1);...
      ones(settings.iir{iirset_ind}.iirorder, 1)*settings.optimset.lb(2); settings.optimset.lb(3)];
   
% optimize
[optimres.x, optimres.resnorm] =...
   lsqnonlin(@optimizer_costfcn, settings.optimset.x0, lb, ub, settings.optimset);

% preliminarily adjust poles/zeros (conjugate complex)
for j = 1: settings.iir{iirset_ind}.iirorder %#ok<FXUP>
   % zeros (create column vector)
   optimres.z(2*j-1, 1) = iir.z(2*j-1) * optimres.x(j);
   optimres.z(2*j  , 1) = conj(optimres.z(2*j-1));
   % poles (create column vector)
   optimres.p(2*j-1 ,1) = iir.p(2*j-1) * optimres.x(j+settings.iir{iirset_ind}.iirorder);
   optimres.p(2*j   ,1) = conj(optimres.p(2*j-1));
end

% preliminary adjust gain
optimres.k = iir.k * optimres.x(end);

% calculate transfer function coefficients
[optimres.b, optimres.a] = zp2tf(optimres.z, optimres.p, optimres.k);



% *******************************************************************************************************
% *******************************************************************************************************
% NESTED: cost function for IIR bandstop optimization
% explanation of settings.optimset:
%    .cf_f:        vector of interesting frequencies (cf. .cf_desired below)
%    .cf_desired:  desired values (complex to match phase too)
%                    1 rho@fmin (should be close to 1)
%                  2:3 first two values of characteristic (left = low frequencies)
%                    4 center value
%                  5:6 last two values of characteristic (right = high frequencies)
%                    7 rho@fs/s (should be close to 1)
%   .cf_w:         weights for [.cf_desired; smoothness of magnitude, |rho|>rhomax]
   function cost = optimizer_costfcn(x)
      % create local copies
      zeros = iir.z;
      poles = iir.p;

      % adjust poles/zeros (conjugate complex)
      for j = 1: settings.iir{iirset_ind}.iirorder %#ok<FXUP>
         % zeros
         zeros(2*j-1) = iir.z(2*j-1) * x(j);
         zeros(2*j)   = conj(zeros(2*j-1));
         % poles
         poles(2*j-1) = iir.p(2*j-1) * x(j+settings.iir{iirset_ind}.iirorder);
         poles(2*j)   = conj(poles(2*j-1));
      end
      
      % adjust gain
      gain = iir.k * x(end);
      
      % keep stable under all circumstances
      if any( abs(zeros) < 0 ) || any( abs(zeros) >= 1 ) ||...
         any( abs(poles) < 0 ) || any( abs(poles) >= 1 )
         cost = ones(size(settings.optimset.cf_w)) / eps; return;
      end
      
      % calculate filter transfer function and frequency response
      [b, a] = zp2tf(zeros, poles, gain);
      %     inside characteristic (fit to borders of curve)
      h1 = abs(freqz(b, a, settings.optimset.cf_f, settings.fs));
      %     outside characteristic (settings.optimset.cf_f may be arbitrary => separate calc)
      h2 = abs(freqz(b, a, f([1:settings.nf_out1+settings.optimset.ind_fmin-1,...
         settings.nf_out1+settings.optimset.ind_fmax+1:end]), settings.fs));

      % cost function
      %     inside characteristic (fit to borders of curve)
      cost1 = ([h1(1); h1(2:3); h1(4); h1(5:6); h1(7)] - settings.optimset.cf_desired) .* settings.optimset.cf_w(1:end-2);
      %     enforce smooth amplitude changes outside characteristic
      if ~isnan(lastopt_abs)
         cost2 = sqrt(sum( (h2 - lastopt_abs([1:settings.nf_out1+settings.optimset.ind_fmin-1,...
            settings.nf_out1+settings.optimset.ind_fmax+1:end])).^2 )) / length(h2) * settings.optimset.cf_w(end-1);
      else
         cost2 = 0;
      end
      %     punishment for |rho| > settings.rhomax outside characteristic
      cost3 = sum( h2(find(abs(h2) > settings.rhomax)) - settings.rhomax) /...
         (settings.nf - settings.optimset.ind_fmax + settings.optimset.ind_fmin - 1) * settings.optimset.cf_w(end); %#ok<FNDSB>
      %     sum
      cost = [cost1; cost2; cost3];
      
%       figure; hold on;
%       plot(f([1:settings.nf_out1+settings.optimset.ind_fmin-1, settings.nf_out1+settings.optimset.ind_fmax+1:end]), h2, 'r.');
%       plot(f, lastopt_abs, 'b.');
%       hold off; grid on;
      
   end

end



% *******************************************************************************************************
% *******************************************************************************************************
% OLD AND RUSTY

% % extrapolation help for Pmin
% disp(sprintf('   WARNING: Providing manual extrapolation help.'));
% meas_tag_picmin.f = [2*meas_tag_picmin.f(1)-meas_tag_picmin.f(2); meas_tag_picmin.f;...
%    2*meas_tag_picmin.f(end)-meas_tag_picmin.f(end-1)];
% meas_tag_picmin.p = [meas_tag_picmin.p(1)-0.67*(meas_tag_picmin.p(2)-meas_tag_picmin.p(1)); meas_tag_picmin.p;
% meas_tag_picmin.p(end)];


%       pav = linspace(max(min(m_pav(:)), settings.pavmin)^(-1/settings.plf),...
%                      min(max(m_pav(:)), settings.pavmax)^(-1/settings.plf), settings.np)' .^
%                      (-settings.plf);
