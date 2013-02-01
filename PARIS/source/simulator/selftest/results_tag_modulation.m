% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_tag_modulation function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
%
% Note that the selftest of tag_modulation uses artificial characteristics that have no physical
% background (for simplicity), i.e., Pav, Pin and Pic are not connected to Za, Zat, and Zic. These
% characteristics are re-used by the selftests of reader_main and tag_main.
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
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
%
%
% ***** Todo *****
%
% *******************************************************************************************************


% *******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% initialize global stuff (assumes that globalinit.m is part of path)
globalinit('silent'); % do not print any messages (some functions will produce warnings here)

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% test function intrand()
disp('Running some quick checks on integer random number generator (INTRAND).');
if ~intrand_isok
   error('Integer random number generator does not seem to be functional. Rerun to verify.');
else
   disp('   ...passed.');
end


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_tag_modulation';

% test setup
%     margin of error for carrier length estimate (added twice)
%     ... tag_modulation has a more sophisticated length calculation; to compensate for the difference
%         (e.g., rounding of fclk to fs, LF to fclk, clock jitter, etc.), this margin (ratio) is added twice
testsettings.lmargin = 0.05;
%     number of magnitude/phase estimator windows per LF period
testsettings.est_nwin = 8; % (good: 4...8; make even)

% partitioning
partitioning.names = {'modulate', 'modulate (forced)', 'filter', 'exception handling'};
partitioning.runs  = [20, 10, 20, 7];

% tag_modulation characteristic settings [min, max] for synthetic characteristics 1 and 2
%     filename
charsettings_mod{1}.filename = 'tagchar_modulator_synth01';
charsettings_mod{2}.filename = 'tagchar_modulator_synth02';
%     size of characteristics (number of points)
charsettings_mod{1}.nf   = 192; % frequency (full range)
charsettings_mod{2}.nf   = 168;
charsettings_mod{1}.nfch = 128; % frequency (nonextrap. range); nf - nfch > 2 !
charsettings_mod{2}.nfch = 129;
charsettings_mod{1}.np   =  96; % power
charsettings_mod{2}.np   = 127;
charsettings_mod{1}.nm   =   2; % modulation impedance; nr = 2 !
charsettings_mod{2}.nm   =   2;
%     available power
charsettings_mod{1}.pav  = [1e-5, 5e-1]; % W
charsettings_mod{2}.pav  = [2e-5, 7e-1]; % W
%     minimum chip input power for functionality; coefficients [a, b] for Pmin[dBm]=a*f[GHz]+b   
charsettings_mod{1}.pmin_coeff = [ -5, -15]; % W
charsettings_mod{2}.pmin_coeff = [-10, -10]; % W
%     frequency; nonextrap. characteristic
charsettings_mod{1}.fch  = [700e6, 1200e6]; % Hz
charsettings_mod{2}.fch  = [750e6, 1000e6]; % Hz
%     frequency; extrap. characteristic (=> sampling frequency)
charsettings_mod{1}.f    = [1e3, 8.2e9]; % Hz ... make {1}.f(2) smaller than {2}.f(2); set > max(fs)/2
charsettings_mod{2}.f    = [1e2,  16e9]; % Hz ... make {2}.f(2) larger  than {1}.f(2); set > max(fs)/2
%     modulation impedance
charsettings_mod{1}.zmod = [complex(150,10), complex(1e4,1e4)]; % Ohm
charsettings_mod{2}.zmod = [complex(200,20), complex(5e4,1e4)]; % Ohm
%     assembly impedance [min, max]
charsettings_mod{1}.rat = [ 10,  500]; % Ohm
charsettings_mod{1}.cat = [ 50, 1000]*1e-15; % F
charsettings_mod{2}.rat = [ 25,  750]; % Ohm 
charsettings_mod{2}.cat = [  5,  750]*1e-15; % F
%     [start, end] values for "nonextrap. char." borders (make sure phase jumps mod-umod are small)
charsettings_mod{1}.val1 =    [0.10, 0.70;  0.6, 0.9]; % magnitude [unmod; mod]
charsettings_mod{1}.val2 = pi*[0.10, 0.90; -0.4, 0.8]; % phase [unmod; mod]
charsettings_mod{2}.val1 =    [0.05, 0.25;  0.3, 0.7]; % real part [unmod; mod]
charsettings_mod{2}.val2 =    [0.15, 0.50;  0.1, 0.7]; % imaginary part [unmod; mod]
%     [start, end] values for chip input impedance (interpolation will be re/im for max variability)
charsettings_mod{1}.val34 = [complex(  5, -23), complex(99, -223)]; % [Ohm]
charsettings_mod{2}.val34 = [complex(127, -77), complex( 7,  -13)]; % [Ohm]
%     maximum reflection coefficient (for "extrapolated" areas)
charsettings_mod{1}.rhomax = complex(0.75, 0.64);
charsettings_mod{2}.rhomax = 0.99;
%     tag is not functional in the range (ratio) nfunc_pav of pav [min, max] for all
%     fch(1:end*nfunc_fch(1)) and fch(1+end*(1-nfunc_fch(2):end)
%     ... use info struct to tune settings (all reasons for a nonfunctional tag should be present)
charsettings_mod{1}.nfunc_pav = [0.25, 0.31];
charsettings_mod{1}.nfunc_fch = [0.05, 0.93];
charsettings_mod{2}.nfunc_pav = [0.69, 0.71];
charsettings_mod{2}.nfunc_fch = [0.06, 0.95];

% tag_power characteristic settings [min, max] for synthetic characteristic
%     filename
charsettings_pwr.filename = 'tagchar_power_synth';
%     size of characteristics (number of points)
charsettings_pwr.n  = 256;
%     available power
charsettings_pwr.pic = [-20, 10]; % dBm
%     vdda (names below in randsettings.pwrchars)
%     ATTENTION: do not include zero (test_tag_modulation checks relative error)
charsettings_pwr.vdda{1} = [  1,  10]; % V at min/max power
charsettings_pwr.vdda{2} = [ 12, 0.3]; % V at min/max power

% standard (nonrandom) setup
%     create a new (base) carrier for each single test?
stdsettings.newcarrier = false;
%     data to modulate; will be truncated randomly below (see randsettings); 
%     (random data might violate protocol and prevent correct timing detection!)
%     .) DO NOT place more than 3 ones/zeros in a row!
%     .) has to start with 1 for t0 detection
%     .) a trailing 1 will be added below to allow for trailing t0 estimation; make sure this does not
%        create protocol violations
stdsettings.data = [1;1;0;1;0;1 ; 1;0;1;1;0;0;0;1;1]; % READ THE NOTES ABOVE
%     carrier setup (frequency will be set below)
stdsettings.reader.oscillator.fstddev = 10e3; % Hz
stdsettings.reader.oscillator.astddev = 1e-7; % amplitude instability (variance)
stdsettings.reader.oscillator.snr     =   80; % dB
%     clock (tag)
stdsettings.tag.clock.fsigma    =    5/3; % percent 
stdsettings.tag.clock.phi0      =     -1; % deg (-1: random)
stdsettings.tag.clock.mode      = 'fclk'; % clock used to generate a signal
%     parameter estimation
stdsettings.other.est_sinusoid.ol =  50; % percent (will be adjusted later)
%     tag_modulation (function under test)
stdsettings.tag.modulation.nfilt       =    32; % samples (channel filters)
stdsettings.tag.modulation.nfft        =  4096; % channels (double-sided spectrum)
stdsettings.tag.modulation.imagtol     =  1e-5; % maximum average imag/real after filterbank synthesis
stdsettings.tag.modulation.ptag_f      =   NaN; % force tag power level (set to NaN for normal operation)

% random setups [min, max]
randsettings.data = [    6, length(stdsettings.data)]; % number of data bits
randsettings.lf   = [ 25e3, 100e3]; % Hz backscatter link frequency {40...640 kHz} [keep small!]
randsettings.fs   = [  3e9,  11e9]; % Hz sampling frequency
randsettings.fc   = [860e6, 960e6]; % Hz carrier frequency (will be set per test, not per run)
randsettings.fclk = [1.5e6, 2.3e6]; % Hz tag clock frequency
randsettings.t0   = [ 1e-5,  5e-5]; % s of unmodulated carrier before and after modulation
randsettings.pav  = [ 9e-6,   1e0]; % W available carrier power (=> for carrier gain)

% random setup: lists
%     power supply characteristic (ATTENTION: WHEN CHANGED, ALSO CHANGE THE SAVED CHAR. AT THE BOTTOM)
randsettings.pwrchars = {'vdda1', 'vdda2'}; 

% test setup for selected tests (creates exceptions)
% Tested: t0=0, smaller nfft (filterbank coeff. recalc), too large filterbank, arbitrary carrier length,
%         forced tag power level, fs > char.fs
% NOT Tested (too complex): extrapolation value for fs > char.fs
testsettings.selected{1} = [0]; % [t0]; t0=0
testsettings.selected{2} = [stdsettings.tag.modulation.nfft/2]; % [nfft]; smaller filterbank
testsettings.selected{3} = [19.5]; % [LF=fclk/?]; LF not a multiple of fclk
testsettings.selected{4} = [8192, 8, 2]; % [nfft, LF=fs/(nfft*?), fclk=LF*?]; too large filterbank
testsettings.selected{5} = [1/4]; % [add. carrier length ?*nfft]; carrier length not multiple of nfft/2
testsettings.selected{6} = [0.5]; % [0...1 index in pav vector]; forced power level
testsettings.selected{7} = [2.1, 320e3, 2]; %[fs=?*char.fs, LF, fclk=LF*?] fs > fs of char.

% other settings
miscsettings.pin_checks = 1e4; % number of random checks for pin matrix


% *******************************************************************************************************
% check some settings and complete/check partitioning 

if length(randsettings.pwrchars) ~= length(charsettings_pwr.vdda)
   error('Lengths of randsettings.pwrchars and charsettings_pwr.vdda have to match.');
end

if length(partitioning.names) ~= length(partitioning.runs)
   error('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create synthetic tag_power characteristics for tests

for i = 1: length(charsettings_pwr.vdda) 
   % copy settings
   tagchar_pwr.settings.n = charsettings_pwr.n;
   tagchar_pwr.pic = linspace(charsettings_pwr.pic(1), charsettings_pwr.pic(2), tagchar_pwr.settings.n)';
   
   % calculate polynomial coefficients a*x+b; [a, b]
   x1 = charsettings_pwr.pic(1);
   x2 = charsettings_pwr.pic(2);
   y1 = charsettings_pwr.vdda{i}(1);
   y2 = charsettings_pwr.vdda{i}(2);
   charsettings_pwr.poly{i} = [(y2-y1)/(x2-x1); (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   
   % create characteristic
   tagchar_pwr.vdda{i} = charsettings_pwr.poly{i}(1) * tagchar_pwr.pic + charsettings_pwr.poly{i}(2);
end


% *******************************************************************************************************
% create synthetic tag_modulation characteristics for tests

% automatic: create/complete settings and "axes" vectors
for i = 1 : length(charsettings_mod)
   % copy/complete size settings
   tagchar_mod{i}.settings.nf      = charsettings_mod{i}.nf;
   tagchar_mod{i}.settings.nfch    = charsettings_mod{i}.nfch;
   tagchar_mod{i}.settings.np      = charsettings_mod{i}.np;
   tagchar_mod{i}.settings.nr      = charsettings_mod{i}.nm;
   tagchar_mod{i}.settings.nf_out1 = round( (charsettings_mod{i}.nf-charsettings_mod{i}.nfch)/2 );
   tagchar_mod{i}.settings.nf_out4 = charsettings_mod{i}.nf - charsettings_mod{i}.nfch - tagchar_mod{i}.settings.nf_out1;
   
   % create "axes" (log spacing for pav; real char. will have different spacing!)
   tagchar_mod{i}.pav = logspace(log10(charsettings_mod{i}.pav(1)), log10(charsettings_mod{i}.pav(2)), charsettings_mod{i}.np)';
   tagchar_mod{i}.zmod = linspace(charsettings_mod{i}.zmod(1), charsettings_mod{i}.zmod(2), charsettings_mod{i}.nm)';  
   tagchar_mod{i}.fch = linspace(charsettings_mod{i}.fch(1), charsettings_mod{i}.fch(2), charsettings_mod{i}.nfch)';
   tagchar_mod{i}.f = [...
      linspace(charsettings_mod{i}.f(1), charsettings_mod{i}.fch(1), tagchar_mod{i}.settings.nf_out1+1)';... % fout1
      tagchar_mod{i}.fch(2:end-1);... % fch (start/end already in fout1, fout4)
      linspace(charsettings_mod{i}.fch(2), charsettings_mod{i}.f(2), tagchar_mod{i}.settings.nf_out4+1)']; % fout4
    
   % calculate polynomial coefficients a*x+b; [a unmod, b unmod; a mod, b mod]
   charsettings_mod{i}.poly1 = zeros(2,2); % magnitude / real part
   charsettings_mod{i}.poly2 = zeros(2,2); % angle / imaginary part
   %    magnitude / real part, unmodulated
   x1 = charsettings_mod{i}.fch(1);
   x2 = charsettings_mod{i}.fch(2);
   y1 = charsettings_mod{i}.val1(1,1);
   y2 = charsettings_mod{i}.val1(1,2);
   charsettings_mod{i}.poly1(1,:) = [(y2-y1)/(x2-x1); (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   %    magnitude / real part, modulated
   y1 = charsettings_mod{i}.val1(2,1);
   y2 = charsettings_mod{i}.val1(2,2);
   charsettings_mod{i}.poly1(2,:) = [(y2-y1)/(x2-x1); (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   %    phase / imaginary, unmodulated
   x1 = log10(charsettings_mod{i}.pav(1));
   x2 = log10(charsettings_mod{i}.pav(2));
   y1 = charsettings_mod{i}.val2(1,1);
   y2 = charsettings_mod{i}.val2(1,2);
   charsettings_mod{i}.poly2(1,:) = [(y2-y1)/(x2-x1); (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   %    phase / imaginary, modulated
   y1 = charsettings_mod{i}.val2(2,1);
   y2 = charsettings_mod{i}.val2(2,2);
   charsettings_mod{i}.poly2(2,:) = [(y2-y1)/(x2-x1); (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   
   % minimum operational power [W]
   tagchar_mod{i}.picopmin = charsettings_mod{i}.pmin_coeff(1) * tagchar_mod{i}.fch/1e9 + charsettings_mod{i}.pmin_coeff(2);
   tagchar_mod{i}.picopmin = 10.^(tagchar_mod{i}.picopmin/10-3); % dBm -> W
   
   % complete settings as far as necessary 
   % (no need to make these characteristics suitable for tagchar_modulation_plots.m)
   tagchar_mod{i}.settings.fs = charsettings_mod{i}.f(2) * 2 * charsettings_mod{i}.nf/(charsettings_mod{i}.nf-1);
   tagchar_mod{i}.settings.rhomax = charsettings_mod{i}.rhomax;
   tagchar_mod{i}.settings.rat = charsettings_mod{i}.rat(1) + rand*(charsettings_mod{i}.rat(2) - charsettings_mod{i}.rat(1));
   tagchar_mod{i}.settings.cat = charsettings_mod{i}.cat(1) + rand*(charsettings_mod{i}.cat(2) - charsettings_mod{i}.cat(1));
end

% create "nonextrap." characteristic for: synth1 (abs linear w. frequ, angle linear w. power)
%     unmodulated
line_f = charsettings_mod{1}.poly1(1,1) *       tagchar_mod{1}.fch  + charsettings_mod{1}.poly1(1,2);
line_p = charsettings_mod{1}.poly2(1,1) * log10(tagchar_mod{1}.pav) + charsettings_mod{1}.poly2(1,2);
tagchar_mod{1}.rho_pav(:,2,:) = repmat(line_f, 1,charsettings_mod{1}.np) .*...
   exp(complex(0,1) * repmat(line_p', charsettings_mod{1}.nfch,1) );
%     modulated
line_f = charsettings_mod{1}.poly1(2,1) *       tagchar_mod{1}.fch  + charsettings_mod{1}.poly1(2,2);
line_p = charsettings_mod{1}.poly2(2,1) * log10(tagchar_mod{1}.pav) + charsettings_mod{1}.poly2(2,2);
tagchar_mod{1}.rho_pav(:,1,:) = repmat(line_f, 1,charsettings_mod{1}.np) .*...
   exp(complex(0,1) * repmat(line_p', charsettings_mod{1}.nfch,1) );

% create "nonextrap." characteristic for: synth2 (re linear w. frequ, im linear w. power)
%     unmodulated
line_f = charsettings_mod{2}.poly1(1,1) *       tagchar_mod{2}.fch  + charsettings_mod{2}.poly1(1,2);
line_p = charsettings_mod{2}.poly2(1,1) * log10(tagchar_mod{2}.pav) + charsettings_mod{2}.poly2(1,2);
tagchar_mod{2}.rho_pav(:,2,:) = repmat(line_f, 1,charsettings_mod{2}.np) +...
   complex(0,1) * repmat(line_p', charsettings_mod{2}.nfch,1);
%     modulated
line_f = charsettings_mod{2}.poly1(2,1) *       tagchar_mod{2}.fch  + charsettings_mod{2}.poly1(2,2);
line_p = charsettings_mod{2}.poly2(2,1) * log10(tagchar_mod{2}.pav) + charsettings_mod{2}.poly2(2,2);
tagchar_mod{2}.rho_pav(:,1,:) = repmat(line_f, 1,charsettings_mod{2}.np) +...
   complex(0,1) * repmat(line_p', charsettings_mod{2}.nfch,1);

% automatic: create "nonfunctional areas", "extrap." characteristic and chip input impedance
for i = 1 : length(charsettings_mod)
   % check if a non-extrapolated characteristic was created
   if ~isfield(tagchar_mod{i}, 'rho_pav')
      error('rho_pav has not been created for all characteristics. Check settings/program.');
   end
   
   % "extrapolated" characteristic
   %     unmodulated
   tagchar_mod{i}.rhoext_pav(:,1,:) = ones(charsettings_mod{i}.nf, charsettings_mod{i}.np) * charsettings_mod{i}.rhomax;
   tagchar_mod{i}.rhoext_pav(1+tagchar_mod{i}.settings.nf_out1:end-tagchar_mod{i}.settings.nf_out4,1,:) = ...
      tagchar_mod{i}.rho_pav(:,1,:);
   %     modulated
   tagchar_mod{i}.rhoext_pav(:,2,:) = ones(charsettings_mod{i}.nf, charsettings_mod{i}.np) * charsettings_mod{i}.rhomax;
   tagchar_mod{i}.rhoext_pav(1+tagchar_mod{i}.settings.nf_out1:end-tagchar_mod{i}.settings.nf_out4,2,:) = ...
      tagchar_mod{i}.rho_pav(:,2,:);
   
   % areas w. nonfunctional tag (also store indices)
   tagchar_mod{i}.indp_nfunc = round(1+charsettings_mod{i}.nfunc_pav*(tagchar_mod{i}.settings.np-1));
   tagchar_mod{i}.indf_nfunc = round(1+charsettings_mod{i}.nfunc_fch*(tagchar_mod{i}.settings.nfch-1));
   tagchar_mod{i}.rho_pav(tagchar_mod{i}.indf_nfunc(1):tagchar_mod{i}.indf_nfunc(2), :,...
      tagchar_mod{i}.indp_nfunc(1):tagchar_mod{i}.indp_nfunc(2)) =...
      nan(tagchar_mod{i}.indf_nfunc(2)-tagchar_mod{i}.indf_nfunc(1)+1, charsettings_mod{i}.nm, ...
      tagchar_mod{i}.indp_nfunc(2)-tagchar_mod{i}.indp_nfunc(1)+1);
   
   % calculate input power (pin) ranges
   % ... spans all possible pin values pav,rho->pin (does not have to be the case for real char)
   pin_min = tagchar_mod{i}.pav(1)   * (1 - max(max(max(1-abs(tagchar_mod{i}.rho_pav)))).^2);
   pin_max = tagchar_mod{i}.pav(end) * (1 - min(min(min(1-abs(tagchar_mod{i}.rho_pav)))).^2);
   
   % create pic -> pin matrix step 1: factors
   % ... = ( |log10(fch)| *  |log10(|zmod|)| )^(1/4)
   % ... ATTENTION: purely artificial; no physical background or connection to Zic, Zat, etc.
   nfch = tagchar_mod{i}.settings.nfch;
   nr   = tagchar_mod{i}.settings.nr;
   np   = tagchar_mod{i}.settings.np;
   tagchar_mod{i}.pin = (...
      abs(log10(repmat(reshape(    tagchar_mod{i}.fch  , [nfch, 1, 1]), [   1, nr, np]))) .*...
      abs(log10(repmat(reshape(abs(tagchar_mod{i}.zmod), [   1,nr, 1]), [nfch,  1, np]))) ).^(1/4);
   %     check for active assembly
   if any(any(any( tagchar_mod{i}.pin <= 1 )))
      error('pic->pin would require an active assembly impedance. Check settings and this script.');
   end
   
   % calculate chip input power vector (pic)
   % ... ATTENTION: This is purely artificial and DOES NOT have any physical background.
   % ... spans all possible pin values pav,rho->pin (does not have to be the case for real char)
   pic_min = pin_min / max(max(max(tagchar_mod{i}.pin)));
   pic_max = pin_max / min(min(min(tagchar_mod{i}.pin)));
   tagchar_mod{i}.pic = logspace(log10(pic_min), log10(pic_max), charsettings_mod{i}.np)';
   
   % create pic -> pin matrix step 2: multiply by pic (this matrix is pin(fch,rmod,pic)
   tagchar_mod{i}.pin = tagchar_mod{i}.pin .*...
      repmat(reshape(tagchar_mod{i}.pic, [1,1,np]), [nfch,nr,1]);
   %     do a few checks
   for dummy = 1 : miscsettings.pin_checks
      % random indices
      ind_fch  = intrand([1, tagchar_mod{i}.settings.nfch]);
      ind_rmod = intrand([1, tagchar_mod{i}.settings.nr  ]);
      ind_p    = intrand([1, tagchar_mod{i}.settings.np  ]);
      % calculate 
      x1 = tagchar_mod{i}.pin(ind_fch, ind_rmod, ind_p);
      x2 = tagchar_mod{i}.pic(ind_p) *...
         ( abs(log10(tagchar_mod{i}.fch(ind_fch))) * abs(log10(abs(tagchar_mod{i}.zmod(ind_rmod)))) )^(1/4);
      % compare
      if (x1 - x2)^2 > eps
         error('Found an error in tagchar_mod{%d}.pin maxtrix (%d,%d,%d).', i, ind_fch, ind_rmod, ind_pic);
      end
   end
      
   % calculate polynomial coefficients a*x+b for Zic [a, b]
   %   ... make sure |Zic| changes for fch and pic (otherwise errors will not be detected)
   charsettings_mod{i}.poly3 = zeros(1,2); % real part
   charsettings_mod{i}.poly4 = zeros(1,2); % imaginary part
   %    real part
   x1 = charsettings_mod{i}.fch(1);
   x2 = charsettings_mod{i}.fch(2);
   y1 = real(charsettings_mod{i}.val34(1));
   y2 = real(charsettings_mod{i}.val34(2));
   charsettings_mod{i}.poly3 = [(y2-y1)/(x2-x1), (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   %    imaginary part
   x1 = log10(pic_min);
   x2 = log10(pic_max);
   y1 = imag(charsettings_mod{i}.val34(1));
   y2 = imag(charsettings_mod{i}.val34(2));
   charsettings_mod{i}.poly4 = [(y2-y1)/(x2-x1), (y2*(x2-x1)-x2*(y2-y1))/(x2-x1)];
   % calculate chip input impedance and capacity
   line_f = charsettings_mod{i}.poly3(1) *       tagchar_mod{i}.fch  + charsettings_mod{i}.poly3(2);
   line_p = charsettings_mod{i}.poly4(1) * log10(tagchar_mod{i}.pic) + charsettings_mod{i}.poly4(2);
   tagchar_mod{i}.zic = complex(repmat(line_f, 1,charsettings_mod{i}.np), repmat(line_p', charsettings_mod{i}.nfch,1));
   tagchar_mod{i}.cic = -1 ./...
      (2*pi*repmat(tagchar_mod{i}.fch, 1, tagchar_mod{i}.settings.np) .* imag(tagchar_mod{i}.zic));
end

% check characteristics
for i = 1 : length(tagchar_mod)
   % check size of characteristics
   %     "non-extrapolated"
   if size(tagchar_mod{i}.rho_pav, 1) ~= length(tagchar_mod{i}.fch)
      error('Size mismatch in characteristic synth%i (rho_pav, fch)', i);
   end
   if size(tagchar_mod{i}.rho_pav, 2) ~= length(tagchar_mod{i}.zmod)
      error('Size mismatch in characteristic synth%i (rho_pav, r)', i);
   end
   if size(tagchar_mod{i}.rho_pav, 3) ~= length(tagchar_mod{i}.pav)
      error('Size mismatch in characteristic synth%i (rho_pav, p)', i);
   end
   %     "extrapolated"
   if size(tagchar_mod{i}.rhoext_pav, 1) ~= length(tagchar_mod{i}.f)
      error('Size mismatch in characteristic synth%i (rho_pav, f)', i);
   end
   if size(tagchar_mod{i}.rhoext_pav, 2) ~= length(tagchar_mod{i}.zmod)
      error('Size mismatch in characteristic synth%i (rho_pav, r)', i);
   end
   if size(tagchar_mod{i}.rhoext_pav, 3) ~= length(tagchar_mod{i}.pav)
      error('Size mismatch in characteristic synth%i (rho_pav, p)', i);
   end
   % check for "active" reflection
   if any(any(any(abs(tagchar_mod{i}.rho_pav) > 1)))
      error('Characteristic synth%i has |rho| > 1. Check settings.', i);
   end
   if any(any(any(abs(tagchar_mod{i}.rho_pav) < 0)))
      error('Characteristic synth%i has |rho| < 0. Check settings and this script.', i);
   end
   % check for short/open chip
   if any(any(abs(tagchar_mod{i}.zic) < eps))
      error('Characteristic synth%i has |Zic| -> 0. Check settings and this script.', i);
   end
   if any(any(abs(tagchar_mod{i}.zic) > 1/eps))
      error('Characteristic synth%i has |Zic| -> inf. Check settings and this script.', i);
   end
end



% *******************************************************************************************************
% create settings
%   ... sampling rate will be taken from characteristics for first test; changing the sampling rate will
%       also cause the filterbank coefficients to be recalculated for all other blocks

% loop variables
index = 1; % current index in partitioning.indices
clen = zeros(size(partitioning.runs)); % needed carrier length in s (will be coarsly estimated)

for i = 1 : partitioning.indices(end)
           
   % start from standard setup
   settings{i} = stdsettings;
   
   % re-initialize carrier/sampling settings only once per block (speed considerations)  
   if i == partitioning.indices(index) + 1 % this is a new block
      index = index + 1;
      fs = randsettings.fs(1) + rand*(randsettings.fs(2) - randsettings.fs(1)); % sampling frequency
      fc = randsettings.fc(1) + rand*(randsettings.fc(2) - randsettings.fc(1)); % carrier frequency
      settings{i}.newcarrier = true; % trigger carrier (re-)creation
   end
   
   % if this is the first block: set sampling rate to characteristic fs (... another test)
   if i == 1
      fs = tagchar_mod{i}.settings.fs;
   end

   % data to modulate (might be changed to [0] below); random length
   settings{i}.data = [settings{i}.data(1 : intrand(randsettings.data)); 1]; % plus trailing one
   
   % test specific settings
   %     test 1: modulate
   if     i > partitioning.indices(1) && i <= partitioning.indices(2)
      settings{i}.mode = 'modulate';
      settings{i}.tag.modulation.force = false;
   %     test 2: forced modulation
   elseif i > partitioning.indices(2) && i <= partitioning.indices(3)
      settings{i}.mode = 'modulate';
      settings{i}.tag.modulation.force = true;
   %     test 3: filter
   elseif i > partitioning.indices(3) && i <= partitioning.indices(4)
      settings{i}.mode = 'filter';
      settings{i}.tag.modulation.force = true;
      settings{i}.data = [];
   %     test 4: special selected tests (exception handling)
   elseif i > partitioning.indices(4) && i <= partitioning.indices(5)
      settings{i}.mode = 'modulate';
      settings{i}.tag.modulation.force = true;
   else
      error('Index out of range. Check partitioning.');
   end  
   
   % complete setup of reader/tag
   % ATTENTION: SOME SETTINGS ARE CHANGED BELOW; CHECK FOR CONFLICTS WHEN MODIFYING THIS BLOCK!
   %     characteristic (select one of the synthetic characteristics ~uniform
   settings{i}.ind_modcharfile = intrand([1,length(charsettings_mod)]);
   settings{i}.tag.modulation.charfile = charsettings_mod{settings{i}.ind_modcharfile}.filename;
   %     "simple" frequencies
   settings{i}.fs = fs; % will be used frequently
   settings{i}.tag.modulation.fs = fs;
   settings{i}.tag.clock.fs = fs;
   settings{i}.reader.oscillator.fs = fs;
   settings{i}.reader.oscillator.fcenter = fc;
   settings{i}.tag.modulation.fc = fc;
   settings{i}.tag.clock.fcenter = randsettings.fclk(1) + rand*(randsettings.fclk(2) - randsettings.fclk(1));
   settings{i}.tag.modulation.fclk = settings{i}.tag.clock.fcenter;
   %     modulation frequency: make sure fclk is multiple of LF
   settings{i}.tag.modulation.lf = randsettings.lf(1) + rand*(randsettings.lf(2)   - randsettings.lf(1)  );
   settings{i}.tag.modulation.lf = settings{i}.tag.modulation.fclk / round(settings{i}.tag.modulation.fclk / settings{i}.tag.modulation.lf);
   %     leading/trailing unmodulated part; will be corrected by tag_modulation during test
   settings{i}.tag.modulation.t0 = randsettings.t0(1) + rand*(randsettings.t0(2)   - randsettings.t0(1)  );
   %     tag clock length: will be set by tag_modulation
   settings{i}.tag.clock.length = NaN;
   %     available (incident) power level; log spacing
   settings{i}.pav = 10^(log10(randsettings.pav(1)) + rand*(log10(randsettings.pav(2)/randsettings.pav(1))));
   %     power supply characteristic
   settings{i}.ind_pwrchar                = intrand([1, length(randsettings.pwrchars)]);
   settings{i}.tag.modulation.pwrcharfile = charsettings_pwr.filename;
   settings{i}.tag.modulation.pwrchar     = randsettings.pwrchars{settings{i}.ind_pwrchar};
   settings{i}.tag.modulation.pfnc_bias   = rand;
   
   % modify settings for exception handling tests
   if i >= partitioning.indices(end-1)+1 && i <= partitioning.indices(end)
      % make sure tags are operational for these tests by setting a rather high power (only NaN in
      % characteristic is left as possible reason for nonfunctionality => manual check at the end)
      settings{i}.pav = tagchar_mod{settings{i}.ind_modcharfile}.pav(end-1);

      % switch tests
      switch i
         % t0 = 0
         case length(settings)-6
            settings{i}.tag.modulation.t0 = testsettings.selected{1}(1);
            
         % smaller filterbank (recalculation of coefficients)
         %     ... should'nt be a problem for other settings
         case length(settings)-5
            settings{i}.tag.modulation.nfft = testsettings.selected{2}(1);
            
         % LF not a multiple of fclk (will be rounded to fclk)
         case length(settings)-4
            settings{i}.tag.modulation.lf = settings{i}.tag.modulation.fclk / testsettings.selected{3}(1);
            
         % too large filterbank (truncation + likely recalculation of coefficients)
         %     increasing NFFT will lead to large group delay (carrier too short)
         %     => increase NFFT only slightly and change LF
         case length(settings)-3
            %     increase FFT size
            settings{i}.tag.modulation.nfft = testsettings.selected{4}(1);
            settings{i}.tag.modulation.lf = settings{i}.fs / ...
               ( settings{i}.tag.modulation.nfft * testsettings.selected{4}(2)) ;
            %     also increase fclk
            settings{i}.tag.clock.fcenter   = settings{i}.tag.modulation.lf * testsettings.selected{4}(3);
            settings{i}.tag.modulation.fclk = settings{i}.tag.modulation.lf * testsettings.selected{4}(3);

         % carrier length not a multiple of length(fft)/2 (will be truncated)
         case length(settings)-2
            settings{i}.addclen_s = settings{i}.tag.modulation.nfft * testsettings.selected{5}(1);

         % forced power level (might cause nonfunc tag => manual check at the end)
         case length(settings)-1
            settings{i}.tag.modulation.ptag_f = ...
               tagchar_mod{ettings{i}.ind_charfile}.pav(1+ round(testsettings.selected{6}(1)*end-1));
            
         % sampling frequency > fs of characteristic (extrapolation)
         % NOTE: The extrapolation value cannot be checked in a simple way; tag_modulation will prevent
         %       carrier frequencies outside characteristic range
         %       => this checks only if an extrapolation is done
         % ATTENTION: THIS HAS TO BE THE LAST TEST OF A BLOCK, AS IT CHANGES FS
         case length(settings)
            % set sampling frequency above 2*fs of characteristic
            fs = tagchar_mod{settings{i}.ind_modcharfile}.settings.fs * testsettings.selected{7}(1);
            settings{i}.fs = fs;
            settings{i}.tag.modulation.fs = fs;
            settings{i}.tag.clock.fs = fs;
            settings{i}.reader.oscillator.fs = fs;
            % cut data to minimum and increase LF to minimize RAM usage
            % ... might also cause filterbank coefficient recalculation
            % ... tag clock has to be a multiple of LF
            settings{i}.data = settings{i}.data(1:randsettings.data(1));
            settings{i}.tag.modulation.lf = testsettings.selected{7}(2); % (not too high: estimators)
            settings{i}.tag.clock.fcenter   = settings{i}.tag.modulation.lf * testsettings.selected{7}(3);
            settings{i}.tag.modulation.fclk = settings{i}.tag.modulation.lf * testsettings.selected{7}(3);
            % and force carrier recreation
            settings{i}.newcarrier = true;
            
         otherwise % (cannot detect all mismatches)
            error('Number of "exception handling" setup does not match partitioning.');
      end
   end
   
   % additional carrier length [samples] (length will be chosen by tag_modulation.m in case of
   % modulation; this length is used in selected tests to check the truncation mechanism)
   settings{i}.addclen_s = 0; % samples

   % phase/amplitude estimator setup
   %     ... window/overlapping size multiples of carrier period 1/f0
   %     ... window size max t_lf/testsettings.est_nwin (samples)
   settings{i}.other.est_sinusoid.f0 = fc / fs;
   settings{i}.other.est_sinusoid.n  =...
      round(ceil( fc/(testsettings.est_nwin*settings{i}.tag.modulation.lf) ) / settings{i}.other.est_sinusoid.f0 );
   settings{i}.other.est_sinusoid.ol =...
      round(round( settings{i}.other.est_sinusoid.n * settings{i}.other.est_sinusoid.ol/100 *...
      settings{i}.other.est_sinusoid.f0 ) / settings{i}.other.est_sinusoid.f0) / settings{i}.other.est_sinusoid.n * 100;
   
   % needed carrier length: (t0 + data + group delay plus margin rounded to multiples of filterbank length)
   % ... duration of one entry in data is 1/2 subcarrier period here (data creates the subcarrier)
   clen_s = ceil( fs * (1+testsettings.lmargin) * (2*settings{i}.tag.modulation.t0 ...
         + (1+testsettings.lmargin) * (length(settings{i}.data) / settings{i}.tag.modulation.lf / 2 ... 
         + settings{i}.tag.modulation.nfft/2/fs * ( settings{i}.tag.modulation.nfilt - 1))) );
   clen_s = clen_s + settings{i}.tag.modulation.nfft/2 - mod(clen_s, settings{i}.tag.modulation.nfft/2);
   % record maximum needed length for this block
   if clen(index-1) < clen_s / fs;
      clen(index-1) = clen_s / fs;
   end
end

% set missing lengths to maximum recorded carrier length (will be verified later)
index = 1; % current index in partitioning.indices
for i = 1 : length(settings)
   % new block?
   if i == partitioning.indices(index) + 1 % this is a new block
      index = index + 1;
   end
   % needed carrier length (again + margin to compensate for uncertainties of clock(end) (see above)
   settings{i}.reader.oscillator.length = clen(index-1) * (1+testsettings.lmargin); % s
end


% *******************************************************************************************************
% calculate results

for i = 1 : length(settings)
   
   % incident power level used by tag_modulation
   if isnan(settings{i}.tag.modulation.ptag_f)
      results{i}.pav = settings{i}.pav;
   else
      results{i}.pav = settings{i}.tag.modulation.ptag_f;
   end
      
   % group delay of filterbank in samples @ fs
   results{i}.grpdel_s = settings{i}.tag.modulation.nfft/2 * ( settings{i}.tag.modulation.nfilt - 1);
   
   % settings modified by tag_modulation (no other possibility on tag)
   %     real LF will be rounded to fclk  
   results{i}.lf = settings{i}.tag.modulation.fclk / round(settings{i}.tag.modulation.fclk / settings{i}.tag.modulation.lf);
   %     t0 will also be rounded to fclk (implicitly) and to filterbank stepsize and saturated to make
   %     sure the modulated part will not be truncated
   results{i}.t0_s = round(settings{i}.fs * ceil(settings{i}.tag.modulation.fclk *...
                           settings{i}.tag.modulation.t0) / settings{i}.tag.modulation.fclk);
   results{i}.t0_s = round(settings{i}.fs*settings{i}.tag.modulation.t0);
   results{i}.t0_s = max(results{i}.t0_s,...
      2*round(settings{i}.fs/settings{i}.tag.modulation.fclk)+settings{i}.tag.modulation.nfft/2);
   %     t0 at end will be calculated by test_tag_modulation to avoid conflicts with LF
   
   % expected reflection coefficient [unmod, mod]
   switch settings{i}.ind_modcharfile
      case 1 % abs/angle
         results{i}.gain(1) = (charsettings_mod{1}.poly1(1,1)*settings{i}.tag.modulation.fc+charsettings_mod{1}.poly1(1,2)) *... % abs
                 exp(complex(0,charsettings_mod{1}.poly2(1,1)*log10(results{i}.pav)+charsettings_mod{1}.poly2(1,2))); % angle
         results{i}.gain(2) = (charsettings_mod{1}.poly1(2,1)*settings{i}.tag.modulation.fc+charsettings_mod{1}.poly1(2,2)) *... % abs
                 exp(complex(0,charsettings_mod{1}.poly2(2,1)*log10(results{i}.pav)+charsettings_mod{1}.poly2(2,2))); % angle
      case 2 % re/im
         results{i}.gain(1) = complex(...
            charsettings_mod{2}.poly1(1,1)*settings{i}.tag.modulation.fc+charsettings_mod{2}.poly1(1,2),... % re
            charsettings_mod{2}.poly2(1,1)*log10(results{i}.pav)+charsettings_mod{2}.poly2(1,2)); % im
         results{i}.gain(2) = complex(...
            charsettings_mod{2}.poly1(2,1)*settings{i}.tag.modulation.fc+charsettings_mod{2}.poly1(2,2),... % re
            charsettings_mod{2}.poly2(2,1)*log10(results{i}.pav)+charsettings_mod{2}.poly2(2,2)); % im
      otherwise
         error('Characteristic index out of range.');
   end
                    
   % is the tag functional ?
   %     shortcut
   ind_c = settings{i}.ind_modcharfile;
   %     input power and chip input power [unmod, mod]
   %     ... see creation of synthetic characteristic above (tagchar_mod{}.pin) for pic -> pin factor
   results{i}.pin_01 = results{i}.pav * (1 - abs(results{i}.gain).^2);
   results{i}.pic_01 = results{i}.pin_01 ./ ( abs(log10(settings{i}.tag.modulation.fc)) *...
      abs(log10( abs(charsettings_mod{settings{i}.ind_modcharfile}.zmod(end:-1:1)) )) ).^(1/4);
   %     average (this is not as accurate as the averaging done by tag_modulation!)
   if ~isempty(settings{i}.data)
      results{i}.pin = results{i}.pin_01(2) * mean(settings{i}.data) + results{i}.pin_01(1) * (1-mean(settings{i}.data));
      results{i}.pic = results{i}.pic_01(2) * mean(settings{i}.data) + results{i}.pic_01(1) * (1-mean(settings{i}.data));
      % bias towards unmodulated (simple supply buffer model)
      results{i}.pin = results{i}.pin_01(1) * settings{i}.tag.modulation.pfnc_bias + results{i}.pin * (1-settings{i}.tag.modulation.pfnc_bias);
      results{i}.pic = results{i}.pic_01(1) * settings{i}.tag.modulation.pfnc_bias + results{i}.pic * (1-settings{i}.tag.modulation.pfnc_bias);      
   else
      results{i}.pin = results{i}.pin_01(1);
      results{i}.pic = results{i}.pic_01(1);
   end
   %     position in characteristic ('extrap' for index 1 and end)
   pos_fch = interp1(tagchar_mod{ind_c}.fch,...
      [1:charsettings_mod{ind_c}.nfch], settings{i}.tag.modulation.fc, 'linear', 'extrap');
   pos_pav = interp1(tagchar_mod{ind_c}.pav,...
      [1:charsettings_mod{ind_c}.np], results{i}.pav, 'linear', 'extrap');
   if pos_pav == charsettings_mod{ind_c}.np; % this is similar to what tag_modulation does
      pos_pav = pos_pav - 1; 
   end
   %     check
   %        out of pav range: can't be functional
   results{i}.func.pav = results{i}.pav > tagchar_mod{ind_c}.pav(1) && results{i}.pav < tagchar_mod{ind_c}.pav(end);
   %        chip input power > threshold? (levels in dBm here)
   results{i}.func.pic = ( charsettings_mod{ind_c}.pmin_coeff(1) * settings{i}.tag.modulation.fc/1e9 +...
      charsettings_mod{ind_c}.pmin_coeff(2) ) <= ( 10*log10(results{i}.pic)+30 );
   %        only nonextrapolated characteristic contains NaNs to indicate a nonfunctional tag
   results{i}.func.nan = results{i}.func.pav && ... % in pav range
         pos_fch >= 1 && pos_fch <= tagchar_mod{ind_c}.settings.nfch && ... % in fch range
         ~any(any(any(isnan(tagchar_mod{ind_c}.rho_pav(floor(pos_fch):floor(pos_fch)+1, :,... % no NaNs...
         floor(pos_pav):floor(pos_pav)+1))))); % ...in vicinity
   %        combine  
   results{i}.func.all = results{i}.func.pav && results{i}.func.pic && results{i}.func.nan;
   
   % expected power supply voltage
   results{i}.vdda = charsettings_pwr.poly{settings{i}.ind_pwrchar}(2) +...
      charsettings_pwr.poly{settings{i}.ind_pwrchar}(1) * (10*log10(results{i}.pic)+30);
   %     saturate to minimum/maximum (tagchar_pwr.vdda are linear curves => min/max are borders)
   if results{i}.vdda > max(tagchar_pwr.vdda{settings{i}.ind_pwrchar})
      results{i}.vdda = max(tagchar_pwr.vdda{settings{i}.ind_pwrchar});
   elseif results{i}.vdda < min(tagchar_pwr.vdda{settings{i}.ind_pwrchar})
      results{i}.vdda = min(tagchar_pwr.vdda{settings{i}.ind_pwrchar});
   end
   
   % expected chip input impedance (re/im)
   %     no support extrapolated support || no (nonext) support or nonfunctional tag and not forced
   %     => gain will be set to zero anyway
   if ~results{i}.func.pav || (~results{i}.func.all && ~settings{i}.tag.modulation.force)
      results{i}.zic        = Inf; % chip input impedance
      results{i}.zin_filter = Inf; % input impedance (including assembly and modulation)
   %     normal calculation (tag_modulation is able to estimate pic)
   else
      results{i}.zic = complex(...
         charsettings_mod{ind_c}.poly3(1)*settings{i}.tag.modulation.fc + charsettings_mod{ind_c}.poly3(2),... % re
         charsettings_mod{ind_c}.poly4(1)*log10(results{i}.pic)         + charsettings_mod{ind_c}.poly4(2)); % im
      zat = complex(tagchar_mod{ind_c}.settings.rat, -1/(2*pi*settings{i}.tag.modulation.fc*tagchar_mod{ind_c}.settings.cat));
      results{i}.zin_filter = 1 / (1/results{i}.zic + 1/zat + 1/charsettings_mod{ind_c}.zmod(end)); % ONLY IN MODE FILTER
   end
   
   % modify reflection coefficient in filter/passive mode
   %     copy data vector, modify below in case no data is modulated
   results{i}.data = settings{i}.data;
   %     in filter mode (transmission coefficient, not reflection coefficient; no modulation, no phase)
   %     ... multiplication by abs(results{i}.zin) handled by test_tag_modulation
   if strcmpi(settings{i}.mode, 'filter') % gain(1) is unmod
      results{i}.gain(1) = sqrt(1 - abs(results{i}.gain(1)).^2); 
      results{i}.gain(2) = results{i}.gain(1); % just to be on the safe side
   end
   %     unmodulated reflection (not functional and not forced)
   if strcmpi(settings{i}.mode, 'modulate') && ~results{i}.func.all && ~settings{i}.tag.modulation.force
      results{i}.gain(2) = results{i}.gain(1);
      results{i}.data = results{i}.data * 0;
   end
   %     zero vector if outside char. support or filter & not functional
   if ~results{i}.func.pav ||...
      strcmpi(settings{i}.mode, 'filter') && ~results{i}.func.all && ~settings{i}.tag.modulation.force
      results{i}.gain(1) = 0;
      results{i}.gain(2) = 0;
      results{i}.data = results{i}.data * 0; % if filter: []*0=[]
   end
   
   % clock length (= length of modulated part [samples @ fclk]); this is likely an upper bound
   results{i}.clk_len = round(length(results{i}.data)/results{i}.lf/2 * settings{i}.tag.clock.fcenter) + 1;
end


% *******************************************************************************************************
% check if tag_modulation thinks it is able to deal with these settings

for i = 1 : length(settings)
   % let tag_modulation complete the settings
   settings_completed = tag_modulation(settings{i}.data, settings{i}.tag.modulation);
   % needed carrier length (can't use critwarn, because logging is switched off)
   if settings_completed.length_s > round(settings{i}.reader.oscillator.length * settings{i}.fs)
      disp('WARNING: tag_modulation reports a longer needed carrier length. Tests will likely fail. Check tag_modulation.m and testsettings.lmargin.');
   end
end


% *******************************************************************************************************
% output some information

% prepare info struct fields
info.nans = zeros(size(tagchar_mod)); % ratio of NaNs in "nonextrapolated" characteristic
%     ratio of functional tags
info.func.pav = zeros(size(partitioning.runs)); % reason: pav outside char. range
info.func.pic = zeros(size(partitioning.runs)); % reason: not enough chip power
info.func.nan = zeros(size(partitioning.runs)); % reason: "nonextrap." char. is NaN
info.func.all = zeros(size(partitioning.runs)); % tag functional?

% ratio of NaNs in "nonextrapolated" characteristic
for i = 1 : length(tagchar_mod)
   info.nans(i) = sum(sum(sum(isnan(tagchar_mod{i}.rho_pav)))) / numel(tagchar_mod{i}.rho_pav);
end

% ratio of functional tags
for i = 1 : length(results)
   % find index for partition struct (partitioning.indices is "end-of-block")
   for j = 1 : length(partitioning.indices)
      if i <= partitioning.indices(j)
         j = j - 1; %#ok<FXSET> "end-of-block" => in-block
         break;
      end
   end
   % tag functional (makes only sense in mode "modulate") ?
   info.func.pav(j) = info.func.pav(j) + results{i}.func.pav;
   info.func.pic(j) = info.func.pic(j) + results{i}.func.pic; 
   info.func.nan(j) = info.func.nan(j) + results{i}.func.nan; 
   info.func.all(j) = info.func.all(j) + results{i}.func.all; 
end
%     average
info.func.pav = info.func.pav ./ partitioning.runs;
info.func.pic = info.func.pic ./ partitioning.runs;
info.func.nan = info.func.nan ./ partitioning.runs;
info.func.all = info.func.all ./ partitioning.runs;

% output ("quick and dirty")
disp(sprintf('\nStatistics:\n***********'));
disp(sprintf('Ratio of NaNs in nonextrap. char 1 through %i.: %.1f%%   %.1f%%', length(tagchar_mod), info.nans*100));
disp('Ratio of nonfunctional tags for reasons:   Pav |   Pic |   NaN |   SUM');
for i = 1 : length(partitioning.runs)
   disp(sprintf('                                         %5.1f | %5.1f | %5.1f | %5.1f [%%]',...
      (1-info.func.pav(i))*100, (1-info.func.pic(i))*100, (1-info.func.nan(i))*100, (1-info.func.all(i))*100));
end
disp(' ');
disp('Please check if these results (statistics) are satisfying:');
disp('   - all reasons for nonfunctional tags have to be present for all blocks except the last');
disp('   - tag should always be functional for the last block');
disp(' ');


% *******************************************************************************************************
% save results

% is there already a .mat-file with this name?
search = dir(strcat(fullfile(directory, filename),'*'));
if ~isempty( search )
   % ask before overwriting
   reply = input('This will overwrite an existing file. Are you REALLY sure y/n [n]?  ', 's');
   if ~strcmpi(reply, 'y')
      return
   end
end

% prepare information header
matfilename    = filename; %#ok<NASGU>
characteristic = 'settings and results for selftest: tag_modulation.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning', 'info');
disp(sprintf('\nSettings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));


% *******************************************************************************************************
% also save synthetic characteristics

% tagchar_power
%     add normal information headers
matfilename    = charsettings_pwr.filename; %#ok<NASGU>
characteristic = 'synthetic (partial) pwr characteristics for test_tag_modulation'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
% extract data
pic = tagchar_pwr.pic;
for i = 1 : length(tagchar_pwr.vdda)
   eval(sprintf('%s = tagchar_pwr.vdda{%d};', randsettings.pwrchars{i}, i));
end
settings = tagchar_pwr.settings;
%     save
save(fullfile(directory, matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'pic',randsettings.pwrchars{:}, 'settings');
disp(sprintf('Synthetic characteristic for tag_power has been saved in %s/%s.mat',...
   pwd, fullfile(directory, matfilename) ));

% tagchar_modulator
for i = 1 : length(tagchar_mod)
   % add normal information headers
   matfilename    = charsettings_mod{i}.filename; %#ok<NASGU>
   characteristic = 'synthetic (partial) mod characteristics for test_tag_modulation'; %#ok<NASGU>
   createdon      = datestr(now, 0); %#ok<NASGU>
   createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
   usedmatfiles   = ''; %#ok<NASGU>
   % extract data
   f          = tagchar_mod{i}.f;
   fch        = tagchar_mod{i}.fch;
   zmod       = tagchar_mod{i}.zmod;
   pav        = tagchar_mod{i}.pav;
   pin        = tagchar_mod{i}.pin;
   rho_pav    = tagchar_mod{i}.rho_pav;
   rhoext_pav = tagchar_mod{i}.rhoext_pav;
   pic        = tagchar_mod{i}.pic;
   zic        = tagchar_mod{i}.zic;
   cic        = tagchar_mod{i}.cic;
   picopmin   = tagchar_mod{i}.picopmin;
   settings   = tagchar_mod{i}.settings;
   % save
   save(fullfile(directory, matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
         'f','fch','zmod','pav','pin', 'rho_pav','rhoext_pav', 'pic','zic','cic', 'picopmin', 'settings');
   disp(sprintf('Synthetic characteristic %i for tag_modulation has been saved in %s/%s.mat',...
      i, pwd, fullfile(directory, matfilename) ));
end


return
% *******************************************************************************************************
% DEBUG

% close all;
% 
% colors = get(0, 'DefaultAxesColorOrder');
% figure; hold on;
% for i = 1: length(charsettings_pwr.vdda) 
%    plot(tagchar_pwr.pic, tagchar_pwr.vdda{i}, 'Color', colors(i,:));
% end
% hold off; grid on; setlegend(randsettings.pwrchars, 'NorthWest');
% title('pwr\_char'); xlabel('Pic [dBm]'); ylabel('Vdda [V]');
% 
% for i = 1 : length(charsettings_mod)
%    figure; hold on;
%    surface(tagchar_mod{i}.fch, tagchar_mod{i}.pav, abs(squeeze(tagchar_mod{i}.rho_pav(:,1,:)))')
%    surface(tagchar_mod{i}.fch, tagchar_mod{i}.pav, abs(squeeze(tagchar_mod{i}.rho_pav(:,2,:)))', 'EdgeColor', 'none')
%    hold off; grid on; set(gca, 'yScale','log');
%    title(sprintf('mod\\_char%d: abs',i)); xlabel('f [Hz]'); ylabel('Pav [W]'); legend('mod', 'umd'); colorbar;
%    
%    figure; hold on;
%    surface(tagchar_mod{i}.fch, tagchar_mod{i}.pav, angle(squeeze(tagchar_mod{i}.rho_pav(:,1,:)))')
%    surface(tagchar_mod{i}.fch, tagchar_mod{i}.pav, angle(squeeze(tagchar_mod{i}.rho_pav(:,2,:)))', 'EdgeColor', 'none')
%    hold off; grid on; set(gca, 'yScale','log');
%    title(sprintf('mod\\_char%d: arg',i)); xlabel('f [Hz]'); ylabel('Pav [W]'); legend('mod', 'umd'); colorbar;
%    
%    figure; hold on;
%    surface(tagchar_mod{i}.fch, tagchar_mod{i}.pic, abs(squeeze(tagchar_mod{i}.zic))')
%    hold off; grid on; set(gca, 'yScale','log'); title(sprintf('mod\\_char%d: |Zic|',i)); xlabel('f [Hz]'); ylabel('Pic [W]'); colorbar;
%    
%    figure; hold on;
%    surface(tagchar_mod{i}.fch, tagchar_mod{i}.pic, angle(squeeze(tagchar_mod{i}.zic))')
%    hold off; grid on; set(gca, 'yScale','log'); title(sprintf('mod\\_char%d: arg(Zic)',i)); xlabel('f [Hz]'); ylabel('Pic [W]'); colorbar;
% end

%    figure; hold on;
%    surface(tagchar_mod{i}.f, tagchar_mod{i}.pav, abs(squeeze(tagchar_mod{i}.rhoext_pav(:,1,:)))')
%    surface(tagchar_mod{i}.f, tagchar_mod{i}.pav, abs(squeeze(tagchar_mod{i}.rhoext_pav(:,2,:)))', 'EdgeColor', 'none')
%    title(sprintf('mod_char%d',i)); xlabel('f [Hz]'); ylabel('Pav [W]'); legend('mod', 'umd'); colorbar;
%    hold off; grid on; legend('mod', 'umd'); colorbar;
   
% figure; hold on;
% surface(tagchar_mod{2}.fch, tagchar_mod{2}.pav, abs(squeeze(tagchar_mod{2}.rho_pav(:,1,:)))')
% surface(tagchar_mod{2}.fch, tagchar_mod{2}.pav, abs(squeeze(tagchar_mod{2}.rho_pav(:,2,:)))', 'EdgeColor', 'none')
% hold off; grid on; set(gca, 'yScale','log'); title('char2: abs'); xlabel('f [Hz]'); ylabel('Pav [W]'); legend('mod', 'umd'); colorbar;
% figure; hold on;
% surface(tagchar_mod{2}.fch, tagchar_mod{2}.pav, angle(squeeze(tagchar_mod{2}.rho_pav(:,1,:)))')
% surface(tagchar_mod{2}.fch, tagchar_mod{2}.pav, angle(squeeze(tagchar_mod{2}.rho_pav(:,2,:)))', 'EdgeColor', 'none')
% hold off; grid on; set(gca, 'yScale','log'); title('char2: arg'); xlabel('f [Hz]'); ylabel('Pav [W]'); legend('mod', 'umd'); colorbar;

% close all;

% figure; hold on;
% surface(tagchar_mod{1}.fch, tagchar_mod{1}.pic, abs(squeeze(tagchar_mod{1}.zic))')
% hold off; grid on; set(gca, 'yScale','log'); title('char1: |Zic|'); xlabel('f [Hz]'); ylabel('Pic [W]'); colorbar;
% figure; hold on;
% surface(tagchar_mod{1}.fch, tagchar_mod{1}.pic, angle(squeeze(tagchar_mod{1}.zic))')
% hold off; grid on; set(gca, 'yScale','log'); title('char1: arg(Zic)'); xlabel('f [Hz]'); ylabel('Pic [W]'); colorbar;
% 
% figure; hold on;
% surface(tagchar_mod{2}.fch, tagchar_mod{2}.pic, abs(squeeze(tagchar_mod{2}.zic))')
% hold off; grid on; set(gca, 'yScale','log'); title('char2: |Zic|'); xlabel('f [Hz]'); ylabel('Pic [W]'); colorbar;
% figure; hold on;
% surface(tagchar_mod{2}.fch, tagchar_mod{2}.pic, angle(squeeze(tagchar_mod{2}.zic))')
% hold off; grid on; set(gca, 'yScale','log'); title('char2: arg(Zic)'); xlabel('f [Hz]'); ylabel('Pic [W]'); colorbar;

% i = 31;
% ind_p = interp1(tagchar_mod{settings{i}.ind_modcharfile}.pic, [1:1:length(tagchar_mod{settings{i}.ind_modcharfile}.pic)], results{i}.pic, 'nearest');
% ind_f = interp1(tagchar_mod{settings{i}.ind_modcharfile}.fch, [1:1:length(tagchar_mod{settings{i}.ind_modcharfile}.fch)], settings{i}.tag.modulation.fc, 'nearest');
% tagchar_mod{settings{i}.ind_modcharfile}.zic(ind_f,ind_p)
% results{i}.zic

for i = 1 : length(settings)
   test.pic(i) = results{i}.pic;
   test.thr(i) = charsettings_mod{ind_c}.pmin_coeff(1) * settings{i}.tag.modulation.fc/1e9 +...
      charsettings_mod{ind_c}.pmin_coeff(2)
   ind_c = settings{i}.ind_modcharfile;
end

figure; plot(10*log10(test.pic))



% *******************************************************************************************************
% old code snippets

% expected chip input impedance (re/im)
%     no nonext support, ext support (otherwise forcing not possible) and forced
%     => tag_modulation uses pic=tagchar_mod.pic(round(pos_pav)) to ensure operability range)
% if ~results{i}.func.nan && results{i}.func.pav && settings{i}.tag.modulation.force
%    results{i}.zic = complex(...
%       charsettings_mod{ind_c}.poly3(1)*settings{i}.tag.modulation.fc                 + charsettings_mod{ind_c}.poly3(2),... % re
%       charsettings_mod{ind_c}.poly4(1)*log10(tagchar_mod{ind_c}.pic(round(pos_pav))) + charsettings_mod{ind_c}.poly4(2)); % im

%    if strcmpi(settings{i}.mode, 'modulate')
%       results{i}.pic = results{i}.pav * (1 - abs(results{i}.gain(1))^2); % gain is reflected power @ avg = unmod
%    elseif strcmpi(settings{i}.mode, 'filter')
%        results{i}.pic = results{i}.pav * abs(results{i}.gain(1))^2; % gain is transmitted power @ unmod
%    else
%       error('Unsupported settings{i}.mode=''%s''.', settings{i}.mode);
%    end    
