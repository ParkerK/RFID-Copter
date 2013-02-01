% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% multi-frequency continuous-wave radar ranging (MFCW)
%    concept THEORY simulator (using only theory values)
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
% ***** Assumptions *****
% = based on linear channel model with linearized tag, AWGN channel (no fading so far)
% - partially compensated systematic errors
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


function [dist_est, exp_var, other, settings] = mfcw_theorysim(settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   dist_est = version;
   return
end


% *******************************************************************************************************
% input parameter check

if nargin ~= 1
   criterr('Wrong number of input arguments.');
end

% % check contents of settings
% %     prepare required data
% expected.name = 'settings';
% expected.reqfields = {'f_i', 'c', 'A_i', 'f_c', 'N_0', 'f_s', 'f_rs', 'h_i', 'n'}; 
% %     check
% errortext = contentcheck(settings, expected);
% %     output
% if ~isempty(errortext)
%    err('Incomplete settings\n%s', errortext);
% end

% check number of secondary carriers
if settings.carriers.N_c < 2
   err('At least two carriers needed. Check settings.N_c.');
end
if length(settings.carriers.f_i) ~= length(settings.carriers.A_i) ||...
   settings.carriers.N_c > length(settings.carriers.f_i)
   err('Lengths in carrier setup do not match.');
end

% check noise correlation coefficient
if abs(settings.awgn.corr_iq) > 1 || abs(settings.awgn.corr_ci) > 1
   err('Correlation coefficients has to be |r|<=1');
end
   

% *******************************************************************************************************
% initialization, random parameters
if ~settings.silent; disp('   Initializing'); end

% ensure proper vector orientation
settings.carriers.f_i = settings.carriers.f_i(:)';
settings.carriers.A_i = settings.carriers.A_i(:)';

% cut length of carrier vectors to settings.carriers.N_c
settings.carriers.f_i = settings.carriers.f_i(1:settings.carriers.N_c);
settings.carriers.A_i = settings.carriers.A_i(1:settings.carriers.N_c);

% set random parameters [mean, std]
%     phase shift of carriers (~normal per carrier) not including demod. mixer delay
settings.carriers.arg_A_i = settings.random.arg_A_i(:,1) + randn(1,settings.carriers.N_c).*settings.random.arg_A_i(:,2);
%     delay of IQ mixer signal compared to sent carrier (~normal; kept positive)
settings.carriers.ms_del = abs( settings.random.ms_del(:,1) + randn * settings.random.ms_del(:,2) );
%     phase shift of tag subcarrier (~uniform)
settings.tag.phi_m = settings.random.phi_m(:,1) + (2*rand-1)*sqrt(3)*settings.random.phi_m(:,2);
%     delay of feedback term (~normal; kept positive)
settings.feedback.delay = abs( settings.random.fb_del(:,1) + randn * settings.random.fb_del(:,2) );

% calculate large-scale delay of channel
settings.channel.delay = settings.dist / settings.c;

% round delays to samples (match mfcw_sim)
% settings.feedback.delay = round(settings.feedback.delay * settings.f_s) / settings.f_s;
% settings.channel.delay  = round(settings.channel.delay  * settings.f_s) / settings.f_s;
%     output
% disp(sprintf('      Rounding delays to samples @ f_s; resolution is %g ns (%g cm).', 1/settings.f_s*1e9, settings.c/settings.f_s*1e2));

% number of samples (estimation window)
settings.N_e = round(settings.t_max * settings.f_rs);


% *******************************************************************************************************
% calculate (complex) gains
if ~settings.silent; disp('   Calculating (complex) gains'); end

% carrier amplitudes (plus random phase shift)
A_i = settings.carriers.A_i .*...
   exp(complex(0, settings.carriers.arg_A_i + 2*pi*settings.carriers.ms_del*settings.carriers.f_c));

% feedback (gain + delay)
K_i = settings.feedback.gain * exp(complex(0, settings.feedback.delay*2*pi*(settings.carriers.f_c+settings.carriers.f_i)));

% channel
%     largescale
H_Mi = settings.c ./ (4*pi*settings.dist*(settings.carriers.f_c+settings.carriers.f_i-settings.tag.f_m)) .*...
   exp(complex(0,-2*pi*(settings.carriers.f_c+settings.carriers.f_i-settings.tag.f_m)*settings.channel.delay)); % H_{n-M}
H_i  = settings.c ./ (4*pi*settings.dist*(settings.carriers.f_c+settings.carriers.f_i                 )) .*...
   exp(complex(0,-2*pi*(settings.carriers.f_c+settings.carriers.f_i                 )*settings.channel.delay)); % H_{n}
H_iM = settings.c ./ (4*pi*settings.dist*(settings.carriers.f_c+settings.carriers.f_i+settings.tag.f_m)) .*...
   exp(complex(0,-2*pi*(settings.carriers.f_c+settings.carriers.f_i+settings.tag.f_m)*settings.channel.delay)); % H_{n+M}

% reader input stage (has to be identical to mfcw_sim)
%     generate baseband filter
[b, a] = cheby2(settings.filter.ord_iq, settings.filter.att_iq, (settings.carriers.f_c-max(settings.carriers.f_i)-settings.tag.f_m)*2/settings.f_s);
%    	calculate complex gains for baseband
G_Mi = freqz(b,a, [settings.carriers.f_i-settings.tag.f_m], settings.f_s); % G_{n-M}
G_i  = freqz(b,a, [settings.carriers.f_i                 ], settings.f_s); % G_{n}
G_iM = freqz(b,a, [settings.carriers.f_i+settings.tag.f_m], settings.f_s); % G_{n+M}

% variable cleanup
clear('b','a');


% *******************************************************************************************************
% obtain reflection coefficient model data
%   ... power calculations are kept as simple as possible, but may be inaccurate (rounded indices, etc)
if ~settings.silent; disp('   Calculating linear reflection coefficient model'); end

% load data
tagchar_mod = loadmat(settings.tag.tagchar_mod);

% calculate reflection coefficient for necessary frequency range and minimum power for support
%     get indices in frequency vector
ind_fmin = floor(interp1(tagchar_mod.f, [1:1:tagchar_mod.settings.nf],...
   min(settings.carriers.f_c, settings.carriers.f_c+min(settings.carriers.f_i))-settings.tag.f_m));
ind_fmax =  ceil(interp1(tagchar_mod.f, [1:1:tagchar_mod.settings.nf],...
   max(settings.carriers.f_c, settings.carriers.f_c+max(settings.carriers.f_i))+settings.tag.f_m));
if isnan(ind_fmin) || isnan(ind_fmax)
   err('Minimum/Maximum frequency outside range %.1f-%.1f MHz.', tagchar_mod.f(1), tagchar_mod.f(2));
end
%     get center value for that frequency range
rho_c = ( squeeze(tagchar_mod.rhoext_pav(ind_fmin:ind_fmax, 1, :)) + squeeze(tagchar_mod.rhoext_pav(ind_fmin:ind_fmax, end, :)) ) / 2;
%     sum over frequency => rho_c(pav)
rho_c = sum(rho_c, 1);
%     find index of first non-NaN value => minimum operational power index
ind = findzeros([1, isnan(rho_c), 1] - 0.1); % -0.1 => bias towards non-NaN
ind_Pmin = ind(1)-1;
ind_Pmax = ind(end)-1;

% calculate power level
if isnan(settings.tag.p_av)
   % ... var(cos) = 1/2; var(a*cos)=a^2/2; var(a+b)=var(a)+var(b)+covar; cos diff frequ: covar=0
   settings.tag.p_av = sum( abs(settings.carriers.A_i).^2 .* abs(H_i).^2 ) / 2;
else
   if ~settings.silent; disp('      Overriding power level!'); end
end
%     output
if settings.verbose && ~settings.silent
   disp(sprintf('      Available power is %.2e W (%.2f dBm).', settings.tag.p_av, 10*log10(settings.tag.p_av)+30));
end

% check power level
settings.tag.p_trunc = false;
if settings.tag.p_av < tagchar_mod.pav(ind_Pmin) || settings.tag.p_av > tagchar_mod.pav(ind_Pmax)
   settings.tag.p_trunc = true;
   if ~settings.silent; disp('      Power level outside operational range. Using last available level.'); end
   if settings.tag.p_av < tagchar_mod.pav(ind_Pmin)
      settings.tag.p_av = tagchar_mod.pav(ind_Pmin) * (1+eps); % make sure interp2 below does not produce NaNs 
   else
      settings.tag.p_av = tagchar_mod.pav(ind_Pmax) * (1-eps);
   end
end
%     output
if settings.verbose && ~settings.silent
   disp(sprintf('      Next characteristic power level: %.2e W (%.2f dBm).', settings.tag.p_av, 10*log10(settings.tag.p_av)+30));
end

% calculate linear reflection coefficient model
[x , y ] = meshgrid(tagchar_mod.pav, tagchar_mod.f);
[xi, yi] = meshgrid(settings.tag.p_av+eps, tagchar_mod.f);
hh0 = interp2(x,y, squeeze(tagchar_mod.rhoext_pav(:, end, :)), xi, yi, 'linear', 1); % interp2c leads to steps in linear model
hh1 = interp2(x,y, squeeze(tagchar_mod.rhoext_pav(:,   1, :)), xi, yi, 'linear', 1);
 rho = ( hh1 + hh0 ) / 2;
drho = rho - hh0;

% figure; hold on;
% plot(tagchar_mod.f, angle(hh0), 'b--');
% plot(tagchar_mod.f, angle(hh1), 'r--');
% plot(tagchar_mod.f, angle(rho), 'b');
% plot(tagchar_mod.f, angle(drho), 'r');
% hold off; grid on; xlim([700,1200]*1e6);

% extract relevant data
%     warnings off for NaN values
warning('off', 'MATLAB:interp1:NaNinY')
warning('off', 'MATLAB:chckxy:IgnoreNaN');
%     linear model at carrier frequencies
 rho_i = interp1c(tagchar_mod.f,  rho, settings.carriers.f_c+settings.carriers.f_i, 'spline');
drho_i = interp1c(tagchar_mod.f, drho, settings.carriers.f_c+settings.carriers.f_i, 'spline');
%     warnings on again
warning('on', 'MATLAB:interp1:NaNinY')
warning('on', 'MATLAB:chckxy:IgnoreNaN');

% double-check for NaNs
if any(isnan(rho_i)) || any(isnan(drho_i))
   error('NaNs found in linear model.');
end

% calulate average chip power level assuming square modulation
%    input power for unmod, mod
if settings.tag.p_trunc
   % use weighted p_av from above
   p_in0 = settings.tag.p_av * abs(settings.carriers.A_i).^2 / sum(abs(settings.carriers.A_i).^2) .* (1-abs(rho_i + drho_i).^2);
   p_in1 = settings.tag.p_av * abs(settings.carriers.A_i).^2 / sum(abs(settings.carriers.A_i).^2) .* (1-abs(rho_i - drho_i).^2);
else
   % calculate p_in directly
   p_in0 = abs(settings.carriers.A_i).^2 .* abs(H_i).^2 .* (1-abs(rho_i - drho_i).^2) / 2;
   p_in1 = abs(settings.carriers.A_i).^2 .* abs(H_i).^2 .* (1-abs(rho_i + drho_i).^2) / 2;
end
%     sum up input power [unmod, square 50% duty, mod]
pin01 = [sum(p_in0), sum(p_in1)];
settings.tag.p_in = [pin01(1), mean(pin01), pin01(2)];
%     chip input power for unmod, mod (simple, but may be inaccurate)
ind_fch  = interp1(tagchar_mod.fch, 1:tagchar_mod.settings.nfch, settings.carriers.f_c+settings.carriers.f_i, 'nearest');
pos_pic = nan(length(settings.carriers.f_i), 2);
for i = 1 : length(settings.carriers.f_i)
   pos_pic(i, 1) = interp1(squeeze(tagchar_mod.pin(ind_fch(i), end, :)), 1:tagchar_mod.settings.np, p_in0(i), 'linear','extrap');
   pos_pic(i, 2) = interp1(squeeze(tagchar_mod.pin(ind_fch(i),   1, :)), 1:tagchar_mod.settings.np, p_in1(i), 'linear','extrap');
end
%     sum up chip input power [unmod, square 50% duty, mod]
pic01 = sum( interp1(1:tagchar_mod.settings.np, tagchar_mod.pic, pos_pic, 'linear','extrap'), 1);
settings.tag.p_ic = [pic01(1), mean(pic01), pic01(2)];

% settings.dist
% 10*log10(settings.tag.p_av)+30
% 10*log10(settings.tag.p_in)+30
% 10*log10(settings.tag.p_ic)+30

% variable cleanup
clear('tagchar_mod', 'hh0','hh1','rho','drho', 'ind_fmin','ind_fmax', 'ind_Pmin','ind_Pmax','ind_Pav',...
   'pin0','pin1','pin01','ind_fch','pos_pic','pic01');


% *******************************************************************************************************
% theory values for components (averages)
if ~settings.silent; disp('   Calculating component theory values'); end

% amplitudes (will be used later)
settings.components.c_Mi = 1/2 * ( A_i .* G_Mi .* H_i .* H_Mi .* drho_i * exp(complex(0, -settings.tag.phi_m)) ); % f_i - f_m
settings.components.c_i  =       ( A_i .* G_i .* (K_i + H_i.^2 .* rho_i) ); % f_i
settings.components.c_iM = 1/2 * ( A_i .* G_iM .* H_i .* H_iM .* drho_i * exp(complex(0, +settings.tag.phi_m)) ); % f_i + f_m

% create vectors
c_Mi = ones(settings.N_e,1) * settings.components.c_Mi;
c_i  = ones(settings.N_e,1) * settings.components.c_i;
c_iM = ones(settings.N_e,1) * settings.components.c_iM;


% *******************************************************************************************************
% add noise
if ~settings.silent; disp('   Adding Gaussian white noise'); end

% calculate theoretical noise variance
% ... bandwidth = f_rs, i.e. the sampling rate of these vectors => N0 is identical to variance
% ... create noise level matching the SNR of main simulator
settings.awgn.var_n = 10^(settings.awgn.N_0(1)/10-3)/settings.awgn.N_0(2);

% output SNR values
% ... noise with var_n will be added to real and imaginary part => 2*var_n in output, i.e. - 3 dB
if settings.verbose && ~settings.silent
   disp(sprintf('      Noise level N0 = %6.2f dBm, Correlation coefficients iq = %g, ci = %g',...
      10*log10(settings.awgn.var_n/1e-3), settings.awgn.corr_iq, settings.awgn.corr_ci));
   disp('      Signal to noise levels [c_{i-M}; c_{i}; c_{i+M}]');
   disp(sprintf('         %6.2f dB', 10*log10( abs(settings.components.c_Mi).^2 / settings.awgn.var_n) - 3));
   disp(sprintf('         %6.2f dB', 10*log10( abs(settings.components.c_i ).^2 / settings.awgn.var_n)));
   disp(sprintf('         %6.2f dB', 10*log10( abs(settings.components.c_iM).^2 / settings.awgn.var_n) - 3));
end

% add noise
%     create correlated part between components
n = sqrt(   abs(settings.awgn.corr_iq) *settings.awgn.var_n)*sign(settings.awgn.corr_iq)*complex(1,1)*randn(settings.N_e,1) + ... % i,q correlated
    sqrt((1-abs(settings.awgn.corr_iq))*settings.awgn.var_n)*complex(randn(settings.N_e,1),randn(settings.N_e,1)); % iq uncorrelated
%     add noise to compontents
c_Mi = c_Mi + ...
   sqrt(   abs(settings.awgn.corr_iq) *(1-abs(settings.awgn.corr_ci))*settings.awgn.var_n)*sign(settings.awgn.corr_iq)*randn(size(c_Mi))*complex(1,1) +... % iq corr, comp uncorr
   sqrt((1-abs(settings.awgn.corr_iq))*(1-abs(settings.awgn.corr_ci))*settings.awgn.var_n)*complex(randn(size(c_Mi)),randn(size(c_Mi))) + ... % iq uncorr, comp uncorr
   sqrt(abs(settings.awgn.corr_ci))*sign(settings.awgn.corr_ci)*repmat(n, 1, settings.carriers.N_c); % iq corr/uncorr, comp corr
c_i = c_i + ...
   sqrt(   abs(settings.awgn.corr_iq) *(1-abs(settings.awgn.corr_ci))*settings.awgn.var_n)*sign(settings.awgn.corr_iq)*randn(size(c_i))*complex(1,1) +...
   sqrt((1-abs(settings.awgn.corr_iq))*(1-abs(settings.awgn.corr_ci))*settings.awgn.var_n)*complex(randn(size(c_i)),randn(size(c_i))) + ...
   sqrt(abs(settings.awgn.corr_ci))*sign(settings.awgn.corr_ci)*repmat(n, 1, settings.carriers.N_c);
c_iM = c_iM + ...
   sqrt(   abs(settings.awgn.corr_iq) *(1-abs(settings.awgn.corr_ci))*settings.awgn.var_n)*sign(settings.awgn.corr_iq)*randn(size(c_iM))*complex(1,1) +...
   sqrt((1-abs(settings.awgn.corr_iq))*(1-abs(settings.awgn.corr_ci))*settings.awgn.var_n)*complex(randn(size(c_iM)),randn(size(c_iM))) + ...
   sqrt(abs(settings.awgn.corr_ci))*sign(settings.awgn.corr_ci)*repmat(n, 1, settings.carriers.N_c);


% *******************************************************************************************************
% receiver: correct known errors
if ~settings.silent; disp('   Removing estimates of systematic influences'); end

% create estimates by adding an error to the true values
%     sent carrier amplitudes (i.e., primarily the phase shifts)
sysest.A_i = ...
   (settings.sysest.A_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* abs(A_i) .* ...
   (settings.sysest.A_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* exp(complex(0,angle(A_i)));
%     reader input stage
sysest.G_Mi = ...
   (settings.sysest.G_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* abs(G_Mi) .* ...
   (settings.sysest.G_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* exp(complex(0,angle(G_Mi)));
sysest.G_i = ...
   (settings.sysest.G_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* abs(G_i) .* ...
   (settings.sysest.G_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* exp(complex(0,angle(G_i)));
sysest.G_iM = ...
   (settings.sysest.G_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* abs(G_iM) .* ...
   (settings.sysest.G_err(1) + randn(1,settings.carriers.N_c)*settings.sysest.G_err(2)) .* exp(complex(0,angle(G_iM)));

% remove estimates from vectors
c_Mi = c_Mi ./ repmat(sysest.A_i./abs(sysest.A_i) .* sysest.G_Mi, settings.N_e, 1);
c_i  = c_i  ./ repmat(sysest.A_i .* G_i , settings.N_e, 1); % also remove abs(A_i)
c_iM = c_iM ./ repmat(sysest.A_i./abs(sysest.A_i) .* sysest.G_iM, settings.N_e, 1);


% *******************************************************************************************************
% receiver: estimate phase
if ~settings.silent; disp('   Estimating phases and variances'); end

% average sidebands
avg_c_Mi       = mean(c_Mi, 1);
other.avg_c_i  = mean(c_i,  1); % directly return this value
avg_c_iM       = mean(c_iM, 1);


% % % disp('theorysim');
% % % [abs(avg_c_Mi), abs(avg_c_iM)]
% % % [angle(avg_c_Mi)*180/pi, angle(avg_c_iM)*180/pi]


% estimate phase differences
%     determine order of phase differences (carrier indices)
if strcmpi(settings.carriers.comb, 'all'); % all combinations
   settings.c_ord2 = [1,2];
   settings.c_ord3 = [2,3; 1,3; 1,2];
   settings.c_ord4 = [3,4; 2,4; 1,4; 2,3; 1,3; 1,2];
   settings.c_ord5 = [4,5; 3,5; 2,5; 1,5; 3,4; 2,4; 1,4; 2,3; 1,3; 1,2];
   settings.c_ord6 = [5,6; 4,6; 3,6; 2,6; 1,6; 4,5; 3,5; 2,5; 1,5; 3,4; 2,4; 1,4; 2,3; 1,3; 1,2];
elseif strcmpi(settings.carriers.comb, 'sim'); % simulator
   settings.c_ord2 = [1,2];
   settings.c_ord3 = [1,2; 1,3]; % some correlation!
   settings.c_ord4 = [1,2; 3,4];
   settings.c_ord5 = [1,2; 3,4; 1,5]; % some correlation!
   settings.c_ord6 = [1,2; 3,4; 5,6];
else
   err('Unsupported mode settings.carriers.comb=''%s'' for order of phase differences.', settings.carriers.comb);
end
%     estimate phase shifts (all sidebands)
if settings.carriers.N_c>=2; D_2 = phase_diff(avg_c_Mi, avg_c_iM, settings.c_ord2); end
if settings.carriers.N_c>=3; D_3 = phase_diff(avg_c_Mi, avg_c_iM, settings.c_ord3); end
if settings.carriers.N_c>=4; D_4 = phase_diff(avg_c_Mi, avg_c_iM, settings.c_ord4); end
if settings.carriers.N_c>=5; D_5 = phase_diff(avg_c_Mi, avg_c_iM, settings.c_ord5); end
if settings.carriers.N_c>=6; D_6 = phase_diff(avg_c_Mi, avg_c_iM, settings.c_ord6); end

% expected variances for simple averaging
if settings.carriers.N_c>=2; exp_var.f2 = dist_var(settings.c_ord2, settings); else exp_var.f2 = NaN; end
if settings.carriers.N_c>=3; exp_var.f3 = dist_var(settings.c_ord3, settings); else exp_var.f3 = NaN; end
if settings.carriers.N_c>=4; exp_var.f4 = dist_var(settings.c_ord4, settings); else exp_var.f4 = NaN; end
if settings.carriers.N_c>=5; exp_var.f5 = dist_var(settings.c_ord5, settings); else exp_var.f5 = NaN; end
if settings.carriers.N_c>=6; exp_var.f6 = dist_var(settings.c_ord6, settings); else exp_var.f6 = NaN; end

% variable cleanup
clear('c_iM', 'c_i', 'c_MI', 'ind');


% *******************************************************************************************************
% estimate distance (different estimators)
if ~settings.silent; disp('   Estimating distance using different models'); end

% extract frequencies from settings (for simpler equations)
w_i = 2*pi*settings.carriers.f_i;

% estimates
%     ideal
dist_est.ref = settings.dist;
%     2-frequency, tag known
domega = -[w_i(1)-w_i(2)];
domega = repmat(domega, 1, 2); % LSB, USB symmetrical
dist_est.f2(1) =  mean( ( D_2 - angle(drho_i(1)) + angle(drho_i(2)) ) * settings.c / (2*domega) );
%     2-frequency (no model)
dist_est.f2(2) = mean( D_2 * settings.c / (2*domega) );
%     3-frequency
if settings.carriers.N_c >= 3
   domega = w_i(settings.c_ord3(:,2))' - w_i(settings.c_ord3(:,1))';
   % no model
   dist_est.f3(1) = mean( D_3 * settings.c ./ (2*repmat(domega', 1, 2)) ); % LSB, USB symmetrical
   % constant model
   b = pinv([sign(domega), 2*domega/settings.c])*(D_3(1:end/2)'+D_3(1+end/2:end)')/2; % [B, tau]
   dist_est.f3(2) = b(2);
   % coarse and detailed estimate
   lambda = 2*pi*settings.c / abs(domega(end));
   d_hat = D_3 * settings.c ./ (2*repmat(domega', 1, 2));
   dist_est.f3(3) = floor(4*mean(d_hat([1,3]))/lambda) * lambda/4 + mean(d_hat([2,4]));
else
   dist_est.f3(1:3) = NaN;
end
%     4-frequency
if settings.carriers.N_c >= 4
   domega = w_i(settings.c_ord4(:,2))' - w_i(settings.c_ord4(:,1))';
   % no model
   dist_est.f4(1) = mean( D_4 * settings.c ./ (2*repmat(domega', 1, 2)) ); % LSB, USB symmetrical
   % constant model
   b = pinv([sign(domega), 2*domega/settings.c])*D_4(1:end/2)'; % [B, tau]
   dist_est.f4(2) = b(2);
   % partial quadratic model with slope correction
   b = pinv([sign(domega).*domega.^2, sign(domega), 2*domega])*(D_4(1:end/2)'+D_4(1+end/2:end)')/2; % [B, C, tau]
   dist_est.f4(3) = b(3) * settings.c;
else
   dist_est.f4(1:3) = NaN;
end
%     5-frequency
if settings.carriers.N_c >= 5
   domega = w_i(settings.c_ord5(:,2))' - w_i(settings.c_ord5(:,1))';
   % no model
   dist_est.f5(1) = mean( D_5 * settings.c ./ (2*repmat(domega', 1, 2)) ); % LSB, USB symmetrical
else
   dist_est.f5(1) = NaN;
end
%     6-frequency
if settings.carriers.N_c >= 6
   domega = w_i(settings.c_ord6(:,2))' - w_i(settings.c_ord6(:,1))';
   % no model
   dist_est.f6(1) = mean( D_6 * settings.c ./ (2*repmat(domega', 1, 2)) ); % LSB, USB symmetrical
else
   dist_est.f6(1) = NaN;
end

% output
if settings.verbose && ~settings.silent
   disp(sprintf('      Estimation window: %.2f us (%.0f samples)', settings.N_e/settings.f_rs*1e6, settings.N_e));
   disp(sprintf('\n                           dist [m] |    err [%%] | avg; noise [%%]'));
              disp('      ------------------------------+------------+-----------------')
   disp(sprintf('      2-f, tag known:       %7.3f | %10.2f |      %10.2f',...
      dist_est.f2(1), 100*(dist_est.f2(1)-settings.dist)/settings.dist, 100*sqrt(exp_var.f2)/settings.dist ));
   disp(sprintf('      2-f (no model):       %7.3f | %10.2f |      %10.2f',...
      dist_est.f2(2), 100*(dist_est.f2(2)-settings.dist)/settings.dist, 100*sqrt(exp_var.f2)/settings.dist ));
   disp(sprintf('      3-f (no model):       %7.3f | %10.2f |      %10.2f',...
      dist_est.f3(1), 100*(dist_est.f3(1)-settings.dist)/settings.dist, 100*sqrt(exp_var.f3)/settings.dist ));
   disp(sprintf('      4-f (no model):       %7.3f | %10.2f |      %10.2f',...
      dist_est.f4(1), 100*(dist_est.f4(1)-settings.dist)/settings.dist, 100*sqrt(exp_var.f4)/settings.dist ));
   disp(sprintf('      5-f (no model):       %7.3f | %10.2f |      %10.2f',...
      dist_est.f5(1), 100*(dist_est.f5(1)-settings.dist)/settings.dist, 100*sqrt(exp_var.f5)/settings.dist ));
   disp(sprintf('      6-f (no model):       %7.3f | %10.2f |      %10.2f',...
      dist_est.f6(1), 100*(dist_est.f6(1)-settings.dist)/settings.dist, 100*sqrt(exp_var.f6)/settings.dist ));
end

end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% phase differences between USB (LSB) of all carriers
function D = phase_diff(avg_c_Mi, avg_c_iM, c_order)

% calculate phase differences
for i = 1 : size(c_order, 1)
   D(i                ) = angle( avg_c_Mi(c_order(i,1)) ./ avg_c_Mi(c_order(i,2)) );
   D(i+size(c_order,1)) = angle( avg_c_iM(c_order(i,1)) ./ avg_c_iM(c_order(i,2)) );
end
end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% variance of estimate using all available sidebands
% ... will only match if each carrier is used only once (i.e. no identical entries in c_order) and for
%     uncorrelated noise (correlation does not work => deactivated)
function dvar = dist_var(c_order, settings)

% extract some data from settings for simpler equations
avg_c_Mi = settings.components.c_Mi; % carrier amplitudes @ w_i - w_m
avg_c_iM = settings.components.c_iM; % carrier amplitudes @ w_i + w_m
w_i      = 2 * pi * settings.carriers.f_i; % carrier offset frequencies
N0       = 10^(settings.awgn.N_0(1)/10-3) / settings.awgn.N_0(2); % (single-sided) noise density

% add up "phase" variances; use identical correlation coefficients; no vectorization for simplicity
% ... already divided by frequency offset!!
phivar = 0;
for i = 1 : size(c_order, 1)
   % lower sideband
   phivar_lsb = 1/abs(avg_c_Mi(c_order(i,1)))^2 + 1/abs(avg_c_Mi(c_order(i,2)))^2;% + ... % uncorr
%           settings.awgn.corr_ci*N0/(settings.N_e*abs(avg_c_Mi(c_order(i,1)))^2*abs(avg_c_Mi(c_order(i,2)))^2 ); % corr
   % upper sideband
   phivar_usb = 1/abs(avg_c_iM(c_order(i,1)))^2 + 1/abs(avg_c_iM(c_order(i,2)))^2;% + ... % uncorr
%           settings.awgn.corr_ci*N0/(settings.N_e*abs(avg_c_iM(c_order(i,1)))^2*abs(avg_c_iM(c_order(i,2)))^2 ); % corr
   % average by frequency
   phivar = phivar + (phivar_lsb + phivar_usb) / ( w_i(c_order(i,1)) - w_i(c_order(i,2)) )^2;
end

% avg. "phase" variance plus -> distance 
dvar = settings.c^2 * N0 / ( 4 * settings.N_e * (2*size(c_order, 1))^2 ) * phivar; % conservative

end





% *******************************************************************************************************
% *******************************************************************************************************
% OLD AND RUSTY

% number of used carriers 
% Ncu = (1+sqrt(8*size(c_order, 1)+1))/2; % all sidebands ... size(c_order,1)=N_c,used*(N_c,used-1)/2
% dvar = settings.c^2 * N0 / ( 4 * settings.N_e * Ncu^2*(Ncu-1)^2 ) * phivar; % all sidebands used

% % expected standard deviations due to noise
% 
% %     2-frequency
% exp_var.f2 = settings.c^2*N0 / (settings.N_e*16*pi^2*(settings.carriers.f_i(2)-settings.carriers.f_i(1))^2) *
% for i = 1 : size(settings.c_ord2,1)
%    D_4(i                        ) = angle( avg_c_Mi(settings.c_ord4(i,1)) ./ avg_c_Mi(settings.c_ord4(i,2)) );
%    D_4(i+size(settings.c_ord4,1)) = angle( avg_c_iM(settings.c_ord4(i,1)) ./ avg_c_iM(settings.c_ord4(i,2)) );
% end
% %     2-frequency (assuming uncorrelated noise)
% %     ... one sideband :     ( c/(2*w_1) )^2 * 1/N * ( 1/SNR1 + 1/SNR2 )
% %     ... two sidebands:     ( c/(2*w_1) )^2 * 1/N * ( 1/SNR1 + 1/SNR2 + 1/SNR3 + 1/SNR4 ) / 4
% %                            i.e., approx ( c/(2*w_1) )^2 * 1/N * ( 1/SNR1 + 1/SNR2 ) / 2

% D_2(1) = angle( avg_c_Mi(1) ./ avg_c_Mi(2) );
% D_2(2) = angle( avg_c_iM(1) ./ avg_c_iM(2) );
% %     three-frequency (all 6 sidebands)
% if settings.carriers.N_c >= 3
%    settings.c_ord3 = [2,3; 1,3; 1,2]; % order of carriers
%    for i = 1 : size(settings.c_ord3,1)
%       D_3(i                        ) = angle( avg_c_Mi(settings.c_ord3(i,1)) ./ avg_c_Mi(settings.c_ord3(i,2)) );
%       D_3(i+size(settings.c_ord3,1)) = angle( avg_c_iM(settings.c_ord3(i,1)) ./ avg_c_iM(settings.c_ord3(i,2)) );
%    end
% end
% %     four-frequency (all 12 sidebands)
% if settings.carriers.N_c >= 4
%    settings.c_ord4 = [3,4; 2,4; 1,4; 2,3; 1,3; 1,2]; % order of carriers
%    for i = 1 : size(settings.c_ord4,1)
%       D_4(i                        ) = angle( avg_c_Mi(settings.c_ord4(i,1)) ./ avg_c_Mi(settings.c_ord4(i,2)) );
%       D_4(i+size(settings.c_ord4,1)) = angle( avg_c_iM(settings.c_ord4(i,1)) ./ avg_c_iM(settings.c_ord4(i,2)) );
%    end
% end


% exp_var.f2 = settings.c^2*N0 / (settings.N_e*16*pi^2*(settings.carriers.f_i(2)-settings.carriers.f_i(1))^2) *...
%    ( 1/abs(settings.components.c_iM(1))^2 + 1/abs(settings.components.c_iM(2))^2); % only one sideband
%     3-frequency

%  rho = ( squeeze(tagchar_mod.rho_pav(:, 1, ind_Pav)) + squeeze(tagchar_mod.rho_pav(:, end, ind_Pav)) ) / 2;
% drho = rho - squeeze(tagchar_mod.rho_pav(:, 1, ind_Pav));

% % index for power level
% ind_Pav = round(interp1(tagchar_mod.pav, [1:1:tagchar_mod.settings.np], settings.tag.p_av, 'linear', 'extrap'));
% %     check
% if ind_Pav < ind_Pmin || ind_Pav > ind_Pmax
%    if ~settings.silent; disp('      Power level outside operational range. Using last available level.'); end
%    if ind_Pav < ind_Pmin
%       settings.tag.p_av = tagchar_mod.pav(ind_Pmin);
%    else
%       settings.tag.p_av = tagchar_mod.pav(ind_Pmax);
%    end
% end
% %     output
% if settings.verbose && ~settings.silent
%    disp(sprintf('      Next characteristic power level: %.2e W (%.2f dBm).', tagchar_mod.pav(ind_Pav), 10*log10(tagchar_mod.pav(ind_Pav))+30));
% end

% % calculate linear reflection coefficient model
%  rho = ( squeeze(tagchar_mod.rho_pav(:, 1, ind_Pav)) + squeeze(tagchar_mod.rho_pav(:, end, ind_Pav)) ) / 2;
% drho = rho - squeeze(tagchar_mod.rho_pav(:, 1, ind_Pav));

% c_i  = c_i  + sqrt((1-settings.ncorr_t)*settings.awgn.var_n)*complex(randn(size(c_i )), randn(size(c_i ))) +...
%    sqrt(settings.awgn.corr)*repmat(n, 1, settings.carriers.N_c);
% c_iM = c_iM + sqrt((1-settings.awgn.corr)*settings.awgn.var_n)*randn(size(c_iM))*complex(1,1) + ...
%    sqrt(settings.awgn.corr)*repmat(n, 1, settings.carriers.N_c);
