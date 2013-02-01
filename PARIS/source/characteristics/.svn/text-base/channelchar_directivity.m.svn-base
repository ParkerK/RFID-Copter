% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: channel - reader/tag antenna directivity
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


% ******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% one last chance to reconsider
disp(sprintf('WARNING: This might overwrite existing files!\n'));
reply = input('Are you REALLY sure y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end
disp(' ');

% paths, etc
%     path to globalinit.m
path = fullfile(fileparts(mfilename('fullpath')), '..');
%     Matlab does not resolve symbolic links; if we're running on Unix, let the system resolve them
if isunix
   [dummy, path] = system(sprintf('readlink -f "%s"', path));
end
%     add this path and call globalinit script
addpath(path); clear('dummy', 'path');
globalinit('init');

% "switch to here"
cd(fileparts(mfilename('fullpath')));



% ******************************************************************************************************
% settings / generate base vectors (angles for azimuth and elevation)

% settings
settings.res        =       1; % number of points per degree
settings.plotfolder = 'plots'; % subfolder for plots

% base vectors for azimuth and elevation
% ... include 0 and 180,360 degree to simplify interpolation and avoid extrapolation (esp. 2D)
settings.az = 0 : 1/settings.res : 360;
settings.el = 0 : 1/settings.res : 180;

return




% ******************************************************************************************************
% UWB tag antenna (measured)

% load data
usedmatfiles = 'meas_uwb-tag_antenna_v2'; 
data         = load(usedmatfiles);

% interpolate to new angles
az_gain = interp1(data.az, data.az_gain, settings.az, 'spline', 'extrap');
el_gain = interp1(data.el, data.el_gain, settings.el, 'spline', 'extrap');

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_uwb-tag-v2'; %#ok<NASGU>
characteristic = 'radiation pattern for uwb tag antenna (avg. 0.5-1.5GHz): gain(az,el) = az_gain(az) * el_gain(el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data
az = settings.az(:); el = settings.el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d  = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure('visible','off');
hold off; polar(el *pi/180, el_gain, 'b.');
hold on;  polar(az *pi/180, az_gain, 'r.');
polar(data.el *pi/180, data.el_gain, 'b-');
polar(data.az *pi/180, data.az_gain, 'r-');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): UWB TAG ANTENNA, 0.5-1.5GHz AVG (UWB)', '', '');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure('visible','off'); hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on;  view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): UWB TAG ANTENNA, 0.5-1.5GHz AVG (UWB)', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');



% ******************************************************************************************************
% Vivaldi 4x1-antenna-array, 20cm spacing, UWB (measured)

% load data
usedmatfiles = 'meas_vivaldi_antenna_4x1-20cm'; 
data         = load(usedmatfiles);

% interpolate to new angles
az_gain = interp1(data.az, data.az_gain, settings.az, 'spline', 'extrap');
el_gain = interp1(data.el, data.el_gain, settings.el, 'spline', 'extrap');

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_vivaldi_4x1-20cm'; %#ok<NASGU>
characteristic = 'radiation pattern for vivaldi 4x1-array (avg. 0.5-1.5GHz): gain(az,el) = az_gain(az) * el_gain(el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data
az = settings.az(:); el = settings.el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d  = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure('visible','off');
hold off; polar(el *pi/180, el_gain, 'b.');
hold on;  polar(az *pi/180, az_gain, 'r.');
polar(data.el *pi/180, data.el_gain, 'b-');
polar(data.az *pi/180, data.az_gain, 'r-');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY, 20cm, 0.5-1.5GHz AVG (UWB)', '', '');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure('visible','off'); hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on;  view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY, 20cm, 0.5-1.5GHz AVG (UWB)', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');



% *******************************************************************************************************
% Vivaldi antenna (custom design) at 0.7:0.1:1.2 GHz

for f = 0.7 : 0.1 : 1.2; % GHz
   
   % load data
   usedmatfiles = sprintf('meas_antenna_vivaldi-v1_%.0fGHz%.0f.mat',floor(f), 10*(f-floor(f)));
   data         = load(usedmatfiles);
   
   % interpolate to new angles
   % ... keep in mind that meshgrid creates transposed matrices (x is column, y is row)
   [ a,  e] = meshgrid(data.az, data.el);
   [ai, ei] = meshgrid(settings.az, settings.el); 
   gain = interp2(a, e, data.gain', ai, ei, 'spline', NaN)';
   %     (re-)normalize
   gain = gain / max(max(gain));
   
   % prepare everything and save
   %     header
   matfilename    = sprintf('channelchar_directivity_vivaldi-v1_%.0fGHz%.0f',floor(f), 10*(f-floor(f)));
   characteristic = sprintf('radiation pattern for custom design vivaldi reader antenna, f=%.1fGHz: gain(az,el) (3D)', f); %#ok<NASGU>
   createdon      = datestr(now, 0); %#ok<NASGU>
   createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
   %     data
   az = settings.az(:); el = settings.el(:);
   full3d = true; % true: gain(az,el);   false: az_gain(az), el_gain(el)
   % save
   save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
      'az','el','gain', 'full3d');
   
   % plots
   %     2D (first hold off to draw grid)
   az_gain = gain(:, find(el==90));
   el_gain = gain(1,            :)';
   figure('visible','off');
   hold off; polar(el *pi/180, el_gain, 'b-');
   hold on;  polar(az *pi/180, az_gain, 'r-');
   hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside');
   setlabels(sprintf('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI CUSTOM DESIGN (%.1f GHz)',f), '', '');
   savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png'); close(gcf);
   %     3D (just for the view)
   [x,y,z] = sph2cart(ai*pi/180, pi/2-ei*pi/180, gain');
   figure('visible','off'); surface(x, y, z, gain', 'EdgeColor','k', 'FaceColor','interp'); 
   grid on; setcolorbar(); view([30, 60]);
   setlabels(sprintf('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI CUSTOM DESIGN (%.1f GHz)',f), 'x', 'y', 'z');
   savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png'); close(gcf);
      
end



% *******************************************************************************************************
% Vivaldi antenna (custom design) at for UWB (0.7-1.2 GHz)

% at first frequency
f = 0.7; % GHz
%     load data
usedmatfiles = sprintf('meas_antenna_vivaldi-v1_%.0fGHz%.0f.mat',floor(f), 10*(f-floor(f)));
data         = load(usedmatfiles);
%     (re-)normalize
data.gain = data.gain / max(max(data.gain));
%     interpolate to new angles
%     ... keep in mind that meshgrid creates transposed matrices (x is column, y is row)
[ a,  e] = meshgrid(data.az, data.el);
[ai, ei] = meshgrid(settings.az, settings.el);
gain = interp2(a, e, data.gain', ai, ei, 'spline', NaN)';

% at all other frequencies
for f = 0.8 : 0.1 : 1.2; % GHz
   % load data
   usedmatfiles = sprintf('%s,\nmeas_antenna_vivaldi-v1_%.0fGHz%.0f.mat',usedmatfiles, floor(f), 10*(f-floor(f)));
   data         = load(sprintf('meas_antenna_vivaldi-v1_%.0fGHz%.0f.mat',floor(f), 10*(f-floor(f))));
   %     (re-)normalize
   data.gain = data.gain / max(max(data.gain));
   %     interpolate to new angles and sum up for averaging
   [ a,  e] = meshgrid(data.az, data.el);
   gain = gain + interp2(a, e, data.gain', ai, ei, 'spline', NaN)';
end

% (re-)normalize (... includes averaging)
gain = gain / max(max(gain));

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_vivaldi-v1_uwb';
characteristic = 'radiation pattern for custom design vivaldi reader antenna, UWB (f=0.7 to 1.2GHz): gain(az,el) (3D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
%     data
az = settings.az(:); el = settings.el(:);
full3d = true; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','el','gain', 'full3d');

% plots
%     2D (first hold off to draw grid)
az_gain = gain(:, find(el==90));
el_gain = gain(1,            :)';
figure('visible','off');
hold off; polar(el *pi/180, el_gain, 'b-');
hold on;  polar(az *pi/180, az_gain, 'r-');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI CUSTOM DESIGN (UWB 0.7-1.2 GHz)', '', '');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'eps'); close(gcf);
%     3D (just for the view)
[x,y,z] = sph2cart(ai*pi/180, pi/2-ei*pi/180, gain');
figure('visible','off'); surface(x, y, z, gain', 'EdgeColor','k', 'FaceColor','interp');
grid on; setcolorbar(); view([30, 60]);
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI CUSTOM DESIGN (UWB 0.7-1.2 GHz)', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'eps'); close(gcf);



% *******************************************************************************************************
% Vivaldi 4x1-antenna-array (10,15,20 cm) at 0.7:0.1:1.2 GHz

for d = 10 : 5 : 20 % cm antenna spacing
   for f = 0.7 : 0.1 : 1.2; % GHz
      
      % load data
      usedmatfiles = sprintf('meas_antenna_vivaldi-4x1-array-v1-%2.0fcm_%.0fGHz%.0f',d, floor(f), 10*(f-floor(f)));
      data         = load(usedmatfiles);
      
      % interpolate to new angles
      % ... keep in mind that meshgrid creates transposed matrices (x is column, y is row)
      [ a,  e] = meshgrid(data.az, data.el);
      [ai, ei] = meshgrid(settings.az, settings.el);
      gain = interp2(a, e, data.gain', ai, ei, 'spline', NaN)';
      %     (re-)normalize
      gain = gain / max(max(gain));
      
      % prepare everything and save
      %     header
      matfilename    = sprintf('channelchar_directivity_vivaldi-4x1-array-v1-%2.0fcm_%.0fGHz%.0f',d, floor(f), 10*(f-floor(f)));
      characteristic = sprintf('radiation pattern for custom design vivaldi 4x1 reader antenna array, f=%.1fGHz, d=%2.0fcm: gain(az,el) (3D)', f,d); %#ok<NASGU>
      createdon      = datestr(now, 0); %#ok<NASGU>
      createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
      %     data
      az = settings.az(:); el = settings.el(:);
      full3d = true; % true: gain(az,el);   false: az_gain(az), el_gain(el)
      % save
      save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
         'az','el','gain', 'full3d');
      
      % plots
      %     2D (first hold off to draw grid)
      az_gain = gain(:, find(el==90));
      el_gain = gain(1,            :)';
      figure('visible','off');
      hold off; polar(el *pi/180, el_gain, 'b-');
      hold on;  polar(az *pi/180, az_gain, 'r-');
      hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside');
      setlabels(sprintf('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY (%.0f cm spacing,   %.1f GHz)',d,f), '', '');
      savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png'); close(gcf);
      %     3D (just for the view)
      [x,y,z] = sph2cart(ai*pi/180, pi/2-ei*pi/180, gain');
      figure('visible','off'); surface(x, y, z, gain', 'EdgeColor','k', 'FaceColor','interp'); 
      grid on; setcolorbar(); view([30, 60]);
      setlabels(sprintf('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY (%.0f cm spacing,   %.1f GHz)',d,f), 'x', 'y', 'z');
      savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png'); close(gcf);
      
   end
end



% *******************************************************************************************************
% Vivaldi 4x1-antenna-array (10,15,20 cm) for UWB (0.7-1.2 GHz)

for d = 10 : 5 : 20 % cm antenna spacing
   
   % at first frequency
   f = 0.7; % GHz
   %     load data
   usedmatfiles = sprintf('meas_antenna_vivaldi-4x1-array-v1-%2.0fcm_%.0fGHz%.0f',d, floor(f), 10*(f-floor(f)));
   data         = load(usedmatfiles);
   %     (re-)normalize
   data.gain = data.gain / max(max(data.gain));
   %     interpolate to new angles
   %     ... keep in mind that meshgrid creates transposed matrices (x is column, y is row)
   [ a,  e] = meshgrid(data.az, data.el);
   [ai, ei] = meshgrid(settings.az, settings.el);
   gain = interp2(a, e, data.gain', ai, ei, 'spline', NaN)';
   
   % at all other frequencies
   for f = 0.8 : 0.1 : 1.2; % GHz
      % load data
      usedmatfiles = sprintf('%s,\nmeas_antenna_vivaldi-4x1-array-v1-%2.0fcm_%.0fGHz%.0f',usedmatfiles, d, floor(f), 10*(f-floor(f)));
      data         = load(sprintf('meas_antenna_vivaldi-4x1-array-v1-%2.0fcm_%.0fGHz%.0f',d, floor(f), 10*(f-floor(f))));
      %     (re-)normalize
      data.gain = data.gain / max(max(data.gain));
      %     interpolate to new angles and sum up for averaging
      [ a,  e] = meshgrid(data.az, data.el);
      gain = gain + interp2(a, e, data.gain', ai, ei, 'spline', NaN)';
   end
   
   % (re-)normalize (... includes averaging)
   gain = gain / max(max(gain));
     
   % prepare everything and save
   %     header
   matfilename    = sprintf('channelchar_directivity_vivaldi-4x1-array-v1-%2.0fcm_uwb',d);
   characteristic = sprintf('radiation pattern for custom design vivaldi 4x1 reader antenna array, d=%2.0fcm, UWB (f=0.7 to 1.2GHz): gain(az,el) (3D)', d); %#ok<NASGU>
   createdon      = datestr(now, 0); %#ok<NASGU>
   createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
   %     data
   az = settings.az(:); el = settings.el(:);
   full3d = true; % true: gain(az,el);   false: az_gain(az), el_gain(el)
   % save
   save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
      'az','el','gain', 'full3d');
   
   % plots
   %     2D (first hold off to draw grid)
   az_gain = gain(:, find(el==90));
   el_gain = gain(1,            :)';
   figure('visible','off');
   hold off; polar(el *pi/180, el_gain, 'b-');
   hold on;  polar(az *pi/180, az_gain, 'r-');
   hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside');
   setlabels(sprintf('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY (%.0f cm spacing, UWB 0.7-1.2 GHz)',d), '', '');
   savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png'); close(gcf);
   %     3D (just for the view)
   [x,y,z] = sph2cart(ai*pi/180, pi/2-ei*pi/180, gain');
   figure('visible','off'); surface(x, y, z, gain', 'EdgeColor','k', 'FaceColor','interp');
   grid on; setcolorbar(); view([30, 60]);
   setlabels(sprintf('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY (%.0f cm spacing, UWB 0.7-1.2 GHz)',d), 'x', 'y', 'z');
   savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png'); close(gcf);
      
end




% *******************************************************************************************************
% Vivaldi 4x1-antenna-array (20 cm) with metal backplane 20cm from array for UWB (0.7-1.2 GHz)
   
% average over all frequencies
az_gain = zeros(size(settings.az));
el_gain = zeros(size(settings.el));
usedmatfiles = '';
for f = 0.7 : 0.1 : 1.2; % GHz
   % load data
   usedmatfiles = sprintf('%s,\nmeas_antenna_vivaldi-4x1-array-v1-20cm_%.0fGHz%.0f_met',usedmatfiles, floor(f), 10*(f-floor(f)));
   data         = load(sprintf('meas_antenna_vivaldi-4x1-array-v1-20cm_%.0fGHz%.0f_met', floor(f), 10*(f-floor(f))));
   %     make sure the data is normalized
   data.az_gain = data.az_gain / max(data.az_gain);
   data.el_gain = data.el_gain / max(data.el_gain);
   
   % interpolate to new angles and sum up for averaging
   az_gain = az_gain + interp1([data.az, 360], [data.az_gain, data.az_gain(1)], settings.az, 'linear');
   el_gain = el_gain + interp1([data.el, 180], [data.el_gain, data.el_gain(1)], settings.el, 'linear');
end

% (re-)normalize (... includes averaging)
az_gain = az_gain / max(az_gain);
el_gain = el_gain / max(el_gain);

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_vivaldi-4x1-array-v1-20cm-met_uwb';
characteristic = 'radiation pattern for custom design vivaldi 4x1 reader antenna array, d=20cm, metal backplane, UWB (f=0.7 to 1.2GHz): gain(az,el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
%     data
az = settings.az(:); el = settings.el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure('visible','off');
hold off; polar(el *pi/180, el_gain, 'b-');
hold on;  polar(az *pi/180, az_gain, 'r-');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY W. METAL BACKPL. (20 cm spacing, UWB 0.7-1.2 GHz)', '', '');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure('visible','off'); hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on;  view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): VIVALDI 4x1 ARRAY W. METAL BACKPL. (20 cm spacing, UWB 0.7-1.2 GHz)', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');
     



% ******************************************************************************************************
% Intermec reader antenna IA39b

% load data
usedmatfiles = 'meas_antenna_intermec-ia39b';
data         = load(usedmatfiles);

% interpolate to new angles
az_gain = interp1(data.az, data.az_gain, settings.az, 'spline', 'extrap');
el_gain = interp1(data.el, data.el_gain, settings.el, 'spline', 'extrap');

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_intermec-ia39b'; %#ok<NASGU>
characteristic = 'radiation pattern for intermec ia39b reader antenna: gain(az,el) = az_gain(az) * el_gain(el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data
az = settings.az(:); el = settings.el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d  = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure('visible','off');
hold off; polar(el *pi/180, el_gain, 'b.');
hold on;  polar(az *pi/180, az_gain, 'r.');
polar(data.el *pi/180, data.el_gain, 'b-');
polar(data.az *pi/180, data.az_gain, 'r-');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): INTERMEC IA39b', '', '');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure('visible','off'); hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on;  view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): INTERMEC IA39b', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');



% ******************************************************************************************************
% ideal lambda/4 dipole 

% dipole antenna pattern (dipole in z-axis)
az = settings.az(:); el = settings.el(:);
az_gain = ones(size(az));
el_gain = ( cos(pi/2 * cos(el*pi/180)) - cos(pi/2) ) ./ sqrt( 1 - cos(el*pi/180).^2 );

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_l4-dipole'; %#ok<NASGU>
characteristic = 'radiation pattern for lambda/4 dipole: gain(az,el) = az_gain(az) * el_gain(el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data
az = az(:); el = el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d  = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure('visible','off');
hold off; polar(el *pi/180, el_gain, 'b');
hold on;  polar(az *pi/180, az_gain, 'r');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): \lambda/4 DIPOLE', '', '');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-2d',matfilename)), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure('visible','off'); hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on; view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): \lambda/4 DIPLOLE', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');



% ******************************************************************************************************
% ideal lambda/2 dipole 

% dipole antenna pattern (dipole in z-axis)
az = settings.az(:); el = settings.el(:);
az_gain = ones(size(az));
el_gain = ( cos(pi * cos(el*pi/180)) - cos(pi) ) ./ sqrt( 1 - cos(el*pi/180).^2 );

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_l2-dipole'; %#ok<NASGU>
characteristic = 'radiation pattern for lambda/2 dipole: gain(az,el) = az_gain(az) * el_gain(el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data
az = az(:); el = el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d  = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure;
hold off; polar(el *pi/180, el_gain, 'b');
hold on;  polar(az *pi/180, az_gain, 'r');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): \lambda/2 DIPOLE', '', '');
savefigure(gcf, fullfile(settings.plotfolder, matfilename), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure; hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on; view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): \lambda/2 DIPLOLE', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');



% ******************************************************************************************************
% approximation for omnidirectional UWB tag antenna based on ideal lambda/4 dipole
% ... plus 2.5, 2.7, 3, 3.4, 3.7 dB main lobe for 0.7:0.1:1.2 GHz => average 3dB

% dipole antenna pattern (dipole in z-axis)
sidelobe = 10^(-3/20); % 3dB sidelobe attenuation (main lobe is at 0dB)
az = settings.az(:); el = settings.el(:);
az_gain = ones(size(az)) .* (1-(1-sidelobe)/2 + (1-sidelobe)/2*cos(2*az*pi/180));
el_gain = ( cos(pi/2 * cos(el*pi/180)) - cos(pi/2) ) ./ sqrt( 1 - cos(el*pi/180).^2 );

% prepare everything and save
%     header
matfilename    = 'channelchar_directivity_uwb-tag-v1'; %#ok<NASGU>
characteristic = 'radiation pattern for custom UWB tag antenna: gain(az,el) = az_gain(az) * el_gain(el) (2x2D)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data
az = az(:); el = el(:);
az_gain(isnan(az_gain)) = 0; az_gain = az_gain(:);
el_gain(isnan(el_gain)) = 0; el_gain = el_gain(:);
full3d  = false; % true: gain(az,el);   false: az_gain(az), el_gain(el)
% save
save(fullfile('channel',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'az','az_gain', 'el','el_gain', 'full3d');

% 2D plots (first hold off to draw grid)
figure;
hold off; polar(el *pi/180, el_gain, 'b');
hold on;  polar(az *pi/180, az_gain, 'r');
hold off; setlegend({'elevation','azimuth'}, 'NorthEastOutside'); 
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): UWB TAG ANTENNA', '', '');
savefigure(gcf, fullfile(settings.plotfolder, matfilename), 'png');

% 3D plot (just for the view)
%     transform from sphere to cartesian coordinates ... gain(az,el) = gain(az) * gain(el)
x = (el_gain .* sin(el*pi/180)) * (az_gain .* cos(az*pi/180))';
y = (el_gain .* sin(el*pi/180)) * (az_gain .* sin(az*pi/180))';
z = (el_gain .* cos(el*pi/180)) * (az_gain                  )';
rho = sqrt( x.^2 + y.^2 + z.^2 );
%     plot
figure; hold on;
surface(x, y, z, rho, 'EdgeColor','none', 'FaceColor','interp');
plot3( az_gain .* cos(az*pi/180), az_gain .* sin(az*pi/180), zeros(size(az)), 'k-'); % horz pattern
plot3( el_gain .* sin(el*pi/180), zeros(size(el)), el_gain .* cos(el*pi/180), 'k-'); % vert pattern
% plot3( cos(az*pi/180), sin(az*pi/180), zeros(size(az)), 'k--'); % unit circe xy
% plot3( sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--'); % unit circle yz
% plot3(-sin(el*pi/180), zeros(size(el)), cos(el*pi/180), 'k--');
hold off; grid on; view([30, 60]); setlegend({'3D=2x2D', 'horizontal, vertical'}, 'NorthWest');
setlabels('ANTENNA RADIATION PATTERN (DIRECTIVITY): UWB TAG ANTENNA', 'x', 'y', 'z');
savefigure(gcf, fullfile(settings.plotfolder, sprintf('%s-3d',matfilename)), 'png');
