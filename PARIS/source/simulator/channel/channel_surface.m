% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - surface (reflection/transmission)
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
% version = channel_surface()
%    Just returns the version number (string).
% gain = channel_surface(settings)
%    Calculates the transmission/reflection coefficient (amplitude gain) of a flat surface using the 
%    Fresnel equations. The incident wave can either be polarized parallel or perpenticular to the 
%    surface, mixed polarizations and depolarization are not considered. The function furthermore assumes
%    a refractive index of n1=1 (air) for the main transmission medium (can be changed via
%    internalsettings within this function). Transmission is modeled as "1-reflection" (no absorbtion).
%     
%    
%
%
% ***** Interface definition *****
% gain = channel_surface(settings)
%    surface   struct containing the surfaces's/wave's parameters
%       .mode      mode for the surface: 0="none/off", 1="transmit", 2="reflect"
%       .n2        refractive index of surface material (typ n2>1)
%       .dim       dimension of the surface's normal vector (e.g. 3 means z in [x,y,z])
%       .pol_dim   dimension of polarization vector (electric field)  (e.g. 3 means z in [x,y,z])
%       .aoi       angle of incidence (0<=aoi<=90) of the wave in degree; can be a multidimensional array
%                  (the aoi is measured between the wave's direction of movement and the surface normal)
%
%    gain        reflection/transmission coefficient (amplitude gain)
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


function gain = channel_surface(settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   gain = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% refractive index of main transmission medium (air: n=1)
internalsettings.n1 = 1; 

% modes for surfaces (have to be identical to settings in CHANNEL_MAIN and CHANNEL_NEWPOS)
internalsettings.surf_off  = 0; % neither reflection, nor transmission => surface inactive
internalsettings.surf_tran = 1; % "reflect" mode
internalsettings.surf_refl = 2; % "transmit" mode


% *******************************************************************************************************
% input parameter checks / prepare parameters

% check contents of surface struct
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'mode', 'n2', 'dim', 'pol_dim', 'aoi'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% neither transmission, nor reflection => return
if settings.mode == internalsettings.surf_off
   gain = 0;
   return
end

% check aoi
if any(settings.aoi(:) < 0) || any(settings.aoi(:) > 90)
   err('Angle of incidence (settings.aoi) has to be 0<=aoi<=90 degree.');
end


% *******************************************************************************************************
% calculate reflection / transmission coefficient

switch settings.mode
   case internalsettings.surf_tran % transmit ... = 1 - reflect for power (...no absorbtion)
      gain = 1 - reflection(settings.aoi*pi/180, internalsettings.n1, settings.n2, settings).^2;
   case internalsettings.surf_refl % reflect
      gain = reflection(settings.aoi*pi/180, internalsettings.n1, settings.n2, settings);
   otherwise
      err('Unsupported surface.mode=%i (0=off, 1=transmit, 2=reflect).', settings.mode);
end

end


% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% reflection coefficient (Fresnel equations)
function gain = reflection(a, n1, n2, settings)
% polarization perpenticular to reflective surface
if settings.dim == settings.pol_dim
   gain = (-n2*cos(a)+n1*sqrt(1-(n1/n2*sin(a)).^2) ) ./ ( n2*cos(a)+n1*sqrt(1-(n1/n2*sin(a)).^2) );
% polarization along (parallel to) reflective surface
else
   gain = ( n1*cos(a)-n2*sqrt(1-(n1/n2*sin(a)).^2) ) ./ ( n1*cos(a)+n2*sqrt(1-(n1/n2*sin(a)).^2) );
end
end

