% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% simple result checker (minimizes lines needed for test routines)
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
% [errors, errortext] = check_result(result, expected, errors, errortext, resultname, method, epsilon)
%    Compares real or complex input variables  (scalar or vector; matrices will be vectorized).
%    The string METHOD defines the type of comparison:
%       mse:    okay if mse(result, expected) <= epsilon(1)
%       range:  okay if all(epsilon(1) <= result-expected <= epsilon(2))
%       relerr: okay if max(|(result-expected)./expected|) <= epsilon(1)
%       equal: (default) okay if result == expected (seperate option for performance reasons)
%    If RESULT and EXPECTED are not equal, the counter ERRORS is incremented and RESULTNAME is added 
%    to ERRORTEXT. If ERRORTEXT is not empty (already reportet errors in prior calls), ERRORTEXT 
%    and RESULTNAME are seperated by a comma.
%    
%
%
% ***** Interface definition *****
% function [errors, errortext] = check_result(result, expected, errors, errortext, resultname,...
%                                             method, epsilon)
%    num          integer that is printed with the error message to help identify the error
%    result       compared to expected (scalar or vector; matrices will be vectorized)
%    expected     expected result (compared to result, scalar or vector; matrices will be vectorized)
%    errortext    cumulative string containing the last error texts
%    resultname   will be added to errortext in case of an error
%    method       (optional) comparison method {'mse', 'range', 'relerr', default: 'equal'}
%    epsilon      (optional) maximum deviation for all methods except 'equal'
%                 .) scalar for method 'mse','relerr' and 'equal' (default: eps)
%                 .) vector [min, max] or scalar (range +/- epsilon) for 'range' (default: +/- eps)
%                 ... eps is the floating point precision of Matlab
%
%    errors       amount of errors including this function call
%    errortext    cumulative string containing the last error texts
%                 including this one (comma seperated)
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

function [errors, errortext] = check_result(num, result, expected, errors, errortext, resultname, method, epsilon)


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input parameters
if nargin < 6
   error('Not enough input arguments.');
elseif nargin == 6
   method  = 'equal';
elseif nargin == 7
   epsilon = eps;
elseif nargin > 8
   error('Too many input arguments.');
end

% -> column vectors
result   = result(:);
expected = expected(:);

% make sure epsilon has the correct format
if strcmpi(method, 'range') % "range" needs a vector
   if isscalar(epsilon) 
      epsilon = [-epsilon, epsilon];
   end
elseif ~strcmpi(method, 'equal') % everything else except 'equal' needs a scalar 
   if ~isscalar(epsilon)
      warn('Parameter epsilon has to be scalar for comparison method "%s".', method);
      epsilon = epsilon(1);
   end
end 

% switch method from 'relerr' to 'equal' if expected contains zeros
if strcmpi(method, 'relerr') && any( expected == 0 )
   warn('Expected value contains zeros. Switching from method "relerr" to "mse".');
   method = 'mse';
end

% range makes no sense for complex numbers (re/imag or abs/phase)?
if strcmpi(method, 'range') && ~(isreal(result) && isreal(expected))
   err('Complex numbers not supported for method ''%s''.', lower(method));
end


% *******************************************************************************************************
% checks

% size
if sum(size(result) ~= size(expected)) %~= 0
   errors = errors + 1;
   if isempty(errortext)
      errortext = sprintf('%i:size mismatch(%s)', num, resultname);
   else
      errortext = sprintf('%s, %i:size mismatch(%s)', errortext, num, resultname);
   end
   return
end
   
% compare to expected results
equal = false;
switch lower(method)
   case 'mse'
      equal = mean(abs((result - expected).^2)) <= epsilon;
   case 'range'
      equal =  epsilon(1) <= min(result - expected) && max(result - expected) <= epsilon(2);
   case 'relerr'
      equal = max(abs((result - expected) ./ expected)) <= epsilon;
   case 'equal'
      equal = all(result == expected);
   otherwise
      err('Unsupported comparison method %s.', method);
end
   
% error text
if ~equal
   errors = errors + 1;    
   if isempty(errortext)
      errortext = sprintf('%i:%s', num, resultname);
   else
      errortext = sprintf('%s, %i:%s', errortext, num, resultname);
   end
end
