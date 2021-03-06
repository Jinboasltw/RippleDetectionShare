function [dat, w] = ft_preproc_denoise(dat, refdat, hilbertflag)

% FT_PREPROC_DENOISE performs a regression of the matrix dat onto refdat, and
% subtracts the projected data. This is for the purpose of removing signals generated
% by coils during continuous head motion tracking, for example.
%
% Use as
%   [dat] = ft_preproc_denoise(dat, refdat, hilbertflag)
% where
%   dat         data matrix (Nchan1 X Ntime)
%   refdat      data matrix (Nchan2 X Ntime)
%   hilbertflag boolean, regress out the real and imaginary parts of the Hilbert 
%                 transformed signal, this is only meaningful for narrow band 
%                 reference data (default = false)
%
% The number of channels of the data and reference data does not have to be the same.
%
% If the data contains NaNs, the output of the affected channel(s) will be all NaN.
%
% See also PREPROC

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<3
  hilbertflag = 0;
end

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:))) || any(isnan(refdat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

n1 = size(dat,2);
n2 = size(refdat,2);
m1 = nanmean(dat,2);
m2 = nanmean(refdat,2);

%remove mean
refdat  = refdat-m2(:,ones(n2,1));
tmpdat  = dat-m1(:,ones(n1,1));

%do hilbert transformation
if hilbertflag>0
  hrefdat = hilbert(refdat')';
  refdat  = [real(hrefdat);imag(hrefdat)];
end

c12 = tmpdat*refdat'; %covariance between signals and references
c1  = refdat*refdat'; %covariance between references and references
w   = (pinv(c1)*c12')'; %regression weights

%subtract
dat = dat-w*refdat;
