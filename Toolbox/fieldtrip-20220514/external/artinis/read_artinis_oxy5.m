function data = read_artinis_oxy55(filename, header, begsample, endsample, chanindx)
% reads Artinix oxy55-files into FieldTrip format
%
% use as
%   header = read_artinis_oxy5(filename)
% or 
%   event  = read_artinis_oxy5(filename, read_event)
% where read_event is a Boolean with the value true, or  
%   data   = read_artinis_oxy55(filename, header, [begsample], [endsample], [chanindx])
% where begsample, endsample and chanindx are optional.
%
% The returned variables will be in FieldTrip style. 
%
% See also FT_READ_HEADER, FT_READ_DATA, READ_ARTINIS_OXU5

% You are using the FieldTrip NIRS toolbox developed and maintained by 
% Artinis Medical Systems (http://www.artinis.com). For more information
% on FieldTrip, see http://www.fieldtriptoolbox.org
% 
% This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% 
% Creative Commons Attribution-ShareAlike 4.0 International License:
% -----------------------------------
% You are free to:
% 
%     Share - copy and redistribute the material in any medium or format
%     Adapt - remix, transform, and build upon the material
%     for any purpose, even commercially.
% 
%     The licensor cannot revoke these freedoms as long as you follow the 
%     license terms.
% 
% Under the following terms:
% 
%     Attribution - You must give appropriate credit, provide a link to 
%                    the license, and indicate if changes were made. You 
%                    may do so in any reasonable manner, but not in any way 
%                    that suggests the licensor endorses you or your use.
% 
%     ShareAlike - If you remix, transform, or build upon the material, 
%                   you must distribute your contributions under the same 
%                   license as the original.
% 
%     No additional restrictions - You may not apply legal terms or 
%                                   technological measures that legally 
%                                   restrict others from doing anything the 
%                                   license permits.
% 
% -----------------------------------
% 
% This toolbox is not to be used for medical or clinical purposes.
% 
% Copyright (c) 2021 by Artinis Medical Systems.
% Contact: askforinfo@artinis.com

if nargin > 2 && islogical(header)
  % header must be defined if begsample and endsample are used
  error('wrong type of header variable.')
end

if nargin == 1 || nargin == 2 && islogical(header) && ~header    
  data = read_oxy5_header(filename);
  if isfield(data, 'opto')
    % ensure that the optode definition is according to the latest standards
    data.opto = ft_datatype_sens(data.opto);
  end
elseif nargin == 2  && islogical(header) && header
  data = read_oxy3_event(filename);
else % nargin > 1 && ~islogical(header)  
  if nargin < 5
    chanindx = 1:header.nChans;
    if nargin < 4
      endsample = header.nSamples;
      if nargin < 3
        begsample = 1;
      end
    end
  end  
  data = read_oxy5_data(filename, 0, header.nSamples); % sample subselection does not work yet
  
  if endsample > size(data, 2)
    warning('Cannot deliver all requested samples, nan''ing %d sample(s)', endsample-size(data, 2));
    data(:, end:endsample) = nan;
  end
  
  data = data(chanindx, begsample:endsample);
end
