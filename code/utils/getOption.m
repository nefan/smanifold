%
%  segframe, Copyright (C) 2009-2012, Stefan Sommer (sommer@diku.dk)
%  https://github.com/nefan/segframe.git
% 
%  This file is part of segframe.
% 
%  segframe is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  segframe is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with segframe.  If not, see <http://www.gnu.org/licenses/>.
%  

function optionOrEmpty = getOption(options,field,varargin)

    hasDefault = false;
    if size(varargin,2) == 1
        hasDefault = true;
        default = varargin{1};
    end

    if isfield(options,field) && getfield(options,field)
        optionOrEmpty = getfield(options,field);
    else if hasDefault
        optionOrEmpty = default;
    else
        optionOrEmpty = [];
    end
    end
end
