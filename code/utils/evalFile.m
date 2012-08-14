%  smanifold, Copyright (C) 2009-2012, Stefan Sommer (sommer@diku.dk)
%  https://github.com/nefan/smanifold.git
% 
%  This file is part of smanifold.
% 
%  smanifold is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  smanifold is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with smanifold.  If not, see <http://www.gnu.org/licenses/>.
%  

%
% evalute contents of file in variable evalFileName
%
% not perfect, but works ok
% 

assert(exist(evalFileName,'file') > 0);
fid = fopen(evalFileName);

tline = fgetl(fid);
cline = [];
while ischar(tline)
    %disp(tline);
    %eval(tline);
    kk = strfind(tline,'%'); % remove comments 
    kkk = strfind(tline,39); % apostrofes
    for jj = 1:length(kk)
        k = kk(jj);
        skip = false;
        if (mod(length(kkk),2) == 0)
            for ii = 1:length(kkk)/2
                if kkk(2*(ii-1)+1) < k && k < kkk(2*(ii-1)+2)
                    skip = true;
                end
            end
        end
        if ~skip && ~isempty(tline)
            tline = tline(1:k-1);
        end
    end
    cline = [cline sprintf('\n') tline];
    tline = fgetl(fid);
end
eval(cline); % eval entire file at once

fclose(fid);
