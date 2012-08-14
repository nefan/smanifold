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

function exit = startmulticore(pmulticoreSettings)


global multicoreSettings;
multicoreSettings = pmulticoreSettings;

exit = false;

if multicoreSettings.processId > 0
    fprintf('slave process %d waiting for jobs...\n',multicoreSettings.processId);
    nrRuns = startmulticoreslave(multicoreSettings.conf.multicoreDir,multicoreSettings.exitFile,multicoreSettings.processId);
    fprintf('slave process %d exiting after %d function evaluations...\n',multicoreSettings.processId,nrRuns);
    
    exit = true;
    return;
end

end
