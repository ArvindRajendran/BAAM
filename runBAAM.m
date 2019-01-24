%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% University of Alberta, Edmonton, Canada
% Laboratory of Advanced Separation Processes (Prof. Dr. Arvind Rajendran)
%
% Project:  ........
% Year:     2019
% MATLAB:   R2018b, Windows 64bit
% Authors:  Vishal Subramanian Balashankar (VS)
%           Ashwin Kumar Rajagopalan (AK)
%          
%
% Purpose:
% ...
%
% Last modified:
% - 2019-01-24, AK: Introduced header for this file
%
% Input arguments:
% - 
%
% Output arguments:
% - 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runBAAM
    % Load information about the simulation and the adsorbent
    [simInfo,adsInfo] = inputParserBAAM;
    
    % Check if the options that are given in the input file is reasonable
    checkOptionsBAAM(simInfo,adsInfo);
    
    % Run the BAAM model
    runBAAM(adsInfo,simInfo);
end
