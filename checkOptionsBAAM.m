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

function checkOptionsBAAM(simInfo,adsInfo)
% Check if mole fraction of component A is empty or if it is above 1 or
% below 0
if isempty(simInfo.molFracFeed_A) || simInfo.molFracFeed_A < 0 || simInfo.molFracFeed_A > 1
    error('BAAM:checkOptionsBAAM:molefractionfeed',...
        'The variable molFracFeed_A is either empty or is below 0 or above 1. Check it!');
end

% Check if molarMass_A of component A is empty
if isempty(simInfo.molarMass_A)
    error('BAAM:checkOptionsBAAM:molarmass',...
        'The variable molarMass_A is empty. This is necessary for BAAM!');
end

% Check if temperature is empty
if isempty(simInfo.temperature)
    error('BAAM:checkOptionsBAAM:molarmass',...
        'The variable temperature is empty. This is necessary for BAAM!');
end

% Check if pressureHigh is empty
if isempty(simInfo.pressureHigh)
    error('BAAM:checkOptionsBAAM:highpressure',...
        'The variable pressureHigh is empty. This is necessary for BAAM!');
end

% Check if pressureLow is empty
if isempty(simInfo.pressureLow)
    error('BAAM:checkOptionsBAAM:lowpressure',...
        'The variable pressureLow is empty. This is necessary for BAAM!');
end

% Check if pressureHigh is greater than pressureLow
if simInfo.pressureHigh < simInfo.pressureLow
    error('BAAM:checkOptionsBAAM:highpressureless',...
        'The variable pressureLow should be smalled than pressureHigh. Check it!');
end

% Check if void fraction is empty or if it is above 1 or below 0
if isempty(simInfo.voidFraction) || simInfo.voidFraction < 0 || simInfo.voidFraction > 1
    error('BAAM:checkOptionsBAAM:voidFraction',...
        'The variable voidFraction is either empty or is below 0 or above 1. Check it!');
end

% Check if pressureHigh is greater than pressureLow
if simInfo.pressureHigh < simInfo.pressureLow
    error('BAAM:checkOptionsBAAM:highpressureless',...
        'The variable pressureLow should be smalled than pressureHigh. Check it!');
end

end