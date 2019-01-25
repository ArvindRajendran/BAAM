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
% Purpose:
% This function performs a sanity check on the inputs provided by the user 
% for the simulation info and the adsorbent info. If there is a serious
% issue this function throws out an error and stops the execution of the
% program. If it is something minor like a variable not entered, certain
% default value will be put into the missing entry and the model would be
% simulated. In such a case there would be a warning that would be thrown
% out by the function
%
% Last modified:
% - 2019-01-25, AK: Added more options that needs to be checked if they
%                   fall within the range or if they are empty
% - 2019-01-24, AK: Introduced header for this file
%
% Input arguments:
% - simInfo:       Structure containing information for the simulation.
%                  These generally do not depend on the adsorbent. This
%                  contains pressure values, temperature, column
%                  properties, etc.,
% - adsInfo:       Structure containing information about the adsorbent.
%                  This contains the adsorbent name, adsorbent density and
%                  DSL parameters for component A and B. In the original
%                  manuscript A is CO2 and B is N2.
%
% Output arguments:
% - simInfo:       Structure containing information for the simulation.
%                  These generally do not depend on the adsorbent. This
%                  contains pressure values, temperature, column
%                  properties, etc.,
% - adsInfo:       Structure containing information about the adsorbent.
%                  This contains the adsorbent name, adsorbent density and
%                  DSL parameters for component A and B. In the original
%                  manuscript A is CO2 and B is N2.
% Note that the output arguments can either be the same as the input
% arguments or could have been altered by this function based on certain
% checks that are put in place. 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [simInfo,adsInfo] = checkOptionsBAAM(simInfo,adsInfo)
% Check if flag for loading adsorbent is empty. If it is empty, then user
% input/default will be used for simulations 
if simInfo.loadAdsorbentInfo == 0
    simInfo.filenameAdsorbentInfo = [];
    warning('BAAM:checkOptionsBAAM:loadadsorbentinfo',...
            'The variable loadAdsorbentInfo is empty. Will load user data or default settings to run the BAAM model!');
end

% Check if mole fraction of component A is empty or if it is above 1 or
% below 0
if isempty(simInfo.molFracFeed_A)
    simInfo.molFracFeed_A = 0.15;
    warning('BAAM:checkOptionsBAAM:molefractionfeed',...
        'The variable molFracFeed_A is empty. Mole fraction of 0.15 has been set for component A in feed!');
end
if simInfo.molFracFeed_A < 0 || simInfo.molFracFeed_A > 1
    error('BAAM:checkOptionsBAAM:molefractionfeed',...
        'The variable voidFraction is either below 0 or above 1. Check it!');
end
    
% Check if molarMass_A of component A is empty
if isempty(simInfo.molarMass_A)
    simInfo.molarMass_A = 44.01;
    warning('BAAM:checkOptionsBAAM:molarmass',...
        'The variable molarMass_A is empty. Molar mass of component A has been set to 44.01 (CO2)!');
end

% Check if temperature is empty
if isempty(simInfo.temperature)
    simInfo.temperature = 298.15;
    warning('BAAM:checkOptionsBAAM:molarmass',...
        'The variable temperature is empty. The feed temperature has been set to 298.15 K!');
end

% Check if pressureHigh is empty
if isempty(simInfo.pressureHigh)
    simInfo.pressureHigh = 1;
    warning('BAAM:checkOptionsBAAM:highpressure',...
        'The variable pressureHigh is empty. The high pressure has been set to 1 bar!');
end

% Check if pressureLow is empty
if isempty(simInfo.pressureLow)
    simInfo.pressureLow = 0.03;
    warning('BAAM:checkOptionsBAAM:lowpressure',...
        'The variable pressureLow is empty. The low pressure has been set to 0.03 bar!');
end

% Check if pressureHigh is greater than pressureLow
if simInfo.pressureHigh < simInfo.pressureLow
    error('BAAM:checkOptionsBAAM:highpressureless',...
        'The variable pressureLow should be smaller than pressureHigh. Check it!');
end

% Check if pressureLowUpperBound is empty
if isempty(simInfo.pressureLowUpperBound)
    simInfo.pressureLowUpperBound = 0.10;
    warning('BAAM:checkOptionsBAAM:lowpressureupperbound',...
        'The variable pressureLowUpperBound is empty. The low pressure upper bound for BAAM has been set to 0.10 bar!');
end
if simInfo.pressureLowUpperBound < simInfo.pressureLow ...
        || simInfo.pressureLowUpperBound > simInfo.pressureHigh
   error('BAAM:checkOptionsBAAM:lowpressureupperbound',...
        'The variable pressureLowUpperBound is either below low pressure or above high pressure. Check it!');
end

% Check if void fraction is empty or if it is above 1 or below 0
if isempty(simInfo.voidFraction)
    simInfo.voidFraction = 0.37;
    warning('BAAM:checkOptionsBAAM:voidfraction',...
        'The variable voidFraction is empty. The void fraction has been set to 0.37!');
end
if simInfo.voidFraction < 0 || simInfo.voidFraction > 1
   error('BAAM:checkOptionsBAAM:voidfraction',...
        'The variable voidFraction is either below 0 or above 1. Check it!');
end

% Check if adiabatic constant of the gas is empty or if it is below 1
if isempty(simInfo.adiabaticConstant)
    simInfo.adiabaticConstant = 1.4;
    warning('BAAM:checkOptionsBAAM:adiabaticconstant',...
        'The variable adiabaticConstant is empty. The adiabatic constant has been set to 1.4!');
end
if simInfo.adiabaticConstant < 1
   error('BAAM:checkOptionsBAAM:adiabaticconstant',...
        'The variable voidFraction is either below 1. For an ideal gas it should be above 1. Check it!');
end

% Check if the pump efficiency is empty or if it is above 1 or below 0
if isempty(simInfo.pumpEfficiency)
    simInfo.pumpEfficiency = 0.72;
    warning('BAAM:checkOptionsBAAM:pumpefficiency',...
        'The variable pumpEfficiency is empty. The vacuum pump efficiency has been set to 0.72!');
end
if simInfo.pumpEfficiency < 0 || simInfo.pumpEfficiency > 1
   error('BAAM:checkOptionsBAAM:pumpefficiency',...
        'The variable pumpEfficiency is either below 0 or above 1. Check it!');
end

% Check if the flag for plotting is empty. If the adsorbent info is input
% by the user then plot the figure at the end of the simulation if not turn
% off the plotting functionality.
if isempty(simInfo.plotFlag)
    if simInfo.loadAdsorbentInfo
        simInfo.plotFlag = false;
        warning('BAAM:checkOptionsBAAM:plotflag',...
            'The variable plotFlag is empty. Plotting has been turned off!');
    else
        simInfo.plotFlag = true;
        warning('BAAM:checkOptionsBAAM:plotflag',...
            'The variable plotFlag is empty. Plotting has been turned off!');
    end
end

% Check if the mole fraction variable for plotting is empty.
if isnan(simInfo.molFracPlotting)
    simInfo.molFracPlotting = [0.15,0.30,0.50,0.75,0.9,1.0];
    warning('BAAM:checkOptionsBAAM:molfracplotting',...
            'The variable molFracPlotting is empty. The final plot would be generated at 6 different default mole fractions!');
end

% Check if the isotherm parameters for component A and site 1 is empty
if isempty(adsInfo.qSaturationSite1_A) || isempty(adsInfo.adsorptionCoefficientSite1_A)...
        || isempty(adsInfo.internalEnergySite1_A)
    error('BAAM:checkOptionsBAAM:isothermcompAsite1',...
        'Some or all the isotherm parameters for component A and site 1 is empty. Check it!');
end

% Check if the isotherm parameters for component B and site 1 is empty
if isempty(adsInfo.qSaturationSite1_B) || isempty(adsInfo.adsorptionCoefficientSite1_B)...
        || isempty(adsInfo.internalEnergySite1_B)
    error('BAAM:checkOptionsBAAM:isothermcompBsite1',...
        'Some or all the isotherm parameters for component B and site 1 is empty. Check it!');
end

% Check if the isotherm parameters for component A and site 2 is empty or
% zero. Set the default value to 0 for all the parameters in this case.
if isempty(adsInfo.qSaturationSite2_A) || isempty(adsInfo.adsorptionCoefficientSite2_A)...
        || isempty(adsInfo.internalEnergySite2_A) || adsInfo.qSaturationSite2_A == 0 ...
        || adsInfo.adsorptionCoefficientSite2_A == 0|| adsInfo.internalEnergySite2_A == 0
    warning('BAAM:checkOptionsBAAM:isothermcompAsite2',...
        'Some or all the isotherm parameters for component A and site 2 is empty or zero. All the parameters set to zero!');
    adsInfo.qSaturationSite2_A = 0;
    adsInfo.adsorptionCoefficientSite2_A = 0;
    adsInfo.internalEnergySite2_A = 0;
end

% Check if the isotherm parameters for component B and site 2 is empty or
% zero. Set the default value to 0 for all the parameters in this case.
if isempty(adsInfo.qSaturationSite2_B) || isempty(adsInfo.adsorptionCoefficientSite2_B)...
        || isempty(adsInfo.internalEnergySite2_B) || adsInfo.qSaturationSite2_B == 0 ...
        || adsInfo.adsorptionCoefficientSite2_B == 0|| adsInfo.internalEnergySite2_B == 0
    warning('BAAM:checkOptionsBAAM:isothermcompBsite2',...
        'Some or all the isotherm parameters for component B and site 2 is empty or zero. All the parameters set to zero!');
    adsInfo.qSaturationSite2_B = 0;
    adsInfo.adsorptionCoefficientSite2_B = 0;
    adsInfo.internalEnergySite2_B = 0;
end
end