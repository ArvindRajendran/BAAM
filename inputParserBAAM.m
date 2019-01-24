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

function [simInfo,adsInfo] = inputParserBAAM(folder)
    % Determine folder (i.e., location of BAAMInfo)
    if nargin == 0
        folder = '.';
    end
    % Parse experiment info file
    fid = fopen([folder,filesep,'BAAMInfo'],'r');
    % Create file if it doesn't exist and ask the user to fill in the
    % fields an restart the code
    if fid <0
        fid = fopen([folder,filesep,'BAAMInfo'],'w');
        % Simulation info
        fprintf(fid,'%% I. SIMULATION INFO \n');
        fprintf(fid,'%% PROCESS INFO\n');
        fprintf(fid,'%% -> Feed Mole Fraction of Component A): molFracFeed_A = [0-1]\n');
        fprintf(fid,'molFracFeed_A = \n');
        fprintf(fid,'%% -> Molar Mass of Component A: molarMass_A = [in g/mol]\n');
        fprintf(fid,'molarMass_A = \n');
        fprintf(fid,'%% -> Feed Temperature: temperature = [in K]\n');
        fprintf(fid,'temperature = \n');
        fprintf(fid,'%% -> High Pressure: pressureHigh = [in bar]\n');
        fprintf(fid,'pressureHigh = \n');
        fprintf(fid,'%% -> Low Pressure: pressureLow = [in bar]\n');
        fprintf(fid,'pressureLow = \n');
        fprintf(fid,'%% -> Bed Void Fraction: voidFraction = [0-1]\n');
        fprintf(fid,'voidFraction = \n');
        fprintf(fid,'%% -> Adiabatic Constant of Gas: adiabaticConstant = [-]\n');
        fprintf(fid,'adiabaticConstant = \n');
        fprintf(fid,'%% -> Efficiency of Vacuum Pump: pumpEfficiency = [0-1]\n');
        fprintf(fid,'pumpEfficiency = \n');
        
        fprintf(fid,'%% PLOTTING INFO\n');
        fprintf(fid,'%% -> Boolean Flag after Simulation Run for Plotting (Yes/No): plotFlag = \n');
        fprintf(fid,'plotFlag = \n');
        fprintf(fid,'%% -> Mole Fractions used for plotting: molFracPlotting = \n');
        fprintf(fid,'%% -> This can be a series of mole fractions that would be used for plottting \n');
        fprintf(fid,'%% -> For example: 0.15, 0.5, 0.75, 0.9 or 0.15; 0.5; 0.75; 0.9 \n');
        fprintf(fid,'molFracPlotting = \n\n');
        
        fprintf(fid,'%% -> Load Adsorbent Info from a File (Yes/No): loadAdsorbentInfo = \n');
        fprintf(fid,'loadAdsorbentInfo = \n\n');
        fprintf(fid,'%% -> Filename to Load Adsorbent Info: filenameAdsorbentInfo = \n');
        fprintf(fid,'%% -> If not loaded from file, go down to enter adsorbentInfo = \n');
        fprintf(fid,'filenameAdsorbentInfo = \n\n');
        
        % Adsorbent info
        fprintf(fid,'%% II. ADSORBENT INFO \n');
        fprintf(fid,'%% *** User Input *** \n');
        fprintf(fid,'%% -> Adsorbent Name: adsorbentName = \n');
        fprintf(fid,'adsorbentName = \n');
        fprintf(fid,'%% -> Adsorbent Density: adsorbentDensity = [in kg/m3]\n');
        fprintf(fid,'adsorbentDensity = \n');        
        % Isotherm parameters
        fprintf(fid,'%% ISOTHERM PARAMETERS\n');
        % Component A
        fprintf(fid,'%% Component A\n');
        fprintf(fid,'%% -> Site 1, Saturation Capacity: qSaturationSite1_A = [in mol/kg]\n');
        fprintf(fid,'qSaturationSite1_A = \n');
        fprintf(fid,'%% -> Site 1, Adsorption Coefficient: adsorptionCoefficientSite1_A = [in m3/mol]\n');
        fprintf(fid,'adsorptionCoefficientSite1_A = \n');
        fprintf(fid,'%% -> Site 1, Internal Energy: internalEnergySite1_A = [in J/mol]\n');
        fprintf(fid,'internalEnergySite1_A = \n');
        fprintf(fid,'%% -> Site 2, Saturation Capacity: qSaturationSite2_A = [in mol/kg]\n');
        fprintf(fid,'qSaturationSite2_A = \n');
        fprintf(fid,'%% -> Site 2, Adsorption Coefficient: adsorptionCoefficientSite2_A = [in m3/mol]\n');
        fprintf(fid,'adsorptionCoefficientSite2_A = \n');
        fprintf(fid,'%% -> Site 2, Internal Energy: internalEnergySite2_A = [in J/mol]\n');
        fprintf(fid,'internalEnergySite2_A = \n\n');
        % Component B
        fprintf(fid,'%% Component B\n');
        fprintf(fid,'%% -> Site 1, Saturation Capacity: qSaturationSite1_B = [in mol/kg]\n');
        fprintf(fid,'qSaturationSite1_B = \n');
        fprintf(fid,'%% -> Site 1, Adsorption Coefficient: adsorptionCoefficientSite1_B = [in m3/mol]\n');
        fprintf(fid,'adsorptionCoefficientSite1_B = \n');
        fprintf(fid,'%% -> Site 1, Internal Energy: internalEnergySite1_B = [in J/mol]\n');
        fprintf(fid,'internalEnergySite1_B = \n');
        fprintf(fid,'%% -> Site 2, Saturation Capacity: qSaturationSite2_B = [in mol/kg]\n');
        fprintf(fid,'qSaturationSite2_B = \n');
        fprintf(fid,'%% -> Site 2, Adsorption Coefficient: adsorptionCoefficientSite2_B = [in m3/mol]\n');
        fprintf(fid,'adsorptionCoefficientSite2_B = \n');
        fprintf(fid,'%% -> Site 2, Internal Energy: internalEnergySite2_B = [in J/mol]\n');
        fprintf(fid,'internalEnergySite2_B = \n');
        
        % Close file and ask the user to fill in the required options
        fclose(fid);
        error('BAAM:BAAMInfoMissing',...
            'Could not find BAAMInfo file in the current folder. An empty file has been created. Please fill in the fields and rerun this script.');

    elseif fid>0 % Parse BAAMInfo file
        % Load file contents
        expInfo = textscan(fid,'%s','delimiter','\n');expInfo = expInfo{1};
        % Parse DB-related fields from BAAMInfo file
        for expLoop = 1:length(expInfo)
            locVal = strsplit(expInfo{expLoop},'=');
            if ~isempty(expInfo{expLoop})
                if ~strcmp(locVal{1,1}(1),'%')
                    % Extract values for all fields except the 'comments'
                    % field
                    locVar = locVal{1};
                    locVal = locVal{2};
                    tempVar = genvarname(locVar,who);
                    if strcmpi(tempVar,'molFracPlotting')
                        % minWavenumber - change the string to a vector of 1xn
                        tempmolFracPlotting = strsplit(locVal,{';',','}); 
                        tempmolFracPlotting = strrep(tempmolFracPlotting,' ','');
                        for i=1:size(tempmolFracPlotting,2)
                            storemolFracPlotting(1,i) = str2double(tempmolFracPlotting {i});
                        end
                        % Store minWavenumber into the BAAMInfo structure
                        fieldValues{expLoop} = storemolFracPlotting;
                    else
                        if ~isnan(str2double(locVal))
                            % Convert value to a double
                            fieldValues{expLoop} = str2double(strrep(locVal,' ',''));
                        else
                            % Leave value as a string
                            fieldValues{expLoop} = strrep(locVal,' ','');
                        end
                    end
                    if strcmpi(tempVar,'molFracFeed_A')...
                            || strcmpi(tempVar,'molarMass_A')...
                            || strcmpi(tempVar,'temperature')...
                            || strcmpi(tempVar,'pressureHigh')...
                            || strcmpi(tempVar,'pressureLow')...
                            || strcmpi(tempVar,'voidFraction')...
                            || strcmpi(tempVar,'adiabaticConstant')...
                            || strcmpi(tempVar,'pumpEfficiency')...
                            || strcmpi(tempVar,'plotFlag')...
                            || strcmpi(tempVar,'molFracPlotting')...
                            || strcmpi(tempVar,'loadAdsorbentInfo')...
                            || strcmpi(tempVar,'filenameAdsorbentInfo')
                        if strcmpi(tempVar,'loadAdsorbentInfo')
                            if strcmpi(fieldValues{expLoop},'yes') || strcmpi(fieldValues{expLoop},'Yes')
                                fieldValues{expLoop} = true;
                                adsInfoStructFlag = false;
                            else
                                fieldValues{expLoop} = false;
                                adsInfoStructFlag = true;
                            end
                        end
                        if strcmpi(tempVar,'plotFlag')
                            if strcmpi(fieldValues{expLoop},'yes') || strcmpi(fieldValues{expLoop},'Yes')
                                fieldValues{expLoop} = true;
                            else
                                fieldValues{expLoop} = false;
                            end
                        end
                        eval(['simInfo.',tempVar, '= fieldValues{expLoop};']);
                    else
                        if adsInfoStructFlag
                            eval(['adsInfo.',tempVar, '= fieldValues{expLoop};']);
                        else
                            adsInfo = [];
                        end
                    end
                end
            end
        end
        % Close file
        fclose(fid);
    end
end
