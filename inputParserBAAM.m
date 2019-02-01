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
% This function is used to parse user data from BAAMInfo file which would
% contain the simulation settings and adsorption properties entered by the
% user. In case the file is not present, the function would generate a new
% file called BAAMInfo with all the necessary fields. In case the user
% needs more info, he/she can refer to the README file. 
%
% Last modified:
% - 2019-02-01, AK: Minor cosmetic changes in the structure of the BAAMInfo
%                   file
% - 2019-01-25, AK: Minor cosmetic changes
% - 2019-01-24, AK: Introduced header for this file
%
% Input arguments:
% - folder:         Parent directory of the BAAMInfo file. If nothing
%                   specified then would look for the file in the current
%                   folder
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
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [simInfo,adsInfo] = inputParserBAAM(folder)
    simInfo = []; adsInfo = [];
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
        fprintf(fid,'%% BATCH ADSORBER ANALOGUE MODEL INFO \n');
        fprintf(fid,'%% Enter the simulation and adsorbent information here. \n');
        fprintf(fid,'%% For more information regarding the fields below check the README file. \n\n');
        % Simulation info
        fprintf(fid,'%% I. SIMULATION INFO \n');
        fprintf(fid,'molFracFeed_A = \n');
        fprintf(fid,'molarMass_A = \n');
        fprintf(fid,'temperature = \n');
        fprintf(fid,'pressureHigh = \n');
        fprintf(fid,'pressureLow = \n');
        fprintf(fid,'pressureLowUpperBound = \n');
        fprintf(fid,'voidFraction = \n');
        fprintf(fid,'adiabaticConstant = \n');
        fprintf(fid,'pumpEfficiency = \n');
        fprintf(fid,'plotFlag = \n');
        fprintf(fid,'molFracPlotting = \n\n');
        
        % Adsorbent info
        fprintf(fid,'%% II. ADSORBENT INFO \n');
        fprintf(fid,'loadAdsorbentInfo = \n');
        fprintf(fid,'filenameAdsorbentInfo = \n\n');
        fprintf(fid,'%% *** User Input. Fill out the fields below if not loading from file. *** \n');
        fprintf(fid,'adsorbentName = \n');
        fprintf(fid,'adsorbentDensity = \n');        
        % Isotherm parameters
        fprintf(fid,'%% ISOTHERM PARAMETERS\n');
        fprintf(fid,'%% Component A\n');
        fprintf(fid,'qSaturationSite1_A = \n');
        fprintf(fid,'adsorptionCoefficientSite1_A = \n');
        fprintf(fid,'internalEnergySite1_A = \n');
        fprintf(fid,'qSaturationSite2_A = \n');
        fprintf(fid,'adsorptionCoefficientSite2_A = \n');
        fprintf(fid,'internalEnergySite2_A = \n\n');
        % Component B
        fprintf(fid,'%% Component B\n');
        fprintf(fid,'qSaturationSite1_B = \n');
        fprintf(fid,'adsorptionCoefficientSite1_B = \n');
        fprintf(fid,'internalEnergySite1_B = \n');
        fprintf(fid,'qSaturationSite2_B = \n');
        fprintf(fid,'adsorptionCoefficientSite2_B = \n');
        fprintf(fid,'internalEnergySite2_B = ');
        
        % Close file and ask the user to fill in the required options
        fclose(fid);
        error('BAAM:BAAMInfoMissing',...
            'Could not find BAAMInfo file in the current folder. An empty file has been created. Please fill in the fields and rerun this script.');

    elseif fid>0 % Parse BAAMInfo file
        % Load file contents
        expInfo = textscan(fid,'%s','delimiter','\n');expInfo = expInfo{1};
        % Parse fields from BAAMInfo file
        for expLoop = 1:length(expInfo)
            locVal = strsplit(expInfo{expLoop},'=');
            if ~isempty(expInfo{expLoop})
                if ~strcmp(locVal{1,1}(1),'%')
                    % Extract values for all fields except the
                    % 'molFracPlotting'
                    locVar = locVal{1};
                    locVal = locVal{2};
                    tempVar = genvarname(locVar,who);
                    if strcmpi(tempVar,'molFracPlotting')
                        % molFracPlotting - change the string to a vector of 1xn
                        tempmolFracPlotting = strsplit(locVal,{';',','}); 
                        tempmolFracPlotting = strrep(tempmolFracPlotting,' ','');
                        for i=1:size(tempmolFracPlotting,2)
                            storemolFracPlotting(1,i) = str2double(tempmolFracPlotting {i});
                        end
                        % Store molFracPlotting into the BAAMInfo structure
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
                            || strcmpi(tempVar,'pressureLowUpperBound')...
                            || strcmpi(tempVar,'voidFraction')...
                            || strcmpi(tempVar,'adiabaticConstant')...
                            || strcmpi(tempVar,'pumpEfficiency')...
                            || strcmpi(tempVar,'plotFlag')...
                            || strcmpi(tempVar,'molFracPlotting')...
                            || strcmpi(tempVar,'loadAdsorbentInfo')...
                            || strcmpi(tempVar,'filenameAdsorbentInfo')
                        % Convert the user input into into boolean for
                        % loadAdsorbentInfo
                        if strcmpi(tempVar,'loadAdsorbentInfo')
                            if strcmpi(fieldValues{expLoop},'yes') || strcmpi(fieldValues{expLoop},'Yes')
                                fieldValues{expLoop} = true;
                                adsInfoStructFlag = false;
                            else
                                fieldValues{expLoop} = false;
                                adsInfoStructFlag = true;
                            end
                        end
                        % Convert the user input into into boolean for
                        % plotFlag
                        if strcmpi(tempVar,'plotFlag')
                            if strcmpi(fieldValues{expLoop},'yes') || strcmpi(fieldValues{expLoop},'Yes')
                                fieldValues{expLoop} = true;
                            else
                                fieldValues{expLoop} = false;
                            end
                        end
                        eval(['simInfo.',tempVar, '= fieldValues{expLoop};']);
                    else
                        % If adsInfo is from user input, save the fields
                        % into adsInfo structure
                        if adsInfoStructFlag
                            eval(['adsInfo.',tempVar, '= fieldValues{expLoop};']);
                        % If adsInfo is from a file, no matter what is
                        % loaded from the BAAMInfo file initialize the
                        % structure to an empty one
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
