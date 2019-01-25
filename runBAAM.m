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
% This is the main wrapper function of the main function BAAM.m that has
% the entire BAAM model. This function calls the other subsidiary functions
% which are necessary to run the model. First, the function calls
% inputParserBAAM to parse out the user data regarding the simulation 
% settings and adsorbent properties. After which the checkOptionsBAAM 
% function is called to perform a sanity check on all the user inputs.
% After which the BAAM model is called for as many adsorbents are requested
% to be simulated. Finally the results are stored in a folder called
% SimulationResults. Refer to the README file for further information
%
% Last modified:
% - 2019-01-25, AK: Added functionalities for loading adsorbent data from a
%                   file, and integrated sanity checks for user input and
%                   the model
% - 2019-01-24, AK: Introduced header for this file
%
% Input arguments:
% - Variable number of arguemnts as input. Maximum number of arguments that
% can be given is 1. This would be a boolean to specify whether the user
% requires detailed command line outputs for each adsorbent being
% simulated. If not this can be left out.
%
% Output arguments:
% - N/A
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runBAAM(varargin)
    simulationStartTime = datetime('now');
    % Get GIT commit ID and save it to the output Structure
    [status,cmdOut] = system('git rev-parse HEAD');
    % If command was successful
    if status == 0
        % Save the first 7 characters of the GIT commit ID
        commitId = cmdOut(1:7);
    end
    % Save the git commit id to the outputStruct
    outputStruct.commitId = commitId;
    
    % Parse information from the BAAMInfo file. If a file does not exist
    % then create a new BAAMInfo file and ask the user to update the
    % information
    [simInfo,adsInfo] = inputParserBAAM;
    
    % Check if the user wants the adsorbent information to come from a file
    % or if the user input the information in the BAAMInfo file
    if simInfo.loadAdsorbentInfo
        % Based on the type of file used, .mat or .xlsx, process the data
        % into the adsInfo structure such that it can be used by the model
        % function for simulation
        [adsInfo] = loadAdsorbentInfoFromFile(simInfo.filenameAdsorbentInfo);
        % Print a warning message saying the information for the adsorbent
        % properties will be loaded from the file
        warning('BAAM:runBAAM:adsorbentinfoload',...
            'All the adsorbent properties will be loaded from the file specified in BAAMInfo!');
    end
    
    % Get the number of adsorbents that would be simulated using BAAM
    numberOfAdsorbents = size(adsInfo,2);
    % Do not print detailed command line outputs if more than 1 adsorbent
    % is simulated. If user wants the output to be printed then print it
    % nevertheless
    if nargin<1
        if numberOfAdsorbents > 1
            silentFlag = true;
        else
            silentFlag = false;
        end
    elseif nargin==1
        silentFlag = varargin{1};
    else
        error('BAAM:runBAAM:numberofarguments',...
            'Too many input arguments to the BAAM wrapper function. Check it!');
    end
    
    % Loop over all the available adsorbents and run BAAM and obtain an
    % output structure and save the output for every individual adsorbent
    for ii = 1:numberOfAdsorbents
        % Some command line output
        if silentFlag
            if (ii==1)
                disp('##########################################################################');
                disp('                  BAAM - Batch Adsorber Analogue Model');
                disp('        Developed by Laboratory of Advanced Seperation Processes');
                disp('                     University Of Alberta, Canada');
                disp('##########################################################################');
                disp(' ');
                fprintf(' --> Number of adsorbents to be simulated using BAAM = %i\n',numberOfAdsorbents);
                disp('- Simulation of BAAM is in progress... ');
                disp(' ');
            end
            fprintf(' --> Adsorbent being simulated = %s\n',adsInfo(ii).adsorbentName);
            fprintf(' --> Number of adsorbents left to be simulated = %i\n',numberOfAdsorbents-ii);
        end
        % Based on the information provided by the user check if all the
        % user inputs are sensible. This is merely for sanity check. If any
        % necessary field is left empty, default values would be used. 
        [simInfoIn,adsInfoIn] = checkOptionsBAAM(simInfo,adsInfo(ii));

        % Execute BAAM after the sanity check and obtain the final output
        [BAAMOutput,~,~] = BAAM(simInfoIn,adsInfoIn,silentFlag);
        % Save the output structure for each adsorbent as ADS_<counterId>
        eval(['outputStruct.ADS_',num2str(ii),' = BAAMOutput;']);
    end
    % Evaluate the total computational time
    simulationStopTime = datetime('now');
    computationalTime = seconds(simulationStopTime-simulationStartTime);
    % Some command line output
    if silentFlag
        disp(' ');
        disp('- Simulation of BAAM for all adsorbents is completed...');
        disp('- Saving of BAAM results in progress...');
    end
    % Save the overall outputStruct into a .mat file. The file is named as
    % BAAMOutput_<git commit id.>_<time of simulation: ddmmyyyy_hhMMss>
    % Create a new folder called SimulationResults if it does not exist
    if exist('SimulationResults') ~= 7
        mkdir('SimulationResults');
    end
    % Save the file in the SimulationResults folder
    saveFileName = ['BAAMOutput_',commitId,'_',datestr(datetime('now'),'ddmmyyyy_hhMMss')];
    save(['SimulationResults',filesep,saveFileName],'outputStruct')
    % Some command line output
    if silentFlag
        fprintf(' --> Results saved in %s.mat\n',saveFileName);
        fprintf(' --> Computational time: %2.1f [s].',computationalTime); fprintf('\n');
        disp('##########################################################################');
    end
end