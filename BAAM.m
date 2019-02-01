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
% The function BAAM simulates the batch adsorber analogue model proposed in
% <Title of the paper, doi>. This code simulates a two component system
% using the LPP cycle with the BAAM model. Component A is the component of
% interest (CO2 in the original manuscript). This function generates two
% outputs, namely a .mat file with the entire simulation result and a
% simple .txt file called BAAM indicating whether the adsorbent satisfied
% the constraints on purity/recovery and energy/working capacity
% information.

% For more information on how to run this function and to understand the
% output from this function, please read the README.md file in the GIT 
% repository.
%
% Last modified:
% - 2019-02-01, AK: Minor structural modifications
% - 2019-01-25, AK: Minor cosmetic changes
% - 2019-01-13, AK: Added additional comment for discretization
% - 2019-01-13, VS: Cosmetic changes incorporated
% - 2019-01-12, AK: Finished major cleaning of the code and output
%                   structure and first functioning verion of the code.
% - 2019-01-09, AK: Cleaned up the adsorption part of the code and the
%                   performance indicators
% - 2019-01-08, AK: Cleaned up the pressurization part of the code
% - 2019-01-07, AK: Cleaned up the variable names, added a new input
%                   structure and cleaned up the blowdown/evacuation part
%                   of the code
% - 2019-01-06, AK: Introduced header for this file
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
% - silentFlag:    Flag for turning off command window output.
%
% Output arguments:
% - BAAMOutput:    Structure containing output from the BAAM model. This
%                  contains performance indicators, and output from each
%                  step that was simulated.
%
% Comments:
% - Hardcoded values need to be changed in case if the pressure ranges are
% different (I have changed it from plow to phigh to be universal can be
% changed back again. Search for harcoded and look at the original code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BAAMOutput]=BAAM(simInfo,adsInfo,silentFlag)
% START OF MAIN FUNCTION
simulationStartTime = datetime('now');

%% COMMAND OUTPUT WINDOW
if ~silentFlag
    disp(' ');
    disp('##########################################################################');
    disp('                  BAAM - Batch Adsorber Analogue Model');
    disp('        Developed by Laboratory of Advanced Seperation Processes');
    disp('                     University Of Alberta, Canada');
    disp('##########################################################################');
    disp(' ');
end

%% STARTUP FORMALITIES OF THE CODE
% ADSORBENT PROPERTIES
adsName = adsInfo.adsorbentName; % Name of the adsorbent being simulated
adsorbentMass = 1; % Mass of the adsorbent [kg]
adsorbentDensity = adsInfo.adsorbentDensity; % Particle Density [kg/m3]

% ADDITIONAL PROPERTIES
voidFraction = simInfo.voidFraction; % Void fraction [-]
columnVolume = adsorbentMass...
    /(adsorbentDensity*(1-voidFraction)); % Volume of the column [m3]
adiabaticConstant = simInfo.adiabaticConstant; % Adiabatic constant [-]
pumpEfficiency = simInfo.pumpEfficiency; % Efficiency of the vacuum pump [-]

% OPERATING CONDITIONS
temperature = simInfo.temperature; % Feed temperature [K]
pressureHigh = simInfo.pressureHigh; % High pressure [bar]
pressureLow = simInfo.pressureLow; % Low pressure [bar]
pressureLowUpperBound = simInfo.pressureLowUpperBound; % Low pressure upper bound for simulation [bar]
molFracFeed_A = simInfo.molFracFeed_A; % Feed mole fraction of component A [-]
molarMass_A = simInfo.molarMass_A; % Molar mass of component A [g/mol]

% UNIVERSAL CONSTANTS
UNIVERSAL_GAS_CONSTANT = 8.314e-5; % Universal gas constant [(m3 bar)/(K mol)]
UNIVERSAL_GAS_CONSTANT_JMK = 8.314; % Universal gas constant [J/mol K]
JOULE_TO_KWH = 2.77778e-7; % Conversion from Joules to kWh

% ISOTHERM PROPERTIES
% Dual-site Langmuir model for a binary system is used in this version of
% BAAM
% Component A
qSaturationSite1_A = adsInfo.qSaturationSite1_A; % Saturation capacity of site 1 [mol/kg]
adsorptionCoefficientSite1_A = adsInfo.adsorptionCoefficientSite1_A; % Adsorption coeffecient of site 1 [m3/mol]
internalEnergySite1_A = adsInfo.internalEnergySite1_A*1000; % Internal energy of site 1 [J/mol]
qSaturationSite2_A = adsInfo.qSaturationSite2_A; % Saturation capacity of site 2 [mol/kg]
adsorptionCoefficientSite2_A = adsInfo.adsorptionCoefficientSite2_A; % Adsorption coeffecient of site 2 [m3/mol]
internalEnergySite2_A = adsInfo.internalEnergySite2_A*1000; % Internal energy of site 2 [J/mol]
adsorptionEqbmConstSite1_A = adsorptionCoefficientSite1_A...
    *exp(-internalEnergySite1_A/(UNIVERSAL_GAS_CONSTANT_JMK*temperature)); % Adsorption equilibrium constant for site 1 [m3/mol]
adsorptionEqbmConstSite2_A = adsorptionCoefficientSite2_A...
    *exp(-internalEnergySite2_A/(UNIVERSAL_GAS_CONSTANT_JMK*temperature)); % Adsorption equilibrium constant for site 2 [m3/mol]

% Component B
qSaturationSite1_B = adsInfo.qSaturationSite1_B; % Saturation capacity of site 1 [mol/kg]
adsorptionCoefficientSite1_B = adsInfo.adsorptionCoefficientSite1_B; % Adsorption coeffecient of site 1 [m3/mol]
internalEnergySite1_B = adsInfo.internalEnergySite1_B*1000; % Internal energy of site 1 [J/mol]
qSaturationSite2_B = adsInfo.qSaturationSite2_B; % Saturation capacity of site 2 [mol/kg]
adsorptionCoefficientSite2_B = adsInfo.adsorptionCoefficientSite2_B; % Adsorption coeffecient of site 2 [m3/mol]
internalEnergySite2_B = adsInfo.internalEnergySite2_B*1000; % Internal energy of site 2 [J/mol]
adsorptionEqbmConstSite1_B = adsorptionCoefficientSite1_B...
    *exp(-internalEnergySite1_B/(UNIVERSAL_GAS_CONSTANT_JMK*temperature));  % Adsorption equilibrium constant for site 1 [m3/mol]
adsorptionEqbmConstSite2_B = adsorptionCoefficientSite2_B...
    *exp(-internalEnergySite2_B/(UNIVERSAL_GAS_CONSTANT_JMK*temperature));  % Adsorption equilibrium constant for site 1 [m3/mol]

% PLOTTING OPTIONS
% Mole fraction vector to obtain the solid phase loading required for
% plotting the plots similar to Fig. 3 in the original manuscript.
molFracPlotting = simInfo.molFracPlotting;

%% SOLVER DEFINITIONS/PRESSURE MATRIX
% ODE SOLVER PROPERTIES
% Have a discretization step size small enough such that the results are
% not influenced by the step size (0.001 seems to be fair value to use). 
deltaPressure = 0.001; % Length of discretized grid for pressure [bar]
% Generate a pressure vector that spans from the high pressure to the low
% pressure with a length of deltaPressure for the ode solver to be used in
% the blowdown step
pressureVector=(pressureHigh:-deltaPressure:pressureLow)';

% PRESSURE MATRIX
% Create a pressure vector for the low pressure that spans a range from
% low pressure (pressureLow) to upper bound of low pressure specifified
% by the user. The pressure vector is in bar.
pressureLowVector = pressureLow:0.001:(pressureLowUpperBound - 0.01);
% Create a pressure vector for the intermediate pressure that spans a range
% from low pressure (pressureLow+0.01 bar) to the high pressure
% (pressureHigh - 0.01 bar).
pressureInterVector = (pressureLow + 0.01):0.01:(pressureHigh - 0.01);
% Generate a matrix with all combinations of low and intermediate pressures
[pressureLowGrid, pressureInterGrid] = meshgrid(pressureLowVector,pressureInterVector);

% Loop over all the entries of the pressure matrix generated above and
% remove entires that have intermediate pressure lower than the low
% pressure
[numRows, numColumns]=size(pressureLowGrid);
% Loop over all the rows
for ii = 1:numRows
    % Loop over all the columns
    for jj = 1:numColumns
        % Check for entries that have pressureInter <= pressureLow
        if pressureInterGrid(ii,jj)<pressureLowGrid(ii,jj)...
                || eq(pressureInterGrid(ii,jj),pressureLowGrid(ii,jj))
            % Reinitialize such entries to NaN
            pressureInterGrid(ii,jj)=NaN;
        end
    end
end

% Initialize the performance indicators with a NaN matrix with 
% length(pressureInterVector) rows and length(pressureLowVector) columns
pressure.low = nan(length(pressureInterVector),length(pressureLowVector)); % Low pressure
pressure.inter = nan(length(pressureInterVector),length(pressureLowVector)); % Intermediate pressure
purity = nan(length(pressureInterVector),length(pressureLowVector)); % Purity
recovery = nan(length(pressureInterVector),length(pressureLowVector)); % Recovery
energyConsumption = nan(length(pressureInterVector),length(pressureLowVector)); % Energy consumption
workingCapacity = nan(length(pressureInterVector),length(pressureLowVector)); % Working capacity
massBalanceError = nan(length(pressureInterVector),length(pressureLowVector)); % Mass balance error
distancePurityRecovery = nan(length(pressureInterVector),length(pressureLowVector)); % Purity/Recovery distance

%% COMMAND OUTPUT WINDOW
if ~silentFlag
    fprintf(' --> Adsorbent = %s\n',adsName);
    fprintf(' --> Temperature = %0.2f [K]\n',temperature);
    fprintf(' --> High Pressure (PH) = %0.2f [bar]\n',pressureHigh);
    fprintf(' --> Low Pressure (PL) = %0.2f [bar]\n',pressureLow);
    fprintf(' --> Feed Mole Fraction (A) = %0.2f [-]\n',molFracFeed_A);
    disp('- Simulation of BAAM is in progress... ');
end

%% RUN THE MODEL
% Simulate the blowdown and the evacuation step using the function
% simulateBlowdownEvacuation and obtain the pressure, molefraction of 
% component A, solid loadings of A and B, total moles desorbed, moles of A 
% and B desorbed and the energy consumption
% Matrix structure: [Pressure MoleFracA SolidLoadingA SolidLoadingB
% MolesDesorbedTotal MolesDesorbedA MolesDesorbedB EnergyConsumption]
outputBlowEvac = simulateBlowdownEvacuation();

%%% TO VISHAL: Why do you even need the matrix from meshgrid? That thing
%%% would slow the code down

% Simulate the pressurization and adsorption steps. Loop over all the low
% and intermediate pressures. 
% Loop over entire low pressure vector
for counterLow = 1:length(pressureLowVector)
    % Find the index (in the output of the simulateBlowdownEvacuation
    % function) of the low pressure that corresponds to the current low
    % pressure value in the loop. Have a tolerance of 10e-6 just in case.
    pressureLowIndex = find(abs(outputBlowEvac(:,1)-pressureLowGrid(1,counterLow))<1e-6);
    
    % Simulate the pressurization step starting from the low pressure with
    % the goal of reaching the high pressure. The input for this function
    % would be the low pressure and mole fraction at the end of the
    % blowdown/evacuation step (outputBlowEvac(:,2) is the mole fraction)
    % outputPress - [MoleFracA SolidLoadingA SolidLoadingB MolesPressTotal]
    outputPress(counterLow,:) = simulatePressurization(pressureLowVector(counterLow),outputBlowEvac(pressureLowIndex,2));

    % Simulate the adsorption step starting from the final state of the
    % pressurization step with the goal of feeding the column with gas at
    % feed conditions to breakthrough the column. The input for this
    % function would be the final state of the pressurization step
    % (outputPress(1))
    % outputAds - [molesAds_Total molesAds_Raffinate SolidLoadingA 
    % SolidLoadingB]
    outputAds(counterLow,:) = simulateAdsorption(outputPress(counterLow,1));

%     % If the current step pressure is equal to the low pressure then save
%     % the solid phase loadings at that condition which would be the final
%     % state of the colum at the end of the pressurization step and save it
%     % to the outputPress structure
%     if (pressureLowVector(counterLow) == pressureLow)
%         solidPhaseLoading_FinalLPP_A = outputPress(counterLow,2);
%         solidPhaseLoading_FinalLPP_B = outputPress(counterLow,3);
%     end
    
    % Loop over entire intermediate pressure vector
    for counterInter = 1:length(pressureInterVector)
        %%% TO VISHAL: Why would the low pressure even be NaN???
        if ~(isnan(pressureInterGrid(counterInter,counterLow)) ...
                || isnan(pressureLowGrid(1,counterLow)))
            % Find the index (in the output of the 
            % simulateBlowdownEvacuation function) of the pressure that 
            % corresponds to the current intermediate pressure value in the
            % loop given the low pressure. Have a tolerance of 10e-6 just 
            % in case.
            pressureInterIndex = find(abs(outputBlowEvac(:,1)-pressureInterGrid(counterInter,counterLow))<1e-6);
            % Get the number of moles desorbed till from the high pressure 
            % to the intermediate pressure (from outputBlowEvac)
            % Component A
            molesHighToInter_A = outputBlowEvac(pressureInterIndex,6);
            % Component B
            molesHighToInter_B = outputBlowEvac(pressureInterIndex,7);
            % Evaluate the number of moles desorbed from the intermediate
            % pressure to the low pressure for component (from
            % outputBlowEvac). This would be the difference of moles
            % adorbed at the low pressure and intermediate pressure
            % Component A
            molesInterToLow_A = outputBlowEvac(pressureLowIndex,6)-outputBlowEvac(pressureInterIndex,6);
            % Component B
            molesInterToLow_B = outputBlowEvac(pressureLowIndex,7)-outputBlowEvac(pressureInterIndex,7);
            
            % Calculate the performance indicators at the given low and
            % intermediate pressure
            % Purity [%] - Moles of component A desorbed from intermediate 
            % to low pressure/Total moles desorbed from intermediate to low
            % pressure
            purityStep = (molesInterToLow_A/(molesInterToLow_A + molesInterToLow_B))*100;
            % Recovery [%] - Moles of component A desorbed from 
            % intermediate to low pressure/Moles fed into the column in the
            % adsorption step of component A (outputAds(1)*molFracFeed_A): 
            % outputAds(1) - molesAds_Total
            recoveryStep = (molesInterToLow_A/(outputAds(counterLow,1)*molFracFeed_A))*100;
            % Energy consumption (Blowdown) [kWh/tonne CO2. evac] - Energy
            % consumption (kWh) in the blowdown step per tonne of component
            % A obtained in the evacuation step (Eq. 17 in the original
            % manuscript)
            % (energy consumption in simulateBlowdownEvacuation is given 
            % in J, and 1 kWh = 2.77778e-7 J, and the molar mass is given 
            % in g/mol and 1e-6 is for g -> tonne)
            energyBlowdownStep = (outputBlowEvac(pressureInterIndex,8)*JOULE_TO_KWH)...
                                        /(molesInterToLow_A*molarMass_A*1e-6);
            % Energy consumption (Evacuation) [kWh/tonne CO2. evac] - 
            % Energy consumption (kWh) in the evacuation step per tonne of 
            % component A obtained in the evacuation step. The energy from
            % the evacuation step is computed as the difference between the
            % energy at low pressure and energy at intermediate pressure
            % (energy consumption in simulateBlowdownEvacuation is given 
            % in J, and 1 kWh = 2.77778e-7 J, and the molar mass is given 
            % in g/mol and 1e-6 is for g -> tonne). (Eq. 17 in the original
            % manuscript)
            energyEvacuationStep = ((outputBlowEvac(pressureLowIndex,8)...
                                    - outputBlowEvac(pressureInterIndex,8))*JOULE_TO_KWH)...
                                    /(molesInterToLow_A*molarMass_A*1e-6);
            % Total Energy consumption [kWh/tonne CO2. evac]
            energyTotalStep = energyBlowdownStep+energyEvacuationStep;
            % Working capacity [mol CO2/m3 of adsorbent] - defined as the
            % number of moles of CO2 desorbed in the evaculation for a
            % m3 of adsorbent
            workingCapacityStep = (molesInterToLow_A)/(columnVolume*(1-voidFraction));
            % Mass balance error for the entire cycle [%]. Moles of gas fed
            % in is through the adsorption step, while the moles of gas out
            % of the cycle is from the blowdown/evacuation step and the
            % difference in the moles of gas fed to the pressurization step
            % using the adsorption outlet.
            % Moles of gas in
            numMolesIn = outputAds(counterLow,1);
            % Moles of gas out
            numMolesOut = (outputAds(counterLow,2)-outputPress(counterLow,4))...
                            + (molesHighToInter_A + molesHighToInter_B)...
                            + (molesInterToLow_A + molesInterToLow_B);
            % Mass balance error
            massBalanceErrorStep = ((numMolesIn - numMolesOut)/numMolesOut)*100;
            % Eucledian distance of purity/recovery from origin (Eq. 22)
            % which would be used to determine if the material can satisfy
            % 95/90 purity recovery which comes from training the
            % simplified model w.r.t the full scale model
            if ~isempty(purityStep) && ~isempty(recoveryStep)
                distancePurityRecoveryStep = sqrt(purityStep^2 + recoveryStep^2);
            else
                % If the purity or recovery is empty, initialize all the
                % performance indicators and elements of the output
                % structure to NaN.
                purityStep = NaN;
                recoveryStep = NaN;
                energyTotalStep = NaN;
                workingCapacityStep = NaN;
                massBalanceErrorStep = NaN;
                distancePurityRecoveryStep = NaN;
            end
        
            % Prepare output structure with pressure, performance 
            % indicators, the mass balance error and the distance of the
            % purity/recovery curve from the origin
            pressure.low(counterInter,counterLow) = pressureLowGrid(counterInter,counterLow); % Low pressure
            pressure.inter(counterInter,counterLow) = pressureInterGrid(counterInter,counterLow); % High pressure
            purity(counterInter,counterLow) = purityStep; % Purity (Eq. 15)
            recovery(counterInter,counterLow) = recoveryStep; % Recovery (Eq. 16)
            energyConsumption(counterInter,counterLow) = energyTotalStep; % Energy consumption (Eq. 17)
            workingCapacity(counterInter,counterLow) = workingCapacityStep; % Working capacity (Eq. 18)
            massBalanceError(counterInter,counterLow) = massBalanceErrorStep; % Mass balance error
            distancePurityRecovery(counterInter,counterLow) = distancePurityRecoveryStep; % Purity/Recovery distance
        end % End of performance indicator condition
    end % End of intermediate pressure loop
end % End of low pressure loop

% Determine if the adsorbent being simulated satisfies the 95% purity and
% 90% recovery constraint based on the model created using the MATLAB
% Statistics and Machine Learning Toolbox (not shown here, but refer to
% Section 5.1 in the original manuscript)
% Find the indices that closely corresponds to the locus of the 95/90 curve
% described in Section 5.1, which would aid in computing the minimum energy
% described in Section 5.2. Even though from the model it is shown that the
% locus would correspond to a distance of 110.25, points that fall between
% 110.2 and 110.3 are taken here.
indexDistance9590 = find(distancePurityRecovery > 110.2...
                            & distancePurityRecovery < 110.3);
% Find the maximum distance from the origin to the purity-recovery pareto
% for a given adsorption
maxDistancePurityRecovery = max(max(distancePurityRecovery));

% Compute the minimum energy if the indexDistance9590 is not empty. Based
% on that compute the minimum energy and the working capacity that would
% correspond to the 95/90% purity/recovery. The minimum energy obtained is
% from BAAM as well as from the full scale model, which in this case is 
% obtained from Eq. 25 (see Section 5.2). Store the results thus obtained 
% from BAAM into a .txt file for future use.
if ~isempty(indexDistance9590)
    % Boolean to save in the output to say whether purity/recovery
    % constraints have been met or not.
    satisfyPurityRecovery = 1;
    % Get the minimum energy consumption from the simplified model and the
    % corresponding index in the energyConsumption matrix.
    [minEnergyBAAM, indexMinimumEnergyBAAM]...
        = min(energyConsumption(indexDistance9590));
    % Based on Eq. 25, compute the energy that could be obtained from the
    % simplified model, but using the energy consumption from the BAAM 
    minEnergyFullModel = 1.1446*minEnergyBAAM + 66.528;
    % Get the working capacity corresponding to the minimum energy
    % consumption from BAAM
    workingCapacity_MinEnergy ...
        = workingCapacity(indexDistance9590(indexMinimumEnergyBAAM));
        
    % Create a new file titled BAAM.txt if it does not exist in the current
    % folder and create the header for the file.
    if (exist('BAAM.txt','file') ~= 2)
        fid = fopen('BAAM.txt', 'at');
        fprintf(fid,'%s %s %s %s %s %s %s\n','Adsorbent','Temperature','MaximumDistance_Pu/Re',...
            'Satisfy_Pu/Re','MinimumEnergy_BAAM', 'MinimumEnergy_FullModel','WorkingCapacity_BAAM_MinEnergy');
        fclose(fid);
    end
    % Save the name of the adsorbent, temperature, distance, whether the
    % purity/recovery constraints have been met or not, minimum energy from
    % BAAM, full model and the working capacity at minimum energy. 
    fid = fopen('BAAM.txt', 'at');
    fprintf(fid,'%s %f %f %f %f %f %f\n',adsName,temperature,maxDistancePurityRecovery,...
        satisfyPurityRecovery,minEnergyBAAM,minEnergyFullModel,workingCapacity_MinEnergy);
    fclose(fid);
else
    % Boolean to save in the output to say whether purity/recovery
    % constraints have been met or not. If it comes into this loop it means
    % that there is no point which has a distance between 110.2 and 110.3
    % in the purity/recovery curve from the origin.
    satisfyPurityRecovery = 0;
    % Create a new file titled BAAM.txt if it does not exist in the current
    % folder and create the header for the file.
    if (exist('BAAM.txt','file') ~= 2)
        fid = fopen('BAAM.txt', 'at');
        fprintf(fid,'%s %s %s %s %s %s %s\n','Adsorbent','Temperature','MaximumDistance_Pu/Re',...
            'Satisfy_Pu/Re','MinimumEnergy_BAAM', 'MinimumEnergy_FullModel','WorkingCapacity_BAAM_MinEnergy');
        fclose(fid);
    end
    % Save the name of the adsorbent, temperature, distance, whether the
    % purity/recovery constraints have been met or not, minimum energy from
    % BAAM, full model and the working capacity at minimum energy. In this
    % case minimum energy from BAAM, full model and the working capacity at
    % minimum energy does not exist.
    fid = fopen('BAAM.txt', 'at');
	fprintf(fid,'%s %f %f %f %s %s %s\n',adsName,temperature,maxDistancePurityRecovery,...
        satisfyPurityRecovery,'-','-','-');
    fclose(fid);
end

% Generate output which is the solid phase loadings at different mole
% fractions for component A and B which would be later used for plotting
% figures like the one shown in Fig. 3 in the original manuscript
pressureForPlotting = 0:0.01:pressureHigh;
% Calculate the solid phase loadings for all the mole fractions given as an
% input by the user
plotting = struct([]);
for ii = 1:length(molFracPlotting)
    plotting(ii).pressure = pressureForPlotting;
    plotting(ii).moleFraction = molFracPlotting(ii);
    [plotting(ii).solidPhaseLoading_A, plotting(ii).solidPhaseLoading_B]...
        = evaluateDSLIsotherm (pressureForPlotting,molFracPlotting(ii));
end

% Prepare the output structure
BAAMOutput.adsInfo = adsInfo; % Adsorbent Info structure
BAAMOutput.simInfo = simInfo; % Simulation Info structure
BAAMOutput.outputBlowEvac = outputBlowEvac; % Output from simulateBlowEvac
BAAMOutput.outputPress = outputPress; % Output from simulatePressurization
BAAMOutput.outputAds = outputAds; % Output from simulateAdsorption
BAAMOutput.pressure = pressure; % Pressure for model
BAAMOutput.purity = purity; % Purity
BAAMOutput.recovery = recovery; % Recovery
BAAMOutput.energyConsumption = energyConsumption; % Energy Consumption
BAAMOutput.workingCapacity = workingCapacity; % Working capacity
BAAMOutput.massBalanceError = massBalanceError; % Mass balance error
BAAMOutput.distancePurityRecovery = distancePurityRecovery; % Distance from origin to Pu/Re curve
BAAMOutput.plotting = plotting; % Plotting data

simulationStopTime = datetime('now');
computationalTime = seconds(simulationStopTime-simulationStartTime);

%% COMMAND OUTPUT WINDOW
if ~silentFlag
    disp('- Model simulation completed.')
    disp('- Summary of the BAAM simulation results:')
    fprintf(' --> Maximum distance (95/90 Pu/Re) = %0.2f [-]\n',maxDistancePurityRecovery);
    fprintf(' --> Satisfy Constraints (1 - Yes, 0 - No) = %i\n',satisfyPurityRecovery);
    if ~isempty(indexDistance9590)
        fprintf(' --> Minimum Energy (BAAM) = %0.2f [kWh/tonne CO2]\n',minEnergyBAAM);
        fprintf(' --> Minimum Energy (Full Model) = %0.2f [kWh/tonne CO2]\n',minEnergyFullModel);
        fprintf(' --> Working Capacity (at minimum energy, BAAM) = %0.2f [mol CO2/m^3 ads.]\n',workingCapacity_MinEnergy);
    end
    fprintf(' --> Computational time: %2.1f [s].',computationalTime); fprintf('\n');
    disp('##########################################################################');
end
% END OF MAIN FUNCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% AUXILLIARY FUNCTIONS
    % These are functions that would be used by the main part of the 
    % simulator that would simulate the different steps of the different 
    % steps of the cycle and would evaluate the DSL isotherm to give the 
    % solid phase loadings of the two components

    %% Function: simulateBlowdownEvacuation
    % simulateBlowdownEvacuation - This function simulates the blodown step
    % and the evacuation step of the LPP cycle. The evolution of the mole
    % fraction of component A w.r.t. P is evaluted in this function using
    % eq. 7 in the original manuscript. The output from this function
    % consists of the pressure, molefraction of component A, solid loadings
    % of A and B, total moles desorbed, moles of A and B desorbed and the
    % energy consumption
    function outputBlowEvac = simulateBlowdownEvacuation()
        % Define the options for the ode solver that would solve eq. 7 in
        % the original manuscript
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        % Solve eq. 7, i.e. evaluate dy/dP for component A using the ode
        % solver, where the pressure spans the range of pressureVector. The
        % initial condition for the ode is the mole fraction of component A
        % in the feed stream. The output from the ode solver is the
        % pressure and the mole fraction of component A. 
        [pressureBlowEvac, molFracBlowEvac] = ode23s(@odeBlowEvac,pressureVector,molFracFeed_A,options);
        
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions obtained from solving
        % eq. 7 in the original manuscript.
        [solidPhaseLoading_BlowEvac_A, solidPhaseLoading_BlowEvac_B]...
                    = evaluateDSLIsotherm(pressureBlowEvac,molFracBlowEvac);
        
        % Initialize the total moles in the column, moles of A and B 
        % desorbed, and the energy consumption to zero vector
        % Total moles desorbed
        molesDesorbed_Total = zeros(size(pressureBlowEvac,1),1);
        % Moles of component A desorbed
        molesDesorbed_A = zeros(size(pressureBlowEvac,1),1);
        % Moles of component B desorbed
        molesDesorbed_B = zeros(size(pressureBlowEvac,1),1);
        % Energy consumption in the desorption step
        energyConsumptionDesorption = zeros(size(pressureBlowEvac,1),1);
        
        % Loop over the pressure range to compute the quantities of
        % interest described above. This would simply translate into
        % computing the quantities of interests for each step taken in the
        % pressure domain
        for counterDesorption = 1:size(pressureBlowEvac,1)-1
            % Initial condition at the step
            % Number of moles in the fluid phase
            % Component A
            fluidMolesDesorption_Init_A = (pressureBlowEvac(counterDesorption)*molFracBlowEvac(counterDesorption)...
                *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
            % Component B
            fluidMolesDesorption_Init_B = (pressureBlowEvac(counterDesorption)*(1-molFracBlowEvac(counterDesorption))...
                *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
            % Number of moles in the solid phase
            % Component A
            solidMolesDesorption_Init_A  = solidPhaseLoading_BlowEvac_A(counterDesorption)*adsorbentMass;
            % Component B
            solidMolesDesorption_Init_B  = solidPhaseLoading_BlowEvac_B(counterDesorption)*adsorbentMass;
            % Total moles in the solid and fluid phase
            % Component A
            totalMolesDesorption_Init_A = fluidMolesDesorption_Init_A + solidMolesDesorption_Init_A;
            % Component B
            totalMolesDesorption_Init_B = fluidMolesDesorption_Init_B + solidMolesDesorption_Init_B;
            
            % Final condition at the step
            % Number of moles in the fluid phase
            % Component A
            fluidMolesDesorption_Final_A = (pressureBlowEvac(counterDesorption+1)*molFracBlowEvac(counterDesorption+1)...
                *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
            % Component B
            fluidMolesDesorption_Final_B = (pressureBlowEvac(counterDesorption+1)*(1-molFracBlowEvac(counterDesorption+1))...
                *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
            % Number of moles in the solid phase
            % Component A
            solidMolesDesorption_Final_A = solidPhaseLoading_BlowEvac_A(counterDesorption+1)*adsorbentMass;
            % Component B
            solidMolesDesorption_Final_B = solidPhaseLoading_BlowEvac_B(counterDesorption+1)*adsorbentMass;
            % Total moles in the solid and fluid phase
            % Component A
            totalMolesDesorption_Final_A = fluidMolesDesorption_Final_A + solidMolesDesorption_Final_A;
            % Component B
            totalMolesDesorption_Final_B = fluidMolesDesorption_Final_B + solidMolesDesorption_Final_B;
            
            % Compute the moles of component A, B and total moles desorbed 
            % at that given step
            % Component A
            molesDesorbedStep_A = totalMolesDesorption_Init_A - totalMolesDesorption_Final_A;
            % Component B
            molesDesorbedStep_B = totalMolesDesorption_Init_B - totalMolesDesorption_Final_B;
            % Total moles
            molesDesorbedStep_Total = molesDesorbedStep_A + molesDesorbedStep_B;

            % Compute the moles of component A, B and the total moles
            % desorbed till that given step
            % Component A
            molesDesorbed_A(counterDesorption+1,1) = molesDesorbed_A(counterDesorption,1) + molesDesorbedStep_A;
            % Component B
            molesDesorbed_B(counterDesorption+1,1) = molesDesorbed_B(counterDesorption,1) + molesDesorbedStep_B;
            % Total moles
            molesDesorbed_Total(counterDesorption+1,1) = molesDesorbed_Total(counterDesorption,1) + molesDesorbedStep_Total;
            
            % Energy consumption till that given step
            energyConsumptionDesorption(counterDesorption+1,1)...
                = energyConsumptionDesorption(counterDesorption,1)...
                + ((adiabaticConstant/(adiabaticConstant-1))...
                *((UNIVERSAL_GAS_CONSTANT_JMK*temperature)/pumpEfficiency)...
                *molesDesorbedStep_Total*((1/pressureBlowEvac(counterDesorption+1))^((adiabaticConstant-1)/adiabaticConstant)-1));
        end
        
        % Prepare output matrix to be used in the successive steps
        % Matrix structure: [Pressure MoleFracA SolidLoadingA SolidLoadingB
        % MolesDesorbedTotal MolesDesorbedA MolesDesorbedB
        % EnergyConsumption]
        outputBlowEvac = [pressureBlowEvac molFracBlowEvac...
            solidPhaseLoading_BlowEvac_A solidPhaseLoading_BlowEvac_B...
                molesDesorbed_Total molesDesorbed_A molesDesorbed_B energyConsumptionDesorption];
        
        % The ode equation corresponding to dy/dP for the blowdown and 
        % evacuation step. The equation is an analytical form obtained from
        % <VISHAL NEEDS TO FILL IN HERE>
        function dydP=odeBlowEvac(P,y)
            dydP(1,1) = (UNIVERSAL_GAS_CONSTANT*temperature*adsorbentMass*y(1,1)*(y(1,1) - 1)*(UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite1_A*qSaturationSite1_A - UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite1_B*qSaturationSite1_B + UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A - UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B + P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B^2*qSaturationSite1_A + P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A - P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B^2*qSaturationSite1_B - P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B^2*qSaturationSite1_A*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) + 2*P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B^2*qSaturationSite1_B*y(1,1) + 2*P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1) + P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A^2*qSaturationSite1_A*y(1,1)^2 + P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B^2*qSaturationSite1_A*y(1,1)^2 + P^2*adsorptionEqbmConstSite1_A^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)^2 + P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A^2*qSaturationSite1_B*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B^2*qSaturationSite1_B*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_A^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)^2 + 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1) + 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_B*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)^2 + 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)^2 - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A*y(1,1)^2 + 2*P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_B*y(1,1)^2 + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite1_B - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*qSaturationSite1_A*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A*y(1,1) + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite1_B*y(1,1) + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite1_B*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1) + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)))/((adsorbentMass*((P*adsorptionEqbmConstSite1_A*qSaturationSite1_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1)) + (P*adsorptionEqbmConstSite2_A*qSaturationSite2_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1)) - (P^2*adsorptionEqbmConstSite1_A*qSaturationSite1_A*y(1,1)*(adsorptionEqbmConstSite1_A - adsorptionEqbmConstSite1_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2 - (P^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)*(adsorptionEqbmConstSite2_A - adsorptionEqbmConstSite2_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2) - adsorbentMass*y(1,1)*((P*adsorptionEqbmConstSite1_A*qSaturationSite1_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1)) - (P*adsorptionEqbmConstSite1_B*qSaturationSite1_B)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1)) + (P*adsorptionEqbmConstSite2_A*qSaturationSite2_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1)) - (P*adsorptionEqbmConstSite2_B*qSaturationSite2_B)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1)) + (P^2*adsorptionEqbmConstSite1_B*qSaturationSite1_B*(adsorptionEqbmConstSite1_A - adsorptionEqbmConstSite1_B)*(y(1,1) - 1))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2 + (P^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*(adsorptionEqbmConstSite2_A - adsorptionEqbmConstSite2_B)*(y(1,1) - 1))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2 - (P^2*adsorptionEqbmConstSite1_A*qSaturationSite1_A*y(1,1)*(adsorptionEqbmConstSite1_A - adsorptionEqbmConstSite1_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2 - (P^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)*(adsorptionEqbmConstSite2_A - adsorptionEqbmConstSite2_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2) + (P*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature))*(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2*(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2);
        end
    end
    
    %% Function: simulatePressurization
	% simulatePressurization - This function is a wrapper function that is
	% used to simulate the pressurization step of the LPP cycle. The system
	% of algaebric equations given by Eq. 11 and 12 are solved
	% simultaneously to determine the moles of gas required to pressurize
	% the column and the mole fraction of component A at the end of the LPP
	% step
    function outputPress = simulatePressurization(pressureLowStep,molFracEvacStep)
        % Define the options for the solver that would solve eq. 11/12 in 
        % the original manuscript
        options = optimoptions('fsolve','Display','off');
        % Solve eq. 11 and 12 to find the number of moles required to
        % pressurize the column back to the high pressure and to determine
        % the mole fraction of the gas that would be needed to pressurize
        % the column. The mole fraction at the end of this step comes from
        % the adsorption step as it is modeled as a plug flow. The output
        % from the solver would be the number of moles and the mole
        % fraction of A.
        % Initial condition for the solver
        % (1) - molFracPress_A; (2) - molesPress_Total
        inSolverPress0 = [0 0];
        % Solve the system of algaebric equations
        outSolverPress = fsolve(@(inSolverPress)...
            solveEquationsPressurization(inSolverPress,pressureLowStep,molFracEvacStep),inSolverPress0,options);
        % Initialize the output of the solver to new variables
        % Mole fraction at the end of the pressurization step
        molFracPress_A = outSolverPress(1);
        % Total moles required to pressurize the column from low pressure
        % to high pressure
        molesPress_Total = outSolverPress(2);
        
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the end of the
        % step
        [solidPhaseLoading_Press_Final_A, solidPhaseLoading_Press_Final_B]...
            = evaluateDSLIsotherm(pressureHigh,molFracPress_A);

        % Prepare output matrix to be used in the successive steps (vector)
        % Matrix structure: [MoleFracA SolidLoadingA SolidLoadingB
        % MolesPressTotal]
        outputPress = [molFracPress_A solidPhaseLoading_Press_Final_A...
                        solidPhaseLoading_Press_Final_B molesPress_Total];
    end

    %% Function: solveEquationsPressurization
	% solveEquationsPressurization - This function contains the equation
	% that is called by the function simulatePressurization which would be
	% used to determine the two quantities of interest for the
	% pressurization step, namely (1) - molFracPress_A; and
    % (2) - molesPress_Total. The output is in the vector outSolverPress
    function outSolverPress ...
            = solveEquationsPressurization(inSolverPress,pressureLowStep,molFracEvacStep)
        % Initial condition at the step
        % Number of moles in the fluid phase
        % Component A
        fluidMolesPressurization_Init_A = (pressureLowStep*molFracEvacStep*...
            columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Component B
        fluidMolesPressurization_Init_B = (pressureLowStep*(1-molFracEvacStep)...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the inital
        % time of the step
        [solidPhaseLoading_Press_Init_A, solidPhaseLoading_Press_Init_B]...
            = evaluateDSLIsotherm(pressureLowStep,molFracEvacStep);
        % Number of moles in the solid phase
        % Component A
        solidMolesPressurization_Init_A = solidPhaseLoading_Press_Init_A*adsorbentMass;
        % Component B
        solidMolesPressurization_Init_B = solidPhaseLoading_Press_Init_B*adsorbentMass;
        % Total moles in the solid and fluid phase
        % Component A
        totalMolesPressurization_Init_A = fluidMolesPressurization_Init_A + solidMolesPressurization_Init_A ;
        % Component B
        totalMolesPressurization_Init_B = fluidMolesPressurization_Init_B + solidMolesPressurization_Init_B;
        
        % Final condition at the step
        % Number of moles in the fluid phase
        % Component A
        fluidMolesPressurization_Final_A = (pressureHigh*inSolverPress(1)...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Component B
        fluidMolesPressurization_Final_B = (pressureHigh*(1-inSolverPress(1))...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the end of the
        % step
        [solidPhaseLoading_Press_Final_A, solidPhaseLoading_Press_Final_B]...
            = evaluateDSLIsotherm(pressureHigh,inSolverPress(1));
        % Number of moles in the solid phase
        % Component A
        solidMolesPressurization_Final_A = solidPhaseLoading_Press_Final_A*adsorbentMass;
        % Component B
        solidMolesPressurization_Final_B = solidPhaseLoading_Press_Final_B*adsorbentMass;
        % Total moles in the solid and fluid phase
        % Component A
        totalMolesPressurization_Final_A = fluidMolesPressurization_Final_A + solidMolesPressurization_Final_A;
        % Component B
        totalMolesPressurization_Final_B = fluidMolesPressurization_Final_B + solidMolesPressurization_Final_B ;
        
        % Equations to be solved by the solver which involves the overall
        % mass balance the component mass balance for component A. These
        % two equations corresponds to Eq. 11 and 12 in the original
        % manuscript
        % Equation 11 in the original manuscript
        outSolverPress(1,1) = (totalMolesPressurization_Init_A+totalMolesPressurization_Init_B)...
            - (totalMolesPressurization_Final_A+totalMolesPressurization_Final_B) + inSolverPress(2);
        % Equation 12 in the original manuscript
        outSolverPress(2,1) = totalMolesPressurization_Init_A - totalMolesPressurization_Final_A...
            + (inSolverPress(2)*inSolverPress(1));
    end

    %% Function: simulateAdsorption
	% simulateAdsorption - This function is a wrapper function that is
	% used to simulate the adsorption step of the LPP cycle. The system
	% of algaebric equations given by Eq. 13 and 14 are solved
	% simultaneously to determine the moles of gas required to be fed 
    % such that the gas breaks through and to determine the number of moles
    % of raffinate that is obtained during the adsorption step
    function outputAds = simulateAdsorption(molFracPressStep)
        % Define the options for the solver that would solve eq. 13/14 in 
        % the original manuscript
        options = optimoptions('fsolve','Display','off');
        % Solve eq. 13 and 14 to find the number of moles required for
        % adsorption step and the number of moles of gas that leaves the
        % column in the raffinate stream. The adsorption step is modeled as
        % a plug that leaves the column with a mole fraction corresponding
        % to the mole fraction at the end of the pressurization step
        % Initial condition for the solver
        % (1) - molesAds_Total; (2) - molesAds_Raffinate
        inSolverAds0 = [500 500];
        % Solve the system of algaebric equations
        outSolverAds = fsolve(@(inSolverAds)...
            solveEquationsAdsorption(inSolverAds,molFracPressStep),inSolverAds0,options);
        % Initialize the output of the solver to new variables
        % Number of moles of gas required in the adsorption step for the
        % breakthrough
        molesAds_Total = outSolverAds(1);
        % Number of moles of gas that leaves as the raffinate product in
        % the adsorption step
        molesAds_Raffinate = outSolverAds(2);

        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the final
        % time of the step
        [solidPhaseLoading_Ads_Final_A, solidPhaseLoading_Ads_Final_B]...
            = evaluateDSLIsotherm(pressureHigh,molFracFeed_A);
        % Prepare output matrix to be used in the successive steps (vector)
        % Matrix structure: [molesAds_Total molesAds_Raffinate 
        % SolidLoadingA SolidLoadingB]
        outputAds = [molesAds_Total molesAds_Raffinate...
            solidPhaseLoading_Ads_Final_A solidPhaseLoading_Ads_Final_B];
    end

    %% Function: solveEquationsAdsorption
	% solveEquationsAdsorption - This function contains the equation
	% that is called by the function simulateAdsorption which would be
	% used to determine the two quantities of interest for the adsorption
	% step, namely (1) - molesAds_Total; and (2) - molesAds_Raffinate. 
    % The output is in the vector outSolverPress
    function outSolverAds = solveEquationsAdsorption(inSolverAds,molFracPressStep)
        % Initial condition at the step
        % Number of moles in the fluid phase
        % Component A
        fluidMolesAdsorption_Init_A = (pressureHigh*molFracPressStep...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Component B
        fluidMolesAdsorption_Init_B = (pressureHigh*(1-molFracPressStep)...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the initial time of
        % the step
        [solidPhaseLoading_Ads_Init_A, solidPhaseLoading_Ads_Init_B]...
            = evaluateDSLIsotherm(pressureHigh,molFracPressStep);
        % Number of moles in the solid phase
        % Component A
        solidMolesAdsorption_Init_A = solidPhaseLoading_Ads_Init_A*adsorbentMass;
        % Component B
        solidMolesAdsorption_Init_B = solidPhaseLoading_Ads_Init_B*adsorbentMass;
        % Total moles in the solid and fluid phase
        % Component A
        totalMolesAdsorption_Init_A = fluidMolesAdsorption_Init_A + solidMolesAdsorption_Init_A;
        % Component B
        totalMolesAdsorption_Init_B = fluidMolesAdsorption_Init_B + solidMolesAdsorption_Init_B;
        
        % Final condition at the step
        % Number of moles in the fluid phase
        % Component A
        fluidMolesAdsorption_Final_A = (pressureHigh*molFracFeed_A...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Component B
        fluidMolesAdsorption_Final_B = (pressureHigh*(1-molFracFeed_A)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the end of the
        % step
        [solidPhaseLoading_Ads_Final_A, solidPhaseLoading_Ads_Final_B]...
            = evaluateDSLIsotherm(pressureHigh,molFracFeed_A);
        % Number of moles in the solid phase
        % Component A
        solidMolesAdsorption_Final_A = solidPhaseLoading_Ads_Final_A*adsorbentMass;
        % Component B
        solidMolesAdsorption_Final_B = solidPhaseLoading_Ads_Final_B*adsorbentMass;
        % Total moles in the solid and fluid phase
        % Component A
        totalMolesAdsorption_Final_A = fluidMolesAdsorption_Final_A + solidMolesAdsorption_Final_A;
        % Component B
        totalMolesAdsorption_Final_B = fluidMolesAdsorption_Final_B + solidMolesAdsorption_Final_B;
        
        % Equations to be solved by the solver which involves the overall
        % mass balance the component mass balance for component A. These
        % two equations corresponds to Eq. 13 and 14 in the original
        % manuscript
        % Equation 13 in the original manuscript  
        outSolverAds(1,1) = (totalMolesAdsorption_Init_A + totalMolesAdsorption_Init_B)...
                   - (totalMolesAdsorption_Final_A + totalMolesAdsorption_Final_B) ...
                   + inSolverAds(1) - inSolverAds(2);
        % Equation 14 in the original manuscript  
        outSolverAds(2,1) = totalMolesAdsorption_Init_A - totalMolesAdsorption_Final_A...
                   + (inSolverAds(1)*molFracFeed_A) - (inSolverAds(2)*molFracPressStep);
    end

    %% Function: evaluateDSLIsotherm
    % evaluateDSLIsotherm - This function evaluates the DSL isotherm and 
    % gives away the solid phase loadings for component A and B at the 
    % given pressure and mole fraction
    function [solidPhaseLoadingLoc_A, solidPhaseLoadingLoc_B]...
                = evaluateDSLIsotherm(pressureIsotherm,molFracIsotherm)
        % Fluid phase concentrations
        % Evaluate the concentration of component A at the given pressure
        % and mole fraction
        concentration_A = (pressureIsotherm.*molFracIsotherm)./(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the concentration of component B at the given pressure
        % and mole fraction
        concentration_B = (pressureIsotherm.*(1-molFracIsotherm))./(UNIVERSAL_GAS_CONSTANT*temperature);
        
        % Solid phase concentrations
        % Evaluate the solid phase loading for component A
        % Site 1 [mol/kg]
        solidPhaseLoadingSite1_A = (qSaturationSite1_A*adsorptionEqbmConstSite1_A*concentration_A)...
            ./(1 + adsorptionEqbmConstSite1_A.*concentration_A + adsorptionEqbmConstSite1_B.*concentration_B);
        % Site 2 [mol/kg]
        solidPhaseLoadingSite2_A = (qSaturationSite2_A*adsorptionEqbmConstSite2_A*concentration_A)...
            ./(1 + adsorptionEqbmConstSite2_A.*concentration_A + adsorptionEqbmConstSite2_B.*concentration_B);
        % Total solid phase loading of component A [mol/kg]
        solidPhaseLoadingLoc_A = solidPhaseLoadingSite1_A + solidPhaseLoadingSite2_A;
        % Evaluate the solid phase loading for component B
        % Site 1 [mol/kg]
        solidPhaseLoadingSite1_B = (qSaturationSite1_B*adsorptionEqbmConstSite1_B*concentration_B)...
            ./(1 + adsorptionEqbmConstSite1_A.*concentration_A + adsorptionEqbmConstSite1_B.*concentration_B);
        % Site 2 [mol/kg]
        solidPhaseLoadingSite2_B = (qSaturationSite2_B*adsorptionEqbmConstSite2_B*concentration_B)...
            ./(1 + adsorptionEqbmConstSite2_A.*concentration_A + adsorptionEqbmConstSite2_B.*concentration_B);
        % Total solid phase loading of component B [mol/kg]
        solidPhaseLoadingLoc_B = solidPhaseLoadingSite1_B + solidPhaseLoadingSite2_B;
    end
end