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
% - 2019-01-08, AK: Cleaned up the pressurization part of the code
% - 2019-01-07, AK: Cleaned up the variable names, added a new input
%                   structure and cleaned up the blowdown/evacuation part
%                   of the code
% - 2019-01-06, AK: Introduced header for this file
%
% Input arguments:
% - 
%
% Output arguments:
% - 
%
% Comments:
% - Hardcoded values need to be changed in case if the pressure ranges are
% different (I have changed it from plow to phigh to be universal can be
% changed back again. Search for harcoded and look at the original code)
% - Changed deltaPressure from 0.0001 to 0.01
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qA1, qA2, qA3, qA4, qA5, qA6,qB1, qB2, qB3, qB4, qB5, qB6, pspand,YALL1,qAfLPPd,qBfLPPd,qAffeed,qBffeed,YOUT]=BAAM(simInfo,adsInfo)
tic

%% STARTUP FORMALITIES OF THE CODE
% ADSORBENT PROPERTIES
adsName = adsInfo.adsorbentName; % Name of the adsorbent being simulated
adsorbentMass = 1; % Mass of the adsorbent [kg]
adsorbentDensity = simInfo.adsorbentDensity; % Particle Density [kg/m3]

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
molFracFeed_A = simInfo.molFracFeed_A; % Feed mole fraction of component A [-]

% UNIVERSAL CONSTANTS
UNIVERSAL_GAS_CONSTANT = 8.314e-5; % Universal gas constant [(m3 bar)/(K mol)]
UNIVERSAL_GAS_CONSTANT_JMK = 8.314; % Universal gas constant [J/mol K]

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

%% SOLVER DEFINITIONS/PRESSURE MATRIX
% ODE SOLVER PROPERTIES
deltaPressure = 0.01; % Length of discretized grid for pressure [bar]
% Generate a pressure vector that spans from the high pressure to the low
% pressure with a length of deltaPressure for the ode solver to be used in
% the blowdown step
pressureVector=(pressureHigh:-deltaPressure:pressureLow)';

% PRESSURE MATRIX
% Create a pressure vector for the low pressure that spans a range from
% low pressure (pressureLow) to 0.1 bar. The upper bound of the low
% pressure vector is HARDCODED. The pressure vector is in bar.
pressureLowVector = pressureLow:0.01:(pressureHigh - 0.01);
% Create a pressure vector for the intermediate pressure that spans a range
% from low pressure (pressureLow+0.01 bar) to the high pressure
% (pressureHigh - 0.01 bar). The bounds for the intermediate pressure
% vector is HARDCODED. The pressure vector is in bar.
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

%%%%%%%%%%%%%%%%%%% DONNO WHAT TO DO %%%%%%%%%%%%%%%%%%% 
YOUT(1,1:9)=0; % For Printing the simulation output
kinc=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% RUN THE MODEL
% Simulate the blowdown and the evacuation step using the function
% simulateBlowdownEvacuation and obtain the pressure, molefraction of 
% component A, solid loadings of A and B, total moles desorbed, moles of A 
% and B desorbed and the energy consumption
[YALL1]=simulateBlowdownEvacuation();

%%% TO VISHAL: Why do you even need the matrix from meshgrid? That thing
%%% would slow the code down

% Simulate the pressurization and adsorption steps. Loop over all the low
% and intermediate pressures. 
% Loop over all entire low pressure vector
for counterLow = 1:length(pressureLowVector)
    % Find the index (in the output of the simulateBlowdownEvacuation
    % function) of the low pressure that corresponds to the current low
    % pressure value in the loop. Have a tolerance of 10e-6 just in case.
    pressureLowIndex = find(abs(YALL1(:,1)-pressureLowGrid(1,counterLow))<1e-6);
    
    % Simulate the pressurization step starting from the low pressure with
    % the goal of reaching the high pressure. The input for this function
    % would be the low pressure and mole fraction at the end of the
    % blowdown/evacuation step (YALL1(:,2) is the mole fraction)
    [yAfLPP,~,qAfLPP, qBfLPP,NLPP] ...
        = simulatePressurization(pressureLowVector(counterLow),YALL1(pressureLowIndex,2));
    
    %%%%%%%%%%%%% DONE TILL HERE %%%%%%%%%%%%%
    
    [Nfeed,Nwaste,~,qAffeed,qBffeed]=mass_feed_VSB(yAfLPP);
    
    for counterInter = 1:length(pressureInterVector)
        if ~(isnan(pressureInterGrid(counterInter,counterLow)) ||  isnan(pressureLowGrid(1,counterLow)))
            
            if (pressureLowGrid(1,counterLow)==pressureLow)
                qAfLPPd=qAfLPP;
                qBfLPPd=qBfLPP;
            end
            
            
            inter=find(abs(YALL1(:,1)-pressureInterGrid(counterInter,counterLow))<10^(-6));
            nAbd=YALL1(inter,6);
            nBbd=YALL1(inter,7);
            nAevc=YALL1(pressureLowIndex,6)-YALL1(inter,6);
            nBevc=YALL1(pressureLowIndex,7)-YALL1(inter,7);
            
            purity=(nAevc/(nAevc+nBevc))*100;
            recovery=(nAevc/(Nfeed*molFracFeed_A))*100;
            enerbd=(YALL1(inter,8)*2.77778e-7)/(nAevc*44*1e-6);
            enerevc=((YALL1(pressureLowIndex,8)-YALL1(inter,8))*2.77778e-7)/(nAevc*44*1e-6);
            enerT=enerbd+enerevc;
            wc=(nAevc)/(columnVolume*(1-voidFraction));
            PurP(counterInter,counterLow)=purity;
            RecP(counterInter,counterLow)=recovery;
            EnerP(counterInter,counterLow)=enerT;
            WcP(counterInter,counterLow)=wc;
            
            ncyclein=Nfeed;
            ncycleout=(Nwaste-NLPP) + (nAbd+nBbd) + (nAevc+nBevc);
            nTotalerr=((ncyclein-ncycleout)/ncycleout)*100;
            YOUT(kinc,:)=[pressureLowGrid(1,counterLow) pressureInterGrid(counterInter,counterLow) purity recovery enerbd enerevc enerT wc nTotalerr];
            kinc=kinc+1;
            
            %fprintf(fid,'%f  %s  %s  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f\n',adsindex,',',Adsorbents,',',pressureLowGrid(1,counterLow),',',pressureInterGrid(counterInter,counterLow),',',purity,',',recovery,',',enerbd,',',enerevc,',',enerT,',',wc,',',nTotalerr);
            
        end
    end
    
end

PurP(PurP==0)=NaN;
PurPmax=PurP(:,1);
RecP(RecP==0)=NaN;
RecPmax=RecP(:,1);
EnerP(EnerP==0)=NaN;
WcP(WcP==0)=NaN;


pspand=0:0.001:1;
pspand=pspand';

[qA1, qB1]=isotherm(pspand,0);
[qA2, qB2]=isotherm(pspand,0.15);
[qA3, qB3]=isotherm(pspand,0.30);
[qA4, qB4]=isotherm(pspand,0.45);
[qA5, qB5]=isotherm(pspand,0.75);
[qA6, qB6]=isotherm(pspand,1);



%fclose(fid);


YOUT(:,10)=sqrt(YOUT(:,3).^2 + YOUT(:,4).^2);

inddis=find(YOUT(:,10)<110.3 & YOUT(:,10)>110.2);
maxdis=max(YOUT(:,10));

if ~isempty(inddis)
    YOUTdis(1:size(inddis,1),:)=YOUT(inddis(1:end,1),:);
    [minenergy, indmine]=min(YOUTdis(:,7));
    Fullenergy=(minenergy*1.1446) + 66.528;
    wcminenergy=YOUTdis(indmine,8);
    
    fid = fopen('BAAM_Results.txt', 'at');
    fprintf(fid,'%s%s%f%s%f%s%f%s%f%s%f\n',Adsorbents,',',temperature,',',maxdis,',',minenergy,',',Fullenergy,',',wcminenergy);
    fclose(fid);
else
    fid = fopen('BAAM_error.txt', 'at');
    fprintf(fid,'%s%s%f%s%f%s%s%s%s\n',Adsorbents,',',temperature,',',maxdis,',','-',',','-');
    fclose(fid);
end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        molesDesorbed_Total = zeros(1,size(pressureBlowEvac,1));
        % Moles of component A desorbed
        molesDesorbed_A = zeros(1,size(pressureBlowEvac,1));
        % Moles of component B desorbed
        molesDesorbed_B = zeros(1,size(pressureBlowEvac,1));
        % Energy consumption in the desorption step
        energyConsumptionDesorption = zeros(1,size(pressureBlowEvac,1));
        
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
    function [molFracPress_A,massBalErrorPress,solidPhaseLoading_Press_Final_A,solidPhaseLoading_Press_Final_B,molesPress_Total]...
            = simulatePressurization(pressureLowStep,molFracEvacStep)
        % Define the options for the ode solver that would solve eq. 11/12 
        % in the original manuscript
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
        
        % Initial condition at the step
        % Number of moles in the fluid phase
        % Component A     
        fluidMolesPressurization_Init_A = (pressureLowStep*molFracEvacStep...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Component B
        fluidMolesPressurization_Init_B = (pressureLowStep*(1-molFracEvacStep)...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the conditions
        % corresponding to the pressure and mole fraction at the inital
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
        totalMolesPressurization_Init_A = fluidMolesPressurization_Init_A + solidMolesPressurization_Init_A;
        % Component B
        totalMolesPressurization_Init_B = fluidMolesPressurization_Init_B + solidMolesPressurization_Init_B;
        
        % Final condition at the step
        % Number of moles in the fluid phase
        % Component A
        fluidMolesPressurization_Final_A = (pressureHigh*molFracPress_A...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Component B
        fluidMolesPressurization_Final_B = (pressureHigh*(1-molFracPress_A)...
            *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions at the conditions
        % corresponding to the pressure and mole fraction at the end of the
        % step
        [solidPhaseLoading_Press_Final_A, solidPhaseLoading_Press_Final_B]...
            = evaluateDSLIsotherm(pressureHigh,molFracPress_A);
        % Number of moles in the solid phase
        % Component A
        solidMolesPressurization_Final_A = solidPhaseLoading_Press_Final_A*adsorbentMass;
        % Component B
        solidMolesPressurization_Final_B = solidPhaseLoading_Press_Final_B*adsorbentMass;
        % Total moles in the solid and fluid phase
        % Component A
        totalMolesPressurization_Final_A = fluidMolesPressurization_Final_A + solidMolesPressurization_Final_A;
        % Component B
        totalMolesPressurization_Final_B = fluidMolesPressurization_Final_B + solidMolesPressurization_Final_B;
        
        % Compute the total moles in the column initially and finally for
        % the pressureization step to compute the mass balance error
        % Total number of moles at the initial time of the step along with
        % the moles for pressurization
        molesPress_Init_Total = totalMolesPressurization_Init_A + (molesPress_Total*molFracPress_A)...
                    + totalMolesPressurization_Init_B + (molesPress_Total*(1-molFracPress_A));
        % Total number of moles at the final time of the step along with
        % the moles for pressurization
        molesPress_Final_Total = totalMolesPressurization_Final_A + totalMolesPressurization_Final_B;
        % Compute the mass balance error
        massBalErrorPress = abs(((molesPress_Init_Total-molesPress_Final_Total)/molesPress_Init_Total)*100);
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
        % corresponding pressures and mole fractions at the conditions
        % corresponding to the pressure and mole fraction at the inital
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
        % corresponding pressures and mole fractions at the conditions
        % corresponding to the pressure and mole fraction at the end of the
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

    %%
    function [Nfeed,Nwaste,mbalfeederror,qAffeed,qBffeed]=mass_feed_VSB(yAfLPP)
        
        
        
        %options=optimset('Algorithm','Levenberg-Marquardt');
        options = optimoptions('fsolve','Display','off');
        
        zfeed=fsolve(@(x) feed_VSB(x,yAfLPP),[500 500],options);
        
        fAifeed=(pressureHigh*yAfLPP*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBifeed=(pressureHigh*(1-yAfLPP)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAifeed, qBifeed]=isotherm(pressureHigh,yAfLPP);
        
        NAifeed=fAifeed+(qAifeed*adsorbentMass);
        NBifeed=fBifeed+(qBifeed*adsorbentMass);
        
        %Final
        
        fAffeed=(pressureHigh*molFracFeed_A*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBffeed=(pressureHigh*(1-molFracFeed_A)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAffeed, qBffeed]=isotherm(pressureHigh,molFracFeed_A);
        
        NAffeed=fAffeed+(qAffeed*adsorbentMass);
        NBffeed=fBffeed+(qBffeed*adsorbentMass);
        
        
        %fprintf('\n\n***Feed STEP***');
        Nfeed=zfeed(1);
        Nwaste=zfeed(2);
        
        Nfeedinitial=NAifeed + (Nfeed*molFracFeed_A) -(Nwaste*yAfLPP) + NBifeed + Nfeed*(1-molFracFeed_A)-Nwaste*(1-yAfLPP);
        Nfeedfinal=NAffeed+NBffeed;
        
        mbalfeederror=abs(((Nfeedinitial-Nfeedfinal)/Nfeedinitial)*100);
        
        
    end

    %%
    function zfeed=feed_VSB(xfeed,yAfLPP)
        
        
        
        %xfeed(1) -> Nfeed
        %xfeed(2) -> NWaste
        
        %Initial
        
        fAifeed=(pressureHigh*yAfLPP*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBifeed=(pressureHigh*(1-yAfLPP)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAifeed, qBifeed]=isotherm(pressureHigh,yAfLPP);
        
        NAiLPP=fAifeed+(qAifeed*adsorbentMass);
        NBiLPP=fBifeed+(qBifeed*adsorbentMass);
        
        %Final
        fAffeed=(pressureHigh*molFracFeed_A*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBfLPP=(pressureHigh*(1-molFracFeed_A)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAffeed, qBffeed]=isotherm(pressureHigh,molFracFeed_A);
        
        NAffeed=fAffeed+(qAffeed*adsorbentMass);
        NBffeed=fBfLPP+(qBffeed*adsorbentMass);
        
        %Overall Balance
        
        zfeed(1,1)=(NAiLPP+NBiLPP)-(NAffeed+NBffeed) +xfeed(1) -xfeed(2);
        zfeed(2,1)=NAiLPP-NAffeed+(xfeed(1)*molFracFeed_A)-(xfeed(2)*yAfLPP);
        
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
toc
end