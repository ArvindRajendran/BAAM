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
feedMolFrac_A = simInfo.feedMolFrac_A; % Feed mole fraction of component A [-]

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
lowPressureVector = pressureLow:0.01:(pressureHigh - 0.01);
% Create a pressure vector for the intermediate pressure that spans a range
% from low pressure (pressureLow+0.01 bar) to the high pressure
% (pressureHigh - 0.01 bar). The bounds for the intermediate pressure
% vector is HARDCODED. The pressure vector is in bar.
interPressureVector = (pressureLow + 0.01):0.01:(pressureHigh - 0.01);
% Generate a matrix with all combinations of low and intermediate pressures
[pressureLowGrid, pressureInterGrid] = meshgrid(lowPressureVector,interPressureVector);

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

%%%%%%%%%%%%% DONE TILL HERE %%%%%%%%%%%%%


for i= 1:size(pressureLowGrid,2)
    
    lowp=find(abs(YALL1(:,1)-pressureLowGrid(1,i))<10^(-6));
    
    [yAfLPP,~,qAfLPP, qBfLPP,NLPP]=mass_LPP_VSB(pressureLowGrid(1,i),YALL1(lowp,2));
    [Nfeed,Nwaste,~,qAffeed,qBffeed]=mass_feed_VSB(yAfLPP);
    
    for j=1:size(pressureLowGrid,1)
        
        if ~(isnan(pressureInterGrid(j,i)) ||  isnan(pressureLowGrid(1,i)))
            
            if (pressureLowGrid(1,i)==pressureLow)
                qAfLPPd=qAfLPP;
                qBfLPPd=qBfLPP;
            end
            
            
            inter=find(abs(YALL1(:,1)-pressureInterGrid(j,i))<10^(-6));
            nAbd=YALL1(inter,6);
            nBbd=YALL1(inter,7);
            nAevc=YALL1(lowp,6)-YALL1(inter,6);
            nBevc=YALL1(lowp,7)-YALL1(inter,7);
            
            purity=(nAevc/(nAevc+nBevc))*100;
            recovery=(nAevc/(Nfeed*feedMolFrac_A))*100;
            enerbd=(YALL1(inter,8)*2.77778e-7)/(nAevc*44*1e-6);
            enerevc=((YALL1(lowp,8)-YALL1(inter,8))*2.77778e-7)/(nAevc*44*1e-6);
            enerT=enerbd+enerevc;
            wc=(nAevc)/(columnVolume*(1-voidFraction));
            PurP(j,i)=purity;
            RecP(j,i)=recovery;
            EnerP(j,i)=enerT;
            WcP(j,i)=wc;
            
            ncyclein=Nfeed;
            ncycleout=(Nwaste-NLPP) + (nAbd+nBbd) + (nAevc+nBevc);
            nTotalerr=((ncyclein-ncycleout)/ncycleout)*100;
            YOUT(kinc,:)=[pressureLowGrid(1,i) pressureInterGrid(j,i) purity recovery enerbd enerevc enerT wc nTotalerr];
            kinc=kinc+1;
            
            %fprintf(fid,'%f  %s  %s  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f\n',adsindex,',',Adsorbents,',',pressureLowGrid(1,i),',',pressureInterGrid(j,i),',',purity,',',recovery,',',enerbd,',',enerevc,',',enerT,',',wc,',',nTotalerr);
            
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
        [pressureBlowEvac, moleFracBlowEvac] = ode23s(@odeBlowEvac,pressureVector,feedMolFrac_A,options);
        
        % Evaluate the solid phase loadings for both the components at the
        % corresponding pressures and mole fractions obtained from solving
        % eq. 7 in the original manuscript.
        [solidPhaseLoading_BlowEvac_A, solidPhaseLoading_BlowEvac_B]...
                    = evaluateDSLIsotherm(pressureBlowEvac,moleFracBlowEvac);
        
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
            fluidMolesDesorption_Init_A = (pressureBlowEvac(counterDesorption)*moleFracBlowEvac(counterDesorption)...
                *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
            % Component B
            fluidMolesDesorption_Init_B = (pressureBlowEvac(counterDesorption)*(1-moleFracBlowEvac(counterDesorption))...
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
            fluidMolesDesorption_Final_A = (pressureBlowEvac(counterDesorption+1)*moleFracBlowEvac(counterDesorption+1)...
                *columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
            % Component B
            fluidMolesDesorption_Final_B = (pressureBlowEvac(counterDesorption+1)*(1-moleFracBlowEvac(counterDesorption+1))...
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
        outputBlowEvac = [pressureBlowEvac moleFracBlowEvac...
            solidPhaseLoading_BlowEvac_A solidPhaseLoading_BlowEvac_B...
                molesDesorbed_Total molesDesorbed_A molesDesorbed_B energyConsumptionDesorption];
        
        % The ode equation corresponding to dy/dP for the blowdown and 
        % evacuation step. The equation is an analytical form obtained from
        % <VISHAL NEEDS TO FILL IN HERE>
        function dydP=odeBlowEvac(P,y)
            dydP(1,1) = (UNIVERSAL_GAS_CONSTANT*temperature*adsorbentMass*y(1,1)*(y(1,1) - 1)*(UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite1_A*qSaturationSite1_A - UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite1_B*qSaturationSite1_B + UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A - UNIVERSAL_GAS_CONSTANT^2*temperature^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B + P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B^2*qSaturationSite1_A + P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A - P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B^2*qSaturationSite1_B - P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B^2*qSaturationSite1_A*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) + 2*P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B^2*qSaturationSite1_B*y(1,1) + 2*P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1) + P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A^2*qSaturationSite1_A*y(1,1)^2 + P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B^2*qSaturationSite1_A*y(1,1)^2 + P^2*adsorptionEqbmConstSite1_A^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)^2 + P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A^2*qSaturationSite1_B*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B^2*qSaturationSite1_B*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_A^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)^2 - P^2*adsorptionEqbmConstSite1_B^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)^2 + 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1) + 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_B*y(1,1) - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)^2 + 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)^2 - 2*P^2*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A*y(1,1)^2 + 2*P^2*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*adsorptionEqbmConstSite2_B*qSaturationSite1_B*y(1,1)^2 + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite1_B - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*qSaturationSite1_A*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B*qSaturationSite1_A*y(1,1) + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_A*qSaturationSite1_B*y(1,1) + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite1_B*y(1,1) - 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_A*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1) + 2*P*UNIVERSAL_GAS_CONSTANT*temperature*adsorptionEqbmConstSite1_B*adsorptionEqbmConstSite2_B*qSaturationSite2_B*y(1,1)))/((adsorbentMass*((P*adsorptionEqbmConstSite1_A*qSaturationSite1_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1)) + (P*adsorptionEqbmConstSite2_A*qSaturationSite2_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1)) - (P^2*adsorptionEqbmConstSite1_A*qSaturationSite1_A*y(1,1)*(adsorptionEqbmConstSite1_A - adsorptionEqbmConstSite1_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2 - (P^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)*(adsorptionEqbmConstSite2_A - adsorptionEqbmConstSite2_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2) - adsorbentMass*y(1,1)*((P*adsorptionEqbmConstSite1_A*qSaturationSite1_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1)) - (P*adsorptionEqbmConstSite1_B*qSaturationSite1_B)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1)) + (P*adsorptionEqbmConstSite2_A*qSaturationSite2_A)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1)) - (P*adsorptionEqbmConstSite2_B*qSaturationSite2_B)/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1)) + (P^2*adsorptionEqbmConstSite1_B*qSaturationSite1_B*(adsorptionEqbmConstSite1_A - adsorptionEqbmConstSite1_B)*(y(1,1) - 1))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2 + (P^2*adsorptionEqbmConstSite2_B*qSaturationSite2_B*(adsorptionEqbmConstSite2_A - adsorptionEqbmConstSite2_B)*(y(1,1) - 1))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2 - (P^2*adsorptionEqbmConstSite1_A*qSaturationSite1_A*y(1,1)*(adsorptionEqbmConstSite1_A - adsorptionEqbmConstSite1_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2 - (P^2*adsorptionEqbmConstSite2_A*qSaturationSite2_A*y(1,1)*(adsorptionEqbmConstSite2_A - adsorptionEqbmConstSite2_B))/(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2) + (P*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature))*(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite1_B + P*adsorptionEqbmConstSite1_A*y(1,1) - P*adsorptionEqbmConstSite1_B*y(1,1))^2*(UNIVERSAL_GAS_CONSTANT*temperature + P*adsorptionEqbmConstSite2_B + P*adsorptionEqbmConstSite2_A*y(1,1) - P*adsorptionEqbmConstSite2_B*y(1,1))^2);
        end
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

    %%
    function [yAfLPP,mbalLPPerror,qAfLPP, qBfLPP,NLPP]=mass_LPP_VSB(Pl,yAfevc)
        
        
        
        %options=optimset('Algorithm','Levenberg-Marquardt');
        options = optimoptions('fsolve','Display','off');
        
        zLPP=fsolve( @(x) LPP_VSB(x,Pl,yAfevc),[0 0],options);
        
        %Initial
        
        fAiLPP=(Pl*yAfevc*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBiLPP=(Pl*(1-yAfevc)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        
        [qAiLPP, qBiLPP]=isotherm(Pl,yAfevc);
        
        NAiLPP=fAiLPP+(qAiLPP*adsorbentMass);
        NBiLPP=fBiLPP+(qBiLPP*adsorbentMass);
        
        %Final
        fAfLPP=(pressureHigh*zLPP(1)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBfLPP=(pressureHigh*(1-zLPP(1))*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAfLPP, qBfLPP]=isotherm(pressureHigh,zLPP(1));
        
        NAfLPP=fAfLPP+(qAfLPP*adsorbentMass);
        NBfLPP=fBfLPP+(qBfLPP*adsorbentMass);
        
        %fprintf('\n***PRESSURISATION STEP***');
        yAfLPP=zLPP(1);
        NLPP=zLPP(2);
        
        NLPPinitial=NAiLPP+(NLPP*yAfLPP)+NBiLPP+(NLPP*(1-yAfLPP));
        NLPPfinal=NAfLPP+NBfLPP;
        mbalLPPerror=abs(((NLPPinitial-NLPPfinal)/NLPPinitial)*100);
        
    end

    function zLPP=LPP_VSB(xLPP,Pl,yAfec)
        
        
        %xLPP(1) -> yAfLPP
        %xLPP(2) -> NLPP
        
        %Initial
        
        fAiLPP=(Pl*yAfec*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBiLPP=(Pl*(1-yAfec)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAiLPP, qBiLPP]=isotherm(Pl,yAfec);
        
        NAiLPP=fAiLPP+(qAiLPP*adsorbentMass);
        NBiLPP=fBiLPP+(qBiLPP*adsorbentMass);
        
        %Final
        fAfLPP=(pressureHigh*xLPP(1)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBfLPP=(pressureHigh*(1-xLPP(1))*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAfLPP, qBfLPP]=isotherm(pressureHigh,xLPP(1));
        
        NAfLPP=fAfLPP+(qAfLPP*adsorbentMass);
        NBfLPP=fBfLPP+(qBfLPP*adsorbentMass);
        
        %Overall Balance
        
        zLPP(1,1)=(NAiLPP+NBiLPP)-(NAfLPP+NBfLPP)+xLPP(2);
        zLPP(2,1)=NAiLPP-NAfLPP+(xLPP(2)*xLPP(1));
        
    end


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
        
        fAffeed=(pressureHigh*feedMolFrac_A*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBffeed=(pressureHigh*(1-feedMolFrac_A)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAffeed, qBffeed]=isotherm(pressureHigh,feedMolFrac_A);
        
        NAffeed=fAffeed+(qAffeed*adsorbentMass);
        NBffeed=fBffeed+(qBffeed*adsorbentMass);
        
        
        %fprintf('\n\n***Feed STEP***');
        Nfeed=zfeed(1);
        Nwaste=zfeed(2);
        
        Nfeedinitial=NAifeed + (Nfeed*feedMolFrac_A) -(Nwaste*yAfLPP) + NBifeed + Nfeed*(1-feedMolFrac_A)-Nwaste*(1-yAfLPP);
        Nfeedfinal=NAffeed+NBffeed;
        
        mbalfeederror=abs(((Nfeedinitial-Nfeedfinal)/Nfeedinitial)*100);
        
        
    end

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
        fAffeed=(pressureHigh*feedMolFrac_A*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        fBfLPP=(pressureHigh*(1-feedMolFrac_A)*columnVolume*voidFraction)/(UNIVERSAL_GAS_CONSTANT*temperature);
        [qAffeed, qBffeed]=isotherm(pressureHigh,feedMolFrac_A);
        
        NAffeed=fAffeed+(qAffeed*adsorbentMass);
        NBffeed=fBfLPP+(qBffeed*adsorbentMass);
        
        %Overall Balance
        
        zfeed(1,1)=(NAiLPP+NBiLPP)-(NAffeed+NBffeed) +xfeed(1) -xfeed(2);
        zfeed(2,1)=NAiLPP-NAffeed+(xfeed(1)*feedMolFrac_A)-(xfeed(2)*yAfLPP);
        
    end



toc
end