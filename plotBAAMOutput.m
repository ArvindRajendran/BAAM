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
% This function plots the model output generated by BAAM. The plots
% generated in this figure are similar to the plots shown in Fig. 3 of the
% original manuscript.
%
% Last modified:
% - 2019-02-02, AK: Removed some plots and fixed the issue with the N2
%                   isotherm
% - 2019-02-01, AK: Introduced header for this file
%
% Input arguments:
% - BAAMOutput:    Structure containing output from the BAAM model. This
%                  contains performance indicators, and output from each
%                  step that was simulated.
%
% Output arguments:
% - N/A
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotBAAMOutput(outputStruct)
% Loop over all the adsorbents in the outputStruct. The number of
% adsorbents is one less than the number of fields in outputStruct. The
% other one has the gitCommit
for adsorbentNumber = 1:length(fields(outputStruct))-1
    % Open a new figure for each adsorbent and turn off the numbering
    figure('units','normalized','Position',[0.1 0.4 0.6 0.4],'Color',[1 1 1],'numbertitle', 'off');
    % Load the BAAM output for each adsorbent and store them in a local
    % variable
    localAdsorbent = eval(['outputStruct.ADS_',num2str(adsorbentNumber)]);
    % Find the number of mole fractions that would be used to plot the
    % isotherm according to Fig. 3 in the original manuscript
    numberOfMoleFractions = length(localAdsorbent.plotting);
    % Get the color code for the isotherm plots
    colorCode = colormap(jet(numberOfMoleFractions));
    % Get the mole fraction from the plotting structure for both component
    % A and component B
    moleFraction_A = [outputStruct.ADS_1.plotting.moleFraction]; moleFraction_B = 1 - moleFraction_A;
    % Convert the mole fractions into strings so that they can be used in
    % the legend of the figure
    legendStr_A = cellstr(num2str(moleFraction_A')); legendStr_B = cellstr(num2str(moleFraction_B'));
    
    % Loop over all the mole fractions to plot the isotherm
    for moleFractionStep = 1:numberOfMoleFractions
        % Plot the isotherm in a semilog plot for component A at the given
        % mole fraction
        subplot(1,2,1); hold on;
        plot(localAdsorbent.plotting(moleFractionStep).pressure,...
            localAdsorbent.plotting(moleFractionStep).solidPhaseLoading_A,...
            'Color',colorCode(moleFractionStep,:),'LineWidth',1);
        
        % Plot the isotherm in a semilog plot for component B at the given
        % mole fraction
        subplot(1,2,2); hold on;
        plot(localAdsorbent.plotting(moleFractionStep).pressure,...
            localAdsorbent.plotting(moleFractionStep).solidPhaseLoading_B,...
            'Color',colorCode(moleFractionStep,:),'LineWidth',1);
    end
    
    % Plot the path taken by the adsorbent for the different steps in the
    % cycle for component A
    subplot(1,2,1); hold on;
    % Blowdown/Evacuation step
    plot(localAdsorbent.outputBlowEvac(:,1),localAdsorbent.outputBlowEvac(:,3),...
        'Color',colorCode(moleFractionStep,:),'Color','k','LineWidth',1.5);
    % Pressurization step
    plot([localAdsorbent.outputBlowEvac(end,1) localAdsorbent.simInfo.pressureHigh],...
        [localAdsorbent.outputBlowEvac(end,3) localAdsorbent.outputPress(1,2)],...
        'Color','k','LineWidth',1.5);
    % Adsorption step
    plot([localAdsorbent.simInfo.pressureHigh localAdsorbent.simInfo.pressureHigh],...
        [localAdsorbent.outputPress(1,2) localAdsorbent.outputAds(end,3)],...
        'Color','k','LineWidth',1.5);
    % Mark the point corresponding to the state alpha in the orignal
    % manuscript and have the state written as text above the point
    scatter(localAdsorbent.simInfo.pressureHigh,localAdsorbent.outputAds(end,3),...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
    text(0.8*localAdsorbent.simInfo.pressureHigh,1.10*localAdsorbent.outputAds(end,3),...
            char(945),'Color','red','FontSize',10,'FontWeight','bold')
    % Mark the point corresponding to the state gamma in the orignal
    % manuscript and have the state written as text above the point    
    scatter(localAdsorbent.outputBlowEvac(end,1),localAdsorbent.outputBlowEvac(end,3),...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
    text(0.8*localAdsorbent.outputBlowEvac(end,1),1.15*localAdsorbent.outputBlowEvac(end,3),...
            char(947),'Color','red','FontSize',10,'FontWeight','bold')
    % Mark the point corresponding to the state delta in the orignal
    % manuscript and have the state written as text above the point
    scatter(localAdsorbent.simInfo.pressureHigh,localAdsorbent.outputPress(1,2),...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
    text(0.8*localAdsorbent.simInfo.pressureHigh,1.2*localAdsorbent.outputPress(1,2),...
            char(948),'Color','red','FontSize',10,'FontWeight','bold')
    % Generate the labels, legends and do some graphic modifications
    xlabel('Total Pressure {\it{P}} [bar]'); ylabel('Solid Phase Loading {\it{q}}_A [mol/kg]'); 
    title('Component A');
    lgd = legend (legendStr_A,'Location','Northwest','NumColumns',2);
    title(lgd,'{\it{y}}_{A}');
    set(gca,'xscale','log'); axH = gca; axH.LineWidth = 1; axH.FontName = 'Helvetica'; axH.FontSize = 10; 
    axH.FontWeight = 'bold'; grid on; box on; xlim([0 inf]); ylim([0 inf]);
    
    % Plot the path taken by the adsorbent for the different steps in the
    % cycle for component B
    subplot(1,2,2); hold on;
    % Blowdown/Evacuation step
    plot(localAdsorbent.outputBlowEvac(:,1),localAdsorbent.outputBlowEvac(:,4),...
        'Color','k','LineWidth',1.5);
    % Pressurization step
    plot([localAdsorbent.outputBlowEvac(end,1) localAdsorbent.simInfo.pressureHigh],...
        [localAdsorbent.outputBlowEvac(end,4) localAdsorbent.outputPress(1,3)],...
        'Color','k','LineWidth',1.5);
    % Adsorption step
    plot([localAdsorbent.simInfo.pressureHigh localAdsorbent.simInfo.pressureHigh],...
        [localAdsorbent.outputPress(1,3) localAdsorbent.outputAds(end,4)],...
        'Color','k','LineWidth',1.5);
    % Mark the point corresponding to the state alpha in the orignal
    % manuscript and have the state written as text above the point
    scatter(localAdsorbent.simInfo.pressureHigh,localAdsorbent.outputAds(end,4),...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
    text(0.8*localAdsorbent.simInfo.pressureHigh,1.05*localAdsorbent.outputAds(end,4),...
            char(945),'Color','red','FontSize',10,'FontWeight','bold')
    % Mark the point corresponding to the state gamma in the orignal
    % manuscript and have the state written as text above the point
    scatter(localAdsorbent.outputBlowEvac(end,1),localAdsorbent.outputBlowEvac(end,4),...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
    text(0.8*localAdsorbent.outputBlowEvac(end,1),1.05*localAdsorbent.outputBlowEvac(end,4),...
            char(947),'Color','red','FontSize',10,'FontWeight','bold')
    % Mark the point corresponding to the state delta in the orignal
    % manuscript and have the state written as text above the point
    scatter(localAdsorbent.simInfo.pressureHigh,localAdsorbent.outputPress(1,3),...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
    text(0.8*localAdsorbent.simInfo.pressureHigh,1.10*localAdsorbent.outputPress(1,3),...
            char(948),'Color','red','FontSize',10,'FontWeight','bold')
    % Generate the labels, legends and do some graphic modifications
    xlabel('Total Pressure {\it{P}} [bar]'); ylabel('Solid Phase Loading {\it{q}}_B [mol/kg]'); 
    title('Component B');
    lgd = legend (legendStr_B,'Location','Northwest','NumColumns',2);
    title(lgd,'{\it{y}}_{B}');
    set(gca,'xscale','log'); axH = gca; axH.LineWidth = 1; axH.FontName = 'Helvetica'; axH.FontSize = 10; 
    axH.FontWeight = 'bold'; grid on; box on; xlim([0 inf]); ylim([0 inf]);
    
    % Get the figure handle and change the title of the figure to the name
    % of the adsorbent
    fig = gcf; fig.Name = eval(['outputStruct.ADS_',num2str(adsorbentNumber),'.adsInfo.adsorbentName']);
end
end