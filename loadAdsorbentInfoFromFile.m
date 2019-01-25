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
% This function is used to load the adsorbent properties from either a .mat
% file or from a .xlsx file. The function then extracts the data and then
% assigns them to the proper variables such that the output from this
% functino can be directly fed to BAAM. Note that this only modifies the
% adsInfo structure and does not modify the simInfo. 
%
% Last modified:
% - 2019-01-25, AK: Introduced header for this file
%
% Input arguments:
% - fileName:       The filename from which the adsorbent properties needs
%                   to be loaded
%
% Output arguments:
% - adsInfo:        Structure containing information about the adsorbent.
%                   This contains the adsorbent name, adsorbent density and
%                   DSL parameters for component A and B. In the original
%                   manuscript A is CO2 and B is N2.
%
% Comments:
% If .xlsx file is to be used and if the user is not aware of the format,
% just run the main script and that would generate the .xslx file. Open the
% file and fill the necessary information, and input the name in the
% BAAMInfo file.
%
% FOR EXCEL:
% Variables should be ordered in the following way: The columns containts
% the different variables and the rows contain the information for each
% adsorbent. Column 1-14: 1. Adsorbent name, 2. adsorbent density, 
% 3. saturation capacity of component A in site 1, 4. adsorption
% coeffecient of component A in site 1, 5. internal energy of component A
% in site 1, 6. saturation capacity of component A in site 2, 7. adsorption
% coefficeent of component A in site 2, 8. internal energy of component A
% in site 2, 9. saturation capacitiy of component B in site 1, 10.
% adsorption coffecient of component B in site 1, 11. internal energy of
% component B in site 1, 12. saturation capacity of component B in site 2,
% 13. adsorption coefficient of component B in site 2, 14. itnernal energy
% of component B in site 2
%
% FOR MATLAB:
% Create a vector of string (numberOfAdsorbent x 1) with the name of the
% adsorbents. Create a matrix in the same order as column 2-14 in the above
% case. The size of this matrix should be (numberOfAdsorbent x 13). The
% vector of adsorbent names should be called adsorbentList and the matrix
% of the adsorbent properties should be called adsorbentProperties
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adsInfo = loadAdsorbentInfoFromFile(fileName)
% Check if the file is an excel file or a .mat file
excelFlag = strfind(fileName,'.xlsx');
% Check if the filename is provided by user nor not or if the file is .xslx
% format or not
if isempty(fileName) || ~isempty(excelFlag)
    % If the excel file exsits then load the data directly
    if exist(fileName) == 2
        [~,~,loadedExcelData] = xlsread(fileName);
	% If the excel file does not exist then create a new excel which the
	% user can use to input the data
    else
        % Column names
        columnNames = {'Adsorbent Name','Adsorbent Density',...
            'Saturation Capacity of Component A (Site 1)',...
            'Adsorption Coefficient of Component A (Site 1)',...
            'Internal Energy of Component A (Site 1)',...
            'Saturation Capacity of Component A (Site 2)',...
            'Adsorption Coefficient of Component A (Site 2)',...
            'Internal Energy of Component A (Site 2)',...
            'Saturation Capacity of Component B (Site 1)',...
            'Adsorption Coefficient of Component B (Site 1)',...
            'Internal Energy of Component B (Site 1)',...
            'Saturation Capacity of Component B (Site 2)',...
            'Adsorption Coefficient of Component B (Site 2)',...
            'Internal Energy of Component B (Site 2)'};
        % Save the file as adsorbentSheet.xlsx and throw an error at the
        % user
        xlswrite('adsorbentSheet.xlsx',columnNames)
        error('BAAM:AdsorbentInfoMissing',...
            'Could not find an adsorbent info file in the current folder. An empty excel file has been created. Please fill in the fields, update BAAMInfo and rerun this script.');
    end
    % Loop through all the available data and convert it into the desirable
    % form
    for numRows = 2:size(loadedExcelData,1)
        adsorbentList(numRows-1,1) = string(loadedExcelData{numRows,1});
        for numColumns = 2:size(loadedExcelData,2)
            adsorbentProperties(numRows-1,numColumns-1) = loadedExcelData{numRows,numColumns};
        end
        % If any of the entries is NaN initalize it to zero
        adsorbentProperties(find(isnan(adsorbentProperties))) = 0;
    end
else
    % If it is a .mat file do a direct load from the fileName provided in
    % the BAAMInfo file
    load(fileName);
end
% Loop over all the available adsorbents and create a structure of arrays
% with all the properties in the form that would be compatible with BAAM
numberOfAdsorbents = size(adsorbentList,1);
for ii = 1:numberOfAdsorbents
    adsInfo(ii).adsorbentName = adsorbentList(ii,1);
    adsInfo(ii).adsorbentDensity = adsorbentProperties(ii,1);
    adsInfo(ii).qSaturationSite1_A = adsorbentProperties(ii,2);
    adsInfo(ii).adsorptionCoefficientSite1_A = adsorbentProperties(ii,3);
    adsInfo(ii).internalEnergySite1_A = adsorbentProperties(ii,4);
    adsInfo(ii).qSaturationSite2_A = adsorbentProperties(ii,5);
    adsInfo(ii).adsorptionCoefficientSite2_A = adsorbentProperties(ii,6);
    adsInfo(ii).internalEnergySite2_A = adsorbentProperties(ii,7);
    adsInfo(ii).qSaturationSite1_B = adsorbentProperties(ii,8);
    adsInfo(ii).adsorptionCoefficientSite1_B = adsorbentProperties(ii,9);
    adsInfo(ii).internalEnergySite1_B = adsorbentProperties(ii,10);
    adsInfo(ii).qSaturationSite2_B = adsorbentProperties(ii,11);
    adsInfo(ii).adsorptionCoefficientSite2_B = adsorbentProperties(ii,12);
    adsInfo(ii).internalEnergySite2_B = adsorbentProperties(ii,13);
end
end