						██████╗  █████╗  █████╗ ███╗   ███╗
						██╔══██╗██╔══██╗██╔══██╗████╗ ████║
						██████╔╝███████║███████║██╔████╔██║
						██╔══██╗██╔══██║██╔══██║██║╚██╔╝██║
						██████╔╝██║  ██║██║  ██║██║ ╚═╝ ██║
						╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝
				      		   BATCH ADSORBER ANALOGUE MODEL
                                   
	
						© 2019 University of Alberta, Canada

					    Laboratory of Advanced Separation Processes
						     Prof. Dr. Arvind Rajendran
							Project: 
					    Readme last updated: 03rd February 2019, AK

## INTRODUCTION
The function BAAM simulates the batch adsorber analogue model proposed in *Analysis of a Batch Adsorber Analogue for Rapid Screening of Adsorbents for Postcombustion CO2 Capture, doi*. This code simulates a two component system using the LPP cycle with the BAAM model. Component A is the component of interest (CO2 in the original article).

## FILE LIST
1. `runBAAM.m` : Main wrapper function 
2. `inputParserBAAM.m` : Parser function to obtain simulation settings and adsorbent properties
3. `loadAdsorbentInfoFromFile.m` : Loads adsobent properties from a .mat or .xslx file
4. `checkOptionsBAAM.m` : Performs a sanity check on the user data
5. `BAAM.m` : Main function that has the BAAM model
6. `plotBAAMOutput`: Plot function that plots Fig. 3 in the original article.
	

## INSTALLATION

### Dependencies

The following dependencies are required for the proper execution of this program.

1. MATLAB 2018b *[required]*


### Installation instructions

1. Clone the full software package from the GitHub server into the preferred installation directory (e.g. Desktop, My Documents). The command is as follows:
```sh
git clone https://github.com/ArvindRajendran/BAAM.git

```
This step would require installation of Git (https://git-scm.com/downloads).
2. Go to folder called BAAM which would be extracted by the above command from GitHub. When using the software for the first time, in MATLAB command window execute the following command
```sh
runBAAM

```
This would throw an error. After which it would create a file called `BAAMInfo`. The details of the simulation settings and adsorbent properties can be filled out in the file. 
3. The user can also load adsorbent properties from a `.mat` or `.xlsx` file. In this case, provide the name of the file which contains the adsorbent properties in the `BAAMInfo` file. Read below for more details. 
4. After filling out the necessary information in the `BAAMInfo` file, execute the command in point 2 again.

## HOW TO USE THE SOFTWARE?
The `BAAMInfo` file consists multiple options. A brief description of each of these options is listed below.

### INPUT: SIMULATION INFO
* `molFracFeed_A`: Mole fraction of the heavy component in the gas phase [0-1]. In the article this was set to 0.15.
* `molarMass_A`: Molar mass of the heavy component [g/mol]. In the article this was set to 44.01.
* `temperature`: Temperature of the feed gas [K]. In the article this was set to 298.15.
* `pressureHigh`: High pressure in the column. This would normally correspond to the pressure in the adsorption step [bar]. In the article this was set to 1.00.
* `pressureLow`: Low pressure in the column. This would normally correspond to the pressure in the evacuation step [bar]. In the article this was set to 0.03.
* `pressureLowUpperBound`: This corresponds to the upper bound for the low pressure range that would be simulated (See Fig. 5 in the original article) [bar]. In the article this was set to 0.10 bar.
* `voidFraction`: Void fraction of the column [0-1]. In the article this was set to 0.37.
* `adiabaticConstant`: Adiabatic constant of the gas to be simulated [>1]. In the article this was set to 1.4.
* `pumpEfficiency`: The efficiency of the vacuum pumps used for the blowdown and evacuation step [0-1].  In the article this was set to 0.72.
* `plotFlag`: A boolen flag which would enable the plotting of necessary figures after each simulation run (Fig. 3) [Yes/No].
* `molFracPlotting`: A vector of mole fractions that would be used for plotting a figure similar to Fig. 3 in the original article. The mole fractions can be separated using a comma or a semicolon. For example: 0.15, 0.5, 0.75, 0.9 or 0.15; 0.5; 0.75; 0.9

### INPUT: ADSORBENT INFO
* `loadAdsorbentInfo`: A boolean flag which would enable adsorbent isotherms and properties from a file [Yes/No]. This would be useful for screening multiple adsorbents.
* `filenameAdsorbentInfo`: Filename of the file that consists the adsorbent properties. The file should be present in the BAAM folder. The file can be a .mat file or a .xlsx file. For further details read additional comments.
* In case `loadAdsorbentInfo` flag is No, then the user can input the adsorbent properties using the below mentioned variables. The isotherm used is a Dual-Site Langmuir.
* `adsorbentName`: Name of the adsorbent that the user which to evaluate using BAAM.
* `adsorbentDensity`: Particle density of the adsorbent [kg/m3].
* `qSaturationSite1_A`: Saturation capacity of component A in site 1 [mol/kg]. Ref. to Eq. 19 (qsb) in the original article.
* `adsorptionCoefficientSite1_A`: Adsorption coeffecient of component A in site 1 [m3/mol]. Ref. to Eq. 19 (b0) in the original article.
* `internalEnergySite1_A`: Internal energy of component A in site 1 [J/mol]. Ref. to Eq. 19 (Ub) in the original article.
* `qSaturationSite2_A`: Saturation capacity of component A in site 2 [mol/kg]. Ref. to Eq. 19 (qsd) in the original article.
* `adsorptionCoefficientSite2_A`: Adsorption coeffecient of component A in site 2 [m3/mol]. Ref. to Eq. 19 (d0) in the original article.
* `internalEnergySite2_A`: Internal energy of component A in site 2 [J/mol]. Ref. to Eq. 19 (Ud) in the original article.
* `qSaturationSite1_B`: Saturation capacity of component B in site 1 [mol/kg]. Ref. to Eq. 19 (qsb) in the original article.
* `adsorptionCoefficientSite1_B`: Adsorption coeffecient of component B in site 1 [m3/mol]. Ref. to Eq. 19 (b0) in the original article.
* `internalEnergySite1_B`: Internal energy of component B in site 1 [J/mol]. Ref. to Eq. 19 (Ub) in the original article.
* `qSaturationSite2_B`: Saturation capacity of component B in site 2 [mol/kg]. Ref. to Eq. 19 (qsd) in the original article.
* `adsorptionCoefficientSite2_B`: Adsorption coeffecient of component B in site 2 [m3/mol]. Ref. to Eq. 19 (d0) in the original article.
* `internalEnergySite2_B`: Internal energy of component B in site 2 [J/mol]. Ref. to Eq. 19 (Ud) in the original article.

### OUTPUT:
* After the successful execution of the BAAM for the desired number of adsorbents, the runBAAM function would create a `.mat` file with the following filename `BAAMOutput_<7 digit GIT commit ID>_<Date of the simulation: ddmmyyyy>_<Time of the simulation: hhMMss>`. The file contains all the necessary output from BAAM. If the user is unaware of the fields in the output, please refer to `BAAM.m` where a detailed explanation is provided for each output field.
* Apart from the `.mat` file, a `.txt` file called `BAAM.txt` is also generated which provides condensed information with information as to the maximum distance of the Pu/Re plot (Eq. 22 in the original article), if the adsorbent satisfied the Pu/Re constraints (Section 5.1 in the original article), energy consumption from BAAM and scaled energy consumption from full model (Eq. 25 in the original article), and the working capacity at minimum energy. 

## ADDITIONAL COMMENTS
* In case, the user is unsure of simulation settings, default values would be imposed on the simulation. These default values are the same as the ones used in the original manuscript. The only requirement in such a case is to load adsorbent properties from a file.
* Two adsorbent property files are provided to guide the user to get used to the software. A .mat file and a .xlsx file with the same adsorbent is provided. The example adsorbents are given in Section 4 of the original article. The name of either of these files can be input in `filenameAdsorbentInfo` along with the file extension.
* In case, the user wants to add more adsorbents, the user may edit the existing .xlsx or .mat file. If the user does not wish to do so, `filenameAdsorbentInfo` can be left empty, with `loadAdsorbentInfo` as Yes. This would create a .xlsx file that contains the necessary fields, and the user can specify as many adsorbents as required.
* This model is written in a general fashion for a two component system using the 4 step VSA cycle with LPP, but has been tested extensively only for the CO2/N2 system.

## REFERENCES
 
Analysis of a Batch Adsorber Analogue for Rapid Screening of Adsorbents for Postcombustion CO2 Capture
Vishal Subramanian Balashankar, Ashwin Kumar Rajagopalan, Ruben de Pauw, Adolfo M. Avila, Arvind Rajendran
<Citation>
<DOI>

## AUTHORS

### Maintainers of the repository
* Vishal Subramanian Balashankar (vishal3@ualberta.ca)

### Project contributors
* Prof. Dr. Arvind Rajendran (arvind.rajendran@ualberta.ca)
* Ashwin Kumar Rajagopalan
* Dr. Ruben de Pauw
* Dr. Adolfo M. Avila

## DISCLAIMER
The software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. The authors do not take any responsibility for correctness or completeness of the software provided.