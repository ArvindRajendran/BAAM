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
					    Readme last updated: 01st February 2019, AK

## FILE LIST
1. `runBAAM.m` : Main wrapper function 
2. `inputParserBAAM.m` : Parser function to obtain simulation settings and adsorbent properties
3. `loadAdsorbentInfoFromFile.m` : Loads adsobent properties from a .mat or .xslx file
4. `checkOptionsBAAM.m` : Performs a sanity check on the user data
5. `BAAM.m` : Main function that has the BAAM model
	

## INSTALLATION

### Dependencies

The following dependencies are required for the proper execution of this program.

1. MATLAB 2018b [required]


### Installation instructions

1. Clone the full software package from the GITHUB server into the preferred installation directory (e.g. Desktop, My Documents). The command is as follows:
```sh
git clone https://github.com/ArvindRajendran/BAAM.git

```
2. When using the software for the first time, in MATLAB command window execute the following command
execute the following command
```sh
runBAAM

```
This would throw an error. After which it would create a file called `BAAMInfo`. The details of the simulation settings and adsorbent properties can be filled out in the file. 
3. The user can also load adsorbent properties from a `.mat` or `.xlsx` file. In this case, provide the name of the file which contains the adsorbent properties in the `BAAMInfo` file. Read below for more details. 
4. After filling out the necessary information `BAAMInfo` file, execute the `runBAAM.m` file again.

## HOW TO USE THE SOFTWARE?
The `BAAMInfo` file consists multiple options. A brief description of each of these options is listed below.

### SIMULATION INFO
* `molFracFeed_A`: Mole fraction of the heavy component in the gas phase [0-1]. In the article this was set to 0.15.
* `molarMass_A`: Molar mass of the heavy component [g/mol]. In the article this was set to 44.01.
* `temperature`: Temperature of the feed gas [K]. In the article this was set to 298.15.
* `pressureHigh`: High pressure in the column. This would normally correspond to the pressure in the adsorption step [bar]. In the article this was set to 1.00.
* `pressureLow`: Low pressure in the column. This would normally correspond to the pressure in the evacuation step [bar]. In the article this was set to 0.03.
* `pressureLowUpperBound`: This corresponds to the upper bound for the low pressure range that would be simulated (See Fig. 5 in the original article) [bar]. In the article this was set to 0.10 bar.
* `voidFraction`: Void fraction of the column [0-1]. In the article this was set to 0.37.
* `adiabaticConstant`: Adiabatic constant of the gas to be simulated [>1]. In the article this was set to 1.4.
* `pumpEfficiency`: The efficiency of the vacuum pumps used for the blowdown and evacuation step [0-1].  In the article this was set to 0.72.
* `plotFlag`: A boolen flag which would enable the plotting of necessary figures after each simulation run (Fig. 3 and 5) [Yes/No].
* `molFracPlotting`: A vector of mole fractions that would be used for plotting a figure similar to Fig. 3 in the original article. The mole fractions can be separated using a comma or a semicolon. For example: 0.15, 0.5, 0.75, 0.9 or 0.15; 0.5; 0.75; 0.9

### ADSORBENT INFO
* `loadAdsorbentInfo`: 

## ADDITIONAL COMMENTS


## AUTHORS

### Maintainers of the repository
* Vishal Subramanian Balashankar (vishal3@ualberta.ca)

### Project contributors
* Prof. Dr. Arvind Rajendran (arvind.rajendran@ualberta.ca)
* Ashwin Kumar Rajagopalan
* Dr. Ruben de Pauw
* Dr. Adolfo M. Avila