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
			    Readme last updated: 25th January 2018, AK

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
2. When using for the first time, execute the `runBAAM.m` file, which would throw an error and would create a `BAAMInfo` file. The details regarding the simulation settings and adsorbent properties can be filled out in the file. If the user wishes to load data from a `.mat` or `.xlsx` file then provide the name of the file in the `BAAMInfo` file. Read below for more details. Then executed again the `runBAAM.m` file.


## Authors

### Maintainers of the repository
* Vishal Subramanian Balashankar (vishal3@ualberta.ca)

### Project contributors
* Prof. Dr. Arvind Rajendran (arvind.rajendran@ualberta.ca)
* Ashwin Kumar Rajagopalan
* Dr. Ruben de Pauw
* Dr. Adolfo M. Avila