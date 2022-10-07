%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Simulation files
Miguel Parada Contzen
Conce, October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The files contained in this folder are for computing the examples in the paper

	Parada Contzen, Miguel. 'Medium voltage level power control using distributed PV generation and hydrogen storage.', 2022.

The data files describe a 116 nodes balanced distribution circuit whose topology remains that of the IEEE 123 nodes benchmark. The nominal voltage is assummed to be V = 15[kV], although it can be changed in the Matlab script files. The original value of 4.16[kV] was discarded for the paper as it implies several voltage droop problems that scape from the scope of the study. Loads and DG is not based on the benchmark information but is rather arbitrary.

Loads are estimated from daily profiles for commercial, industrial, and residential classes, and a nominal arbitrary maximum power for summer consumption. The ZIP parameters of the load are those stated in Kundur, 'Power System Stability and Control', 1994, Ch. 7.

Photo-voltaic injection at the nodes is simulated from an arbitrary nominal power, and typical daily profiles of global irradiation for every month. The profiles are for the Carriel Sur Airport (CCP), Concepción, Chile (latitude -36.78055, longitude -73.05083).

The following is the detailed description of the files contained.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data Files:
	- nodes.csv: Information of the 116 nodes of the system. Separator=semicolon. Each row corresponds to a node. The columns correspond to:
		- n: number identifaction of the node.
		- x: x-position for drawing the circuit.
		- y: y-position for drawing the circuit.
		- Pmax[kW]: maximum total load at the load in summer in kW.
		- res pc: residential share in per cent for the load mix at the node.
		- com pc: commercial share in per cent for the load mix at the node.
		- ind pc: industrial share in per cent for the load mix at the node.
		- Ppv[kW]: Maximum power by Photo-voltaic DG in kW at the node.
		- Peol[kW]: Maximum power by Wind power DG in kW at the node.
		- Shunt[kVAr]: Shunt capacitors compensation in kVAr at any moment.
		
	- branches.csv: Information of the branches in the circuit. Separator=semicolon. Each row corresponds to a branch. The columns correspond to:
		- n: origin node of the branch.
		- m: destiny node of the branch.
		- R: 3-phase equivalent balanced resistance of the branch in Ohm.
		- X: 3-phase equivalent balanced inductance of the branch in Ohm at 50Hz.
		
	- LoadProfiles.csv: Daily load profiles for different classes. Separator=semicolon. Each row corresponds to an hour. The columns correspond to:
		- hr: Hour of the day (0 to 23).
		- Com: commercial load in per unit of the maximum class load at the time.
		- Res: residential load in per unit of the maximum class load at the time.
		- Ind: industrial load in per unit of the maximum class load at the time.

	- PerfilEolCarrielSur.csv: Daily wind profiles. Separator=comma. 
		The raw data was obtained from https://climatologia.meteochile.gob.cl/application/informacion/fichaDeEstacion/360019 in november 2020, and comprehends data from 1966 to the date. It was afterwards proccessed to present it this file. 
		Each row corresponds to an hour of a month, in such a way that the first 24 rows correspond to hours 0 to 23 of January, rows 25 to 48 correspond to February, and so on. The columns correspond to:
		- First column: Mean wind direction (0 to 360 degrees) for the hour in the month.
		- Second column: Mean wind speed in [m/s] for the hour in the month.
		- Third column: Standard deviation of wind speed in [m/s] with respect to mean value in the previous column.

	- PerfilSolCarrielSur.csv: Daily irradiation profiles. Separator=comma.
		The raw data was taken from the Explorador Solar website (http://solar.minenergia.cl/inicio) of the Chilean Ministery of Energy and comprehends simulated solar irradiation data from 2004 to 2016. The data was afterwards processed to present in this file.
		Each row corresponds to an hour of a month, in such a way that the first 24 rows correspond to hours 0 to 23 of January, rows 25 to 48 correspond to hours 0 to 23 of February, and so on. The columns correspond to:
		- First column: Mean global radiation in [W/m^2] for the hour in the month.
		- Second column: Mean direct radiation in [W/m^2] for the hour in the month.
		- Third column: Mean difuse radiation in [W/m^2] for the hour in the month. 
		- Fourthcolumn: Mean difuse reflected radiation in [W/m^2] for the hour in the month.
		- Fifth column: Mean horizontal global radiation in [W/m^2] for the hour in the month.
		- Sixth column: Mean normal direct radiation in [W/m^2] for the hour in the month.
		- Seventh column: Mean normal difuse radiation in [W/m^2] for the hour in the month.
		- Eighth column: Mean normal direct radiation in [W/m^2] for the hour in the month.
		- Ninthcolumn: Mean Temperature in [C degree]
		- Tenth column: Mean Wind speed at 2[m] in [m/s]
		- Eleventh column: Geo. shadows probability for the hour in the month.
		- Twelfth column: Clouds probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Matlab simulation scripts:
	- Examplespaper.m: Closed and open loop simulation of described system and power control algorithm.