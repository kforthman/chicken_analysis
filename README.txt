Chicken task analysis code

Written by Katie Forthman, 4/10/19


Step 1: Formatting Raw Data
If there are raw data that need to be formatted for analysis, start by running the script 'formatRaw.m'. If needed, the path from the file where the data are located and the path to the file where the formatted data will be stored may be altered. The formatting will create 12 '.mat' files for each participant, one containing the data for each block and each trial. If the # of blocks/trials are changed in future versions of the task, this must also be corrected in the 'formatRaw.m' file.

Step 2: Running the Analysis
Start by opening 'getDataInfo.m' and making sure the path to the directory is correct. 
Next, run 'fitAdaptivityModel.m'. This is the primary analysis file. This script will output in the pathTo directory estimates for 
	- H_subjective :=  subjective hazard rate
	- noise_in_DV  :=  noise in the decision variable
	- lapse_rate   :=  lapses in attention
These estimates are calculated for each participant for each block.
(To see how noise_in DV and lapse_rate affect belief, please see choicehatProbability.nb)

The optimization function is 'fitAdaptivityModel_err.m'. This is the function that fmincon will attempt to reduce.

Step 3: Creating the Figures
Figures are created by the R file 'createFigures.R' in the main directory. The resulting figures are saved into the folder 'Figures'.

Step 4: Examining Stability of the Model
By the nature of the optimization function, there is the possibility that the estimation is not the ideal estimation. We needed to therefore test the stability of the model based on different 'seed' values, or initial values. The estimations for different initial values were attained by the script 'Models/fitAdaptivityModel_stabilityTest.m'. The results of this test are examined with the file 'H_stability.R'.