Chicken task analysis code

Written by Katie Forthman, 4/10/19

**********************************************************************

Step 1: Formatting Raw Data

If there are raw data that need to be formatted for analysis, start by running the script 'formatRaw.m'. If needed, the path from the file where the data are located and the path to the file where the formatted data will be stored may be altered. If the # of blocks/trials are changed in future versions of the task, this must also be corrected in the 'formatRaw.m' file. The formatted data will be saved to the folder Data/Raw, one mat file for each participant. Each mat file will contain 12 data structures named dataTXBY (where X is the time-point and Y is the block number), one containing the data for each block and each time-point.

The resulting raw data files will contain the following information:

	ID 	:= Participant ID

	trial 	:= The trial, either T1 or T2. Heretofore referred to as time-points.

	block 	:= The block number. There are 6 possible

	N 	:= The number of trials (number of instances the participant is shown the egg in a different position and are asked to make an estimation)

	sigma 	:= The distribution of eggs around the chicken

	muall 	:=  raw [x,y] coordinates of the two chickens (left and right). Center of the screen = left chicken X + ((right chicken X - left chicken X)/2)

	X 	:= [x,y] coordinates of each egg on each trial

	pred 	:= the chicken the participant selected

	muinds 	:= The actual chicken that generated the egg (1 for left, 2 for right)

	rt 	:= reaction time

	H 	:= current true hazard rate

	signaled := whether the correct was revealed or not at the end of each trial

	cp 	:= change-point. 1 if the correct chicken is different from the previous trial. 0 otherwise.

	r 	:= number of trials since last change-point

	Hset 	:= Lists the hazard rate of each block for this particular pattern

Aside: It is a bit confusing, but 'trial' refers to two things. First, trial refers to the time-point at which the participant did the task, either T1 or T2. These time-points are spaced 24 hrs apart. To avoid confusion, I try to refer to T1 and T2 as time-points rather than trials. Second, trial also refers to each instance the participant is shown the egg in a different position and are asked to make an estimation. You will have to go by context to decide the meaning. Use the definitions above for clarification.

**********************************************************************

Step 2: Running the Analysis

Start by opening 'getDataInfo.m' and making sure the path to the directory is correct. 
Next, run 'fitAdaptivityModel.m'. This is the primary analysis file. This script will output in the pathTo directory both a csv and mat file that contain the same table. This table lists ID, trial, block, sigma, and H_true (same as defined above) and the model estimates for 
	H_subjective	:=  subjective hazard rate
	noise_in_DV	:=  noise in the decision variable
	lapse_rate	:=  lapses in attention
These estimates are calculated for each participant for each block.
(To see how noise_in DV and lapse_rate affect belief, please see choicehatProbability.nb)
The table also includes
	pct		:= percent correct
	

The optimization function is 'fitAdaptivityModel_err.m'. This is the function that fmincon will attempt to reduce.

**********************************************************************

Step 3: Creating the Figures

Figures are created by the R file 'createFigures.R' in the main directory. The resulting figures are saved into the folder 'Figures'.

**********************************************************************

Step 4: Examining Stability of the Model

By the nature of the optimization function, there is the possibility that the estimation is not the ideal estimation. We needed to therefore test the stability of the model based on different 'seed' values, or initial values. The estimations for different initial values were attained by the script 'Models/fitAdaptivityModel_stabilityTest.m'. The results of this test are examined with the file 'H_stability.R'.