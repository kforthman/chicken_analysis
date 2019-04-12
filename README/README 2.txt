# RAW DATA FOLDER

Folder with raw participant data files

As used in Glaze et al. 2018 - A bias-variance trade-off governs individual differences in on-line learning in an unpredictable environment.

Questions should be directed to Alex Filipowicz, alsfilip@gmail.com

Each data file , data_XX.mat, contains a structure for each subject, where XX is the subject ID. The names of the variables aren't super intuitive so here's a key to each subject's structure:

Sub structures:
- dataX (where X can range from 1-6): Data structure with performance/stimuli in each session, where X indicates the sessions number. Inside each dataX structure you'll find the following variables:
	- N: total trials for the session
	- signaled: whether the true mean was revealed or not at the end of each trial
	- tau:
	- sigma: 
	- n: number of trials (redundant)
	- muall: raw [x,y] coordinates of the two triangle means (left and right). Center of the screen = left triangle X + ((right triangle X - left triangle X)/2)
	- X: [x,y] coordinates of each star on each trial
	- pred: The actual triangle that generated the star (1 for left, 2 for right)
	- cp: whether or not the star appears to have come from the same triangle as the previous trial. Note that this is not the ACTUAL change-point of the triangle
	- r: trial since last change-point as defined above
	- muinds: the triangle the star appears to have come from (not which one it actually came from)
	- rt: reaction times
	- H: current hazard rate
	- Hset: Set of hazard rates for the current session (either two hazard rates - 1000 trials each, or 5 hazard rates - 400 trials each)
	
	-NOTE that in a few subjects a bug in the code generating task parameters stored the inverse of task hazard rate in data.H. All analysis code here checks for that and corrects if present.
	
- generate_condition_ids: script to generate subject trials and conditions

	