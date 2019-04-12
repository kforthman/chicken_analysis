# ANALYSIS FOLDER

Required to run fitting algorithms for subject data - these variables are called by "getDataInfo.m"

As used in Glaze et al. 2018 - A bias-variance trade-off governs individual differences in on-line learning in an unpredictable environment.

Questions should be directed to Alex Filipowicz, alsfilip@gmail.com

This folder contains data structures with subject and session information used for fitting.

- data_info: Structure containing the following matrices:
	- fileList_: names of all of the subject data files
	- Lgood: index for subjects that should be included in the analysis (1) or excluded (0) (based on the sessions they completed - any subject with fewer than 2 sessions were excluded)
	- sessionsPerSubject: Number of sessions completed by each subject

