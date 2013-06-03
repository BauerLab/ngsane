NGSANE is a framework for advanced production informatics of Next Generation Sequencing libraries.

Version: 0.0.1

#########################################################
# Setup:

Define the environment variable NGSANE_BASE and point it to the directory that
contains the NGSANE framework.
You may want to set up $NGSANE_BASE in a cluster model or in your .bash_rc
You can also specify it in the project-based config file.

#########################################################
# Structure of NGSANE:

ngsane
- bin/ 
	contains the trigger.sh script that is used to launch any and all 
	jobs by supplying an appropriate config file

- conf/ 
	contains the config file for the cluster environment at hand. 
	rename or copy the sample_header.sh into header.sh and populate the 
	variables contained within header.sh with the specific values on your system
	
- doc/
	contains documents describing the use of NGSANE
	
- mods/ 
	contains all modules currently available from withing NGSANE
	as well as helper scripts to submit jobs on the cluster such as 
	prepareJobSubmission.sh (formerly known as pbsTemp.sh)
	
- sampleConfigs/
	contains a sample config file for various pipelines/modules than can be 
	triggered in the NGSANE framework
	copy the corresponding file(s) to your data and customize as approriate
	
- tools/
	contains various helper scripts used within the modules, 
	mostly tapping into python and R

#########################################################
# How to create a new project:

NGSANE requires a config file for each of your projects and a (set of) fastq 
files that need to be placed in a simple but predefined structure.

Make a new folder for your project and create the following folder structure:
<PROJECT>
- fastq
- - <EXPERIMENT1>
- - - <LIBRARY1_READ1>.fq[.gz]
- - - <LIBRARY1_READ2>.fq[.gz] (if paired library)
- - - <LIBRARY2_READ1>.fq[.gz]
- - - <LIBRARY2_READ2>.fq[.gz] (if paired library)
- config.txt

Here <PROJECT> is your project name. The fastq folder contains all raw data
grouped by EXPERIMENT (or any other meaningful way). Each EXPERIMENT folder
can contain multiple sequencing libraries. You probably want them gzipped to
save space. You can add additional EXPERIMENT folders if you like. However, 
all libraries should have the same ending as well as the same read-pair 
identification (e.g. have the suffix _R1.fq.gz for read1 of a zipped fastq).
for examples of config files for different pipelines look in the 
ngsane/sampleConfigs/.

#########################################################
# How to run NGSANE jobs:

The main access point is trigger.sh, which can be invoked to run in different
modes
>trigger.sh <CONFIG> <TASK>

options for TASK:
  "empty"    start dry-run: make dirs, delete old files, print what will be
done
  fetchdata  get data from remote server (via smbclient)
  pushresult puts results to remote server (via smbclient)
  armed      submit tasks to the queue
  direct     run task directly (e.g. on node after qrsh)
  verify     check the pbs logfiles for errors
  html       check the pbs logfiles for errors and and make summary HTML page
  clean      clean up

Following the folder structure example from above, you can submit jobs
using a properly configured config.txt form within the PROJECT folder as 
follows:

PROJECT>trigger.sh ./config.txt armed

#########################################################	
Contact:
	Denis.Bauer - at - csiro.au
	f.buske - at - garvan.org.au
