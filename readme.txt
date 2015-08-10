NGSANE is a framework for advanced production informatics of Next Generation 
Sequencing libraries.

Version: v0.5.2.0


################################################################################
# Associated publication

"NGSANE: A Lightweight Production Informatics Framework for High Throughput 
Data Analysis"
Buske FA, French HJ, Smith MA, Clark SJ, Bauer DC 
Bioinformatics 2014 May 15;30(10):1471-2. doi: 10.1093/bioinformatics/btu036. 
Epub 2014 Jan 26.

################################################################################
# Setup:

Define the environment variable NGSANE_BASE and point it to the directory that
contains the NGSANE framework. You may want to set up $NGSANE_BASE environment 
variable on the cluster or in your .bash_rc or .profile You can also specify it 
in the project-based config file.

For more information check out the wiki at
https://github.com/BauerLab/ngsane/wiki

################################################################################
# Structure of NGSANE:

ngsane
- bin/ 
	contains the trigger.sh script that is used to launch any and all 
	jobs by supplying an appropriate config file

- conf/ 
	contains the config file for the cluster environment at hand. 
	Rename the sample_header.sh template into header.sh and populate the 
	variables contained within with the values specific to your system

- core/ 
	contains the core scripts handling job submission, logging and report 
        generation etc.
	
- doc/
	contains documents/presentations/publication describing the framework
	
- mods/ 
	contains all modules currently available from within NGSANE
	
- sampleConfigs/
	contains a sample config file for various pipelines/modules than can be 
	triggered in the NGSANE framework
	copy the corresponding file(s) to your data and customise as appropriate
	
- tools/
	contains various helper scripts used within the modules, 
	mostly tapping into Python and R

################################################################################
# How to create a new project:

NGSANE requires a config file for each of your projects and a (set of) fastq 
files that need to be placed in a simple but predefined structure.

Make a new folder for your project and create the following folder structure:
<PROJECT>
- fastq
- - <EXPERIMENT1>
- - - <LIBRARY1_READ1>.fastq[.gz]
- - - <LIBRARY1_READ2>.fastq[.gz] (if paired library)
- - - <LIBRARY2_READ1>.fastq[.gz]
- - - <LIBRARY2_READ2>.fastq[.gz] (if paired library)
- config.txt

Here <PROJECT> is your project name. The fastq folder contains all raw data
grouped by EXPERIMENT (or any other meaningful way). Each EXPERIMENT folder
can contain multiple sequencing libraries. You probably want them gzipped to
save space. You can add additional EXPERIMENT folders if you like. However, 
all libraries should have the same ending as well as the same read-pair 
identification (e.g. have the suffix _R1.fastq.gz for read1 of a zipped fastq).
for examples of config files for different pipelines look in the 
ngsane/sampleConfigs/.

################################################################################
# How to run NGSANE jobs:

The main access point is trigger.sh, which can be invoked to run in different
modes
>trigger.sh <CONFIG> <TASK>

options for TASK:
  [empty]    start a dry-run: create folders,prints what will be done
  new        detect new data that has not been processed yet.
  fetchdata  get data from remote server (via smbclient)
  pushresult puts results to remote server (via smbclient)
  armed      submit tasks to the queue
  direct     run task directly (e.g. on node after qrsh)
  postonly   run only the post analysis steps of a task (if available)
  verify     check the pbs logfiles for errors
  recover    pick up unfinished business (interrupted jobs)
  html       check the pbs logfiles for errors and and make summary HTML page
  report     alias for html
  trackhubs      generate trackhubs
  clean          removes all dummy files

Following the folder structure example from above, you can submit jobs
using a properly configured config.txt form within the PROJECT folder as 
follows:

PROJECT>trigger.sh ./config.txt armed

################################################################################	
Contact:
	Denis.Bauer - at - csiro.au
	f.buske - at - garvan.org.au
