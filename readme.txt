NGSANE is a framework for advanced production informatics of Next Generation Sequencing libraries.

Version: 0.0.1

NGSANE is structured as follows:

ngsane
- bin/ 
	contains the trigger.sh script that is used to launch any and all jobs by supplying an appropriate config file

- conf/ 
	contains the config file for the cluster environment at hand. 
	rename or copy the sample_header.sh into header.sh and populate the variables contained within header.sh
	with the specific values on your system
	
- doc/
	contains documents describing the use of NGSANE
	
- mods/ 
	contains all modules currently available from withing NGSANE
	as well as helper scripts to submit jobs on the cluster such as 
	prepareJobSubmission.sh (formerly known as pbsTemp.sh)
	
- sampleConfigs/
	contains a sample config file for various pipelines/modules than can be triggered in the NGSANE framework
	copy the corresponding file(s) to your data and customize as approriate
	
- tools/
	contains various helper scripts used within the modules, mostly tapping into python and R
	
	
Contact:
	d.bauer@csiro.org
	f.buske@garvan.org.au
	