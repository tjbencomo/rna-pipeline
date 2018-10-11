# RNA-Seq Pipeline Guide
This file contains Tomas's notes from building the pipeline as well as general install tips. 
## Dependencies
Before using the pipeline, several programs must be properly configured.
### Languages
Python and Perl must be available to run. Load Python and Perl with the following commands:
```
ml python/3.6.1
ml perl
ml R
```
This loads Python and Perl from the sherlock shared software, placing them in your $PATH variable. It is recommended you place these commands in `/.profile` so that they autoload when your session starts.
### Programs
#### RSEM
RSEM is already installed in the `$PI_HOME` directory. All shared software packages for the lab can be found in `$PI_HOME/software`. Add the following command to `/.profile` so that RSEM is accessible from startup:
```
export PATH=$PATH:/home/groups/carilee/software/RSEM/
```
#### STAR
STAR is available from Sherlock's shared software. Use the following command to load STAR:
```
ml biology star
```
It is recommended that you also add this command to `/.profile` so STAR is available on startup.
#### HTSeq-Count
HTSeq-Count is a Python program and must be installed from pip. 
```
pip3 install --user HTSeq
```
This installs `HTSeq` into the users local python3 library. To make `HTSeq` accessible from anywhere:
```
export PATH=$PATH:/home/users/$USER/.local/bin
```
`$USER` is your sherlock username
### Reference Files
`rna-pipeline` has two options to handle the reference files. The user can manually specify reference files at runtime through command line flags or they can use a `config.ini` file to set the necessary reference file paths once and allow the program to automatically find them.  
To make `config.ini` first create the file inside the `/rna-pipeline` directory:
```
cd rna-pipeline/
touch config.ini
```
Then edit the file so that it looks similar to this:
```
[python]
genomeDirectory=GENOME DIRECTORY HERE
rsemDirectory=RSEM DIRECTORY HERE
htseqGFFFFile=GFF FILE HERE
```

