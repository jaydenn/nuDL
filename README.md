# nuDL v1.1
> Dec 2022

nuDL is a tool for the calculation of coherent neutrino-nucleus scattering rates and derivation of discovery limits on BSM physics.

## Install
Dependancies:
GSL (the GNU Science Library)
If you do not have the GSL package installed, install it by typing (on Debian systems):
```
sudo apt-get install libgsl0-dev
```
or otherwise try:
```
sudo yum install gsl-devel
```
Optional dependancy:
git

Either download the code from GitHub and extract it from the archive, or clone the git repository by typing:
```
git clone https://github.com/jaydenn/nuDL.git
```
future updates can be obtained by typing:
```
git pull origin master
```
Once GSL is installed you can build nuDL by typing 'make' in the main directory.

## Running nuDL
There are two ways to run nuDL:
1. Run in interactive mode by typing
```
./nuDL -i
```
In this mode the program will ask you questions about what you would like it to compute.
    
2. Run in config file mode by typing
```
./nuDL
```
or
```
./nuDL -c configFile.dat
```
In this mode the program will compute based on information specified in the configuration file. Typing "./nuDL" runs with the default config.dat file. Alternatively you can use the -c flag to specify your own configuration file (useful for running batch jobs), see the config.dat file for the required format.

## Output files
The output files will be placed in the location specified by the root variable, the default is `./results/CN_`. All output files will then by found in the results directory with filenames prepended by CN_

Output filenames are created dynamically and will contain information about the type of computation completed. The files will also contain a header line(s) with more information. 

## Detectors
The detectors are specified in the file 'detectors.ini', feel free to modify these or add your own detector (just give it a unique name), it can then be refered to by name in the config.dat file.
