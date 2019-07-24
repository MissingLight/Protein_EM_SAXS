pdbstruct = getpdb('5CYT', 'ToFile', 'electron_transport.pdb');
%pdbstruct_Model2 = pdbread('nicotinic_receptor.pdb', 'ModelNum', 2);
molviewer(pdbstruct)
%molviewer('http://www.rcsb.org/pdb/files/2FK0.pdb.gz')

%% -- example usage: plot the atoms of 3IJU.pdb -- 

PDBdata = pdb2mat('3IJU.pdb'); % read in data from PDB file 
plot3(PDBdata.X, PDBdata.Y, PDBdata.Z, '.'); % make a 3D plot of data 

%% -- example usage: transl                                        te the atoms of 3IJU.pdb by 10 angstroms in x direction -- 

PDBdata = pdb2mat('3IJU.pdb'); % read in data from PDB file 
PDBdata.X = PDBdata.X + 10; % translate coordinates 
PDBdata.outfile = '3IJU_tran10angXdir.pdb'; % update file name 
mat2pdb(PDBdata); % output data in PDB format

%% -- example usage: make a PDB with 30 atoms in random places within a 10 angstrom box -- 
data.X = rand(1,20)*10; 
data.Y = rand(1,20)*10; 
data.Z = rand(1,20)*10; 
mat2pdb(data)

%% -- pdb2mat.m -- 
% 
% It's function is simple: give it a PDB file and out comes a 
% matlab-friendly data structure. In cumbersomely large PDB's (such as those that 
% include solvent), this can shave off a good amount of time relative to 
% many programs. Unfortunately there is no easy way to hasten the slowest 
% step, which is turning strings into doubles. 
%

%% -- mat2PDB.m -- 
% 
% this function creates a PDB from coordinate data. Represent all inputs as 
% a structure field for it to be read. The output format is as given in 
% online documentation (as of July 2012 when writing this program) 
% http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

%To speed up the program, comment out any lines for data not being used. 
%Commenting one line that converts numeric data speeds the program up by roughly 7-8%.

%example usage : this plots the atoms of 3IJU.pdb

atoms = fastPDBRead('3IJU.pdb') 
plot3(atoms.X, atoms.Y, atoms.Z, '.');