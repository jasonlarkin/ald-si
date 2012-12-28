% Anharmonic Lattice Dynamics Code
% LJ potential
% FCC crystal
% 7/20/2012
% Yusuke Masao

% real dimention

format long

 sigma_Ar = 3.4E-10;    %[m]
 epsilon_Ar = 1.67E-21; %[J]
 mass_Ar = 6.6326E-26;  %[kg]
 kb = 1.3806E-23;
 hbar = 1.054E-34;
 ev = 1.60217646E-19; % [J]

% non-dimention
pi = atan(1) * 4;
sigma = sigma_Ar/sigma_Ar;
epsilon = epsilon_Ar/epsilon_Ar;
cutoff = 2.5*sigma;
lattice_const = [1.5636,1.5636,1.5636]*sigma;

drj = 1.0E-6*sigma;
drk = 1.0E-6*sigma;


%--------- simulation system---------------
cell = [4,4,4]; % (x,y,z)
system = [cell(1,1) * lattice_const(1,1),cell(1,2) * lattice_const(1,2),cell(1,3) * lattice_const(1,3)];
trans_vec  =...
	    [ 0.0 0.0 0.0;
	      0.0 0.5 0.5;
	      0.5 0.0 0.5;
	      0.5 0.5 0.0];
reci_vec=...
	    [-1.0  1.0  1.0;
	      1.0 -1.0  1.0;
	      1.0  1.0 -1.0];
%-----------------------------------------------------------
%-------------------  FCs CALCULATION ----------------------
%-----------------------------------------------------------

tic

i=1;
%(1) Build calculation system
Atom = load('./data/ir_atom.dat');

[Natom comp] = size(Atom);
for ndata = 1:Natom
  Atom(ndata,:) = Atom(ndata,:) .* lattice_const;
end 

[ATOM_DATA] = all_system(Atom);
[Natom comp] = size(ATOM_DATA);


%[Fi_eq] = force_on_i(i,ATOM_DATA,system,cutoff,epsilon,sigma)


% (2) CALCULATION ABOUT 2nd-order and 3rd order FCs 
[FC2_a] = Analytical_FC2(i,ATOM_DATA,system,cutoff,epsilon,sigma);
[FC2_n] = FC2_FEM(i,ATOM_DATA,system,cutoff,epsilon,sigma,drj);

[FC3_a] = Analytical_FC3(i,ATOM_DATA,system,cutoff,epsilon,sigma);
[FC3_n] = FC3_FEM(i,ATOM_DATA,system,cutoff,epsilon,sigma,drj,drk);

toc
%-----------------------------------------------------------
%-------------------     SAVE DATA    ----------------------
%-----------------------------------------------------------
save('./data/atom.dat','-ascii','-double','ATOM_DATA');
save('./data/FC2_a.dat','-ascii','-double','FC2_a');
save('./data/FC2_n.dat','-ascii','-double','FC2_n');
save('./data/FC3_a.dat','-ascii','-double','FC3_a');
save('./data/FC3_n.dat','-ascii','-double','FC3_n');
