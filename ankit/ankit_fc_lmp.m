
%--------------------------------------------------------------------------
%-------------------------INPUT--------------------------------------------
%--------------------------------------------------------------------------
LD.constant.kb = 1.3806E-23;                    %aJ/k (1.3806E-23 J/K)
LD.constant.hbar = 1.054E-34;                %J/s
LD.constant.i = sqrt(-1);
LD.constant.c = 29979245800.00019;      %cm/s
LD.constant.s2ps = 10E-13;
LD.constant.ang2m = 10E-11;

format long

%--------------------------------------------------------------------------
%-------------------------READ---------------------------------------------
%--------------------------------------------------------------------------

LD.FC2.id = load('./2x/FULLorder2.dat');
LD.FC3.id = load('./2x/FULLorder3.dat');
%LD.FC4.id = load('./2nn/FULLorder4.dat');

LD.pos = load('./2x/FULLpos.dat');

%--------------------------------------------------------------------------
%LD Input
%--------------------------------------------------------------------------
%define supercell integers
LD.Nx = 2;                LD.Ny = 2;                LD.Nz = 2;

%cell parameters
LD.cell = [5.430 5.430 5.430 90 90 90 0 0 0 0 0 0];
%lattice vectors
LD.latvec =...
        [0.0 0.5 0.5;
		0.5 0.0 0.5;
		0.5 0.5 0.0];
LD.latvec = LD.latvec*LD.cell(1);
LD.latvec_rec = [-1.0 1.0 1.0;
			1.0 -1.0 1.0;
			1.0 1.0 -1.0];

LD.x.ucell.frac =...
    [0.1250000000000000  0.1250000000000000  0.1250000000000000
    0.8750000000000000  0.8750000000000000  0.8750000000000000];
                
LD.x.ucell.gulp = [0.8750000000000000  0.8750000000000000  0.8750000000000000 0 1 1 1
                   0.1250000000000000  0.1250000000000000  0.1250000000000000 0 1 1 1];
%--------------------------------------------------------------------------
%LD Input
%--------------------------------------------------------------------------

%change to Ang
LD.pos(:,2:4) = LD.pos(:,2:4)*LD.cell(1); 

plot3( LD.pos(:,2),LD.pos(:,3),LD.pos(:,4),'.' )

%--------------------------------------------------------------------------
%Prepare LAMMPS Data
%--------------------------------------------------------------------------
lammps.alpha = (180/pi)* ...
    atan2(norm(cross(LD.latvec(1,:),LD.latvec(2,:))),...
    dot(LD.latvec(1,:),LD.latvec(2,:)));
lammps.beta = lammps.alpha;
lammps.gamma = lammps.beta;

lammps.lx = LD.cell(1)*LD.Nx;
lammps.xy = LD.cell(2)*cos(LD.cell(6));
lammps.xz = LD.cell(3)*cos(LD.cell(5));
lammps.ly = sqrt(LD.cell(2)^2 - lammps.xy^2);
lammps.yz = ( LD.cell(2)*LD.cell(3)*cos(LD.cell(4)) - ...
    lammps.xy*lammps.xz ) / lammps.ly;
lammps.lz = sqrt( LD.cell(3)^2 - lammps.xz^2 - lammps.yz^2 );

%set the lammps ids
    LD.id(1:size(LD.pos,1),1) = 1:size(LD.pos,1); 
    LD.type(1:size(LD.pos,1),1) = 1;
%initialize  
    LD.dx = 0.00000001;
    LD.precision = int2str(log10(ceil(1/LD.dx))+1);
    
    LD.pos_tmp = LD.pos(:,:);    
    LAMMPS.force0 = ankit_lmp_force(LD);
    
    
%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
tic
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FC2
%--------------------------------------------------------------------------

for iFC = 1:1:size(LD.FC2.id,1)
%initialize    
    LD.pos_tmp = LD.pos(:,:);
%increment position for iFC,2
    LD.pos_tmp(LD.FC2.id(iFC,2)+1, LD.FC2.id(iFC,4) +1 ) = ...
        LD.pos_tmp(LD.FC2.id(iFC,2)+1, LD.FC2.id(iFC,4) +1 ) + LD.dx;
%calc changed force    
    LAMMPS.force = ankit_lmp_force(LD);
%change in force on iFC,1
    LD.FC2.phi(iFC) =...
        -(LAMMPS.force(LD.FC2.id(iFC,1)+1, LD.FC2.id(iFC,3) ) - ...
        LAMMPS.force0(LD.FC2.id(iFC,1)+1, LD.FC2.id(iFC,3) ) ) / LD.dx ;
%--------------------------------------------------------------------------
%pause  
%--------------------------------------------------------------------------  
end

plot(LD.FC2.phi)
    
%--------------------------------------------------------------------------
%pause  
%--------------------------------------------------------------------------     
    
%--------------------------------------------------------------------------
%FC3
%--------------------------------------------------------------------------

for iFC = 1:1:size(LD.FC3.id,1)
    
%initialize    
    LD.pos_tmp = LD.pos(:,:);
%increment position beta
    LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) - LD.dx;
%increment position gamma
    LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) - LD.dx;
%calc changed force
    LAMMPS.force1 = ankit_lmp_force(LD);

%initialize    
    LD.pos_tmp = LD.pos(:,:);
%increment position beta
    LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) - LD.dx;
%increment position gamma
    LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) + LD.dx;
%calc changed force
    LAMMPS.force2 = ankit_lmp_force(LD);
    
%initialize    
    LD.pos_tmp = LD.pos(:,:);
%increment position beta
    LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) + LD.dx;
%increment position gamma
    LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) - LD.dx;
%calc changed force
    LAMMPS.force3 = ankit_lmp_force(LD);
    
%initialize    
    LD.pos_tmp = LD.pos(:,:);
%increment position beta
    LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,2) + 1 , LD.FC3.id(iFC,5) + 1 ) + LD.dx;
%increment position gamma
    LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) = ...
        LD.pos_tmp(LD.FC3.id(iFC,3) + 1 , LD.FC3.id(iFC,6) + 1 ) + LD.dx;
%calc changed force
    LAMMPS.force4 = ankit_lmp_force(LD);

    
%change in force on alpha
    LD.FC3.phi(iFC) =...
        -...
        (...
        LAMMPS.force4(LD.FC3.id(iFC,1) + 1 , LD.FC3.id(iFC,4) ) - ...
        LAMMPS.force3(LD.FC3.id(iFC,1) + 1 , LD.FC3.id(iFC,4) ) - ...
        LAMMPS.force2(LD.FC3.id(iFC,1) + 1 , LD.FC3.id(iFC,4) ) + ...
        LAMMPS.force1(LD.FC3.id(iFC,1) + 1 , LD.FC3.id(iFC,4) ) ...
        )...
        / (4*LD.dx^2);
%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 

end

plot(LD.FC3.phi)

%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
%FC4
%--------------------------------------------------------------------------

%set disp in Ang
%     LD.dx = 0.00001; 
% for iFC = 1:1:size(LD.FC4.id,1)
% %--------------------------------------------------------------------------
% % 1 - dr
% %--------------------------------------------------------------------------    
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force1 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force2 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------    
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force3 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------
% % 2 - dr
% %--------------------------------------------------------------------------    
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) + LD.dx;
%     LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force4 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------    
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) + LD.dx;
%     LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force5 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------    
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) + LD.dx;
%     LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force6 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------
% % 3 - dr
% %--------------------------------------------------------------------------    
% %initialize    
%     LD.pos_tmp = LD.pos(:,:);
% %increment position
%     LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,2)+1, LD.FC4.id(iFC,6) +1 ) + LD.dx;
%     LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,3)+1, LD.FC4.id(iFC,7) +1 ) + LD.dx;
%     LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) = ...
%         LD.pos_tmp(LD.FC4.id(iFC,4)+1, LD.FC4.id(iFC,8) +1 ) + LD.dx;
% %calc changed force
%     LAMMPS.force7 = ankit_lmp_force(LD);
% %--------------------------------------------------------------------------
%  
%     LD.FC4.phi(iFC) =...
%         -...
%         (...
%         LAMMPS.force7(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) - ...
%         LAMMPS.force6(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) - ...
%         LAMMPS.force5(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) - ...
%         LAMMPS.force4(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) + ...
%         LAMMPS.force3(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) + ...
%         LAMMPS.force2(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) + ...
%         LAMMPS.force1(LD.FC4.id(iFC,1)+1, LD.FC4.id(iFC,5) ) ...
%         )...
%         / (LD.dx^3);
% end
% 
% plot(LD.FC4.phi)
% 
% %--------------------------------------------------------------------------
% toc
% %--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%pause  
%-------------------------------------------------------------------------- 


dlmwrite('./PHI2.dat',LD.FC2.phi','delimiter',' ','precision', '%15.10f');
dlmwrite('./PHI3.dat',LD.FC3.phi','delimiter',' ','precision', '%15.10f');
dlmwrite('./PHI4.dat',LD.FC4.phi','delimiter',' ','precision', '%15.10f');


