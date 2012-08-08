function [ATOM] = all_system(atom_ir)

[Natom_ir,comp] = size(atom_ir);

% expand the irreducible list into whole system
ATOM = repmat(0,79,3);
ndata = 1;
ATOM(1,:) = [atom_ir(1,:)];
ndata = ndata +1;

%atom 2 (12 types) 3 permutations and 4 rotations
j = 2;

x = 1;
y = 2;
z = 3;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 2;
y = 1;
z = 3;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 3;
y = 2;
z = 1;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);

%atom 3 (6 types) 3 permutations and 2 rotations
j = 3;

x = 1;
y = 2;
z = 3;
[ATOM,ndata] = calc2(ATOM,atom_ir,j,ndata,x,y,z);
x = 2;
y = 1;
z = 3;
[ATOM,ndata] = calc2(ATOM,atom_ir,j,ndata,x,y,z);
x = 3;
y = 2;
z = 1;
[ATOM,ndata] = calc2(ATOM,atom_ir,j,ndata,x,y,z);

%atom 4 (24 types) 3 permutations and 8 rotations
j = 4;

x = 1;
y = 2;
z = 3;
[ATOM,ndata] = calc3(ATOM,atom_ir,j,ndata,x,y,z);
x = 2;
y = 1;
z = 3;
[ATOM,ndata] = calc3(ATOM,atom_ir,j,ndata,x,y,z);
x = 3;
y = 2;
z = 1;
[ATOM,ndata] = calc3(ATOM,atom_ir,j,ndata,x,y,z);

%atom 5 (12 types) 3 permutations and 4 rotations
j = 5;

x = 1;
y = 2;
z = 3;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 2;
y = 1;
z = 3;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 3;
y = 2;
z = 1;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);

%atom 6 (6 types) 4 permutations and 4 rotations
j = 6;

x = 1;
y = 2;
z = 3;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 1;
y = 3;
z = 2;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 2;
y = 1;
z = 3;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 2;
y = 3;
z = 1;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 3;
y = 1;
z = 2;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z);
x = 3;
y = 2;
z = 1;
[ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z); 



%-----------------functions-----------------

function [ATOM,ndata] = calc1(ATOM,atom_ir,j,ndata,x,y,z)
  for rot_1 = 1:2
  for rot_2 = 1:2
    test_pos = [atom_ir(j,x),atom_ir(j,y),atom_ir(j,z)];

    if x == 1
      alpha = 2;
      beta = 3;
      elseif y == 1
	alpha = 1;
	beta = 3;
      else
	alpha = 1;
	beta = 2;
    end

    if rot_1 == 2
      test_pos(1,alpha) = -test_pos(1,alpha);
    end
    if rot_2 == 2
      test_pos(1,beta) = -test_pos(1,beta);
    end
    ATOM(ndata,:) = [test_pos];
    ndata = ndata + 1;
  end
  end

end

%--------------------------------------------

function [ATOM,ndata] = calc2(ATOM,atom_ir,j,ndata,x,y,z)
  for rot_1 = 1:2
    test_pos = [atom_ir(j,x),atom_ir(j,y),atom_ir(j,z)];
    if rot_1 == 2
      test_pos(1,x) = -test_pos(1,x);
    end
    ATOM(ndata,:) = [test_pos];
    ndata = ndata + 1;
  end
end

%---------------------------------------------

function [ATOM,ndata] = calc3(ATOM,atom_ir,j,ndata,x,y,z)
  for rot_x = 1:2
  for rot_y = 1:2
  for rot_z = 1:2
    test_pos = [atom_ir(j,x),atom_ir(j,y),atom_ir(j,z)];
    if rot_x == 2
      test_pos(1,1) = -test_pos(1,1);
    end
    if rot_y == 2
      test_pos(1,2) = -test_pos(1,2);
    end
    if rot_z == 2
      test_pos(1,3) = -test_pos(1,3);
    end
    ATOM(ndata,:) = [test_pos];
    ndata = ndata +1;
  end
  end
  end
end

%----------------------------------------------

end




