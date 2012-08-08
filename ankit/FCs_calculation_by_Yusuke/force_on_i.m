function [Fi] = force_on_i(i,Atom,system,cutoff,epsilon,sigma)

  % initialize
  [Num_atom,comp] = size(Atom);
  Rij = repmat(0,Num_atom,comp+1); %(x,y,z,distsnce)

  % distance
  for j = 1:Num_atom
    Rij(j,1:3) = Atom(j,:) - Atom(i,:);
    % periodic boundary
    % x
    if abs(Rij(j,1) + system(1,1)) < abs(Rij(j,1))
      Rij(j,1) = Rij(j,1) + system(1,1);
      elseif abs(Rij(j,1) - system(1,1)) < abs(Rij(j,1))
	Rij(j,1) = Rij(j,1) - system(1,1);
    end
    % y
    if abs(Rij(j,2) + system(1,2)) < abs(Rij(j,2))
      Rij(j,2) = Rij(j,2) + system(1,2);
      elseif abs(Rij(j,2) - system(1,2)) < abs(Rij(j,2))
	Rij(j,2) = Rij(j,2) - system(1,2);
    end
    % z
    if abs(Rij(j,3) + system(1,3)) < abs(Rij(j,3))
      Rij(j,3) = Rij(j,3) + system(1,3);
      elseif abs(Rij(j,3) - system(1,3)) < abs(Rij(j,3))
	Rij(j,3) = Rij(j,3) - system(1,3);
    end
    Rij(j,4) = sqrt(sum(Rij(j,1:3).*Rij(j,1:3)));
  end

  % all force on i
  Fi = [0;0;0];
  for j = 1:Num_atom	
    if j == i
      continue;
    end
    % calculate potential, its 1st order derivation, and its 2nd order
    if Rij(j,4)<cutoff & Rij(j,4)>0
      [phi1] = calc_phi1(Rij(j,:),cutoff,epsilon,sigma);
      else
      phi1 = [0;0;0];
    end
  Fi = Fi + phi1; 
  end

%-------------------- sub function ----------------------
function [phi1] = calc_phi1(Rij,cutoff,epsilon,sigma)
  phi1 = [0;0;0];
  if Rij(1,4) < cutoff
    for alpha = 1:3
      phi1(alpha,1) = 4*epsilon * ...
	(-12*sigma^12*Rij(1,4)^(-13) + 6*sigma^6*Rij(1,4)^(-7)) * Rij(1,alpha)/Rij(1,4);
    end
  end
end
%-------------------- sub function ----------------------

end