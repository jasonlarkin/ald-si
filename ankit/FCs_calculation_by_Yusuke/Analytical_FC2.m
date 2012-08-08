function [FC2] = Analytical_FC2(i,Atom,system,cutoff,epsilon,sigma)

  % initialize
  [Num_atom,comp] = size(Atom);
  Rij = repmat(0,Num_atom,comp+1); %(x,y,z,distsnce)
  FC2 = repmat(0,3*Num_atom,3);

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
  for j = 1:Num_atom	
    if j == i
      continue;
    end
    % calculate potential, its 1st order derivation, and its 2nd order
    if Rij(j,4)<cutoff & Rij(j,4)>0
      PHI = repmat(0,3,3);
      [PHI] = calc_PHI(Rij(j,:),cutoff,epsilon,sigma);
      FC2(3*(j-1)+1:3*j,:) = - PHI;
      FC2(3*(i-1)+1:3*i,:) = FC2(3*(i-1)+1:3*i,:) + PHI;
    end
  end

%------------------- sub functions ----------------------
  function [PHI] = calc_PHI(Rij,cutoff,epsilon,sigma)
    PHI = repmat(0,3,3);
      for alpha = 1:3
      for beta = 1:3
	[phi1] = calc_phi1(Rij(1,4),cutoff,epsilon,sigma);
	[phi2] = calc_phi2(Rij(1,4),cutoff,epsilon,sigma);
	PHI(alpha,beta) = Rij(1,alpha)*Rij(1,beta)/(Rij(1,4)^2) * (phi2 - phi1/Rij(1,4));
	if alpha == beta
	  PHI(alpha,beta) = PHI(alpha,beta) + phi1/Rij(1,4);
	end
      end
      end
  end
%------------------- sub functions ----------------------
%first order energy derivative
  function [phi1] = calc_phi1(r,cutoff,epsilon,sigma)
    if r < cutoff
      phi1 = 4*epsilon * (-12*sigma^12*r^(-13) + 6*sigma^6*r^(-7));
      else
	phi1 = 0;
    end
  end
%second order enrgy derivative
  function [phi2] = calc_phi2(r,cutoff,epsilon,sigma)
    if r < cutoff
      phi2 = 4*epsilon * (12*13*sigma^12*r^(-14) - 6*7*sigma^6*r^(-8));
      else
	phi2 = 0;
    end
  end
%-------------------------------------------------------


end
