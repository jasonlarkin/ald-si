function [FC3] = Analytical_FC3(i,Atom,system,cutoff,epsilon,sigma)
drj = 1.0E-6*sigma;

  % initialize
  [Num_atom,comp] = size(Atom);
  Rij = repmat(0,Num_atom,comp+1); %(x,y,z,distsnce)
  FC3 = repmat(0,3*Num_atom,3*3);

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

  for j = i:Num_atom
    if j == i
      continue;
    end
    PHI = repmat(0,3,3*3);
    for alpha = 1:3
    for beta  = 1:3
    for gamma = 1:3
      [phi1] = calc_phi1(Rij(j,4),cutoff,epsilon,sigma);
      [phi2] = calc_phi2(Rij(j,4),cutoff,epsilon,sigma);
      [phi3] = calc_phi3(Rij(j,4),cutoff,epsilon,sigma);

      PHI(alpha,3*(beta-1)+gamma) = Rij(j,alpha)*Rij(j,beta)*Rij(j,gamma)/(Rij(j,4)^3) * ( phi3 - 3*phi2/Rij(j,4) + 3*phi1/(Rij(j,4)^2) );
      if alpha == beta
	r1 = Rij(j,gamma);
	else
	r1 = 0;
      end
      if beta == gamma
	r2 = Rij(j,alpha);
	else
	r2 = 0;
      end
      if gamma == alpha
	r3 = Rij(j,beta);
	else
	r3 = 0;
      end
      PHI(alpha,3*(beta-1)+gamma) = PHI(alpha,3*(beta-1)+gamma) + (r1 + r2 + r3)/Rij(j,4) * (phi2/Rij(j,4) - phi1/(Rij(j,4)^2));
    end
    end
    end
    FC3(3*(j-1)+1:3*j,:) = - PHI;
    % FC3_ijj = - FC3_iij
    % FC3_iii = - sum(FC3_iij) for all j
    FC3(3*(i-1)+1:3*i,:) = FC3(3*(i-1)+1:3*i,:) - PHI;
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
%third order energy derivative
  function [phi3] = calc_phi3(r,cutoff,epsilon,sigma)
    if r < cutoff
      phi3 = 4*epsilon * (-12*13*14*sigma^12*r^(-15) + 6*7*8*sigma^6*r^(-9));
      else
	phi3 = 0;
    end
  end
%-------------------------------------------------------
end