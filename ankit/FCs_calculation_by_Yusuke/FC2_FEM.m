function [FC2] = FC2_FEM(i,ATOM_DATA,system,cutoff,epsilon,sigma,drj)

  [Natom, comp] = size(ATOM_DATA);

  FC2 = repmat(0,3*Natom,3);
  for j = 1:Natom
    for beta = 1:3
      Atom_j_moved_f = ATOM_DATA;
      Atom_j_moved_f(j,beta) = Atom_j_moved_f(j,beta) + drj;
      [Fi_f] = force_on_i(i,Atom_j_moved_f,system,cutoff,epsilon,sigma);

      Atom_j_moved_b = ATOM_DATA;
      Atom_j_moved_b(j,beta) = Atom_j_moved_b(j,beta) - drj;
      [Fi_b] = force_on_i(i,Atom_j_moved_b,system,cutoff,epsilon,sigma);

      FC2(3*(j-1)+1:3*j,beta) = - (Fi_f - Fi_b)/(2*drj);
    end
  end

end