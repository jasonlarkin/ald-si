function [FC3] = FC3_FEM(i,ATOM_DATA,system,cutoff,epsilon,sigma,drj,drk)

% Using central difference

[Natom,comp] = size(ATOM_DATA);

% calc Phi_ijj
% Phi_iij = -Phi_ijj
FC3_ijj = repmat(0,3*Natom,9);
for j = 1:Natom
  for beta = 1:3
    test_Atom_j_f = ATOM_DATA;
    test_Atom_j_f(j,beta) = test_Atom_j_f(j,beta) + drj;

    test_Atom_j_b = ATOM_DATA;
    test_Atom_j_b(j,beta) = test_Atom_j_b(j,beta) - drj;

    for gamma = 1:3
      test_Atom_jk_ff = test_Atom_j_f;
      test_Atom_jk_ff(j,gamma) = test_Atom_jk_ff(j,gamma) + drk;
      [Fi_ff] = force_on_i(i,test_Atom_jk_ff,system,cutoff,epsilon,sigma);

      test_Atom_jk_fb = test_Atom_j_f;
      test_Atom_jk_fb(j,gamma) = test_Atom_jk_fb(j,gamma) - drk;
      [Fi_fb] = force_on_i(i,test_Atom_jk_fb,system,cutoff,epsilon,sigma);

      test_Atom_jk_bf = test_Atom_j_b;
      test_Atom_jk_bf(j,gamma) = test_Atom_jk_bf(j,gamma) + drk;
      [Fi_bf] = force_on_i(i,test_Atom_jk_bf,system,cutoff,epsilon,sigma);

      test_Atom_jk_bb = test_Atom_j_b;
      test_Atom_jk_bb(j,gamma) = test_Atom_jk_bb(j,gamma) - drk;
      [Fi_bb] = force_on_i(i,test_Atom_jk_bb,system,cutoff,epsilon,sigma);

      dF = (Fi_ff - Fi_fb - Fi_bf + Fi_bb);

      FC3_f_ijj(3*(j-1)+1:3*j,3*(gamma-1)+beta) = - dF/(4*drj*drk);
    end
  end
  % average
  [FC3_f_ijj(3*(j-1)+1:3*j,:)] = average(FC3_f_ijj(3*(j-1)+1:3*j,:));
end

  FC3 = FC3_f_ijj;


%----------------------------------- subfunction ---------------------------------
  function  [A] = average(A); % A = 3*3*3 matrix (3*9)
    B = A;

    for alpha = 1:3
    for beta = 1:3
    for gamma = 1:3
      A(alpha,3*(gamma-1)+beta)=...
	( B(alpha,3*(gamma-1)+beta) + B(alpha,3*(beta-1)+gamma) ...  
	+ B(beta,3*(alpha-1)+gamma) + B(beta,3*(gamma-1)+alpha) ...
	+ B(gamma,3*(beta-1)+alpha) + B(gamma,3*(alpha-1)+beta) )/6;
    end
    end
    end
  end
%----------------------------------- subfunction ---------------------------------



end