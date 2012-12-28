7/20/2012
written by Yusuke Masao





The code "yusuke_FCs.m" is a script-m-file to run the 2nd and 3rd order FCs calculation.
As one example I carry out a calculation on LJ FCC perfect crystal.
If necessary, Please improve this code when you apply this for other systems.

Contents of this code are like this:

(1) Bulding calculation system
    Reading the position of irreducible atom, function "all_system.m" extend it into whole system.
    This function-m-file is only valid for FCC crystal whose cutoff radious is 2.5*sigma.
    irreducible atom data should be written as like this:

    pos_x1   pos_y1   pos_z1
    pos_x2   pos_y2   pos_z2
	.
	.
	.
    pos_xn   pos_yn   pos_zn



(2) FCs calculation
    In both calculations to obtain 2nd and 3rd order FCs, I have prepared two methods, an analytical approach and a numerical one (FEM). But be careful, "ANALYTICAL APPROACH IS ONLY VALID FOR LJ POTENTIAL!"
    As for a numerical approach, especially 3rd order FCs, this calculation is expensive so that it is very important to choose proper triplets of atoms (i,j,k) which give non-zero FCs.
    For example, because LJ potential is two bady potential, FCs shold be zero when (i,j,k) are totally different (See my document "7.13.pdf"). The code "FC3_FEM.m" only calculate FCs on (i,j,j).
    Generally, FCs on (i,i,j) can be calculated with summation rule, so that FC(i,i,j) = -FC(i,j,j), because 3rd FCs should be zero when k is not equal to j. Moreover, self term FC(i,i,i) should be zero in LJ FCC crystal (it's easy to calculate analytically).

    A function-m-file "force_on_i.m" calculate total force on atom i. Its result is output as a 3*1 matrix. 
