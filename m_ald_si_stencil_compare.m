LD(1).FC2.phi = load('./ankit/15/0.00001/PHI2.dat');
LD(1).FC3.phi = load('./ankit/15/0.00001/PHI3.dat');

LD(2).FC2.phi = load('./ankit/stencil/0.00001/PHI2.dat');
LD(2).FC3.phi = load('./ankit/stencil/0.00001/PHI3.dat');


plot(...
    1E6*LD(1).FC2.phi./LD(2).FC2.phi,'.'...
    )
pause
plot(...
    LD(1).FC3.phi./LD(2).FC3.phi,'.'...
    )