1;


maxIter   = 100;
tol       = 1e-9;

matrica_Susjedstava_G = [0,1,1,0,0;
                      0,0,1,0,1;
                      0,0,0,1,1;
                      0,1,0,0,0;
                      0,0,0,0,0];

matrica_Susjedstava_S = [0,1,0;
                        0,0,1;
                        0,0,0];

Z0 = ones(size(matrica_Susjedstava_G,1),size(matrica_Susjedstava_S,1));

Z = blondel_slicnost(matrica_Susjedstava_S,matrica_Susjedstava_G,Z0,maxIter,tol);

Z
