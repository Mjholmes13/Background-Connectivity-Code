Y_G = readmatrix("C:\Users\Owner\Documents\Graduate Documents\Research\HCP\Data\Y_G.csv")
DesignMat = readmatrix("C:\Users\Owner\Documents\Graduate Documents\Research\HCP\Data\DesignMat.csv")
DesignMatT = DesignMat'

wavespecs0.ndim=1;        %%% applied to rows
wavespecs0.boundary='per'; 
wavespecs0.wavelet='db1';
wavespecs0.compress=1;
[D1,wavespecs]=DWT_rows(Y_G,wavespecs0)
D = D1'
[X, wavespecsX] = DWT_rows(DesignMatT, wavespecs)
X = X'

%wavespecs0.wavelet='db4';
wavespecs.T = 284
n=wavespecs.T;
boundary=wavespecs.boundary;
nlevels=8;
extend=1;
[W,Kj]=Get_DWT("db1",wavespecs.T,boundary,1,8);
Wtrans=W'; %W is T^w x T

writematrix(W,"C:\Users\Owner\Documents\Graduate Documents\Research\HCP\Data\W.csv")

writestruct(wavespecs, "wavespecs_YG.xml")