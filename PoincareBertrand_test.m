



echnks = [chnkr_l,chnkr_r];
test_FP
echnks(3) = chnkr;
make_geom
echnks(4) = chnkr;
close all
%%
nnode = 62;
ts = linspace(-pi/d_l,pi/d_l,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d_l);
nkappa = length(kappa);

l=2; N = 40; a = 15; M = 1e4;
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,kappa,l);


%%

src.r = [-0.6;4];
errs = zeros(nkappa,length(echnks));
for i = 1:length(echnks)

chnkr = echnks(i);

ising = 0;
nsub = 1;
double = @(s,t) chnk.lap2dquas.kern(s,t,'d',kappa,d,s0_l,sn_l,l,ising,nsub);
hilbert = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa,d,s0_l,sn_l,l,ising,nsub);
opts = [];
opts.sing = 'smooth';
opts.quad = 'native';

opts2 = [];
opts2.sing = 'smooth';
opts.quad = 'native';


% building system matrix

start = tic;
D = chunkermat(chnkr, double, opts2);
H = chunkermat(chnkr, hilbert, opts2);     

D = reshape(D,nkappa,chnkr.npt,chnkr.npt);
H = reshape(H,nkappa,chnkr.npt,chnkr.npt);


fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate',nu);
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

D_0 = chunkermat(chnkr, double, opts);
H_0 = chunkermat(chnkr, hilbert, opts2); 

D_0 = reshape(D_0,1,chnkr.npt,chnkr.npt);
H_0 = reshape(H_0,1,chnkr.npt,chnkr.npt);

alpha = reshape(exp(1i*kappa(:)*d),[nkappa,1,1]); 
for ii = -nsub:nsub
    if ii ~= 0 
        Hsub = chunkerkernevalmat(chnkr + ii*[d;0],hilbert,chnkr);
        Dsub = chunkerkernevalmat(chnkr + ii*[d;0],double,chnkr);
  
        Hsub = reshape(Hsub,[1,chnkr.npt,chnkr.npt]);
        Dsub = reshape(Dsub,[1,chnkr.npt,chnkr.npt]);

        H_0 = H_0 + alpha.^(ii).*Hsub;
        D_0 = D_0 + alpha.^(ii).*Dsub;

    end
end

dmat = D + D_0; hmat = H + H_0;


rhs = chnk.lap2dquas.kern(src,chnkr,'s',kappa,d,s0_l,sn_l,l);
for k = 1:nkappa
    H = squeeze(hmat(k,:,:));
    D = squeeze(dmat(k,:,:));
    R = rhs(k:nkappa:end);
a = 0.25*H * (H * R);
b = -0.25*R + D * (D * R);

errs(k,i) = norm(a-b) / norm(a);
end
end
max(errs,[],1)