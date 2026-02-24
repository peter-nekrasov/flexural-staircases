



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

ising = 1;
double = kernel(@(s,t) chnk.lap2dquas.kern(s,t,'d',kappa,d,s0_l,sn_l,l,ising));
double.sing = 'smooth';
hilbert = kernel(@(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa,d,s0_l,sn_l,l,ising));
hilbert.sing = 'pv';

dmat = chunkermat(chnkr,double);
hmat = chunkermat(chnkr,hilbert);

rhs = chnk.lap2dquas.kern(src,chnkr,'s',kappa,d,s0_l,sn_l,l);
for k = 1:nkappa
    H = hmat(k:nkappa:end,:);
    D = dmat(k:nkappa:end,:);
    R = rhs(k:nkappa:end);
a = 0.25*H * (H * R);
b = -0.25*R + D * (D * R);

errs(k,i) = norm(a-b) / norm(a);
end
end
max(errs,[],1)