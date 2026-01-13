function [sys,sn,l] = clamped_mat(chnkr,zk,nu,kappa,d)
nkappa = length(kappa);

l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);

%%
ising = 1;
fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate',kappa,d,sn,[],[],l,ising);

curv = signed_curvature(chnkr);
curv = curv(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = reshape(sys,nkappa,2*chnkr.npt,2*chnkr.npt);

sys = sys - reshape(0.5*eye(2*chnkr.npt),1,2*chnkr.npt,2*chnkr.npt);
sys(:,2:2:end,1:2:end) = sys(:,2:2:end,1:2:end) + reshape(curv.*eye(chnkr.npt),1,chnkr.npt,chnkr.npt);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)


end