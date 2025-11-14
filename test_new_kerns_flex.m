zk = 1;
zk = 0.3;

d = 1;

kappa = linspace(-pi/d,pi/d,10);
kappa = kappa - 0.3*1i*sin(kappa*d);
% kappa = kappa(10);
kappa =0*kappa;
kappa = 0;

ht = 1.02*d; hb = -1.02*d;
l = 2;
l = 1;
npxy = 100;

% if abs(zk) > 1e-10
%     skern = kernel('h','s',zk);
%     s2trkern = kernel([kernel('h','s',zk);kernel('h','sp',zk)]);
% else
%     skern = kernel('l','s');
%     s2trkern = kernel([kernel('l','s');kernel('l','sp')]);
% end


tic;
[pxys, cs] = build_flex_pxys(zk,kappa,d,ht,hb,l,npxy);
toc;


%%
src = [0;0];

nplot = 100;
XX = linspace(-3*d/2,3*d/2,nplot);
% YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ = [X(:).';Y(:).'];

tic;
[val, grad, hess,third,fourth] = qflex_green(src,targ,zk,kappa,d,pxys,cs,l,1);
toc;
val = reshape(val, length(kappa), []);
grad = reshape(grad, length(kappa), [],2);
hess = reshape(hess, length(kappa), [],3);
third = reshape(third, length(kappa), [],4);
fourth = reshape(fourth, length(kappa), [],5);


figure(1);clf
% h = pcolor(X,Y,reshape(real(val(10,:)),size(X))); h.EdgeColor = 'None';
h = pcolor(X,Y,reshape(real(grad(1,:,1)),size(X))); h.EdgeColor = 'None';
colorbar

%%
src = [0;0];

nplot = 100;
XX2 = linspace(0.05,3*d/2,nplot);
% YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ = [X(:).';Y(:).'];
tic;
[val_true, grad_true, hess_true, third_true, fourth_true] =  quasi_flex_dual_sum(X(:).',Y(:).',zk,kappa,d);
[val, grad, hess,third,fourth] = qflex_green([0;0],targ,zk,kappa,d,pxys,cs,l,1);
% [val, grad, hess] = new_green([0;0],targ,zk,kappa,d,pxys,cs,l,1);
toc;
val = reshape(val, length(kappa), []);
grad = reshape(grad, length(kappa), [],2);
hess = reshape(hess, length(kappa), [],3);
third = reshape(third, length(kappa), [],4);
fourth = reshape(fourth, length(kappa), [],5);


norm(val - val_true,'fro')/norm(val_true,'fro')
norm(grad - grad_true,'fro')/norm(grad_true,'fro')
norm(hess - hess_true,'fro')/norm(hess_true,'fro')
norm(third - third_true,'fro')/norm(third_true,'fro')
norm(fourth - fourth_true,'fro')/norm(fourth_true,'fro')




% %%
figure(2);clf
i = 1;
subplot(1,3,1)
h = pcolor(X,Y,reshape(real(grad(i,:,2)),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,2)
h = pcolor(X,Y,reshape(real(grad_true(i,:,2)),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,3)
h = pcolor(X,Y,reshape(abs(grad(i,:,2)-grad_true(i,:,2)),size(X))); h.EdgeColor = 'None';
colorbar

%%
figure(5);clf
i = 1;
subplot(1,3,1)
h = pcolor(X,Y,reshape(real(fourth(i,:,2)),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,2)
h = pcolor(X,Y,reshape(real(fourth_true(i,:,2)),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,3)
h = pcolor(X,Y,reshape(log10(abs(fourth(i,:,2)-fourth_true(i,:,2))),size(X))); h.EdgeColor = 'None';
colorbar



figure(3);clf
i = 1;
subplot(1,3,1)
h = pcolor(X,Y,reshape(real(val(i,:)),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,2)
h = pcolor(X,Y,reshape(real(val_true(i,:)),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,3)
h = pcolor(X,Y,reshape(abs(val(i,:)-val_true(i,:)),size(X))); h.EdgeColor = 'None';
colorbar


%%

% tic;
% [pxys, cs] = build_pxys(zk,kappa,d,ht,hb,skern,s2trkern,l,npxy);
% toc;
% tic;
% val = new_green([0;0],targ,zk,kappa,d,pxys,cs,l,1);
% toc;
% 
% tic;
% skern_hq = kernel('hq','s',zk,kappa,d);
% toc;
% tic;
% val_old = skern_hq.eval(struct('r',src),struct('r',targ));
% toc;
% 
