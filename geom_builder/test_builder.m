d = 1.2;

rs = [(-1:1/10:0.2).' 0*(-1:1/10:0.2).'+1; 0.2 0.8; -0.1 0.7; -0.5 0.6; -0.5 0; -2/3 -1/2; -0.3 -0.7; 0 -1; 0.3 -0.6; 2/3 -1/2; 0.6 0; 2/3 1/2; 1/2 0.8; (0.5:0.1:1).' (0.5:0.1:1).'*0+1;].';
rs = [(-1:0.3:-0.1).' 0*(-1:0.3:-0.1).'+1; -0.1 0.8; -0.5 0.6; -2/3 -1/2; 0 -1; 2/3 -1/2; 2/3 1/2; 1/2 0.8; (0.5:0.1:1).' (0.5:0.1:1).'*0+1;].';
% rs(1,:) = d/2*rs(1,:);
coefs = get_splines(rs);

ts = linspace(-1,1,length(rs));
tsub = linspace(-1,1,300);

[rsub,dsub,d2sub] = geom_eval(tsub,coefs);

figure(1); clf
plot(rsub(1,:),rsub(2,:)); hold on
plot(rs(1,:),rs(2,:),'x')

figure(2); clf
t = tiledlayout(1,2);
nexttile
plot(tsub,rsub(1,:),'x-')
hold on
plot(ts,rs(1,:),'x')
nexttile
plot(tsub,rsub(2,:),'x-')
hold on
plot(ts,rs(2,:),'x')

figure(3); clf
t = tiledlayout(1,2);
nexttile
plot(tsub,dsub(1,:),'x-')
nexttile
plot(tsub,dsub(2,:),'x-')

figure(4); clf
t = tiledlayout(1,2);
nexttile
plot(tsub,d2sub(1,:),'x-')
nexttile
plot(tsub,d2sub(2,:),'x-')

cparams = []; cparams.ta = -1; cparams.tb = 1;
cparams.ifclosed = 0;cparams.eps = 1e-6;
chnkr = chunkerfunc(@(t) geom_eval(t,coefs),cparams);
chnkr = reverse(chnkr);
chnkr.npt

figure(5); clf;
plot(chnkr)
hold on
quiver(chnkr)




