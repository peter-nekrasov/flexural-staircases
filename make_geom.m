d = 2;
zk = 1;



cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.maxchunklen = 2/zk;cparams.ifclosed = 0;
nch = 20; A = -1;
% chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
chnkr0 = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr1 = reverse(chnkr0);


chnkr2 = chunkerfunc(@starfish,struct('eps',1e-10)); 
chnkr2 = move(chnkr2,[], [0;2.5],0.3,0.5);
src = []; src.r = [0;0]; src.n = [1;0];


chnkr = merge([chnkr1, chnkr2]);
figure(1);clf
plot(chnkr,'linewidth',2)
axis equal


rend = chunkends(chnkr,[1,chnkr1.nch]);
rend = rend(:,[2,3]);

verts = [rend, rend - [0;10]];
edge2verts = [[1;2], [4;3], [3;1], [2;4],[NaN;NaN]];
fchnk = cell(1,5);
fchnk{1} = chnkr0;
fchnk{5} = chnkr2;
cgrph = chunkgraph(verts,edge2verts,fchnk);


figure(2);clf
% plot(cgrph)
plot_regions(cgrph)
hold on
plot(chnkr,'.','linewidth',2)
hold off
ylim([-1.5,inf])
% axis equal
