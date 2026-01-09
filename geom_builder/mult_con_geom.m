d = 2; 
A = -0.6;

cparams = []; cparams.ta = -d/2; cparams.tb = d/2; 
cparams.maxchunklen = 1/(pi/d);

chnkr1 = chunkerfunc(@(t) cos_func(t,d,A),cparams);
cparams.ta = 0; cparams.tb = 2*pi;  cparams.maxchunklen = 2;
chnkr2 = chunkerfunc(@(t)starfish(t,3),cparams);
chnkr2 = 0.4*chnkr2 + [0;4*abs(A)];


% rend = chunkends(chnkr1,[1,chnkr1.nch]);
% rend = rend(:,[1,end]);
rend = [[-d/2;-A],[d/2;-A]];
verts = [rend, rend - [0;10]];
edge2verts = [[1;2], [4;3], [3;1], [2;4], [NaN;NaN]];
fchnk = cell(1,5);
fchnk{1} = chnkr1;
fchnk{5} = chnkr2;
cgrph = chunkgraph(verts,edge2verts,fchnk);


chnkr1 = reverse(chnkr1);
chnkr = merge([chnkr1,chnkr2]);



figure(6);clf
plot(cgrph,'.')
hold on
plot(chnkr,'.')
hold off
axis equal

figure(7);clf
plot_regions(cgrph)

