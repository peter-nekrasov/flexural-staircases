
nkappa = length(kappas);
npoles = 4;
poles_vec = NaN*zeros(npoles,nkappa);

for i = 1:nkappa
    for j = 1:npoles
    if length(poles{i}) > j-1
        poles_vec(j,i) = poles{i}(j);
    end
    end
end

ibad = find(all(isnan(poles_vec),1));
poles_vec(:,ibad) = [];
kappas(ibad) = [];

figure(1);clf
plot(kappas, poles_vec,'o-','LineWidth',2)
hold on
plot(kappas, kappas,'k--','LineWidth',1.5)
hold off
xlabel('$\xi$','Interpreter','latex')
ylabel('$k$','Interpreter','latex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')


% exportgraphics(gcf,'free_disp.pdf','Resolution',200)