function plotMAP(x1,x2,w1,w2,t)

figure(11)
hold on
title(['t=' num2str(t-1)]);
plot(x1,x2,'ob')
plot(w1,w2,'r','linewidth',1)
plot(w1',w2','r','linewidth',1)
plot(w1,w2,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
hold off

end