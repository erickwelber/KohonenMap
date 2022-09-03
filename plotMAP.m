function plotMAP(x1,x2,w1,w2,t)

figure(2)
plot(x1,x2,'ob')
hold on
plot(w1,w2,'r','linewidth',2)
plot(w1',w2','r','linewidth',2)
plot(w1,w2,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
hold off
title(['t=' num2str(t)]);
drawnow
end