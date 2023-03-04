function plotITERA(x1,x2,w1,w2,t,max_radius,j1_c,j2_c,N,Amplitude)

if(t==1)
    figure(2)
    %     subplot(3,3,1)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(20)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');

    title(['t=' num2str(t)])
    hold off

end

if(t==10)
    figure(3)
    %     subplot(3,3,2)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(30)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==20)
    figure(4)
    %     subplot(3,3,3)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(40)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==50)
    figure(5)
    %     subplot(3,3,4)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(50)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==100)
    figure(6)
    %     subplot(3,3,5)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(60)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==150)
    figure(7)
    %     subplot(3,3,6)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(70)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==200)
    figure(8)
    %     subplot(3,3,7)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(80)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==250)
    figure(9)
    %     subplot(3,3,8)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)])
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(90)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

if(t==300)
    figure(10)
    %     subplot(3,3,9)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',1)
    plot(w1',w2','r','linewidth',1)
    plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    title(['t=' num2str(t)]);
    hold off

    % neurônios
    TAM = length(w1);
    for i=1:TAM
        for j=1:TAM
            figure(100)
            if(i==j1_c && j==j2_c)
                stem3(w1(i,j),w2(i,j),0,'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
            else
                stem3(w1(i,j),w2(i,j),0,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            end
            hold on
        end
    end

    % entradas
    TAM = length(x1);
    for i=1:TAM
        stem3(x1(i),x2(i),0,'bo');
    end

    % sinal dos vizinhos
    for neighbour_radius=1:1:max_radius
        jj1=j1_c - neighbour_radius;
        jj2=j2_c;
        if (jj1>=1)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c + neighbour_radius;
        jj2=j2_c;
        if (jj1<=N)
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c - neighbour_radius;
        if (jj2>=1) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
        jj1=j1_c;
        jj2=j2_c + neighbour_radius;
        if (jj2<=N) % to stay in the matrix
            % plot do sinal Gaussiano do vizinho
            stem3(w1(jj1,jj2),w2(jj1,jj2),Amplitude(jj1,jj2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        end
    end

    % sinal do vencedor
    stem3(w1(j1_c,j2_c),w2(j1_c,j2_c),Amplitude(j1_c,j2_c),'ko','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    title(['t=' num2str(t)])
    hold off
end

end