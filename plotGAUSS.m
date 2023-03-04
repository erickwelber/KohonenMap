function plotGAUSS(x1,x2,w1,w2,max_radius,j1_c,j2_c,N,Amplitude)

% neurônios
TAM = length(w1);
for i=1:TAM
    for j=1:TAM
        figure(19)
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

hold off

end