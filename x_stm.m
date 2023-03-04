function dxj = x_stm(t,xj,MDETN,yi,Bj,Dij,aj,x_past,N)

TERM1 = 0;
dxj = 0;

for i=1:N^2
    TERM1  = TERM1 + Dij(i)*x_past(i);
end   

dxj = -aj*xj + TERM1 + Bj*sum(MDETN.*yi);

end