function der=fderivad(sinal)


% Função que determina a derivada de um sinal através de delta Y. der=fderivada(sinal).
%	Autor:	João Paulo Ramos Teixeira

fim=length(sinal)-1;
der=zeros(fim,1);
i=1:fim;
der(i)=sinal(i+1)-sinal(i);
end
