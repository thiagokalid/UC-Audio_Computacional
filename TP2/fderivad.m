function der=fderivad(sinal)


% Fun��o que determina a derivada de um sinal atrav�s de delta Y. der=fderivada(sinal).
%	Autor:	Jo�o Paulo Ramos Teixeira

fim=length(sinal)-1;
der=zeros(fim,1);
i=1:fim;
der(i)=sinal(i+1)-sinal(i);
end
