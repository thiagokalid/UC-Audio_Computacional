function energia=fenerg2(sinal,espacamento)

% Fun��o que determina a energia de 1 amostra do sinal espa�adas de espa�amento. energia=fenergia2(sinal,espa�amento).
%	Autor: Jo�o Paulo Teixeira

clear energia;
i=1:espacamento:length(sinal);
energia=sinal(i).^2;
end

