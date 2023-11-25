function energia=fenerg2(sinal,espacamento)

% Função que determina a energia de 1 amostra do sinal espaçadas de espaçamento. energia=fenergia2(sinal,espaçamento).
%	Autor: João Paulo Teixeira

clear energia;
i=1:espacamento:length(sinal);
energia=sinal(i).^2;
end

