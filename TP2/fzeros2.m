function M = fzeros2(sinal1,N,espacamento)

% Função M=fzeros(sinal1,N,espacamento) que determina a taxa de passagem por zero do sinal com janela (rectangular) de comprimento N e espaçamento dos elementos de saida, espaçamento, a escolher.
%	Autor:	João Paulo Teixeira

n=length(sinal1);
clear M;
sinal2=zeros(N,1);
M=zeros(round(n/espacamento),1);
fim=(n-N)/espacamento;
j=1:N-1;
for i=1:fim,
  sinal2=sinal1(i*espacamento+1:i*espacamento+N);
  taxa=(abs(sign(sinal2(j+1))-sign(sinal2(j))));
  M(i)=sum(taxa);
end
%clg;
%subplot(211);plot(sinal1);
%title('Taxa de passagem por zero com N=32 com espaçamento da saida de 10 com o sinal cons_p');
%subplot(212);plot(M);
end