function media=fmedia(sinal1,janela,espacamento)

% Função que determina a média deslizante do sinal com numero de elementos de saida a escolher.
%	Autor:	João Paulo Teixeira

n=length(sinal1);
media=zeros(fix(n/espacamento),1);
jan2=round(janela/2);
for i=jan2:espacamento:n-janela/2,
  media(round(i/espacamento)+1)=sum(sinal1(i-jan2+1:i+jan2))/janela;
end;
for j=1:round(janela/2/espacamento),
  media(j)=sum(sinal1(1:j+jan2))/(j+jan2);
end;
v=round(i/espacamento+2);
f=(fix(n/espacamento));
for j=v:f,
  media(j)=sum(sinal1((j-1)*espacamento-jan2+1:n))/(n+1-((j-1)*espacamento-jan2+1));
end
end

