function [sinal,restos]=fgerimp(Fs,f0,janela,resto)

%FUnç~ao geradora de impulsos glotais `a frequ^encia fundamental f0 num
%segmento de comprimento janela a uma frequ^encia de amostragem Fs. O
%sinal resto e restos permitem que a funç~ao comece a gerar o primeiro
%impulo do segmento no ponto onde terminou o ´ultimo impulso do
%segmento anterior.

salto=round(Fs/f0);
a=0.90;
B=[0 -a*exp(1)*log(a)];
A=[1 -2*a a^2];
comp_res=length(resto);
sinal=zeros(1,janela-comp_res);
for i=1:salto:(janela-comp_res);sinal(i)=1;end
%plot(sinal)
%sinal
%pause
sinal=[sinal, zeros(1,i+salto-(janela-comp_res)-1)];
%i+salto-(janela-comp_res)
%plot(sinal)
%pause
sinal=filter(B,A,sinal);
restos=sinal(janela-comp_res:length(sinal));
sinal=[resto,sinal(1:janela-comp_res)];
%plot(sinal)
%pause
