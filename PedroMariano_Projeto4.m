#****************************************************
#                                                   *
#       Oscilador Harmonico Quantico                *
#                                                   *
#****************************************************

# Programa que calcula a evolucao de um estado quantico num oscilador
# harmonico unidimensional.

#***********PARTE QUE PODE SER ALTERADA: ************************************************

# Parametros iniciais:
evolu_temp = 'Sim';  # Sim para adicionar um momento ao pacote de onda
k0 = 20;             # numero de onda
m = 1;               # massa
e = 0;               # estado desejado, 0 a 3
w = 50;              # frequencia angular

r = 0.075;           # parametro, r = dt/(2*dx^2)
dt = 10^-4;          # definicao temporal
T = 2*pi/w;          # periodo
t = 0:dt:T;          # intervalo temporal

#***********FIM**************************************************************************


# Outros parametros iniciais:
dx = sqrt(dt/(2*r)); # definicao espacial
x = -1:dx:1;         # intervalo espacial
x = x';              # transforma em vetor coluna
a = m*w;             # definicao, escalar
y = sqrt(a)*x;       # definicao, vetor coluna

# Partes real e imaginaria R, I de psi:
R = zeros(length(x),length(t)); # parte real
I = zeros(length(x),length(t)); # parte imaginaria
# Vetores para a normalizacao N, potencial V, e energia E:
N = zeros(length(t),1);         # vetor para verificar a normalizacao
V = zeros(length(x),1);         # potencial
E = zeros(length(t),1);         # vetor para as energias
# Vetores para os valores esperados:
X = zeros(length(t),1);         # vetor para o valor esperado da posicao
P = zeros(length(t),1);         # vetor para o valor esperado do momento
grad_V = zeros(length(t),1);    # vetor para o valor esperado do grad(V)
# Polinomios hermitianos (vetor coluna para cada polinomio):
hermit = [ones(length(x),1) sqrt(2)*y 1/sqrt(2)*(2*y.^2-1) 1/sqrt(3)*(2*y.^3-3*y)];

# Condicoes iniciais:
V(:,1) = m/2 * w^2 * x.^2;          # potencial do OH, vetor coluna
psi0 = (a/pi)^(1/4) * exp(-y.^2/2); # funcao de onda do estado fundamental, vetor coluna
psi = hermit(:,e+1).*psi0;          # funcao de onda do estado desejado, vetor coluna

# Adiciona a evolucao temporal se desejado:
if evolu_temp=='Sim'
  psi = psi.*exp(i*k0*x);
endif

# Partes real e imaginaria de psi: 
R(:,1) = real(psi);   # vetor coluna para cada t
I(:,1) = imag(psi);   # vetor coluna para cada t


# Energias esperadas:
for n=0:5
  E_esperada(n+1) = w*(n+1/2);
endfor

#***********CALCULO: ********************************************************************

# Calculo das partes real e imaginaria de psi por RK2:
for n=1:length(t)-1; # para cada t
  
  Rmeio = zeros(length(x),1); # vetor
  Imeio = zeros(length(x),1); # vetor
  
  for k=2:length(x)-1; # para cada x
    
    # k1:
    k1_R = -r*(I(k+1,n)-2*I(k,n)+I(k-1,n)) + V(k)*I(k,n)*dt;
    k1_I = +r*(R(k+1,n)-2*R(k,n)+R(k-1,n)) - V(k)*R(k,n)*dt;
    
    # k 1/2:
    Rmeio(k) = R(k,n) + k1_R/2;
    Imeio(k) = I(k,n) + k1_I/2;
    
  endfor

  # agora tenho os vetores Imeio(:,n) e Rmeio(:,n) inteiramente calculados 
  # posso calcular os k2:
  
  for k=2:length(x)-1; #para cada x novamente
    
    #k2:
    k2_R = -r*(Imeio(k+1)-2*Imeio(k)+Imeio(k-1)) + V(k)*Imeio(k)*dt;
    k2_I = +r*(Rmeio(k+1)-2*Rmeio(k)+Rmeio(k-1)) - V(k)*Rmeio(k)*dt;
    
    #Novos valores:
    R(k,n+1) = R(k,n) + k2_R;
    I(k,n+1) = I(k,n) + k2_I;
    
  endfor
  
endfor

# A funcao psi:
psi = complex(R,I); # junto as partes real e imaginaria
psi2 = abs(psi).^2; # o modulo dessa funcao ao quadrado

#***********NORMALIZACAO: ***************************************************************

# Calcular a normalizacao usando a funcao trapz:
for n=1:length(t);
  N(n) = trapz(psi2(:,n))*dx; # N para cada t
endfor

# Grafico da normalizacao:
figure
args = {"linewidth",2};
plot(t,N,args{:});hold on;     # plota N x t
title('Norma N(t)');           # coloca um título
xlabel('t');                   # coloca uma legenda no eixo x
ylabel('N(t)');                # coloca uma legenda no eixo y
ylim([0 1.3]);                 # define os limites no eixo y
set(gca,'FontSize',20)         # tamanho dos numeros nos eixos

#***********ENERGIA: ********************************************************************

# Funcao que calcula o integrando do valor esperado da energia:
function f = integrando_energia(psi,j,n,dt,r,V)
  # retorna o integrando para calculo do valor esperado de energia:
  # psi*(j,n) H psi(j,n), onde psi* reprensenta o complexo conj. de psi
  f = conj(psi(j,n))*((-r/dt)*(psi(j+1,n)-2*psi(j,n)+psi(j-1,n)) + V(j).*psi(j,n));
endfunction

# Calcular a energia:
for n=1:length(t); 
  for j=2:length(x)-1
    E(n) = E(n) + integrando_energia(psi,j,n,dt,r,V)*dx; # E para cada t
  endfor
endfor

# Grafico da energia:
figure
plot(t,real(E),'b',args{:},t,imag(E),'r',args{:}); # plota partes Real e Imag de E x t
for n=1:length(E_esperada)                         # plota os valores esperados de energia:
  line([t(1) t(length(t))],[E_esperada(n) E_esperada(n)],"linestyle",'--')
  texto = strcat('n=',num2str(n-1));                         # gera uma string
  text(t(length(t))*1.03,E_esperada(n),texto,'FontSize',16)  # coloca um texto
endfor
titulo = strcat('Energia E(t). Estado=',num2str(e));# gera uma string
title(titulo)                          # coloca um título
xlabel('t');                           # coloca uma legenda no eixo x
ylabel('E(t)');                        # coloca uma legenda no eixo y
xlim([0 t(length(t))*1.6]);            # define os limites no eixo x
ylim([-50 300]);                       # define os limites no eixo y
legenda = legend('Re E(t)','Im E(t)'); # coloca uma legenda
set(legenda,'FontSize',20)             # tamanho da legenda
set(gca,'FontSize',20)                 # tamanho dos numeros nos eixos

#***********A FORMA DE PSI (t=0): *******************************************************

figure
p1=plotyy(x,psi2(:,1),x,V);               # plota psi^2 e V
titulo = strcat('Psi e V. t=',num2str(t(1)),'. Estado=',num2str(e));# gera uma string
title(titulo);                            # coloca um titulo
xlabel('x');                              # coloca uma legenda no eixo x
ylabel(p1(1),'|psi(x,t)|**2');            # coloca uma legenda no eixo y da esquerda
ylabel(p1(2),'V(x)');                     # coloca uma legenda no eixo y da direita
ylim([0 6]);                              # mantem a escala do grafico constante
#xlim([-1 1]);                            # mantem a escala do grafico constante
set(p1,'FontSize',20)                     # tamanho dos numeros nos eixos


#***********VALOR ESPERADO DE X: ********************************************************

if evolu_temp=='Sim'

  # Funcao que calcula o integrando do valor esperado da posicao:
  function f = integrando_posicao(psi,j,n,x)
    # retorna o integrando para calculo do valor esperado da posicao:
    # psi*(j,n) X psi(j,n), onde psi* reprensenta o complexo conj. de psi
    f = conj(psi(j,n)) .* x(j) .* psi(j,n);
  endfunction

  # Calcular o valor esperado da posicao:
  for n=1:length(t); 
    for j=2:length(x)-1
      X(n) = X(n) + integrando_posicao(psi,j,n,x)*dx; # <x> para cada t
    endfor
  endfor

  # Grafico da posicao esperada:
  figure
  plot(t,real(X),args{:},t,imag(X),args{:});# plota partes Real e Imag de <x> x t
  titulo = strcat('<x>(t). Estado=',num2str(e));    # gera uma string
  title(titulo)                             # coloca um título
  xlabel('t');                              # coloca uma legenda no eixo x
  ylabel('<x>(t)');                         # coloca uma legenda no eixo y
  legenda = legend('Re <x>(t)','Im <x>(t)');# coloca uma legenda
  set(legenda,'FontSize',22)                # tamanho da legenda
  set(gca,'FontSize',20)                    # tamanho dos numeros nos eixos

  # Derivada de X em relacao a t:
  for n=2:length(t)-1
    dXdt(n) = (X(n+1)-X(n-1))/(2*dt);
  endfor
  dXdt(n+1) = dXdt(n);


  #***********VALOR ESPERADO DE P: ********************************************************

  # Funcao que calcula o integrando do valor esperado do momento:
  function f = integrando_momento(psi,j,n,x,dx)
    # retorna o integrando para calculo do valor esperado do momento:
    # psi*(j,n) P psi(j,n), onde psi* reprensenta o complexo conj. de psi
    f = conj(psi(j,n)) .* -i .* (psi(j+1,n)-psi(j-1,n))/(2*dx);
  endfunction

  # Calcular o valor esperado do momento:
  for n=1:length(t); 
    for j=2:length(x)-1
      P(n) = P(n) + integrando_momento(psi,j,n,x,dx)*dx; # <p> para cada t
    endfor
  endfor

  # Grafico de dXdt por <p> (teo. de Ehrenfest):
  figure
  p1 = plotyy(t,dXdt,t,P);     # plota d<x>/dt e <p> x t
  titulo = strcat('d<x>/dt e <p>. Estado=',num2str(e));# gera uma string
  title(titulo)                             # coloca um título
  xlabel('t');                              # coloca uma legenda no eixo x
  ylabel(p1(1),'d<x>/dt');                  # coloca uma legenda no eixo y a esquerda
  ylabel(p1(2),'<p>');                      # coloca uma legenda no eixo y a direita
  legenda = legend('d<x>/dt','<p>','location','north');# coloca uma legenda
  set(legenda,'FontSize',22)                # tamanho da legenda
  set(p1,'FontSize',20)                     # tamanho dos numeros nos eixos

  # Derivada de P em relacao a t:
  for n=2:length(t)-1
    dPdt(n) = (P(n+1)-P(n-1))/(2*dt);
  endfor
  dPdt(n+1) = dPdt(n);


  #***********VALOR ESPERADO GRAD V: ******************************************************

  # Funcao que calcula o integrando do valor esperado de -grad(V):
  function f = integrando_gradV(psi,j,n,x,V,dx)
    # retorna o integrando para calculo do valor esperado de -grad(V):
    # psi*(j,n) <-grad(V)> psi(j,n), onde psi* reprensenta o complexo conj. de psi
    f = conj(psi(j,n)) .* (V(j-1)-V(j+1))/(2*dx) .* psi(j,n);
  endfunction

  # Calcular o valor esperado de -grad(V):
  for n=1:length(t); 
    for j=2:length(x)-1
      grad_V(n) = grad_V(n) + integrando_gradV(psi,j,n,x,V,dx)*dx; # <-grad(V)> para cada t
    endfor
  endfor

  # Grafico de dPdt por -grad(V) (teo. de Ehrenfest pt. 2):
  figure
  p1 = plotyy(t,dPdt,t,grad_V);   # plota d<p>/dt e -grad(V) x t
  titulo = strcat('d<p>/dt e <-grad(V)>. Estado=',num2str(e));# gera uma string
  title(titulo)                             # coloca um título
  xlabel('t');                              # coloca uma legenda no eixo x
  ylabel(p1(1),'d<p>/dt');                  # coloca uma legenda no eixo y a esquerda
  ylabel(p1(2),'<-grad(V)>');                      # coloca uma legenda no eixo y a direita
  legenda = legend('d<p>/dt','<-grad(V)>','location','northwest');# coloca uma legenda
  set(legenda,'FontSize',22)                # tamanho da legenda
  set(p1,'FontSize',20)                     # tamanho dos numeros nos eixos



  #***********EVOLUCAO TEMPORAL DE PSI: ***************************************************

  # Animacao da evolucao temporal de psi:
  figure1 = figure(10);
  set(figure1,'color','white');
  winsize = get(figure1,'Position');
  winsize(1:2) = [0 0];
  clear M;
  count=1;
  p1=plotyy(x,psi2(:,1),x,V);               # plota psi^2 e V
  titulo = strcat('t=',num2str(t(1)));      # gera uma string
  title(titulo);                            # coloca um titulo
  xlabel('x');                              # coloca uma legenda no eixo x
  ylabel(p1(1),'|psi(x,t)|**2');            # coloca uma legenda no eixo y da esquerda
  ylabel(p1(2),'V(x)');                     # coloca uma legenda no eixo y da direita
  ylim([0 6]);                              # mantem a escala do grafico constante
  #xlim([-1 1]);                            # mantem a escala do grafico constante
  set(p1,'FontSize',20)                     # tamanho dos numeros nos eixos
  pause(10^-4);                             # tempo entre cada imagem
  M(count)=getframe(figure1);               # armazena a frame em 'M' - Octave
  delete(p1);                               # deleta o gráfico (ja armazenou em 'M')
  count=count+1;                            # atualiza o contador  
  for j=1:20:length(t)                        # coluna = tempo
    p1=plotyy(x,psi2(:,j),x,V);               # plota psi^2 e V
    titulo = strcat('t=',num2str(t(j)));      # gera uma string
    title(titulo);                            # coloca um titulo
    xlabel('x');                              # coloca uma legenda no eixo x
    ylabel(p1(1),'|psi(x,t)|**2');            # coloca uma legenda no eixo y da esquerda
    ylabel(p1(2),'V(x)');                     # coloca uma legenda no eixo y da direita
    ylim([0 6]);                              # mantem a escala do grafico constante
    #xlim([-1 1]);                            # mantem a escala do grafico constante
    set(p1,'FontSize',20)                     # tamanho dos numeros nos eixos
    pause(10^-4);                             # tempo entre cada imagem
    M(count)=getframe(figure1);               # armazena a frame em 'M' - Octave
    delete(p1);                               # deleta o gráfico (ja armazenou em 'M')
    count=count+1;                            # atualiza o contador  
  endfor

endif

#***********FIM *************************************************************************