import csv
import matplotlib.pyplot as plt
import math ## script comparacao dos dados do mesmo sensor em leituras diferentes


## vetores da leitura 2
temp1 = []
vel   = []
altura = [] 
vol = []
temp2 = []
#########################################################################################################################################
##                        VARIAVEL CÁLCULO VEL SAÍDA DA ÁGUA , TEMPO TOTAL DE SAÍDA E TAXA DE SAÍDA dm/dt     
#########################################################################################################################################
v0  = (2-0.765) ## volume de ar zinicial e final na garrafa
vf  = 2
h   = 0.0001 # variacao de temp usada na integracao 
pa  = 99900 # pressao atmosférica
psi = 65.2
p0  = 6894.7572932 * psi  #pressao inicial dentro da garrafa
lvol = int((vf-v0)/h) #quantidade de interacoes l para o calculo do tempo esperado para saída completa da agua
A   = 0.000615752 # area da sessao do pico da garrafa para calculo da variação da pressao
j   = 1.4  # expoente expansao adiábitica
rho = 1000 # densidade da água
k1  = 0
k2  = 0
k3  = 0 
k4  = 0
cont = 0
tfa  = 0 # tempo de exaustao da agua do foguete
def f_vol(t, y):
			p = p0*((v0/y))
			return (y + h*A*((2*(p-pa)/rho)**j))
def calcK1_vol(h,xn,yn):
  return h*(f_vol(xn,yn))
  
def calcK2_vol(h, xn, yn):
	return h*(f_vol(xn + 0.5*h, yn + 0.5*k1))

def calcK3_vol(h, xn, yn):
	return h*f_vol(xn + 0.5*h, yn + 0.5*k2)

def calcK4_vol(h, xn, yn):
	return h*f_vol(xn + h, yn + k3)
	
	
for cont in range(lvol+1):
	if cont == 0:
		k1 = calcK1_vol(h,0, v0)
		k2 = calcK2_vol(h,0, v0)	
		k3 = calcK3_vol(h,0, v0)
		k4 = calcK4_vol(h,0, v0)
		ynMais1 = v0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
		temp2.append(cont*h)
		vol.append(ynMais1)
	else:
		k1 = calcK1_vol(h,cont*h, yn)
		k2 = calcK2_vol(h,cont*h, yn)	
		k3 = calcK3_vol(h,cont*h, yn)
		k4 = calcK4_vol(h,cont*h, yn)
		ynMais1 = yn + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
		temp2.append(cont*h)
		vol.append(ynMais1)
		if  1.9999999 < ynMais1  :
			tfa = temp2[cont]
			print(tfa)
			print(ynMais1)
			break
	cont = cont +1
	yn = ynMais1
		
			
#########################################################################################################################################
somaAlt=0 # variavel armazena altura
u = rho # densidade da agua 
pi = p0 # pressao inicial dentro da garrafa
pf = pa # pressao final dentro da garrafa
m0 = 1.04# massa do foguete + massa da agua
mf =0.275 # massa do foguete # tempo total de descarga da agua 
tf = 4  # tempo do lançamento considerado
R = (m0-mf)/tfa # taxa de descarga de agua dm/dt
g = 9.81 # acel da gravidade 
coluAgua = 0.15 # altura da coluna de agua dentro da garrafa
ve = ((2*(pi - pf)/u - 2*g*coluAgua)**(1/2)) # calculo da velocidade de saída da água de dentro do foguete
print(ve)
k = 0.00055724 # constante da força de arrasto 
l= int(tf/h) #quantidade de interacoes l no intervalo de tempo considerado [t0, tf]

y0 = 0 # velocidade inicial 
yn = 0  # velocidade no instante n 
k1 = 0
k2 = 0
k3=0 
k4 = 0
cont_alt = 0## contadores pra saltar linhas com texto

def f_alt(t, y):
		if t <= tfa:
			return (R*ve/(m0-R*t) -g - k*y*y/(m0-R*t))
		else:
			return ( -g -k*y*y/mf)
			
def calcK1_alt(h,xn,yn):
  return h*(f_alt(xn,yn))
  
def calcK2_alt(h, xn, yn):
	return h*(f_alt(xn + 0.5*h, yn + 0.5*k1))

def calcK3_alt(h, xn, yn):
	return h*f_alt(xn + 0.5*h, yn + 0.5*k2)

def calcK4_alt(h, xn, yn):
	return h*f_alt(xn + h, yn + k3)
	
	
for cont_alt in range(l+1):
	if cont_alt == 0:
		k1 = calcK1_alt(h,0, y0)
		k2 = calcK2_alt(h,0, y0)	
		k3 = calcK3_alt(h,0, y0)
		k4 = calcK4_alt(h,0, y0)
		ynMais1 = y0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*0.984807753
		temp1.append(cont_alt*h)
		vel.append(ynMais1)
		somaAlt = (ynMais1-y0)*0.5*h
		altura.append(somaAlt)
	else:
		k1 = calcK1_alt(h,cont_alt*h, yn)
		k2 = calcK2_alt(h,cont_alt*h, yn)	
		k3 = calcK3_alt(h,cont_alt*h, yn)
		k4 = calcK4_alt(h,cont_alt*h, yn)
		 
		ynMais1 = yn + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*0.984807753
		temp1.append(cont_alt*h)
		vel.append(ynMais1)
		somaAlt+= (yn*h + (ynMais1-yn)*h*0.5) # integra velocidade e armazena em somaAlt
		altura.append(somaAlt)
	cont_alt = cont_alt +1
	yn = ynMais1
	
	
plt.subplot(2, 1, 1)
plt.plot(temp1, altura)
plt.title('altura x tempo ')
plt.ylabel('altura/aceleracao')
plt.xlabel('time ')
plt.grid()
plt.show()
