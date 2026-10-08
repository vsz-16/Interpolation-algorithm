import numpy as np
import sympy as sp #serve pra fazer as coisas de calculo
import math #pra agilizar com o fatorial

def Entrada():
    n = int(input("Qual a ordem de seu polinômio?\n"))
    #Recebe já no formato np. Split() separa os valores da string
    x = np.array(input(f"Digite os {n+1} valores de X:\n").split(), dtype = float)
    y = np.array(input(f"\nDigite os {n+1} valores de Y:\n").split(), dtype = float)
    z = float(input("\nDigite o valor que quer interpolar:\n"))
    return x, y, z

def Polinomio(x, y, z, expressao=None):
    n = len(x) #Pra saber a ordem do polinomio sem precisar perguntar
    l = np.zeros(n) #Esse vetor guarda os resultados do produtório
    for i in range(0, n):
        Xj = np.delete(x, i) #Tira o Xi do vetor x
        numerador = z - Xj
        denominador = x[i] - Xj
        l[i] = np.prod(numerador / denominador)
    resultado = np.sum(y * l)

    if expressao:
        Truncamento(x, z, expressao)
    return resultado

def DispositivoL(x, y, z, expressao=None):
    n = len(x)
    g = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j: #Pra preencher a diagonal principal
                g[i, j] = z - x[i]
            else:
                g[i, j] = x[i] - x[j]
    gi = np.prod(g, axis = 1)
    gd = np.prod(np.diag(g)) #O np.diag() cria um vetor unidimensional da diagonal da matriz, assim dá pra fazer o produtóriod dela, como no outro método
    parcela = np.sum(y/gi)
    resultado = gd * parcela

    if expressao:
        Truncamento(x, z, expressao)
    return resultado

def Newton(x, y, z, expressao=None):
    n = len(x)

    t = np.zeros((n, n)) #Matriz que vai servir como tabela pra guardar as dd
    t[:, 0] = y
    for j in range(1, n): #Percorre as colunas que representam as ordens. Ficou ao contrário por conta da recurssividade da conta
        for i in range(n - j): # n-j pois a matriz/tabela é triangular
            t[i, j] = (t[i+1, j-1] - t[i, j-1]) / (x[i+j] - x[i]) #Conta das diferenças divididas

    resultado = t[0, 0]
    produtox = 1.0

    for i in range(1, n):
    # Acumula o produtório (X - X0)(X - X1)...
        produtox = produtox * (z - x[i-1])
        # Adiciona o novo termo à soma total
        resultado = resultado + (t[0, i] * produtox)

    if expressao:
        Truncamento(x, z, expressao)

    return resultado

def Gregory(x, y, z, expressao = None):
    n = len(x)
    h = x[1] - x[0]
    
    # verifica se todos os X estão realmente espaçados por h
    if not np.allclose(np.diff(x), h): #np.diff(X) do NumPy subtrai todos os vizinhos do vetor de uma vez. O np.allclose verifica se todas essas subtrações resultaram num valor idêntico ao h
        print("Erro: O método de Gregory-Newton exige valores de X igualmente espaçados.")
        return None
    
    # Cria a Tabela de diferenças
    d = np.zeros((n, n))
    d[:, 0] = y
    for j in range(1, n): # Constrói a tabela fazendo APENAS subtrações
        for i in range(n - j):
            d[i, j] = d[i+1, j-1] - d[i, j-1]

    u = (z - x[0]) / h # Aplica a mudança de variável u_x = (X - X0) / h
    
    # Calcula o polinômio
    # Pn(x) = Y0 + sum( (Delta^i Y_0 / i!) * prod(u - j) )
    resultado = d[0, 0]
    produtou = 1.0
    fatorial = 1.0
    
    for i in range(1, n):
        produtou = produtou * (u - (i - 1))
        fatorial = fatorial * i
        resultado = resultado + ((d[0, i] / fatorial) * produtou) # Adiciona a parcela à soma total

    if expressao:
        Truncamento(x, z, expressao)
    return resultado

def Truncamento(X, z, expressao_funcao):
    m = len(X) - 1 # O grau do polinômio é o número de pontos menos 1
    
    # Define 'x' como uma variável algébrica
    x = sp.Symbol('x')
    f = sp.sympify(expressao_funcao)
    
    # 1. Calcula a derivada contínua de ordem (m+1)
    derivada = sp.diff(f, x, m + 1)
    
    # 2. Encontra o maior valor em módulo da derivada no intervalo [X_0, X_m]
    # O np.linspace varre o intervalo testando pontos para achar o pico máximo
    pontos = np.linspace(X[0], X[-1], 200)
    max_derivada = max([abs(float(derivada.subs(x, p))) for p in pontos])
    
    # 3. Calcula o produtório das distâncias: (x - X0)(x - X1)...(x - Xm)
    produtorio = 1.0
    for xi in X:
        produtorio *= (z - xi)
        
    # 4. Aplica a fórmula completa do erro
    erro = (max_derivada / math.factorial(m + 1)) * abs(produtorio)
    
    print(f"Erro de Truncamento T_{m}({z}): {erro:.4f}")
    
    return erro

def Main():
    while True:
        print("=== MENU PRINCIPAL ===")
        print("1- Polinômio de Lagrange")
        print("2- Dispositivo prático de Lagrange")
        print("3- Método de Newton")
        print("4- Método de Gregory-Newton")
        print("0- Sair do programa")
        op = input("Escolha a opção desejada\n")

        if op in ["1", "2", "3", "4"]:
            calc_erro = input("\nVocê tem a função original para calcular o erro de truncamento? (s/n)\n").strip().lower()
            expressao = None
            if calc_erro == 's':
                expressao = input("Digite a expressão da função:\n")
        if op == "1":
            x, y, z = Entrada()
            print(f"O resultado da interpolação é: {Polinomio(x, y, z, expressao)}")
        if op == "2":
            x, y, z = Entrada()
            print(f"O resultado da interpolação é: {DispositivoL(x , y, z, expressao)}")
        if op == "3":
            x, y, z = Entrada()
            print(f"O resultado da interpolação é: {Newton(x, y, z, expressao)}")
        if op == "4":
            x, y, z = Entrada()
            print(f"O resultado da interpolação é: {Gregory(x, y, z, expressao)}")
        elif op == "0":
            print("Encerrando programa...")
            break
Main()