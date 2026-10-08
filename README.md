# Interpolação Polinomial em Python

Programa de linha de comando que realiza **interpolação polinomial** a partir de um conjunto de pontos `(x, y)` e estima o valor da função em um ponto desejado. Implementa quatro métodos clássicos de Cálculo Numérico e, opcionalmente, calcula o **erro de truncamento** quando a função original é conhecida.

## Métodos implementados

| # | Método | Função | Exige espaçamento uniforme entre os dados? |
|---|--------|--------|-----------------------------|
| 1 | Polinômio de Lagrange | `Polinomio` | Não |
| 2 | Dispositivo prático de Lagrange | `DispositivoL` | Não |
| 3 | Método de Newton (diferenças divididas) | `Newton` | Não |
| 4 | Método de Gregory-Newton (diferenças finitas) | `Gregory` | **Sim** |

## Requisitos

- Python 3.8+
- [NumPy](https://numpy.org/)
- [SymPy](https://www.sympy.org/)

```bash
pip install numpy sympy
```

O módulo `math` já faz parte da biblioteca padrão do Python.

## Como executar

```bash
python interpolacao.py
```

(substitua `interpolacao.py` pelo nome do seu arquivo)

O programa exibe um menu:

```
=== MENU PRINCIPAL ===
1- Polinômio de Lagrange
2- Dispositivo prático de Lagrange
3- Método de Newton
4- Método de Gregory-Newton
0- Sair do programa
```

Depois de escolher um método (1 a 4), o programa pergunta se você possui a função original para calcular o erro de truncamento. Em seguida, solicita:

1. **A ordem do polinômio** `n` (serão necessários `n + 1` pontos);
2. Os **`n + 1` valores de X**, separados por espaço;
3. Os **`n + 1` valores de Y**, separados por espaço;
4. O **valor a interpolar** `z`.

### Exemplo de uso

Interpolando `f(x) = x²` com três pontos, no ponto `z = 2.5`:

```
=== MENU PRINCIPAL ===
...
Escolha a opção desejada
1

Você tem a função original para calcular o erro de truncamento? (s/n)
s
Digite a expressão da função:
x**2
Qual a ordem de seu polinômio?
2
Digite os 3 valores de X:
1 2 3
Digite os 3 valores de Y:
1 4 9

Digite o valor que quer interpolar:
2.5
Erro de Truncamento T_2(2.5): 0.0000
O resultado da interpolação é: 6.25
```

> A expressão da função deve seguir a sintaxe do Python/SymPy (por exemplo `x**2`, `sin(x)`, `exp(x)`).

## Descrição dos métodos

### 1. Polinômio de Lagrange (`Polinomio`)

Calcula cada polinômio base de Lagrange:

```
L_i(z) = ∏_{j≠i} (z - x_j) / (x_i - x_j)
```

e monta o polinômio interpolador:

```
P(z) = Σ y_i · L_i(z)
```

No código, `np.delete` remove o `x_i` do vetor, e `np.prod` calcula o produtório de forma vetorizada.

### 2. Dispositivo prático de Lagrange (`DispositivoL`)

Variação organizada em forma de matriz. Constrói uma matriz `g` em que:

- a **diagonal principal** guarda `z - x_i`;
- as demais posições guardam `x_i - x_j`.

O resultado é obtido por:

```
P(z) = [∏ (z - x_i)] · Σ ( y_i / ∏_j g_ij )
```

onde o produtório das linhas (`gi`) forma os denominadores e o produtório da diagonal (`gd`) é o fator comum a todos os termos.

### 3. Método de Newton (`Newton`)

Usa a **tabela de diferenças divididas**, construída por:

```
f[x_i, ..., x_{i+j}] = (f[x_{i+1}, ..., x_{i+j}] - f[x_i, ..., x_{i+j-1}]) / (x_{i+j} - x_i)
```

O polinômio é escrito na forma de Newton:

```
P(z) = f[x0] + f[x0,x1](z - x0) + f[x0,x1,x2](z - x0)(z - x1) + ...
```

Os coeficientes são a primeira linha da tabela triangular. Não exige pontos igualmente espaçados.

### 4. Método de Gregory-Newton (`Gregory`)

Versão para pontos **igualmente espaçados** (passo `h`). Primeiro verifica o espaçamento com `np.diff` e `np.allclose`; se os pontos não forem uniformes, exibe uma mensagem de erro e retorna `None`.

Constrói a tabela de **diferenças finitas** (apenas subtrações) e aplica a mudança de variável:

```
u = (z - x0) / h
```

O polinômio é calculado por:

```
P(z) = y0 + Σ ( Δ^i y0 / i! ) · u(u - 1)...(u - i + 1)
```

## Erro de truncamento (`Truncamento`)

Quando a função original é informada, o programa estima o limitante do erro de truncamento:

```
|E(z)| ≤ ( M / (m+1)! ) · |(z - x0)(z - x1)...(z - xm)|
```

em que:

- `m` é o grau do polinômio (número de pontos − 1);
- `M` é o máximo de `|f^(m+1)(x)|` no intervalo `[x0, xm]`.

Para isso, o **SymPy** calcula simbolicamente a derivada de ordem `m + 1`, e o `np.linspace` avalia essa derivada em 200 pontos do intervalo para encontrar o maior valor em módulo. O resultado é exibido no formato `T_m(z)`.

## Estrutura do código

| Função | Responsabilidade |
|--------|------------------|
| `Entrada()` | Lê a ordem, os vetores X e Y e o ponto `z` a interpolar |
| `Polinomio()` | Interpolação de Lagrange |
| `DispositivoL()` | Dispositivo prático de Lagrange |
| `Newton()` | Interpolação por diferenças divididas |
| `Gregory()` | Interpolação por diferenças finitas |
| `Truncamento()` | Cálculo do erro de truncamento |
| `Main()` | Menu principal e laço de execução |

## Observações

- Os valores de X devem ser **distintos** entre si (caso contrário há divisão por zero).
- O método de Gregory-Newton só funciona com X igualmente espaçados.
- O erro de truncamento é um **limitante** estimado por amostragem de 200 pontos; não é o erro exato.
- Para interpolar, `z` deve estar preferencialmente dentro do intervalo `[x0, xm]`; fora dele o processo vira extrapolação e o erro tende a crescer.
