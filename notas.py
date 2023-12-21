#!/usr/bin/env python3

# Número de notas a serem inseridas
num_notas = int(input("Quantas notas você deseja inserir? "))

# Inicialização de variáveis
total = 0

# Loop para receber as notas
for i in range(num_notas):
    try:
        # Entrada da nota do usuário
        nota = float(input(f"Insira a nota {i + 1}: "))
                                    
        # Verifica se a nota está dentro do intervalo válido
        if 0 <= nota <= 10:
        # Soma a nota ao total
            total += nota
        else:
            print("Nota inválida. Deve estar entre 0 e 10.")
            i -= 1  

    except ValueError:
        print("Entrada inválida. Insira um número válido.")
        i -= 1  

# Cálculo da média
media_disciplina = total / num_notas
           
# Impressão da média na tela
print(f"A média da disciplina é: {media_disciplina}")
