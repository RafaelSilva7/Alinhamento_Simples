import numpy as np
from matriz_score.matrizScore import matriz_score


def alignment_local(sequence_a, sequence_b, gap_penalty, matriz):
    """
    Realiza o alinhamento local entre duas sequencias de nucleotideos
    :param sequence_a: primeira sequencia
    :param sequence_b: segunda sequencia
    :param gap_penalty: valor do gap penalty
    :param matriz: nome da matriz de score a ser utilizado
    :return: retorna uma string com o melhor alinhamento possivel
    """
    i = len(sequence_a) + 1
    j = len(sequence_b) + 1

    f = np.zeros((i, j), dtype=np.int)

    # Preenche a primeira linha e primeira coluna com os valores base
    for i in range(0, len(sequence_a)+1):
        f[i, 0] = 0
    for j in range(0, len(sequence_b)+1):
        f[0, j] = 0

    # Realiza o preenchimento dos valores de score da matriz pontuacoes
    for i in range(1, len(sequence_a)+1):
        for j in range(1, len(sequence_b)+1):
            value1 = f[i-1, j-1] + matriz_score(sequence_a[i-1], sequence_b[j-1], matriz)
            value2 = f[i-1, j] - gap_penalty
            value3 = f[i, j-1] - gap_penalty
            f[i, j] = max(value1, value2, value3, 0)

    # Processo de anlinhamento das sequencias (traceback)
    alignment_a = ""
    alignment_b = ""
    str_aux = ""
    num_score = 0

    # Encontra o maior valor
    greater = 0
    for i in range(1, len(sequence_a) + 1):
        for j in range(1, len(sequence_b) + 1):
            if greater < f[i, j]:
                greater = f[i, j]
                x = i
                k = j

    while x > 0 or k > 0:
        if f[x, k] == 0:
            break

        # Caso I, ir para a diagonal
        if (k > 0 and x > 0) and (f[x, k] == f[x-1, k-1] + matriz_score(sequence_a[x-1], sequence_b[k-1], matriz)):
            alignment_b = sequence_b[k-1] + alignment_b
            alignment_a = sequence_a[x-1] + alignment_a
            if sequence_a[x-1] == sequence_b[k-1]:
                str_aux = "|" + str_aux
                num_score += 1
            else:
                str_aux = " " + str_aux
            x -= 1
            k -= 1
        # Caso II, ir para cima
        elif (x > 0 and k == 0) or (f[x][k] == (f[x-1, k] - gap_penalty) and x > 0):
            alignment_b = "-" + alignment_b
            alignment_a = sequence_a[x-1] + alignment_a
            str_aux = " " + str_aux
            x -= 1
        # Caso III, ir para a esquerda
        else:
            alignment_b = sequence_b[k - 1] + alignment_b
            alignment_a = "-" + alignment_a
            str_aux = " " + str_aux
            k = k - 1

    str_aux = alignment_a + "\n" + str_aux + "    score: " + str(num_score) + "\n" + alignment_b
    return "Melhor Alinhamento:\n" + str_aux
