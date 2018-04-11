import numpy as np
from matriz_score.matrizScore import matriz_score


def alignment_global(sequence_a, sequence_b, gap_penalty, matriz):
    """
    Realiza o alinhamento global de duas sequencias de aminoacidos
    :param sequence_a: primeria sequencia
    :param sequence_b: segunda sequencia
    :param gap_penalty: valor do gap penalty
    :param matriz: nome da matriz de score a ser utilizada
    :return: retorna uma string com o melhor alinhamento possivel
    """
    i = len(sequence_a) + 1
    j = len(sequence_b) + 1

    f = np.zeros((i, j), dtype=np.int)

    # Preenche a primeira linha e primeira coluna valor base
    for i in range(0, len(sequence_a) + 1):
        f[i, 0] = (gap_penalty * i) * -1
    for j in range(0, len(sequence_b) + 1):
        f[0, j] = (gap_penalty * j) * -1

    # Realiza o preenchimento dos valores de score da matriz pontuacoes
    for i in range(1, len(sequence_a) + 1):
        for j in range(1, len(sequence_b) + 1):
            value1 = f[i-1, j-1] + matriz_score(sequence_a[i-1], sequence_b[j-1], matriz)
            value2 = f[i-1, j] - gap_penalty
            value3 = f[i, j-1] - gap_penalty
            f[i, j] = max(value1, value2, value3)

    # Processo de anlinhamento das sequencias (traceback)
    alignment_a = ""
    alignment_b = ""
    str_aux = ""
    num_score = 0

    while i > 0 or j > 0:
        # Caso I, ir para a diagonal
        if (j > 0 and i > 0) and (f[i, j] == f[i-1, j-1] + matriz_score(sequence_a[i-1], sequence_b[j-1], matriz)):
            alignment_b = sequence_b[j-1] + alignment_b
            alignment_a = sequence_a[i-1] + alignment_a
            if sequence_a[i-1] == sequence_b[j - 1]:
                str_aux = "|" + str_aux
                num_score += 1
            else:
                str_aux = " " + str_aux
            i -= 1
            j -= 1
        # Caso II, ir para cima
        elif (i > 0 and j == 0) or (f[i][j] == (f[i-1, j] - gap_penalty) and i > 0):
            alignment_b = "-" + alignment_b
            alignment_a = sequence_a[i-1] + alignment_a
            str_aux = " " + str_aux
            i -= 1
        # Caso III, ir para a esquerda
        else:
            alignment_b = sequence_b[j-1] + alignment_b
            alignment_a = "-" + alignment_a
            str_aux = " " + str_aux
            j = j - 1

    str_aux = alignment_a + "\n" + str_aux + "    score: " + str(num_score) + "\n" + alignment_b
    return "Melhor Alinhamento:\n" + str_aux
