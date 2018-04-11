def matriz_score(x, y, matriz):
    """
    Calcula o valor de score (match/mismatch) do alinhamento entre os dois aminoacidos
    com base no valores da Blosum62
    :param x: primeiro aminoacido (linha)
    :param y: segundo aminoacido (coluna)
    :param matriz: nome da tabela de score a ser utilizado pela funcao
    :return: retorna o valor do score (na forma de inteiro)
    """
    if matriz == "global":
        return score_global(x, y)

    local = "./matrizes/" + matriz + ".txt"
    with open(local) as file:

        amino_acids = " ARNDCQEGHILKMFPSTWYVBZX*"
        i = 0
        while amino_acids[i] != x:
            i += 1

        j = 0
        while amino_acids[j] != y:
            j += 1

        j *= 3
        for k in range(0, i + 1):
            amino_acids = file.readline()

        score = amino_acids[j - 1] + amino_acids[j]
    return int(score)


def score_global(x, y):
    """
    Verifica o valor de scoring dos aminoacidos informados
    :param x: primeiro aminoacido
    :param y: segundo aminoacido
    :return: valor de scoring
    """
    if x == y:
        return 5
    return -4
