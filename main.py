from needleman.needleman import alignment_global
from waterman.waterman import alignment_local
import time

def main():
    """
    Funcao main
    :return: void
    """
    for tam in range(40,105,10):
        with open("./sequence_protein1.txt", 'r') as file:
            sequence_a = file.readline()
            sequence_a = file.read(tam)
            sequence_a = sequence_a.replace("\n","")
            file.close()

        with open("./sequence_protein2.txt", 'r') as file:
            sequence_b = file.readline()
            sequence_b = file.read(tam)
            sequence_b = sequence_b.replace("\n","")
            file.close()

        file = open("./saida/alignment_protein" + str(tam) + ".txt", 'w')
        file.write("Entrada (tamanho[" + str(tam) + "]):\n1-  " + sequence_a + "\n2-  " + sequence_b + "\n\n\n")
        matriz1 = "global"
        matriz2 = "blosum62"
        
        time_start = time.time()
        alignment = alignment_global(sequence_a, sequence_b, 2, matriz1)
        time_end = time.time()

        header = "Algoritmo: needleman-wunsch (alinhamento global)\n" \
                 "Match/Mismatch: " + matriz1 + "    Gap: -2\n" \
                 "Tempo execucao: " + str(time_end - time_start) + "s\n"
        file.write(header + alignment)

        time_start = time.time()
        alignment = alignment_global(sequence_a, sequence_b, 2, matriz2)
        time_end = time.time()
        
        header = "\n\n\nAlgoritmo: Smith-Waterman (alinhamento local)\n" \
                 "Match/Mismatch: " + matriz2 + "    Gap: -2\n" \
                 "Tempo execucao: " + str(time_end - time_start) + "s\n"
        file.write(header + alignment)

        file.close()


if __name__ == '__main__':
    main()
