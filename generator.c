#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define N 17
#define N1 7
#define N2 10
#define K 42
#define KOD_LEN 136

// Mapowanie krawędzi bloku bipartite (7x10) na indeks w kodzie 136-bitowym
int get_kod_idx(int r, int c) {
    int i = r;
    int j = c + N1;
    return (i * N - (i * (i + 1)) / 2) + (j - i - 1);
}

void generate_to_file(const char *filename, long count) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Nie mozna otworzyc pliku do zapisu");
        return;
    }

    int total_slots = N1 * N2; // 70 mozliwych krawedzi w grafie dwudzielnym 7x10
    int *edges = malloc(total_slots * sizeof(int));
    char kod[KOD_LEN + 1];

    printf("Generowanie %ld losowych grafow do pliku %s...\n", count, filename);

    for (long g = 0; g < count; g++) {
        // 1. Przygotuj tablice z dokladnie K krawedziami
        for (int i = 0; i < total_slots; i++) {
            edges[i] = (i < K) ? 1 : 0;
        }

        // 2. Potasuj krawedzie (Fisher-Yates shuffle)
        for (int i = total_slots - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            int temp = edges[i];
            edges[i] = edges[j];
            edges[j] = temp;
        }

        // 3. Buduj kod 136-znakowy
        memset(kod, '0', KOD_LEN);
        kod[KOD_LEN] = '\0';

        for (int i = 0; i < total_slots; i++) {
            if (edges[i] == 1) {
                int r = i / N2;
                int c = i % N2;
                kod[get_kod_idx(r, c)] = '1';
            }
        }

        // 4. Zapisz do pliku
        fprintf(f, "%s\n", kod);

        // Wyswietlaj postep co 10 000 grafow
        if (g % 10000 == 0 && g > 0) {
            printf("Postep: %ld/%ld...\n", g, count);
        }
    }

    fclose(f);
    free(edges);
    printf("Gotowe! Zapisano %ld grafow.\n", count);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Uzycie: %s <nazwa_pliku> <liczba_grafow>\n", argv[0]);
        printf("Przyklad: %s grafy.txt 100000\n", argv[0]);
        return 1;
    }

    srand(time(NULL));
    const char *filename = argv[1];
    long count = atol(argv[2]);

    generate_to_file(filename, count);

    return 0;
}