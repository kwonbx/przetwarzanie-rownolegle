#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define N 17
#define K 42
#define KOD_LEN 136 // n*(n-1)/2

// Funkcja sprawdzająca spójność grafu za pomocą BFS
int is_connected(int adj[N][N]) {
    int visited[N] = {0};
    int queue[N];
    int head = 0, tail = 0;

    // Startujemy od wierzchołka 0
    visited[0] = 1;
    queue[tail++] = 0;

    while (head < tail) {
        int v = queue[head++];
        for (int i = 0; i < N; i++) {
            if (adj[v][i] && !visited[i]) {
                visited[i] = 1;
                queue[tail++] = i;
            }
        }
    }

    // Sprawdzamy czy odwiedzono wszystkie N wierzchołków
    return (tail == N);
}

void generate_connected_to_file(const char *filename, long count) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Błąd zapisu");
        return;
    }

    int edges[KOD_LEN];
    int adj[N][N];
    char kod[KOD_LEN + 1];

    printf("Generowanie %ld spójnych grafów do pliku %s...\n", count, filename);

    long found = 0;
    while (found < count) {
        // 1. Inicjalizacja krawędzi (dokładnie K jedynek)
        for (int i = 0; i < KOD_LEN; i++) edges[i] = (i < K) ? 1 : 0;

        // 2. Tasowanie krawędzi (Fisher-Yates)
        for (int i = KOD_LEN - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            int temp = edges[i];
            edges[i] = edges[j];
            edges[j] = temp;
        }

        // 3. Budowa macierzy sąsiedztwa do sprawdzenia spójności
        memset(adj, 0, sizeof(adj));
        int idx = 0;
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                if (edges[idx] == 1) {
                    adj[i][j] = adj[j][i] = 1;
                }
                idx++;
            }
        }

        // 4. Jeśli graf jest spójny, zapisz go w formacie "Kod"
        if (is_connected(adj)) {
            for (int i = 0; i < KOD_LEN; i++) {
                kod[i] = edges[i] ? '1' : '0';
            }
            kod[KOD_LEN] = '\0';
            fprintf(f, "%s\n", kod);
            found++;

            if (found % 10000 == 0) {
                printf("Postęp: %ld/%ld...\n", found, count);
            }
        }
    }

    fclose(f);
    printf("Gotowe! Zapisano %ld grafów spójnych.\n", count);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Użycie: %s <nazwa_pliku> <liczba_grafow>\n", argv[0]);
        return 1;
    }

    srand(time(NULL));
    generate_connected_to_file(argv[1], atol(argv[2]));
    return 0;
}