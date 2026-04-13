#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N 17
#define KOD_LEN 136
#define EPS 1e-6  // Poluzowana precyzja dla liczb całkowitych

void get_eigenvalues(double m[N][N], double evals[N]) {
    double v[N][N];
    for(int i=0; i<N; i++) for(int j=0; j<N; j++) v[i][j] = (i==j);

    // Zwiększona liczba iteracji dla lepszej zbieżności
    for(int it=0; it<500; it++) {
        double max_v = 0;
        int p=0, q=1;
        for(int i=0; i<N; i++) {
            for(int j=i+1; j<N; j++) {
                if(fabs(m[i][j]) > max_v) { 
                    max_v = fabs(m[i][j]); 
                    p=i; q=j; 
                }
            }
        }
        if(max_v < 1e-12) break; // Macierz jest już niemal diagonalna

        double theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
        double t;
        if (fabs(theta) < 1e10) {
            t = (theta >= 0) ? 1.0/(theta + sqrt(1+theta*theta)) : -1.0/(-theta + sqrt(1+theta*theta));
        } else {
            t = 0.5 / theta;
        }
        
        double c = 1.0 / sqrt(1+t*t);
        double s = t*c;
        double tau = s/(1.0+c);
        double apq = m[p][q];
        
        m[p][p] -= t*apq; 
        m[q][q] += t*apq; 
        m[p][q] = m[q][p] = 0;
        
        for(int i=0; i<N; i++) {
            if(i!=p && i!=q) {
                double aip = m[i][p], aiq = m[i][q];
                m[i][p] = m[p][i] = aip - s*(aiq + aip*tau);
                m[i][q] = m[q][i] = aiq + s*(aip - aiq*tau);
            }
        }
    }
    for(int i=0; i<N; i++) evals[i] = m[i][i];
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Użycie: %s <plik_z_grafami>\n", argv[0]);
        return 1;
    }

    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("Błąd pliku"); return 1; }

    char kod[1024]; // Bufor większy niż KOD_LEN dla bezpieczeństwa
    double adj[N][N], evals[N];
    int id = 1;

    while (fscanf(f, "%136s", kod) == 1) {
        if (strlen(kod) != KOD_LEN) continue;

        // Czyszczenie i budowa macierzy
        for(int i=0; i<N; i++) for(int j=0; j<N; j++) adj[i][j] = 0.0;
        
        int char_idx = 0;
        for(int i=0; i<N; i++) {
            for(int j=i+1; j<N; j++) {
                if(kod[char_idx++] == '1') {
                    adj[i][j] = adj[j][i] = 1.0;
                }
            }
        }

        get_eigenvalues(adj, evals);

        int is_integral = 1;
        for(int i=0; i<N; i++) {
            double rounded = round(evals[i]);
            if(fabs(evals[i] - rounded) > EPS) {
                is_integral = 0;
                break;
            }
        }

        if(is_integral) {
            // Wypisujemy dokładnie w formacie, o który prosiłeś
            printf("{\"id\":\"%d\",\"Graph6\":\"P???????????\",\"Kod\":\"%s\"},\n", id++, kod);
            fflush(stdout);
        }
    }

    fclose(f);
    return 0;
}