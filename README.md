#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 3 // number of equations

// Define the system of ODEs
void system(double t, double y[], double dydt[], double a, double b, double c, double A, double C) {
    dydt[0] = y[1];                               // dy1/dt = y2
    dydt[1] = y[2];                               // dy2/dt = y3
    dydt[2] = -a * y[2] - b * y[1] - c * y[0] + C + A * t; // dy3/dt
}

// Runge-Kutta 4th Order Method
void rk4(double t, double y[], double h, double a, double b, double c, double A, double C) {
    double k1[N], k2[N], k3[N], k4[N], y_temp[N];
    int i;

    system(t, y, k1, a, b, c, A, C);
    for (i = 0; i < N; i++) y_temp[i] = y[i] + 0.5 * h * k1[i];
    system(t + 0.5 * h, y_temp, k2, a, b, c, A, C);
    for (i = 0; i < N; i++) y_temp[i] = y[i] + 0.5 * h * k2[i];
    system(t + 0.5 * h, y_temp, k3, a, b, c, A, C);
    for (i = 0; i < N; i++) y_temp[i] = y[i] + h * k3[i];
    system(t + h, y_temp, k4, a, b, c, A, C);

    for (i = 0; i < N; i++)
        y[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

void plot_results() {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set terminal pngcairo size 800,600\n");
        fprintf(gnuplotPipe, "set output 'computed_curves.png'\n");
        fprintf(gnuplotPipe, "set multiplot layout 3,1 title 'Computed Curves' font ',14'\n");
        
        fprintf(gnuplotPipe, "set title 'Position vs Time'\n");
        fprintf(gnuplotPipe, "set xlabel 'Time'\n");
        fprintf(gnuplotPipe, "set ylabel 'Position'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' using 1:2 with lines title 'Position'\n");
        
        fprintf(gnuplotPipe, "set title 'Velocity vs Time'\n");
        fprintf(gnuplotPipe, "set xlabel 'Time'\n");
        fprintf(gnuplotPipe, "set ylabel 'Velocity'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' using 1:3 with lines title 'Velocity'\n");
        
        fprintf(gnuplotPipe, "set title 'Acceleration vs Time'\n");
        fprintf(gnuplotPipe, "set xlabel 'Time'\n");
        fprintf(gnuplotPipe, "set ylabel 'Acceleration'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' using 1:4 with lines title 'Acceleration'\n");
        
        fprintf(gnuplotPipe, "unset multiplot\n");
        
        // Phase portraits
        fprintf(gnuplotPipe, "set terminal pngcairo size 800,600\n");
        fprintf(gnuplotPipe, "set output 'phase_portraits.png'\n");
        fprintf(gnuplotPipe, "set multiplot layout 3,1 title 'Phase Portraits' font ',14'\n");
        
        fprintf(gnuplotPipe, "set title 'Position vs Velocity'\n");
        fprintf(gnuplotPipe, "set xlabel 'Position'\n");
        fprintf(gnuplotPipe, "set ylabel 'Velocity'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' using 2:3 with lines title 'x vs x'\n");
        
        fprintf(gnuplotPipe, "set title 'Position vs Acceleration'\n");
        fprintf(gnuplotPipe, "set xlabel 'Position'\n");
        fprintf(gnuplotPipe, "set ylabel 'Acceleration'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' using 2:4 with lines title 'x vs x'''\n");
        
        fprintf(gnuplotPipe, "set title 'Velocity vs Acceleration'\n");
        fprintf(gnuplotPipe, "set xlabel 'Velocity'\n");
        fprintf(gnuplotPipe, "set ylabel 'Acceleration'\n");
        fprintf(gnuplotPipe, "plot 'data.txt' using 3:4 with lines title 'x' vs x'''\n");
        
        fprintf(gnuplotPipe, "unset multiplot\n");
        pclose(gnuplotPipe);
    } else {
        printf("Error opening pipe to GNU plot.\n");
    }
}

int main() {
    double a_values[] = {1.0, 0.5, 1.5};
    double b_values[] = {0.5, 1.0, 0.0};
    double c_values[] = {0.5, 0.5, 1.0};
    double A_values[] = {1.0, 0.5, 1.5};
    double C = 0.5;
    
    double t0 = 0.0, t_end = 20.0, h = 0.01; // time variables
    int steps = (int)((t_end - t0) / h);
    
    double y[N] = {1.0, 0.0, 0.0}; // Initial conditions: {x(0), x'(0), x''(0)}

    FILE *fp = fopen("data.txt", "w");
    if (!fp) {
        printf("Error opening file!\n");
        return -1;
    }

    for (int param_set = 0; param_set < 3; param_set++) {
        double a = a_values[param_set];
        double b = b_values[param_set];
        double c = c_values[param_set];
        double A = A_values[param_set];

        for (int i = 0; i <= steps; i++) {
            double t = t0 + i * h;
            fprintf(fp, "%lf %lf %lf %lf\n", t, y[0], y[1], y[2]);
            rk4(t, y, h, a, b, c, A, C);
        }
        fprintf(fp, "\n\n"); // Separate data sets with blank lines
    }

    fclose(fp);
    printf("Data written to data.txt\n");

    plot_results();

    return 0;
}
