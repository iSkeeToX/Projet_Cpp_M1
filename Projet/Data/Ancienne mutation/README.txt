int nx = 30, ny = 30;
int Nparts = (nx*ny)/9;

int N_Temps = 100;
int N_Steps = 10 * Nparts;

int N_Stat = 10 * Nparts;

int pop = 100, NbrGenes = 21;

double mean = 0, stdev = 20;

double acceptation = 1;

Darwin D = Darwin(pop, NbrGenes, mean, stdev, nx, ny);

for(int i=0; i<100; i++){
    D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "Fitness_Sortie_Du_Cul.txt", acceptation);
}


Fitness_function : -10*pow(Size - 6, 6) - 10*pow(SizeHoles - 1, 6) - pow(Vol - 7, 4) - 10*pow(porosity - 1./7, 2) - 10*pow(std::abs(Sphericity - sqrt(M_PI*sqrt(3))), 3);
