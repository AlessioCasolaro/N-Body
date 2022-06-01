#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>


#define SOFTENING 1e-9f

typedef struct {
  float x, y, z, vx, vy, vz;//Struttura body: Posizione x,y,z | Velocità x,y,z
} Body;

void randomizeBodies(Body *body, int n_bodies,int seed);
void saveBodies(Body *body,int n_bodies,int n_iter,int seed, char name[80]);

void body_force(Body *body, float dt, int n_bodies, int displs, int sendcounts);
void updateBodies(Body *body, float dt, int displs, int sendcounts);

int main(int argc, char** argv){

  	const float dt = 0.01f;

	MPI_Status status;

	int world_rank;					//Rank del processo 
	int world_size;					//Numero di processi 
  	
	int *sendcounts;				//Numero elementi per processo
	int *displs;					//Displacement per la posizione
	
	int n_iter;						//Numero di iterazioni
	int n_bodies;					//Numero di bodies
	int seed;						//Seme per la srand

	double start_time, end_time;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 					
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
	MPI_Barrier(MPI_COMM_WORLD); /* tutti i processi sono inizializzati */
	start_time = MPI_Wtime();

	//Tipo dati della struttura
	MPI_Datatype body_type, oldtype[1]; //body_type sara' il nuovo tipo 
	oldtype[0] = MPI_FLOAT;
	int blockcounts[1];
	blockcounts[0] = 6;//per i 6 campi della struttura Body
	MPI_Aint offset[1];
	offset[0] = 0;

	MPI_Type_create_struct(1, blockcounts, offset, oldtype, &body_type);
  	MPI_Type_commit(&body_type);

	//Input
	if (argc == 4 ){
			n_bodies = atoi(argv[1]);//Numero di Bodies
			n_iter = atoi(argv[2]);//Numero di iterazioni
			seed = atoi(argv[3]);//Seme randomico per la srand
	}
	else{
		if(world_rank==0)
			printf("Too few arguments.\nNeeded Number of bodies, Number of iteraction and Seed for Srand.\n");
		MPI_Finalize();
        return 0;
		
	}

	int bytes = n_bodies * sizeof(Body);//Grandezza della struttura
	float* buffer = (float*)malloc(bytes); //Buffer per la struttura
  	Body *body = (Body*)buffer;
	

	randomizeBodies(body,n_bodies,seed);//Randomizza i bodies 
	saveBodies(body,n_bodies,n_iter,seed,"input.txt");//Salva i bodies nel file input.txt

  	sendcounts = malloc(sizeof(int)*world_size);
  	displs = malloc(sizeof(int)*world_size);

	//Calcolo displs | sendcounts
  	int rest = n_bodies % world_size;
	int lenght = n_bodies / world_size;	
	int index = 0;
	for(int i=0; i<world_size; i++){
		sendcounts[i] = lenght;
		if(rest > 0){
			rest--;
			sendcounts[i] = sendcounts[i] + 1;
		}
		displs[i] = index;
		index += sendcounts[i];
	}

	for(int i = 0; i < n_iter; i++){

		//Aggiorno la velocita' dei bodies
	    body_force(body, dt, n_bodies, displs[world_rank], sendcounts[world_rank]);

		//Aggiorno la posizione dei bodies usando le velocita' aggiornate
		updateBodies(body, dt, displs[world_rank], sendcounts[world_rank]);
		
		MPI_Allgatherv(MPI_IN_PLACE, sendcounts[world_rank], body_type, body, sendcounts, displs, body_type, MPI_COMM_WORLD);		
	}

	MPI_Barrier(MPI_COMM_WORLD);
	end_time = MPI_Wtime();

	if(world_rank == 0){
		printf("Time spent %.5f\n", end_time-start_time);
		saveBodies(body,n_bodies,n_iter,seed,"output.txt");//Salva i bodies aggiornati nel file output.txt
  	}

	MPI_Type_free(&body_type);
  	free(sendcounts);
  	free(displs);
  	free(buffer);

	MPI_Finalize();

	return EXIT_SUCCESS;
}


void randomizeBodies(Body *body,int n_bodies,int seed){
	srand(seed);//Seme per la generazione random
	for(int i = 0; i < n_bodies; i++){
		body[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		body[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		body[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		body[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		body[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		body[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
	}
}


void body_force(Body *body, float dt, int n_bodies, int displs, int sendcounts){
    for(int i = displs; i < displs + sendcounts; i++){//Inizio e fine partizione
      float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

      for(int j = 0; j < n_bodies; j++){
        float dx = body[j].x - body[i].x;
        float dy = body[j].y - body[i].y;
        float dz = body[j].z - body[i].z;
        float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * invDist3; 
		Fy += dy * invDist3; 
		Fz += dz * invDist3;
      }
      body[i].vx += dt*Fx;
      body[i].vy += dt*Fy;
      body[i].vz += dt*Fz;
    }
}

void updateBodies(Body *body, float dt, int displs, int sendcounts){
	for(int i = displs; i < displs + sendcounts; i++){
		//Aggiorno le posizioni
		body[i].x += body[i].vx*dt;
		body[i].y += body[i].vy*dt;
		body[i].z += body[i].vz*dt;
	}
}
void saveBodies(Body *body,int n_bodies,int n_iter,int seed,char name[80]){
	FILE *fp;
	fp = fopen(name, "w");
	fprintf(fp, "Numero bodies: %d\tNumero iterazioni: %d\tSeed: %d\n",n_bodies,n_iter,seed);
	for (int i = 0; i < n_bodies; i++)
		fprintf(fp, "%f %f %f %f %f %f\n", body[i].x, body[i].y, body[i].z, body[i].vx, body[i].vy, body[i].vz);
	fclose(fp);
	fflush(fp);
}

