#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>




#define faux 0
#define ndims 2
#define NB_VOISINS 4
#define N 0
#define E 1
#define S 2
#define W 3
#define min(a,b) (a<=b?a:b)

/* rang dans le communicateur initial */
int rang;
/* nombre de processus */
int nb_procs;

/* nombre de lignes*/
/*METTRE LE NOMBRE DE LIGNES DE L'IMAGE*/
   int Nlig=512; /*direction des x*/
/* nombre de colonnes*/

/*METTRE LE NOMBRE DE LIGNES DE L'IMAGE*/
   int Mcol=512; /*direction des y*/

 /*********Nombre d'étapes dans la décomposition multiéchelle*******/

   int decomposition_step_number=9;

  /*********Pas spatial dans la discrétisation différences finies********/
  /*********En image, la distance entre 2 pixels est supposée égale à 1**/
  /*********Le pas spatial est donc égal à 1*****************************/
   double  h=1.0;

   /********Paramètre epsilon pour éviter les divisions par 0************/
   double epsilon=0.000001;


   /**********nombre d'itérations dans la méthode de point fixe*/
   int Itermax=200;

void initialisation_mpi(int, char**);
void finalisation_mpi();
void domaine(MPI_Comm ,int ,int *,int,int ,int * ,int *,int * ,int ,int );
void voisinage(MPI_Comm,int *,int *, int *);
double **allocarray(int,int);
void printarr(double **, int,int, char *);
void ecrire_mpi(double *,int,int,int *,MPI_Comm, char*);
int malloc3ddouble(double ****, int , int , int );

int main(int argc, char *argv[]) {

  /* coordonnées dans la grille */
  int coords[ndims];
  /* tableau des dimensions dans la grille */
  int dims[ndims];
  /* communicateur topologie cartésienne */
  MPI_Comm comm2d;
  int periods[ndims];
  const int reorganisation=faux;
  /* tableau contenant les voisins du sous-domaine courant (haut,bas,gauche,droite)*/
  int voisin[NB_VOISINS];

  /*nombre total de points intérieurs dans la direction x et la direction y*/
  int ntx, nty;
  /* ntx --> direction des lignes -->nombre de lignes total-2 (première ligne et dernière ligne puisqu'on ne fait
  pas de calcul sur ces lignes du fait des conditions au bord de type Neumann homogènes)*/

  /* nty --> direction des colonnes -->nombre de colonnes total-2 (première colonne et dernière colonne puisqu'on ne fait
  pas de calcul sur ces colonnes du fait des conditions au bord de type Neumann homogènes)*/


  double t1, t2;

  double lambda0=0.0005;




  /*
        \
 ------- y                      coords[1]/dims[1]
 |      /                     /|\
 |                             |
 |                             |
\ /                            |
 x                             |       \
                               --------- coords[0]/dims[0]
                                       /
*/



  /* Initialisation de MPI */

  initialisation_mpi(argc,argv);

  /* Creation de la topologie cartesienne 2D */
  /*Le nombre de points dans la direction x correspond au nombre de lignes-2 (uniquement les points intérieurs) */
  /*Le nombre de points dans la direction y correspond au nombre de colonnes-2 (uniquement les points intérieurs) */
  ntx=Nlig-2;
  nty=Mcol-2;

  /* Connaître le nombre de processus selon x et le nombre de processus
     selon y en fonction du nombre total de processus */
  if(argc >= 2) {
	   dims[0] = (int)  *argv[1];
  } else {
    dims[0] = 0;
  }
  dims[1] = 0;
  MPI_Dims_create(nb_procs,ndims,dims);

  /* Creation de la grille de processus 2D sans periodicite */

  periods[0] = periods[1] = faux;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorganisation, &comm2d);
  MPI_Comm_rank(comm2d,&rang);

  if(rang == 0) {
    printf("Execution code crack_detection avec %d processus MPI\n"
	   "Taille du domaine : ntx=%d nty=%d\n"
	   "Dimension de la topologie : %d suivant y (colonnes), %d suivant x (lignes)\n"
	   "-----------------------------------------\n",
	   nb_procs, ntx, nty, dims[0], dims[1]);
  }

  int tab_bounds[4]; /*sx,ex,sy,ey*/
  /* Determinination des indices de chaque sous domaine */

  domaine(comm2d,rang,coords,ntx,nty,dims,tab_bounds,periods,reorganisation,nb_procs);
  printf("Je suis le rang %d\n"
	   "Ma coord 0 dans la direction des colonnes est: coord0=%d\n"
	   "Ma coord 1 dans la direction opposee a la direction des lignes est: coord1=%d\n"
	   "-----------------------------------------\n",
	   rang, coords[0], coords[1]);
  //if(rang==3){
  printf("Je suis le rang %d \n"
  "coordonnee en ligne du coin superieur gauche sx=%d\n"
  "coordonnee en colonne du coin superieur gauche sy=%d\n"
  "coordonnee en ligne du coin inferieur droit ex=%d\n"
  "coordonnee en colonne du coin inferieur droit ey=%d\n"
  "-----------------------------------------\n",
	   rang,tab_bounds[0],tab_bounds[2],tab_bounds[1],tab_bounds[3]);
//  }





  /* Recherche de ses 4 voisins pour chaque processus */
  voisinage(comm2d,voisin,coords,dims);


  int code;
  MPI_File descripteur;
  MPI_Offset deplacement_initial;
  MPI_Status statut;
  /* Ouverture du fichier "image_initiale_f.bin" en lecture */
  code = MPI_File_open(comm2d, "image_barbara.bin", MPI_MODE_RDONLY,MPI_INFO_NULL, &descripteur);
  /* Test pour savoir si ouverture du fichier est correcte */
  if (code != MPI_SUCCESS) {
    fprintf(stderr, "ATTENTION erreur lors ouverture du fichier");
    MPI_Abort(comm2d, 2);
  }

  MPI_Datatype mysubarray;/****type sous-matrice****/


  double * f_local=malloc((tab_bounds[1]-tab_bounds[0]+3)*(tab_bounds[3]-tab_bounds[2]+3)*sizeof(double));
  double ** f_local_mat=allocarray((tab_bounds[1]-tab_bounds[0]+3),(tab_bounds[3]-tab_bounds[2]+3));

  if(f_local==NULL)
    printf("Erreur dans l'allocation mémoire de f_local-- \n");

  int starts[2] = {tab_bounds[0]-1,tab_bounds[2]-1};
  int subsizes[2]  = {tab_bounds[1]-tab_bounds[0]+3,tab_bounds[3]-tab_bounds[2]+3};
  int bigsizes[2]  = {Nlig, Mcol};


  MPI_Type_create_subarray(2,bigsizes, subsizes, starts,
                                 MPI_ORDER_C, MPI_DOUBLE, &mysubarray);
  MPI_Type_commit(&mysubarray);


  deplacement_initial=0;

  MPI_File_set_view(descripteur,deplacement_initial,MPI_DOUBLE,mysubarray,"native",MPI_INFO_NULL);
  MPI_File_read(descripteur,f_local,(tab_bounds[1]-tab_bounds[0]+3)*(tab_bounds[3]-tab_bounds[2]+3),MPI_DOUBLE,&statut);

  MPI_File_close(&descripteur);


  for(int i=0;i<(tab_bounds[1]-tab_bounds[0]+3);i++){

     f_local_mat[i]=&(f_local[i*(tab_bounds[3]-tab_bounds[2]+3)]);

  }

  // test
  ecrire_mpi(f_local, ntx, nty, tab_bounds, comm2d, "image_initiale.dat");

   double c01,c02,c03,c04;
   double c0;
   double t;
   double ux,uy,Gradu;
   /* Mesure du temps en seconde dans la boucle en temps */
   t1 = MPI_Wtime();
   int size_x = (tab_bounds[1]-tab_bounds[0]+3);
   int size_y = (tab_bounds[3]-tab_bounds[2]+3);
   double ***u_current;
   malloc3ddouble(&u_current,size_x,size_y,decomposition_step_number);
   double ***u_new;
   malloc3ddouble(&u_new,size_x,size_y,decomposition_step_number);

   MPI_Datatype type_column;
   MPI_Type_vector(size_x , 1 , size_y*decomposition_step_number ,MPI_DOUBLE, &type_column);
   MPI_Type_commit(&type_column);
   MPI_Datatype type_row;
   MPI_Type_vector(size_y , 1 , decomposition_step_number ,MPI_DOUBLE, &type_row);
   MPI_Type_commit(&type_row);

   double *right_col;
   double *left_col;
   double *up_row;
   double *down_row;

   if (voisin[N] != MPI_PROC_NULL) {
     up_row = malloc(size_y*sizeof(double));
   }
   if (voisin[S] != MPI_PROC_NULL) {
     down_row = malloc(size_y*sizeof(double));
   }
   if (voisin[W] != MPI_PROC_NULL) {
     left_col = malloc(size_x*sizeof(double));
   }
   if (voisin[E] != MPI_PROC_NULL) {
     right_col = malloc(size_x*sizeof(double));
   }


   /*BOUCLE PRINCIPALE*/
  for(int l=0;l<decomposition_step_number;l++){
      /*Initialisation de u_l*/
      for(int i=0; i< size_x; i++) {
         for(int j=0; j< size_y; j++){
           u_current[i][j][l] = f_local_mat[i][j];
         }
      }

      /*algo*/
      for(int k=0;k<Itermax;k++){

          for(int i=1;i<size_x-1; i++){
               for(int j=1;j<size_y-1; j++){

               ux=(u_current[i+1][j][l]-u_current[i][j][l]);
               uy=(u_current[i][j+1][l]-u_current[i][j-1][l])/2.0;
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c01=1.0/Gradu;


               ux=(u_current[i][j][l]-u_current[i-1][j][l]);
               uy=(u_current[i-1][j+1][l]-u_current[i-1][j-1][l])/2.0;
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c02=1.0/Gradu;



               ux=(u_current[i+1][j][l]-u_current[i-1][j][l])/2.0;
               uy=(u_current[i][j+1][l]-u_current[i][j][l]);
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c03=1.0/Gradu;



               ux=(u_current[i+1][j-1][l]-u_current[i-1][j-1][l])/2.0;
               uy=(u_current[i][j][l]-u_current[i][j-1][l]);
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c04=1.0/Gradu;



               c0=2.0*h*h*lambda0+c01+c02+c03+c04;

               t=2.0*h*h*lambda0*f_local_mat[i][j]+c01*u_current[i+1][j][l]+c02*u_current[i-1][j][l]+c03*u_current[i][j+1][l]+c04*u_current[i][j-1][l];

               u_new[i][j][l]=(1.0/c0)*t;

             }



          }

        MPI_Request req_recv[NB_VOISINS] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        MPI_Request req_send[NB_VOISINS];
        MPI_Request req_recv_corners[3] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
        MPI_Request req_send_corners[3];
        double up_left, up_right, down_left;

        /*Communications pour les bords*/

        if (voisin[N] != MPI_PROC_NULL) {
          // MPI_Isend(&(u_current[0][0][l]), size_y, MPI_DOUBLE, voisin[N], N, comm2d, &(req_send[N]));
          MPI_Isend(&(u_current[0][0][l]), 1, type_row, voisin[N], N, comm2d, &(req_send[N]));
          MPI_Irecv(up_row, size_y, MPI_DOUBLE, voisin[N], MPI_ANY_TAG, comm2d, &(req_recv[N]));
        }
        if (voisin[S] != MPI_PROC_NULL) {
          // MPI_Isend(&(u_current[size_x-1][0][l]), size_y, MPI_DOUBLE, voisin[S], S, comm2d, &(req_send[S]));
          MPI_Isend(&(u_current[size_x-1][0][l]), 1, type_row, voisin[S], S, comm2d, &(req_send[S]));
          MPI_Irecv(down_row, size_y, MPI_DOUBLE, voisin[S], MPI_ANY_TAG, comm2d, &(req_recv[S]));
        }
        if (voisin[W] != MPI_PROC_NULL) {
          MPI_Isend(&(u_current[0][0][l]), 1, type_column, voisin[W], W, comm2d, &(req_send[W]));
          MPI_Irecv(left_col, size_x, MPI_DOUBLE, voisin[W], MPI_ANY_TAG, comm2d, &(req_recv[W]));
        }
        if (voisin[E] != MPI_PROC_NULL) {
          MPI_Isend(&(u_current[0][size_y-1][l]), 1, type_column, voisin[E], E, comm2d, &(req_send[E]));
          MPI_Irecv(right_col, size_x, MPI_DOUBLE, voisin[E], MPI_ANY_TAG, comm2d, &(req_recv[E]));
        }


        // attente des bords
        MPI_Waitall(NB_VOISINS, req_recv, MPI_STATUSES_IGNORE);

        /*Communications pour les coins*/
        if (voisin[N]!= MPI_PROC_NULL) {
          if (voisin[W]!= MPI_PROC_NULL) {
              MPI_Irecv(&up_left, 1, MPI_DOUBLE, voisin[N], MPI_ANY_TAG, comm2d, &(req_recv_corners[0]));
              //MPI_Wait(&(req_recv[W]),MPI_STATUS_IGNORE);
              MPI_Isend(&(left_col[0]), 1, MPI_DOUBLE, voisin[N], N+W, comm2d, &(req_send_corners[2]));
          }
          if (voisin[E]!= MPI_PROC_NULL) {
              MPI_Irecv(&up_right, 1, MPI_DOUBLE, voisin[N], MPI_ANY_TAG, comm2d, &(req_recv_corners[1]));
          }
        }

        if (voisin[S]!= MPI_PROC_NULL) {
          if (voisin[W]!= MPI_PROC_NULL) {
              //MPI_Wait(&(req_recv[W]),MPI_STATUS_IGNORE);
              MPI_Isend(&(left_col[size_x-1]), 1, MPI_DOUBLE, voisin[S], S+W, comm2d, &(req_send_corners[0]));
              MPI_Irecv(&down_left, 1, MPI_DOUBLE, voisin[S], MPI_ANY_TAG, comm2d, &(req_recv_corners[2]));

          }
          if (voisin[E]!= MPI_PROC_NULL) {
              //MPI_Wait(&(req_recv[E]),MPI_STATUS_IGNORE);
              MPI_Isend(&(right_col[size_x-1]), 1, MPI_DOUBLE, voisin[S], S+E, comm2d, &(req_send_corners[1]));
          }
        }

        MPI_Waitall(3, req_recv_corners, MPI_STATUSES_IGNORE);

        /*GESTION DES BORDS*/
        if (voisin[N] == MPI_PROC_NULL) { // si on est sur le bord haut
          for(int j=1; j<size_y-1; j++){
            u_new[0][j][l]=u_new[1][j][l];
          }
          if (voisin[E] == MPI_PROC_NULL ) { // si on est sur le coin haut-droit
            u_new[0][size_y-1][l]=u_new[1][size_y-2][l];
          }
          if (voisin[W] == MPI_PROC_NULL ) { // si on est sur le coin haut-gauche
            u_new[0][0][l]=u_new[1][1][l];
          }
        } else {
              /*gestion du bord haut (cas général)*/
              for(int j=1;j<size_y-1; j++){

                ux=(u_current[1][j][l]-u_current[0][j][l]);
                uy=(u_current[0][j+1][l]-u_current[0][j-1][l])/2.0;
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c01=1.0/Gradu;


                ux=(u_current[0][j][l]-up_row[j]);
                uy=(up_row[j+1]-up_row[j-1])/2.0;
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c02=1.0/Gradu;



                ux=(u_current[1][j][l]-up_row[j])/2.0;
                uy=(u_current[0][j+1][l]-u_current[0][j][l]);
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c03=1.0/Gradu;



                ux=(u_current[1][j-1][l]-up_row[j-1])/2.0;
                uy=(u_current[0][j][l]-u_current[0][j-1][l]);
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c04=1.0/Gradu;



                c0=2.0*h*h*lambda0+c01+c02+c03+c04;

                t=2.0*h*h*lambda0*f_local_mat[0][j]+c01*u_current[1][j][l]+c02*up_row[j]+c03*u_current[0][j+1][l]+c04*u_current[0][j-1][l];

                u_new[0][j][l]=(1.0/c0)*t;

            }
            /*coin haut-gauche*/
            if (voisin[W] != MPI_PROC_NULL) {
                ux=(u_current[1][0][l]-u_current[0][0][l]);
                uy=(u_current[0][1][l]-left_col[0])/2.0;
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c01=1.0/Gradu;


                ux=(u_current[0][0][l]-up_row[0]);
                uy=(up_row[1]-up_left)/2.0;
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c02=1.0/Gradu;



                ux=(u_current[1][0][l]-up_row[0])/2.0;
                uy=(u_current[0][1][l]-u_current[0][0][l]);
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c03=1.0/Gradu;



                ux=(left_col[1]-up_left)/2.0;
                uy=(u_current[0][0][l]-left_col[0]);
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c04=1.0/Gradu;



                c0=2.0*h*h*lambda0+c01+c02+c03+c04;

                t=2.0*h*h*lambda0*f_local_mat[0][0]+c01*u_current[1][0][l]+c02*up_row[0]+c03*u_current[0][1][l]+c04*left_col[0];

                u_new[0][0][l]=(1.0/c0)*t;
            }
            /*cas haut-droit*/
            if (voisin[E] != MPI_PROC_NULL) {
                  ux=(u_current[1][size_y-1][l]-u_current[0][size_y-1][l]);
                  uy=(right_col[0]-u_current[0][size_y-2][l])/2.0;
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c01=1.0/Gradu;


                  ux=(u_current[0][size_y-1][l]-up_row[size_y-1]);
                  uy=(up_right-up_row[size_y-2])/2.0;
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c02=1.0/Gradu;



                  ux=(u_current[1][size_y-1][l]-up_row[size_y-1])/2.0;
                  uy=(right_col[0]-u_current[0][size_y-1][l]);
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c03=1.0/Gradu;



                  ux=(u_current[1][size_y-2][l]-up_row[size_y-2])/2.0;
                  uy=(u_current[0][size_y-1][l]-u_current[0][size_y-2][l]);
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c04=1.0/Gradu;



                  c0=2.0*h*h*lambda0+c01+c02+c03+c04;

                  t=2.0*h*h*lambda0*f_local_mat[0][size_y-1]+c01*u_current[1][size_y-1][l]+c02*up_row[size_y-1]+c03*right_col[0]+c04*u_current[0][size_y-2][l];

                  u_new[0][size_y-1][l]=(1.0/c0)*t;
            }


        }

        if (voisin[S] == MPI_PROC_NULL ) { // si on est sur le bord bas
          for(int j=1; j<size_y-1; j++){
            u_new[size_x-1][j][l]=u_new[size_x-2][j][l];
          }
          if (voisin[E] == MPI_PROC_NULL ) { // si on est sur le coin bas-droit
            u_new[size_x-1][size_y-1][l]=u_new[size_x-2][size_y-2][l];
          }
          if (voisin[W] == MPI_PROC_NULL ) { // si on est sur le coin bas-gauche
            u_new[size_x-1][0][l]=u_new[size_x-2][1][l];
          }
        } else {
            /*Gestion du bord sud*/
            for(int j=1;j<size_y-1; j++){

                ux=(down_row[j]-u_current[size_x-1][j][l]);
                uy=(u_current[size_x-1][j+1][l]-u_current[size_x-1][j-1][l])/2.0;
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c01=1.0/Gradu;


                ux=(u_current[size_x-1][j][l]-u_current[size_x-2][j][l]);
                uy=(u_current[size_x-2][j+1][l]-u_current[size_x-2][j-1][l])/2.0;
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c02=1.0/Gradu;



                ux=(down_row[j]-u_current[size_x-2][j][l])/2.0;
                uy=(u_current[size_x-1][j+1][l]-u_current[size_x-1][j][l]);
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c03=1.0/Gradu;



                ux=(down_row[j-1]-u_current[size_x-2][j-1][l])/2.0;
                uy=(u_current[size_x-1][j][l]-u_current[size_x-1][j-1][l]);
                Gradu=sqrt(epsilon+ux*ux+uy*uy);
                c04=1.0/Gradu;



                c0=2.0*h*h*lambda0+c01+c02+c03+c04;

                t=2.0*h*h*lambda0*f_local_mat[size_x-1][j]+c01*down_row[j]+c02*u_current[size_x-2][j][l]+c03*u_current[size_x-1][j+1][l]+c04*u_current[size_x-1][j-1][l];

                u_new[size_x-1][j][l]=(1.0/c0)*t;

              }
            /*coin bas-gauche*/
            if (voisin[W] != MPI_PROC_NULL) {
              ux=(down_row[0]-u_current[size_x-1][0][l]);
              uy=(u_current[size_x-1][1][l]-left_col[size_x-1])/2.0;
              Gradu=sqrt(epsilon+ux*ux+uy*uy);
              c01=1.0/Gradu;


              ux=(u_current[size_x-1][0][l]-u_current[size_x-2][0][l]);
              uy=(u_current[size_x-2][1][l]-left_col[size_x-2])/2.0;
              Gradu=sqrt(epsilon+ux*ux+uy*uy);
              c02=1.0/Gradu;



              ux=(down_row[0]-u_current[size_x-2][0][l])/2.0;
              uy=(u_current[size_x-1][1][l]-u_current[size_x-1][0][l]);
              Gradu=sqrt(epsilon+ux*ux+uy*uy);
              c03=1.0/Gradu;



              ux=(down_left-left_col[size_x-2])/2.0;
              uy=(u_current[size_x-1][0][l]-left_col[size_x-1]);
              Gradu=sqrt(epsilon+ux*ux+uy*uy);
              c04=1.0/Gradu;



              c0=2.0*h*h*lambda0+c01+c02+c03+c04;

              t=2.0*h*h*lambda0*f_local_mat[size_x-1][0]+c01*down_row[0]+c02*u_current[size_x-2][0][l]+c03*u_current[size_x-1][1][l]+c04*left_col[size_x-1];

              u_new[size_x-1][0][l]=(1.0/c0)*t;
            }
            /*coin bas-droit*/
            if (voisin[E] != MPI_PROC_NULL) {
                  ux=(down_row[size_y - 1]-u_current[size_x-1][size_y - 1][l]);
                  uy=(right_col[size_x-1] - u_current[size_x-1][size_y - 2][l])/2.0;
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c01=1.0/Gradu;


                  ux=(u_current[size_x-1][size_y - 1][l]-u_current[size_x-2][size_y - 1][l]);
                  uy=(right_col[size_x-2]-u_current[size_x-2][size_y - 2][l])/2.0;
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c02=1.0/Gradu;



                  ux=(down_row[size_y - 1]-u_current[size_x-2][size_y - 1][l])/2.0;
                  uy=(right_col[size_x-1]-u_current[size_x-1][size_y - 1][l]);
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c03=1.0/Gradu;



                  ux=(down_row[size_y - 2]-u_current[size_x-2][size_y - 2][l])/2.0;
                  uy=(u_current[size_x-1][size_y - 1][l]-u_current[size_x-1][size_y - 2][l]);
                  Gradu=sqrt(epsilon+ux*ux+uy*uy);
                  c04=1.0/Gradu;

                  c0=2.0*h*h*lambda0+c01+c02+c03+c04;

                  t=2.0*h*h*lambda0*f_local_mat[size_x-1][size_y - 1]+c01*down_row[size_y - 1]+c02*u_current[size_x-2][size_y - 1][l]+c03*right_col[size_x-1]+c04*u_current[size_x-1][size_y - 2][l];

                  u_new[size_x-1][size_y - 1][l]=(1.0/c0)*t;
            }
        }

        if (voisin[E] == MPI_PROC_NULL ) { // si on est sur le bord droit
          for(int i=1;i<size_x-1; i++){
            u_new[i][size_y-1][l]=u_new[i][size_y-2][l];
          }
      } else {
          for(int i=1;i<size_x-1; i++){

               ux=(u_current[i+1][size_y-1][l]-u_current[i][size_y-1][l]);
               uy=(right_col[i]-u_current[i][size_y-2][l])/2.0;
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c01=1.0/Gradu;


               ux=(u_current[i][size_y-1][l]-u_current[i-1][size_y-1][l]);
               uy=(right_col[i-1]-u_current[i-1][size_y-2][l])/2.0;
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c02=1.0/Gradu;



               ux=(u_current[i+1][size_y-1][l]-u_current[i-1][size_y-1][l])/2.0;
               uy=(right_col[i]-u_current[i][size_y-1][l]);
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c03=1.0/Gradu;



               ux=(u_current[i+1][size_y-2][l]-u_current[i-1][size_y-2][l])/2.0;
               uy=(u_current[i][size_y-1][l]-u_current[i][size_y-2][l]);
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c04=1.0/Gradu;



               c0=2.0*h*h*lambda0+c01+c02+c03+c04;

               t=2.0*h*h*lambda0*f_local_mat[i][size_y-1]+c01*u_current[i+1][size_y-1][l]+c02*u_current[i-1][size_y-1][l]+c03*right_col[i]+c04*u_current[i][size_y-2][l];

               u_new[i][size_y-1][l]=(1.0/c0)*t;

          }
      }

        if (voisin[W] == MPI_PROC_NULL ) { // si on est sur le bord gauche
          for(int i=1;i<size_x-1; i++){
            u_new[i][0][l]=u_new[i][1][l];
          }
      } else {
          for(int i=1;i<size_x-1; i++){

               ux=(u_current[i+1][0][l]-u_current[i][0][l]);
               uy=(u_current[i][1][l]-left_col[i])/2.0;
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c01=1.0/Gradu;


               ux=(u_current[i][0][l]-u_current[i-1][0][l]);
               uy=(u_current[i-1][1][l]-left_col[i-1])/2.0;
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c02=1.0/Gradu;



               ux=(u_current[i+1][0][l]-u_current[i-1][0][l])/2.0;
               uy=(u_current[i][1][l]-u_current[i][0][l]);
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c03=1.0/Gradu;



               ux=(left_col[i+1]-left_col[i-1])/2.0;
               uy=(u_current[i][0][l]-left_col[i]);
               Gradu=sqrt(epsilon+ux*ux+uy*uy);
               c04=1.0/Gradu;



               c0=2.0*h*h*lambda0+c01+c02+c03+c04;

               t=2.0*h*h*lambda0*f_local_mat[i][0]+c01*u_current[i+1][0][l]+c02*u_current[i-1][0][l]+c03*u_current[i][1][l]+c04*left_col[i];

               u_new[i][0][l]=(1.0/c0)*t;

          }
      }

        for(int i=0;i<size_x; i++){
           for(int j=0;j<size_y; j++){
               u_current[i][j][l]=u_new[i][j][l];
           }
        }


      }/*Fin de la boucle en k */

     /*********Actualisation du param�tre lambda0 et de la fonction f*********/

     lambda0=2.0*lambda0;
     for(int i=0;i<size_x; i++){
        for(int j=0;j<size_y; j++){
            f_local_mat[i][j]=f_local_mat[i][j]-u_new[i][j][l];
        }
     }


  }

  /* Mesure du temps a la sortie de la boucle */
   t2 = MPI_Wtime();



  /* Ecriture des resultats pour chaque processus */
  /* FONCTION ecriture_mpi */
  double* vect = malloc(sizeof(double)*size_x*size_y);
  for(int l=0;l<decomposition_step_number;l++){
    for(int i=0;i<size_x; i++){
       for(int j=0;j<size_y; j++){
          vect[j+i*size_y] = u_new[i][j][l];
       }
     }
     char filename[14];
     sprintf(filename, "final_u_%d.dat", l);
     ecrire_mpi(vect, ntx, nty, tab_bounds, comm2d, filename);
  }


  /* Affichage du temps de convergence par le processus 3 */
  if (rang == 0) {
    printf("Convergence en %f secs\n", t2-t1);
  }
  /****Libération mémoire****/
  MPI_Type_free(&type_column);
  MPI_Type_free(&type_row);
  finalisation_mpi();

  return 0;
}



  /**************************************************************************************************/
  /*********************************INITIALISATION***************************************************/
  /**************************************************************************************************/
  void initialisation_mpi(int argc, char* argv[]) {
  /* Initialisation de MPI */
  MPI_Init(&argc, &argv);

  /* Savoir quel processus je suis */
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);

  /* Connaitre le nombre total de processus */
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
 }

  /**************************************************************************************************/
  /*********************************FINALISATION***************************************************/
  /**************************************************************************************************/

  void finalisation_mpi() {
  /* Desactivation de MPI */
  MPI_Finalize();
  }


  /**************************************************************************************************/
  /*********************************CREATION DOMAINE*************************************************/
  /**************************************************************************************************/

  void domaine(MPI_Comm comm2d,int rang,int * coords,int ntx,int nty,int * dims,int * tab_bounds,int * periods,int reorganisation,int nb_procs) {



  int nx,ny,positionx,positiony,restex,restey;
  int sx,ex,sy,ey;

  /* Connaître mes coordonnees dans la topologie */
  MPI_Cart_coords(comm2d,rang,ndims,coords);


  /* Calcul pour chaque processus de ses indices de debut et de fin suivant x */

  /*Nombre de points dans la direction x*/

  nx=ntx/dims[1];
  restex=ntx % dims[1];
  positionx=coords[1];
  sx=1+positionx*nx+min(restex,positionx);/*Indice de depart dans la direction des x*/
  if(positionx<restex){
  nx=nx+1;
  }


  ex=sx+nx-1; /*Indice de fin dans la direction des x*/
  ex=min(ex,ntx+1);/*correction si besoin pour le dernier bloc*/


  /**********************************************************/
  /*Pour renumeroter selon la direction opposee a dims[1]   */
  /**********************************************************/



  sx=ntx-ex+1;
  ex=sx+nx-1;

  /*Nombre de points dans la direction y*/

  ny=nty/dims[0];
  restey=nty % dims[0];
  positiony=coords[0];
  sy=1+positiony*ny+min(restey,positiony);/*Indice de depart dans la direction des y*/
  if(positiony<restey){
  ny=ny+1;
  }

  ey=sy+ny-1;/*Indice de fin*/
  ey=min(ey,nty+1);
 /*correction si besoin pour le dernier bloc*/

  tab_bounds[0]=sx;
  tab_bounds[1]=ex;
  tab_bounds[2]=sy;
  tab_bounds[3]=ey;


  }

  /**************************************************************************************************/
  /*********************************VOISINAGE********************************************************/
  /**************************************************************************************************/

  void voisinage(MPI_Comm comm2d,int * voisin,int * coords, int * dims) {


  /* Recherche des voisins Nord et Sud */
  MPI_Cart_shift(comm2d, 0, 1, &(voisin[W]), &(voisin[E]));

  /* Recherche des voisins Ouest et Est */
  MPI_Cart_shift(comm2d, 1, 1, &(voisin[S]), &(voisin[N]));



  }



   double **allocarray(int Nlig,int Mcol) {
        double **array2 = malloc( Nlig* sizeof( double * ) );
        int i;

        if( array2 != NULL ){
            array2[0] = malloc(Nlig * Mcol * sizeof( double ) );
            if( array2[ 0 ] != NULL ) {
                for( i = 1; i < Nlig; i++ )
                    array2[i] = array2[0] + i * Mcol;
            }

            else {
                free(array2);
                array2 = NULL;
                printf("Erreur dans l'allocation mémoire -- phase 2\n");
            }
        }
        return array2;
    }

   void printarr(double **data, int nlig,int mcol, char *str) {
    printf("-- %s --\n", str);
    for (int i=0; i<nlig; i++) {
        for (int j=0; j<mcol; j++) {
            printf("%3f ", data[i][j]);
        }
        printf("\n");
    }
   }



 void ecrire_mpi(double *v_local_vect,int ntx,int nty,int * tab_bounds,MPI_Comm comm2d, char* nom){

  int code;
  MPI_File descripteur;
  int profil_tab[ndims], profil_sous_tab[ndims], coord_debut[ndims];
  MPI_Datatype type_sous_tab, type_sous_tab_vue;
  int profil_tab_vue[ndims], profil_sous_tab_vue[ndims], coord_debut_vue[ndims];
  MPI_Offset deplacement_initial;
  MPI_Status statut;

  /* Ouverture du fichier "final_v.dat" en écriture */
  code = MPI_File_open(comm2d, nom, MPI_MODE_WRONLY+MPI_MODE_CREATE,
		MPI_INFO_NULL, &descripteur);

 /* Test pour savoir si ouverture du fichier est correcte */
  if (code != MPI_SUCCESS) {
    fprintf(stderr, "ATTENTION erreur lors ouverture du fichier");
    MPI_Abort(comm2d, 2);
  }

  /* Creation du type derive type_sous_tab qui definit la matrice
   * sans les cellules fantomes */
  profil_tab[0] = tab_bounds[1]-tab_bounds[0]+3;
  profil_tab[1] = tab_bounds[3]-tab_bounds[2]+3;

  /* Profil du sous tableau */
  profil_sous_tab[0] =tab_bounds[1]-tab_bounds[0]+1;
  profil_sous_tab[1] =tab_bounds[3]-tab_bounds[2]+1;

  /* Coordonnees de depart du sous tableau */
  coord_debut[0] = 1;
  coord_debut[1] = 1;

  /* Creation du type_derive type_sous_tab */
  MPI_Type_create_subarray(ndims, profil_tab, profil_sous_tab, coord_debut,
			    MPI_ORDER_C, MPI_DOUBLE, &type_sous_tab);

  /* Validation du type_derive type_sous_tab */
  MPI_Type_commit(&type_sous_tab);

  /* Creation du type type_sous_tab_vue  pour la vue sur le fichier */
  /* Profil du tableau global */
  /*On ne tient pas compte des bords de l'image pour plus de simplicité*/
  profil_tab_vue[0] = ntx;
  profil_tab_vue[1] = nty;

  /* Profil du sous tableau */
  profil_sous_tab_vue[0] = tab_bounds[1]-tab_bounds[0]+1;
  profil_sous_tab_vue[1] = tab_bounds[3]-tab_bounds[2]+1;

  /* Coordonnees de depart du sous tableau */
  coord_debut_vue[0] = tab_bounds[0]-1;
  coord_debut_vue[1] = tab_bounds[2]-1;

  /* Creation du type_derive type_sous_tab_vue */
  MPI_Type_create_subarray(ndims, profil_tab_vue, profil_sous_tab_vue, coord_debut_vue,
			   MPI_ORDER_C, MPI_DOUBLE, &type_sous_tab_vue);

  /* Validation du type_derive type_sous_tab_vue */
  MPI_Type_commit(&type_sous_tab_vue);

  /* Définition de la vue sur le fichier a partir du debut */
  deplacement_initial = 0;
  MPI_File_set_view(descripteur, deplacement_initial, MPI_DOUBLE,
		    type_sous_tab_vue, "native", MPI_INFO_NULL);

  /* Ecriture du tableau u par tous les processus avec la vue */
  MPI_File_write_all(descripteur, v_local_vect, 1, type_sous_tab, &statut);

  /* Fermeture du fichier */
  MPI_File_close(&descripteur);
}

int malloc3ddouble(double ****array, int n, int m, int k) {

    /* allocate the n*m contiguous items */
    double *p = (double *)malloc(n*m*k*sizeof(double));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (double ***)malloc(n*sizeof(double**));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n; i++){
       (*array)[i] = (double **)malloc(sizeof(double *)*m);
    	  for(int j=0;j<m;j++)
	       (*array)[i][j]=&(p[i*m*k+j*k]);
    }
    return 0;
}

/************Utilisation de cette dernière fonction****************/
 /*double *** ndarray=NULL;
  malloc3ddouble(&ndarray,Nlig,Mcol,decomposition_step_number);*/
