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

// pour itérer facilement sur un tableau
#define foreach(item, array) \
    for(int keep = 1, \
            count = 0,\
            size = sizeof (array) / sizeof *(array); \
        keep && count != size; \
        keep = !keep, count++) \
      for(item = (array) + count; keep; keep = !keep)

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

  
  int i,j,k,l;
   double c01,c02,c03,c04;
   double c0;
   double t;
   double ux,uy,Gradu;

/*******Parametre lambda initial****************/
double lambda0=0.0005;

 
  

  void initialisation_mpi(int, char**);
  void finalisation_mpi();
  void domaine(MPI_Comm ,int ,int *,int,int ,int * ,int *,int * ,int ,int );
  void voisinage(MPI_Comm,int *,int *, int *);
  double **allocarray(int,int); 
  void printarr(double **, int,int, char *); 
  void ecrire_mpi(double *,int,int,int *,MPI_Comm);
  int malloc3ddouble(double ****, int , int , int );
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
  if(argc >= 2){
	dims[0] = (int)  *argv[1];
  }else{ 
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
  if(rang==3){
  printf("Je suis le rang %d \n"
  "coordonnee en ligne du coin superieur gauche sx=%d\n"
  "coordonnee en colonne du coin superieur gauche sy=%d\n"
  "coordonnee en ligne du coin inferieur droit ex=%d\n"
  "coordonnee en colonne du coin inferieur droit ey=%d\n"
  "-----------------------------------------\n", 
	   rang,tab_bounds[0],tab_bounds[2],tab_bounds[1],tab_bounds[3]);
  }


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

  // allocation locale de mémoire pour chaque processus

  double ***ndarray = NULL;
  malloc3ddouble(&ndarray,subsizes[0],subsizes[1],decomposition_step_number);

 
  

   /* Mesure du temps en seconde dans la boucle en temps */
   t1 = MPI_Wtime();

   /*BOUCLE PRINCIPALE*/
   MPI_Barrier(comm2d);
    printf("Debut de la boucle principale pour le processus %d\n", rang);
   for (l=0; l<decomposition_step_number;l++){
       // Initialisation
        for (i=0; i < subsizes[0]; i++){
            for (j=0; j < subsizes[1]; j++){
                ndarray[i][j][l] = f_local_mat[i][j];
            }
        }

        for (k=0; k < Itermax; k++){
            for (i=1; i < subsizes[0]-1; i++){
                for (j=1; j < subsizes[1]-1; j++){

                    ux = ndarray[i+1][j][l] - ndarray[i][j][l];
                    uy = (ndarray[i][j+1][l] - ndarray[i][j-1][l])/2.0;
                    Gradu=sqrt(epsilon+ux*ux+uy*uy);
                    c01=1.0/Gradu;


                    ux=(ndarray[i][j][l]-ndarray[i-1][j][l]);
                    uy=(ndarray[i-1][j+1][l]-ndarray[i-1][j-1][l])/2.0;
                    Gradu=sqrt(epsilon+ux*ux+uy*uy);
                    c02=1.0/Gradu;



                    ux=(ndarray[i+1][j][l]-ndarray[i-1][j][l])/2.0;
                    uy=(ndarray[i][j+1][l]-ndarray[i][j][l]);
                    Gradu=sqrt(epsilon+ux*ux+uy*uy);
                    c03=1.0/Gradu;



                    ux=(ndarray[i+1][j-1][l]-ndarray[i-1][j-1][l])/2.0;
                    uy=(ndarray[i][j][l]-ndarray[i][j-1][l]);
                    Gradu=sqrt(epsilon+ux*ux+uy*uy);
                    c04=1.0/Gradu;



                    c0=2.0*h*h*lambda0+c01+c02+c03+c04;

                    t=2.0*h*h*lambda0*f_local_mat[i][j]+c01*ndarray[i+1][j][l]+c02*ndarray[i-1][j][l]+c03*ndarray[i][j+1][l]+c04*ndarray[i][j-1][l];

                    ndarray[i][j][l]=(1.0/c0)*t;
                }
            }

            // Gestion aux bords

            // déclaration des tableaux contenant les bords
            double *top = malloc((subsizes[0]-1)*sizeof(double));
            double *bottom = malloc((subsizes[0]-1)*sizeof(double));
            double *right = malloc((subsizes[0]-1)*sizeof(double));
            double *left = malloc((subsizes[0]-1)*sizeof(double));

            for (i = 0; i < subsizes[0]-1; i++){
                top[i] = ndarray[i][1][l];
                bottom[i] = ndarray[i][subsizes[0]-1][l];
            }

            for (j = 0; j < subsizes[1]-1; j++){
                left[i] = ndarray[i][1][l];
                bottom[i] = ndarray[i][subsizes[0]-1][l];
            }

            // on itère pour chaque voisin présent
            foreach(double *v, voisin){

                // voisin du dessus
                if (*v == 0){
                    for (j = 1; j < subsizes[0]; j++){
                        MPI_Send(&ndarray[0][j][l], 1, MPI_DOUBLE, voisin[N], 10,  comm2d);
                        MPI_Recv(&bottom[j], 1, MPI_DOUBLE, voisin[S], 10, comm2d, &statut);
                    }
                }

                // voisin à droite
                else if (*v == 1){
                    for (i = 1; j < subsizes[1]; i++){
                        MPI_Send(&ndarray[i][subsizes[1]][l], 1, MPI_DOUBLE, voisin[E], 10, comm2d);
                        MPI_Recv(&left[i], 1, MPI_DOUBLE, voisin[W], 10, comm2d, &statut);
                    }
                    
                }

                // voisin du dessous
                else if (*v == 2){
                   for (j = 1; j < subsizes[0]; j++){
                       MPI_Recv(&top[j], 1, MPI_DOUBLE, voisin[N], 20, comm2d, &statut);
                       MPI_Send(&ndarray[subsizes[0]][j][l], 1, MPI_DOUBLE, voisin[S], 10, comm2d);
                   }
                }

                // voisin à gauche
                else if (*v == 3){
                    for (i = 1; i < subsizes[1]; i++){
                        MPI_Send(&ndarray[0][i][l], 1, MPI_DOUBLE, voisin[W], 20, comm2d);
                        MPI_Recv(&right[i], 1, MPI_DOUBLE, voisin[E], 20, comm2d, &statut);
                    }
                }
            }
            
        }
   }


  /* Mesure du temps a la sortie de la boucle */
   t2 = MPI_Wtime();

  

  /* Ecriture des resultats pour chaque processus */
  /* FONCTION ecriture_mpi */ 
   
  /* Affichage du temps de convergence par le processus 3 */
  if (rang == 3) {
   
  
  printf("Convergence en %f secs\n", t2-t1);
  }
  /****Libération mémoire****/

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

  

 void ecrire_mpi(double *v_local_vect,int ntx,int nty,int * tab_bounds,MPI_Comm comm2d){

  int code;
  MPI_File descripteur;
  int profil_tab[ndims], profil_sous_tab[ndims], coord_debut[ndims];
  MPI_Datatype type_sous_tab, type_sous_tab_vue;
  int profil_tab_vue[ndims], profil_sous_tab_vue[ndims], coord_debut_vue[ndims];
  MPI_Offset deplacement_initial;
  MPI_Status statut;

  /* Ouverture du fichier "final_v.dat" en écriture */
  code = MPI_File_open(comm2d, "final_v.dat", MPI_MODE_WRONLY+MPI_MODE_CREATE, 
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
