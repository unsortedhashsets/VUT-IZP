/***********************************************************
*     Soubor:  proj3.c                                     *
*     Datum:   3.12.2018                                   *
*     Autor:   Abramov Mikhail, xabram00@stud.fit.vutbr.cz *
*     Projekt: IZP c. 3 - Prace s datovymi strukturami     *
***********************************************************/
/**
  * Kostra programu pro 3. projekt IZP 2018/19
  *
  * Jednoducha shlukova analyza: 2D nejblizsi soused.
  * Single linkage
  */
 #include <stdio.h>
 #include <stdlib.h>
 #include <assert.h>
 #include <math.h> // sqrtf
 #include <limits.h> // INT_MAX
 
 /*****************************************************************
  * Ladici makra. Vypnout jejich efekt lze definici makra
  * NDEBUG, napr.:
  *   a) pri prekladu argumentem prekladaci -DNDEBUG
  *   b) v souboru (na radek pred #include <assert.h>
  *      #define NDEBUG
  */
 #ifdef NDEBUG
 #define debug(s)
 #define dfmt(s, ...)
 #define dint(i)
 #define dfloat(f)
 #else
 
 // vypise ladici retezec
 #define debug(s) printf("- %s\n", s)
 
 // vypise formatovany ladici vystup - pouziti podobne jako printf
 #define dfmt(s, ...) printf(" - "__FILE__":%u: "s"\n",__LINE__,__VA_ARGS__)
 
 // vypise ladici informaci o promenne - pouziti dint(identifikator_promenne)
 #define dint(i) printf(" - " __FILE__ ":%u: " #i " = %d\n", __LINE__, i)
 
 // vypise ladici informaci o promenne typu float - pouziti
 // dfloat(identifikator_promenne)
 #define dfloat(f) printf(" - " __FILE__ ":%u: " #f " = %g\n", __LINE__, f)
 
 #endif
 
 /*****************************************************************
  * Deklarace potrebnych datovych typu:
  *
  * TYTO DEKLARACE NEMENTE
  *
  *   struct obj_t - struktura objektu: identifikator a souradnice
  *   struct cluster_t - shluk objektu:
  *      pocet objektu ve shluku,
  *      kapacita shluku (pocet objektu, pro ktere je rezervovano
  *          misto v poli),
  *      ukazatel na pole shluku.
  */
 
 struct obj_t {
     int id;
     float x;
     float y;
 };
 
 struct cluster_t {
     int size;
     int capacity;
     struct obj_t *obj;
 };
 
 /*****************************************************************
  * Deklarace potrebnych funkci.
  *
  * PROTOTYPY FUNKCI NEMENTE
  *
  * IMPLEMENTUJTE POUZE FUNKCE NA MISTECH OZNACENYCH 'TODO'
  *
  */
 
 /*
  Inicializace shluku 'c'. Alokuje pamet pro cap objektu (kapacitu).
  Ukazatel NULL u pole objektu znamena kapacitu 0.
 */
 void init_cluster(struct cluster_t *c, int cap)
 {
 assert(c != NULL);
 assert(cap >= 0);
     c->size = 0;    
     if (cap > 0) {
         if ((c->obj = malloc(cap*sizeof(struct obj_t)))) { //memory allocation for object
             c->capacity = cap;
             return;
         }
     }
     c->capacity = 0;
     c->obj = NULL;
     return;
 }
 
 /*
  Odstraneni vsech objektu shluku a inicializace na prazdny shluk.
  */
 void clear_cluster(struct cluster_t *c)
 {
     free(c->obj);
     c->size = 0;
     c->capacity = 0;
     c->obj = NULL;
 }
 
 /// Chunk of cluster objects. Value recommended for reallocation.
 const int CLUSTER_CHUNK = 10;
 
 /*
  Zmena kapacity shluku 'c' na kapacitu 'new_cap'.
  */
 struct cluster_t *resize_cluster(struct cluster_t *c, int new_cap)
 {
     // TUTO FUNKCI NEMENTE
     assert(c);
     assert(c->capacity >= 0);
     assert(new_cap >= 0);
 
     if (c->capacity >= new_cap)
         return c;
 
     size_t size = sizeof(struct obj_t) * new_cap;
 
     void *arr = realloc(c->obj, size);
     if (arr == NULL)
         return NULL;
 
     c->obj = (struct obj_t*)arr;
 c->capacity = new_cap;
 return c;
 }
 
 /*
  Prida objekt 'obj' na konec shluku 'c'. Rozsiri shluk, pokud se do nej objekt
  nevejde.
  */
 void append_cluster(struct cluster_t *c, struct obj_t obj)
 {
     if (c->capacity <= c->size)
     {
         resize_cluster(c, c->capacity+CLUSTER_CHUNK); //resize on CHUNK volume
     }
     c->obj[c->size++] = obj;     //object append
 }
 
 /*
  Seradi objekty ve shluku 'c' vzestupne podle jejich identifikacniho cisla.
  */
 void sort_cluster(struct cluster_t *c);
 
 /*
  Do shluku 'c1' prida objekty 'c2'. Shluk 'c1' bude v pripade nutnosti rozsiren.
  Objekty ve shluku 'c1' budou serazeny vzestupne podle identifikacniho cisla.
  Shluk 'c2' bude nezmenen.
  */
 void merge_clusters(struct cluster_t *c1, struct cluster_t *c2)
 {
 assert(c1 != NULL);
 assert(c2 != NULL);
 int OLDsize = c1->size;
 for(int i=0; i < c2->size; i++){
     append_cluster(c1, c2->obj[i]);
 }
 if(c2->size > 0 && c1->size == OLDsize + c2->size) {
         sort_cluster(c1);
     }
 }
 
 /**********************************************************************/
 /* Prace s polem shluku */
 
 /*
  Odstrani shluk z pole shluku 'carr'. Pole shluku obsahuje 'narr' polozek
  (shluku). Shluk pro odstraneni se nachazi na indexu 'idx'. Funkce vraci novy
  pocet shluku v poli.
 */
 int remove_cluster(struct cluster_t *carr, int narr, int idx)
 {
     assert(idx < narr);
     assert(narr > 0);
     int NEWnarr=narr - 1;
     clear_cluster(&carr[idx]);
     for(int i=idx;i<NEWnarr;i++) {
         carr[i] = carr[i + 1];
     }
     return NEWnarr;
 }
 
 /*
  Pocita Euklidovskou vzdalenost mezi dvema objekty.
  */
 float obj_distance(struct obj_t *o1, struct obj_t *o2)
 {
     assert(o1 != NULL);
     assert(o2 != NULL);
     float D_EUKLID = sqrtf((o1->x-o2->x)*(o1->x-o2->x)+(o1->y-o2->y)*(o1->y-o2->y)); //EUKLID distance for nearest neighbrouh
     return D_EUKLID;
 }
 
 /*
  Pocita vzdalenost dvou shluku.
 */
 float cluster_distance(struct cluster_t *c1, struct cluster_t *c2)
 {
     assert(c1 != NULL);
     assert(c1->size > 0);
     assert(c2 != NULL);
     assert(c2->size > 0);
     float MIN_DIST;
     float DIST;
     MIN_DIST = obj_distance(&c1->obj[0], &c2->obj[0]);
     for (int i = 1; i < c1->size; i++) {
         for (int j = i; j < c2->size; j++) {
             DIST = obj_distance(&c1->obj[i], &c2->obj[j]);
             if (DIST<MIN_DIST){
                 MIN_DIST=DIST;
             }
         }
     }
     return MIN_DIST;
 }
 
 /*
  Funkce najde dva nejblizsi shluky. V poli shluku 'carr' o velikosti 'narr'
  hleda dva nejblizsi shluky. Nalezene shluky identifikuje jejich indexy v poli
  'carr'. Funkce nalezene shluky (indexy do pole 'carr') uklada do pameti na
  adresu 'c1' resp. 'c2'.
 */
 void find_neighbours(struct cluster_t *carr, int narr, int *c1, int *c2)
 {
     assert(narr > 0);
     if(narr==1){    //clusters number control
         *c1=*c2=0;
         return;
     }
     float MINdist=-1;	//base distance for MIN searching
     float dist;	//variable distance
     for(int i=0;i<narr;i++){  //matrix of min distance between clusters creation 
         for(int j=i+1;j<narr;j++){
             dist=cluster_distance(&carr[i],&carr[j]);	
             if(MINdist==-1||dist<MINdist){
                 MINdist=dist;
                 *c1=i;
                 *c2=j;
             }
         }
     }
 }
 
 // pomocna funkce pro razeni shluku
 static int obj_sort_compar(const void *a, const void *b)
 {
     // TUTO FUNKCI NEMENTE
     const struct obj_t *o1 = (const struct obj_t *)a;
     const struct obj_t *o2 = (const struct obj_t *)b;
     if (o1->id < o2->id) return -1;
     if (o1->id > o2->id) return 1;
     return 0;
 }
 
 /*
  Razeni objektu ve shluku vzestupne podle jejich identifikatoru.
 */
 void sort_cluster(struct cluster_t *c)
 {
     // TUTO FUNKCI NEMENTE
     qsort(c->obj, c->size, sizeof(struct obj_t), &obj_sort_compar);
 }
 
 /*
  Tisk shluku 'c' na stdout.
 */
 void print_cluster(struct cluster_t *c)
 {
     // TUTO FUNKCI NEMENTE
     for (int i = 0; i < c->size; i++)
     {
         if (i) putchar(' ');
         printf("%d[%g,%g]", c->obj[i].id, c->obj[i].x, c->obj[i].y);
     }
     putchar('\n');
 }
 
 /*
  Ze souboru 'filename' nacte objekty. Pro kazdy objekt vytvori shluk a ulozi
  jej do pole shluku. Alokuje prostor pro pole vsech shluku a ukazatel na prvni
  polozku pole (ukalazatel na prvni shluk v alokovanem poli) ulozi do pameti,
  kam se odkazuje parametr 'arr'. Funkce vraci pocet nactenych objektu (shluku).
  V pripade nejake chyby uklada do pameti, kam se odkazuje 'arr', hodnotu NULL.
 */
 int load_clusters(char *filename, struct cluster_t **arr) {
     assert(arr != NULL);
     FILE *file = fopen(filename, "r");	    //file opening
     if (!file) {			    //ERROR if file do not exist	
         printf("reading ERROR (1)\n");
         return -1;
     }
     struct obj_t obj;			//objects declaration
     struct cluster_t *cluster;		//clusters declaration
     int RealFILEstring_n = 0;		//strings count from file
     int FILEstring_n = 0;		//strings from head of file			
     char string[50];			//standard string -> 50 symbols
     char error;			//an array for validation
     int O_ID = 0;			//object ID
     int O_X = 0;			//object X
     int O_Y = 0;			//object Y
     while (fgets(string, 12, file)) {  //reading symbols from FILE
         RealFILEstring_n++;		//strings count from file = strings count from file + 1
         if (RealFILEstring_n == 1) {   //on first line control of the file head
             if (sscanf(string, "count=%d%[^\n]", &FILEstring_n, &error) != 1) { //head-form control
                 printf("form error in the %d string\n", FILEstring_n); 
                 fclose(file);
                 return -1;
             }
             if (FILEstring_n <= 0) { 						 //head-volume control
                 printf("volume error N=%d\n (1)", FILEstring_n);
                 fclose(file);
                 return -1;
             }
             if ( !(*arr = malloc(FILEstring_n*sizeof(struct cluster_t)))) {	//allocation control
               fclose(file);  
	       return -1;
             }
             for (int i=0;i<FILEstring_n;i++){					//clusters inicialisation
                 init_cluster(&(*arr)[i], 0);
             }
             continue;
         }
         if (sscanf(string, "%d %d %d%[^\n]", &O_ID, &O_X, &O_Y, &error) != 3 //objects validation
	   || O_X < 0 || O_X > 1000||O_Y < 0 || O_Y> 1000 || O_ID < 1) {
             printf("volume error N=%d (2)\n ", FILEstring_n);
	     for (int i=0;i<(FILEstring_n);i++){
	        clear_cluster(&(*arr)[i]);
	     }
	     fclose(file);
	     free(*arr);
	     return -1;
         }
         cluster = &(*arr)[RealFILEstring_n - 2];	//CLUSTER inicialisation
         obj.x = O_X;					
         obj.y = O_Y;
         obj.id = O_ID;
         append_cluster(cluster, obj);		//object to cluster
     }
 fclose(file);						//the end of file reading
 return (RealFILEstring_n-1);				//return of OBJECTS count
 } 
 /*
  Tisk pole shluku. Parametr 'carr' je ukazatel na prvni polozku (shluk).
  Tiskne se prvnich 'narr' shluku.
 */
 void print_clusters(struct cluster_t *carr, int narr)
 {
     printf("Clusters:\n");
     for (int i = 0; i < narr; i++)
     {
         printf("cluster %d: ", i);
         print_cluster(&carr[i]);
     }
 }
 int main(int argc, char *argv[]) {
     struct cluster_t *clusters;              // clusters structure declaration
     int NUMBER_CLUSTERS_COMMAND = 1;         // minimum number of clusters
     int NUMBER_CLUSTERS_FILE = 1;            // minimum number of clusters
     if (argc > 1 && argc <= 3) {	      // arguments validation
         if (argc == 3) {
             char *err = NULL;
             NUMBER_CLUSTERS_COMMAND = (int) strtol(argv[2], &err, 10); //NUMBER_CLUSTERS_COMMAND validation
             if (*err || NUMBER_CLUSTERS_COMMAND < 1) {			//filtration of error entrance
                 printf("N must be whole number greater than ZERO\n");
                 return -1;
             }
         }
         NUMBER_CLUSTERS_FILE = load_clusters(argv[1], &clusters); //function LOAD calling  
	 if (NUMBER_CLUSTERS_FILE==-1) return -1;
         if(NUMBER_CLUSTERS_FILE<NUMBER_CLUSTERS_COMMAND){		//entrance clusters validation
	   printf("Clusters number (%d) in file is lower than Clusters number in command(%d)\n",NUMBER_CLUSTERS_FILE,NUMBER_CLUSTERS_COMMAND);
	   for (int i=0;i<NUMBER_CLUSTERS_FILE;i++){			//clear clusters memmory free
	      clear_cluster(&clusters[i]);
	   }
	   free(clusters);
	   return -1;
	 }
         while(NUMBER_CLUSTERS_FILE>NUMBER_CLUSTERS_COMMAND) {		//inicialisation of clusterisation process
	     int c1_idx, c2_idx;
             find_neighbours(clusters, NUMBER_CLUSTERS_FILE, &c1_idx, &c2_idx);			//function find calling
             merge_clusters(&clusters[c1_idx], &clusters[c2_idx]);				//function merge calling
             NUMBER_CLUSTERS_FILE = remove_cluster(clusters, NUMBER_CLUSTERS_FILE, c2_idx);	//function remove calling
         }
	 sort_cluster(clusters);
         print_clusters(clusters, NUMBER_CLUSTERS_FILE);		//function print calling
         for (int i=0;i<NUMBER_CLUSTERS_FILE;i++){			//clear clusters memmory free
	 clear_cluster(&clusters[i]);
	 }
	 free(clusters);
         return 0;
     }
     else{
       printf("The number of arguments must be 2-3 (filename) or (filename N)\n");
       return -1;
    }
 }     
