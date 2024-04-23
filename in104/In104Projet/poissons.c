#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.>

#define IMAGE_WIDTH 1200
#define IMAGE_HEIGHT 800
#define PI 3.14159265358979323846
#define NB_POISSONS 100
struct vecteur{
    double i;
    double j;
};

struct poisson{
    double x;
    double y;
    struct vecteur v;
};

double norm2(struct vecteur v){
    double n = sqrt(v.i*v.i + v.j*v.j);
    return n;
}

double distance(struct vecteur v1, struct vecteur v2){
    return sqrt((v2.i - v1.i) * (v2.i - v1.i) + (v2.j - v1.j) * (v2.j - v1.j));
}

struct vecteur r(struct poisson pi, struct poisson pj){
    struct vecteur rij = {pj.x-pi.x, pj.y-pi.y};
    struct vecteur *r = &rij;
    double n = norm2(rij);
    r->i = r->i/n;
    r->j = r->j/n;
    return rij;
}

//Fonction qui renvoie la résultante de deux vecteurs 
struct vecteur somme_vecteurs(struct vecteur vi, struct vecteur vj)
{
    struct vecteur resultante;
    resultante.i = vi.i+vj.i;
    resultante.j = vi.j+vj.j;
    return resultante;
}

//Fonction qui détermine la direction privilégiée d_i

struct vecteur dir_priv_tau(struct poisson p, struct poisson* zor, unsigned int nr, struct poisson* zoo, unsigned int no, struct poisson* zoa, unsigned int na){
    struct vecteur d_i;
    struct vecteur d_r = {0, 0};
    struct vecteur d_o = {0, 0};
    struct vecteur d_a = {0, 0};
    
    if (nr ==0 && na ==0 && no==0) //le poisson n'a aucun voisins
    {
        d_i = p.v;
    }

    if (nr > 0) //des poissons voisins se trouvent dans la zone de répulsion
    {
        for(int j=0; j<nr;++j){
            d_r = somme_vecteurs(d_r, r(p,zor[j]));
        }
        struct vecteur *dr = &d_r;
        dr->i = -dr->i;
        dr->j = -dr->j;
    }
    else 
    {
        for(int j=0; j<no;++j){
            d_o = somme_vecteurs(d_o, zoo[j].v);
        }
        struct vecteur *d0 = &d_o;
        d0->i = -d0->i;
        d0->j = -d0->j;
        for(int j=0; j<na;++j){
            d_a = somme_vecteurs(d_a, r(p,zoa[j]));
        }
        if (no>0 && na==0) //des poissons voisins se trouvent uniquement dans la zone d'orientation
        {
            d_i = d_o;
        }

        if (na>0 && no ==0) //des poissons voisins se trouvent uniquement dans la zone d'attraction
        {
            d_i = d_a;
        }

        if (na>0 && no>0) //des poissons voisins se trouvent dans la zone d'attraction ET d'orientation
        {
            if (somme_vecteurs(d_a,d_o).i==0 && somme_vecteurs(d_a,d_o).j==0)//Si les forces résultantes sont égales à 0
            {
                d_i=p.v;
            }
            else 
            {
                d_i=somme_vecteurs(d_a,d_o);
                struct vecteur *di = &d_i;
                di->i = 0.5 * di->i;
                di->j = 0.5 * di->j;
            }
        }
    }
    return d_i;
}

double angle_entre_vecteurs(struct vecteur v1, struct vecteur v2) {
    double produit_scalaire = v1.i * v2.i + v1.j * v2.j;
    double norme_vecteur1 = norm2(v1);
    double norme_vecteur2 = norm2(v2);
    double cos_angle = produit_scalaire / (norme_vecteur1 * norme_vecteur2);
    
    // Utiliser la fonction acos() pour obtenir l'angle en radians
    double angle_radians = acos(cos_angle);
    
    return angle_radians;
}

//Simulation du mouvement des poissons

void simulation(struct poisson* poissons,double tau, double theta)
{     
    //Parcourons l'ensemble des poissons
    for (int i = 0; i<NB_POISSONS;i++) 
    {
        struct vecteur v1 = {poissons[i].x,poissons[i].y}; //vecteur position du poisson i 


        /////// DETERMINONS LES TABLEAUX VOISINS POUR LE POISSON I ////////////
        
        //Initialisation des tableaux des voisins pour chaque zone 
        struct poisson* zor = malloc (NB_POISSONS*sizeof(NB_POISSONS));
        struct poisson* zoa = malloc (NB_POISSONS*sizeof(NB_POISSONS));
        struct poisson* zoo = malloc (NB_POISSONS*sizeof(NB_POISSONS));

        int nr; //Compteurs du nb de poissons dans chaque zone
        int no;
        int na;

        //Remplissons ces tableaux pour le poisson i 
        for(int j = 0; j<NB_POISSONS; j++)  //on parcourt les voisins du poisson i pour leur assigner chacun une zone 
        {
            if (j!=i)
            {
                struct vecteur v2 = {poissons[j].x,poissons[j].y};
                unsigned int dist = distance(v1,v2);

                if (dist <=1) //Si le poisson est dans la zone de répulsion
                {
                    zor[nr]=poissons[j];
                    nr++;
                }

                if (1<dist && dist<=16)
                {
                    zoo[no]=poissons[j];
                    no++;
                }

                if (dist>16 && dist<=31)
                {
                    zoa[na]=poissons[j];
                    na++;
                }
            }
            free(zoo);
            free(zoa);
            free(zor);
        }
        struct vecteur d_i = dir_priv_tau(poissons[i], zor, nr, zoo, no, zoa, na);
        double angle = angle_entre_vecteurs(d_i, poissons[i].v);
    }
}

// Window
void render (SDL_Renderer *renderer , SDL_Texture **texture ) {
    SDL_SetRenderDrawColor ( renderer , 255 , 255 , 255 , 255) ;
    SDL_RenderClear ( renderer ) ;

    SDL_SetRenderDrawColor ( renderer , 0 , 0 , 255 , 255) ;
    SDL_Rect rect = { 400 , 400 , 10 , 10 } ;
    SDL_RenderFillRect ( renderer , &rect ) ;
    SDL_RenderCopy ( renderer , *texture , NULL, &rect ) ;
    // SDL RenderCopyEx ( renderer , ∗texture , NULL, &destRect , angle , NULL, SDL_FLIP_NONE);

    SDL_RenderPresent ( renderer ) ;
}


void loadTexture ( SDL_Renderer *renderer , SDL_Texture **texture ) {
    SDL_Surface *surface = IMG_Load( "fish_bitmap.jpg") ;
    if ( surface == NULL) {
        fprintf( stderr , "Failed to load image : %s \n" , IMG_GetError( ) ) ;
        SDL_Quit( ) ;
        exit( 1 ) ;
    }
    *texture = SDL_CreateTextureFromSurface ( renderer , surface) ;

    SDL_FreeSurface ( surface ) ;

    if (*texture == NULL) {
        fprintf( stderr , "Failed to create texture : %s \n" , SDL_GetError( ) ) ;
        SDL_Quit ( ) ;
        exit( 1 ) ;
    }
}

int main()
{
    if ( SDL_Init (SDL_INIT_VIDEO) < 0) {
        fprintf ( stderr , "SDL initialization failed : %s \n" , SDL_GetError( ) ) ;
        return 1 ;
    }

    SDL_Window *window = SDL_CreateWindow( "N-Body Simulation ", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, IMAGE_WIDTH, IMAGE_HEIGHT, SDL_WINDOW_SHOWN) ;
    if (window == NULL) {
        fprintf ( stderr , "Window creation failed : %s \n" , SDL_GetError( ) ) ;
        SDL_Quit ( ) ;
        return 1 ;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer (window , -1, SDL_RENDERER_ACCELERATED) ;
    if ( renderer == NULL) {
        fprintf( stderr , "Renderer creation failed : %s \n" , SDL_GetError( ) ) ;
        SDL_DestroyWindow(window) ;
        SDL_Quit ( ) ;
        return 1 ;
    }
    SDL_Texture *texture ;
    loadTexture(renderer, &texture);

    

    //Création des poissons 
    struct poisson poissons[NB_POISSONS];  //tableau contenant tous les poissons
    for (int i =0;i<NB_POISSONS; i++)  //On parcourt tous les poissons
    {
        //on donne à chaque poisson une position aléatoire dans la fenetre
        poissons[i].x = (double)(rand()% IMAGE_WIDTH);
        poissons[i].y = (double)(rand() %IMAGE_HEIGHT);

        //On donne à chaque poisson un vecteur direction initial nul
        poissons[i].v.i = 0.0;
        poissons[i].v.j = 0.0;
    }

    SDL_Event event ;
    int quit = 0 ;
    while ( ! quit ) {
        while ( SDL_PollEvent(&event ) != 0) {
            if ( event.type == SDL_QUIT) {
                quit = 1 ;
            }
        }

        //Simulation du mouvement des poissons
        simulation(poissons,)

        // Render the updated positions
        render ( renderer , &texture ) ;

        // Delay to control the frame rate
        // SDL Delay ( 1 ) ;
    }

    

    SDL_DestroyRenderer ( renderer) ;
    SDL_DestroyWindow(window) ;
    SDL_Quit ( ) ;

    return 0 ;

}