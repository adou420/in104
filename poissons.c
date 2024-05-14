#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#define WINDOW_WIDTH 1200
#define WINDOW_HEIGHT 800

#define FISH_HEIGHT 20
#define FISH_RATIO 1.549
#define FISH_WIDTH FISH_RATIO * FISH_HEIGHT


#define PI 3.14159265358979323846
#define NB_POISSONS 100
#define S 5.0
#define THETA 0.5

#define V_I_INIT 5;
#define V_J_INIT 6;

struct vecteur{
    double i;
    double j;
};

struct poisson{
    double x;
    double y;
    struct vecteur v;
};

//fonction qui renvoie la norme d'un vecteur 
double norm2(struct vecteur v){
    return sqrt(v.i*v.i + v.j*v.j);
}

//fonction qui renvoie la distance entre deux positions
double distance(struct vecteur v1, struct vecteur v2){
    return sqrt((v2.i - v1.i) * (v2.i - v1.i) + (v2.j - v1.j) * (v2.j - v1.j));
}

//fonction qui renvoie le vecteur unitaire liant deux poissons
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

double angle_entre_vecteurs(struct vecteur v1, struct vecteur v2) {
    return atan2(v2.j * v1.i - v2.i * v1.j, v2.i * v1.i + v2.j * v1.j);  //renvoie un angle entre -pi et pi

}

//Fonction qui détermine la direction privilégiée d_i du poisson qu'on considère
struct vecteur dir_priv_tau(struct poisson p, struct poisson* zor, int nr, struct poisson* zoo, int no, struct poisson* zoa, int na){
   
    struct vecteur d_r = {0, 0};
    struct vecteur d_o = {0, 0};
    struct vecteur d_a = {0, 0};
    
    if (nr ==0 && na ==0 && no==0) //le poisson n'a aucun voisins
    {
        return p.v;
    }

    if (nr > 0) //des poissons voisins se trouvent dans la zone de répulsion
    {
        for(int j = 0; j < nr; ++j) {
            d_r = somme_vecteurs(d_r, r(p, zor[j]));
        }

        d_r.i *= -1 ;
        d_r.j *= -1 ;

        return d_r;
    }
    
    //on est mnt dans le cas où il n'y a aucun poissons voisins dans la zone de répulsion
    for(int j=0; j<no;++j){
        d_o = somme_vecteurs(d_o, zoo[j].v);
    }
    d_o.i *= -1;
    d_o.j *= -1;

    for(int j=0; j<na;++j){
        d_a = somme_vecteurs(d_a, r(p, zoa[j]));
    }
    

    if (no>0 && na==0) //des poissons voisins se trouvent uniquement dans la zone d'orientation
    {
        return d_o;
    }

    else if (na>0 && no ==0) //des poissons voisins se trouvent uniquement dans la zone d'attraction
    {
        return d_a;
    }

    else //(na>0 && no>0) //des poissons voisins se trouvent dans la zone d'attraction ET d'orientation
    {
        if (somme_vecteurs(d_a, d_o).i==0 && somme_vecteurs(d_a, d_o).j==0)//Si les forces résultantes sont égales à 0
        {
            return p.v;
        }
        else 
        {
            struct vecteur d_i=somme_vecteurs(d_a, d_o);
            d_i.i=0.5*d_i.i;
            d_i.j=0.5*d_i.j;
            return d_i;
        }
    }
}

// On utilise l'algorithme de Box-Muller pour générer un nombre aléatoire gaussien
double generate_random_noise(double mean, double stddev) {
   
    double u1 = ((double) rand() / RAND_MAX);
    double u2 = ((double) rand() / RAND_MAX);

    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    return mean + stddev * z0;
}

// Faisons en sorte que les poissons ne se superposent pas

void resolve_collision(struct poisson* p1, struct poisson* p2) {
    double distance = sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));  // redondant avec la fonction distance, il va falloir arranger ça
    double min_distance = FISH_WIDTH;

    if (distance < min_distance) {
        // Calculons le vecteur de séparation
        double dx = p2->x - p1->x;
        double dy = p2->y - p1->y;
        double angle = atan2(dy, dx);
        double separation = min_distance - distance;

        // Déplacer les poissons pour les séparer
        p1->x -= cos(angle) * separation / 2;
        p1->y -= sin(angle) * separation / 2;
        p2->x += cos(angle) * separation / 2;
        p2->y += sin(angle) * separation / 2;
    }
    return;
}

void check_collisions(struct poisson* poissons) {
    for (int i = 0; i < NB_POISSONS - 1; i++) {
        for (int j = i + 1; j < NB_POISSONS; j++) {
            resolve_collision(&poissons[i], &poissons[j]);
        }
    }
    return;
}





//Simulation du mouvement des poissons

void simulation(struct poisson* poissons, double tau, double alpha)
{        
    //Parcourons l'ensemble des poissons
    for (int i = 0; i < NB_POISSONS;i++) 
    {
        struct vecteur v1 = {poissons[i].x, poissons[i].y}; //vecteur position du poisson i 


        /////// DETERMINONS LES TABLEAUX VOISINS POUR LE POISSON I ////////////
        
        //Initialisation des tableaux des voisins pour chaque zone 
        struct poisson* zor = malloc (NB_POISSONS*sizeof(struct poisson));
        struct poisson* zoa = malloc (NB_POISSONS*sizeof(struct poisson));
        struct poisson* zoo = malloc (NB_POISSONS*sizeof(struct poisson));

        int nr=0; //Compteurs du nb de poissons dans chaque zone
        int no=0;
        int na=0;

        //Remplissons ces tableaux pour le poisson i 
        for(int j = 0; j < NB_POISSONS; j++)  //on parcourt les voisins du poisson i pour leur assigner chacun une zone 
        {
            struct vecteur v2 = {poissons[j].x, poissons[j].y}; 
            struct vecteur v1_oppose = {-v1.i, -v1.j};
            
            if (j != i && angle_entre_vecteurs(poissons[i].v, somme_vecteurs(v2, v1_oppose)) < alpha/2)
            {
                int dist = distance(v1, v2);

                if (dist <= 1) //Si le poisson est dans la zone de répulsion
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

            
        }


        //Calculons di+tau
        struct vecteur d_i = dir_priv_tau(poissons[i], zor, nr, zoo, no, zoa, na);
        
        // Ajout de bruit gaussien à la direction préférée
        double noise_x = generate_random_noise(0.0,0.1);
        double noise_y = generate_random_noise(0.0,0.1);
        struct vecteur noisy_direction = {d_i.i + noise_x, d_i.j + noise_y};

        poissons[i].v = noisy_direction;    

        //Calculons les nouvelles positions du poisson i
        double nv_x = poissons[i].x + poissons[i].v.i*tau*S;
        double nv_y = poissons[i].y + poissons[i].v.j*tau*S;
        

        // Vérifions si le poisson atteint le bord de la fenêtre et inversons sa direction si nécessaire
        if (nv_x < 0 || nv_x + FISH_WIDTH > WINDOW_WIDTH) {
            poissons[i].v.i *= -1;
        }
        if (nv_y < 0 || nv_y + FISH_HEIGHT > WINDOW_HEIGHT) {
            poissons[i].v.j *= -1;
        }

        // Mettons à jour les positions des poissons
        poissons[i].x = nv_x;
        poissons[i].y = nv_y;

        free (zoa);
        free (zor);
        free (zoo);
    }

    // Faisons en sorte que les poissons ne se superposent pas
    //check_collisions(poissons);
    
}


void loadTexture(SDL_Renderer *renderer, SDL_Texture **texture) {
    SDL_Surface *surface = IMG_Load("fish.png");
    if(surface == NULL) {
        fprintf(stderr, "Failed to load image : %s \n", IMG_GetError());
        SDL_Quit();
        exit(1);
    }
    *texture = SDL_CreateTextureFromSurface(renderer, surface);

    SDL_FreeSurface(surface);

    if (*texture == NULL) {
        fprintf(stderr, "Failed to create texture : %s \n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }
}

// Window
void render(SDL_Renderer *renderer, SDL_Texture **texture, struct poisson* p) {
    // if (p->v.j < 0 && p->y + FISH_HEIGHT < 0)
    // {
    //     p->y = WINDOW_HEIGHT + FISH_HEIGHT;
    // }

    // if (p->v.j > 0 && p->y > WINDOW_HEIGHT)
    // {
    //     p->y = -FISH_HEIGHT;
    // }

    // if (p->v.i < 0 && p->x + FISH_WIDTH < 0)
    // {
    //     p->x = WINDOW_WIDTH + FISH_WIDTH;
    // }

    // if (p->v.i > 0 && p->x > WINDOW_WIDTH)
    // {
    //     p->x = -FISH_WIDTH;
    // }

    
    SDL_Rect rect = {(int)p->x, (int)p->y, FISH_WIDTH, FISH_HEIGHT };
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    // SDL_RenderCopy(renderer, *texture, NULL, &rect);

    // Rotation
    struct vecteur const vecteur_horizontal = {1,0};
    double alpha = angle_entre_vecteurs(p->v,vecteur_horizontal);
    if (alpha < 0) 
    {
        alpha *= -1;
    } 
    else
    {
        alpha = 2*PI - alpha;
    }
    SDL_RenderCopyEx(renderer, *texture, NULL, &rect, (180/PI)*alpha, NULL, SDL_FLIP_NONE); //rotation de l'image ATTENTION, l'angle doit être en degrés

    int const x0 = (int)p->x + FISH_WIDTH / 2;
    int const y0 = (int)p->y + FISH_HEIGHT / 2;
    SDL_RenderDrawLine(renderer, x0, y0, x0 + p->v.i * 5, y0 + p->v.j * 5);
}

int main()
{
    if (SDL_Init (SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL initialization failed : %s \n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow("N-Body Simulation ",
                                          SDL_WINDOWPOS_UNDEFINED, 
                                          SDL_WINDOWPOS_UNDEFINED,
                                          WINDOW_WIDTH,
                                          WINDOW_HEIGHT,
                                          SDL_WINDOW_SHOWN);
    if (window == NULL) {
        fprintf(stderr, "Window creation failed : %s \n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer (window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == NULL) {
        fprintf(stderr, "Renderer creation failed : %s \n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    SDL_Texture *texture;
    loadTexture(renderer, &texture);

    

    //Création des poissons 
    struct poisson* poissons = malloc(NB_POISSONS * sizeof(struct poisson));  //tableau contenant tous les poissons
    for (int i = 0; i < NB_POISSONS; i++)  //On parcourt tous les poissons
    {
        //on donne à chaque poisson une position aléatoire dans la fenetre
        poissons[i].x = (double)(rand() % WINDOW_WIDTH);
        poissons[i].y = (double)(rand() % WINDOW_HEIGHT);

        //On donne à chaque poisson un vecteur direction initial 
        poissons[i].v.i = V_I_INIT;
        poissons[i].v.j = V_J_INIT;
    }

    //On crée un tableau pour garder en mémoire les anciennes positions + directions des poissons
    struct poisson * ancien_poissons = malloc(NB_POISSONS*sizeof(struct poisson));


    SDL_Event event;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }
        }

        // On garde en mémoire l'ancien tableau de poissons
        memcpy(ancien_poissons,poissons,NB_POISSONS*sizeof(struct poisson));


        //Simulation du mouvement des poissons : on met a jour le tableau des poissons
        simulation(poissons, 0.1, 4.36);

        // Render the updated positions

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);
        // SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);

        for (int i = 0; i < NB_POISSONS; i++) {
            render(renderer, &texture, poissons + i);
        }

        SDL_RenderPresent(renderer);

        // Delay to control the frame rate
        // SDL_Delay(30);
    }

    

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    free(poissons);
    free(ancien_poissons);
    return 0;
}