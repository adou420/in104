#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>

#define WINDOW_WIDTH 1200
#define WINDOW_HEIGHT 800

#define FISH_HEIGHT 20
#define FISH_RATIO 1.549
#define FISH_WIDTH FISH_RATIO * FISH_HEIGHT

#define PRED_HEIGHT 20
#define PRED_RATIO 1.549
#define PRED_WIDTH FISH_RATIO * FISH_HEIGHT

#define PI 3.14159265358979323846

#define NB_POISSONS 100
#define V_I_INIT sqrt(2)/2
#define V_J_INIT sqrt(2)/2
#define V_PREDA 20.0 
#define RAYON_CHASSE 100


double S = 30.0;
double ALPHA = 4.36;
double RAYON_REPULSION = 0.5;
double RAYON_ALIGN = 16;
double RAYON_ATTRAC = 31;

struct vecteur{
    double i;
    double j;
};

struct poisson{
    double x;
    double y;
    struct vecteur v;
    int type; //0 pour les poissons et 1 pour le prédateur
};

//fonction qui renvoie la norme d'un vecteur 
double norm2(struct vecteur v){
    return sqrt(v.i*v.i + v.j*v.j);
}

//fonction qui renvoie le vecteur liant deux poissons
struct vecteur r(struct poisson pi, struct poisson pj){
    struct vecteur rij = {pj.x-pi.x, pj.y-pi.y};
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
struct vecteur dir_priv_tau (struct poisson p, struct poisson predateur, struct poisson* zor, int nr, struct poisson* zoo, int no, struct poisson* zoa, int na, bool np){
   
    struct vecteur d_r = {0, 0};
    struct vecteur d_o = {0, 0};
    struct vecteur d_a = {0, 0};
    
    if (np)
    {
        struct vecteur d_p = r(predateur,p);
        d_p.i *= 1/norm2(d_p) ;
        d_p.j *= 1/norm2(d_p) ;
        return d_p;
    }

    if (nr == 0 && na == 0 && no == 0) //le poisson n'a aucun voisins
    {
        return p.v;
    }

    if (nr > 0) //des poissons voisins se trouvent dans la zone de répulsion
    {
        for(int j = 0; j < nr; ++j) {
            d_r = somme_vecteurs(d_r, r(p, zor[j]));
        }

        d_r.i *= -1/norm2(d_r) ;
        d_r.j *= -1/norm2(d_r) ;

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
        d_o.i *= 1 / norm2(d_o);
        d_o.j *= 1 / norm2(d_o);
        return d_o;
    }

    else if (na>0 && no ==0) //des poissons voisins se trouvent uniquement dans la zone d'attraction
    {
        d_a.i *= 1 / norm2(d_a);
        d_a.j *= 1 / norm2(d_a);
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
            struct vecteur d_i = somme_vecteurs(d_a, d_o);
            d_i.i *= 0.5;
            d_i.j *= 0.5;

            d_i.i *= 1 / norm2(d_i);
            d_i.j *= 1 / norm2(d_i);
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

void separateur_2poissons(struct poisson* p1, struct poisson* p2) {
    double distance = sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));  
    double min_distance = FISH_WIDTH;

    if (distance < min_distance) {
        // Calculons le vecteur de séparation
        double dx = p2->x - p1->x;
        double dy = p2->y - p1->y;
        double angle = atan2(dy, dx);
        double separation = min_distance - distance;

        // Déplaçons les poissons pour les séparer
        p1->x -= cos(angle) * separation / 2;
        p1->y -= sin(angle) * separation / 2;
        p2->x += cos(angle) * separation / 2;
        p2->y += sin(angle) * separation / 2;
    }
    return;
}

void separation_poissons(struct poisson* poissons) {
    for (int i = 0; i < NB_POISSONS - 1; i++) {
        for (int j = i + 1; j < NB_POISSONS; j++) {
            separateur_2poissons(&poissons[i], &poissons[j]);
        }
    }
    return;
}

//Fonction pour générer le mouvement du prédateur
void mvt_predateur(struct poisson* predateur, struct poisson* poissons, double tau) {
    
    //Déterminons le poisson le plus proche 
    struct poisson poisson_le_plus_proche = poissons[0];
    double min_distance = sqrt((predateur->x - poissons[0].x) * (predateur->x - poissons[0].x) + (predateur->y - poissons[0].y) * (predateur->y - poissons[0].y));
    for (int i=1; i < NB_POISSONS; i++){
        double distance = sqrt((predateur->x - poissons[i].x) * (predateur->x - poissons[i].x) + (predateur->y - poissons[i].y) * (predateur->y - poissons[i].y));
        if (distance < min_distance){
            poisson_le_plus_proche = poissons[i];
            min_distance = distance;
        }
    }

    //Mettre à jour la direction du prédateur pour qu'il aille vers le poisson le plus proche
    predateur->v.i = r(*predateur, poisson_le_plus_proche).i/norm2(r(*predateur, poisson_le_plus_proche));
    predateur->v.j = r(*predateur, poisson_le_plus_proche).j/norm2(r(*predateur, poisson_le_plus_proche));


    // Déplace le prédateur dans la nouvelle direction 
    predateur->x += predateur->v.i * tau * V_PREDA;
    predateur->y += predateur->v.j * tau * V_PREDA;
}



//Simulation du mouvement des poissons

void simulation(struct poisson* poissons, double tau, struct poisson* predateur)
{        
    //Parcourons l'ensemble des poissons
    for (int i = 0; i < NB_POISSONS;i++) 
    {
        struct vecteur v1 = {poissons[i].x, poissons[i].y}; //vecteur position du poisson i 

        bool np = false;
        double dist_pred = sqrt((v1.i - predateur->x) * (v1.i - predateur->x) + (v1.j - predateur->y) * (v1.j - predateur->y)); 
        if (dist_pred < RAYON_CHASSE)
        {
            np = true;
        }

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

            //Vérifions si le poisson voisin est dans l'angle mort
            bool hors_angle_mort = true;
            double angle = angle_entre_vecteurs(poissons[i].v,r(poissons[i],poissons[j]));
            if (fabs(angle) > PI - ALPHA/2)
            {hors_angle_mort = false;}
            
            if ( j != i && hors_angle_mort)
            {
                int dist = sqrt((v1.i - v2.i) * (v1.i - v2.i) + (v1.j - v2.j) * (v1.j - v2.j)); 

                if (dist <= RAYON_REPULSION) //Si le poisson est dans la zone de répulsion
                {
                    zor[nr] = poissons[j];
                    nr++;
                }

                if (RAYON_REPULSION < dist && dist <= RAYON_ALIGN)
                {
                    zoo[no] = poissons[j];
                    no++;
                }

                if (dist > RAYON_ALIGN && dist <= RAYON_ATTRAC)
                {
                    zoa[na] = poissons[j];
                    na++;
                }
            }

            
        }


        //Calculons di+tau
        struct vecteur d_i = dir_priv_tau(poissons[i], *predateur, zor, nr, zoo, no, zoa, na, np);   //d_i est unitaire
        
        // Ajout de bruit gaussien à la direction préférée
        double noise_x = generate_random_noise(0.0,0.1);
        double noise_y = generate_random_noise(0.0,0.1);
        struct vecteur noisy_direction = {d_i.i + noise_x, d_i.j + noise_y};

        //on rend le vecteur unitaire
        noisy_direction.i *= 1 / norm2(noisy_direction);
        noisy_direction.j *= 1 / norm2(noisy_direction);

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
    //separation_poissons(poissons);

    mvt_predateur(predateur, poissons, tau);
    
}


void loadTexture(SDL_Renderer *renderer, SDL_Texture **texture, const char* image) {
    //Texture des poissons
    SDL_Surface *surface = IMG_Load(image);
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

void travers_bords (struct poisson* p) {

    int HEIGHT = FISH_HEIGHT;
    int WIDTH = FISH_WIDTH;
    if (p->type == 1){
        HEIGHT = PRED_HEIGHT;
        WIDTH = PRED_WIDTH;
    }
    if (p->v.j < 0 && p->y + HEIGHT < 0)
    {
        p->y = WINDOW_HEIGHT + HEIGHT;
    }

    if (p->v.j > 0 && p->y > WINDOW_HEIGHT)
    {
        p->y = -HEIGHT;
    }

    if (p->v.i < 0 && p->x + WIDTH < 0)
    {
        p->x = WINDOW_WIDTH + WIDTH;
    }

    if (p->v.i > 0 && p->x > WINDOW_WIDTH)
    {
        p->x = -WIDTH;
    }
    return;
}

//Window
void render(SDL_Renderer *renderer, SDL_Texture **texture, SDL_Texture **pred_texture, struct poisson* p) {
    
    //travers_bords(p);
    
    SDL_Rect rect = {(int)p->x, (int)p->y, FISH_WIDTH, FISH_HEIGHT };
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    // Rotation
    struct vecteur const vecteur_horizontal = {1,0};
    double theta = angle_entre_vecteurs(p->v,vecteur_horizontal);
    if (theta < 0) 
    {
        theta *= -1;
    } 
    else
    {
        theta = 2*PI - theta;
    }
    
    
    if (p->type == 0)
    {
        SDL_RenderCopyEx(renderer, *texture, NULL, &rect, (180/PI)*theta, NULL, SDL_FLIP_NONE); //rotation de l'image ATTENTION, l'angle doit être en degrés
    }

    else if (p->type == 1)
    {
        SDL_RenderCopyEx(renderer, *pred_texture, NULL, &rect, (180 / PI) * theta, NULL, SDL_FLIP_NONE);
    }
    

    int const x0 = (int)p->x + FISH_WIDTH / 2;
    int const y0 = (int)p->y + FISH_HEIGHT / 2;
    SDL_RenderDrawLine(renderer, x0, y0, x0 + p->v.i * 5, y0 + p->v.j * 5);

}

#define SLIDER_WIDTH 200
#define SLIDER_HEIGHT 20
#define SLIDER_X 50
#define SLIDER_Y 700

typedef struct {
    SDL_Rect rect;
    double value;
    double min;
    double max;
}Slider;

bool sliderActive = false;
Slider *activeSlider = NULL;

void drawSlider(SDL_Renderer *renderer, Slider *slider) {
    // Draw the background of the slider
    SDL_SetRenderDrawColor(renderer, 200, 200, 200, 255);
    SDL_RenderFillRect(renderer, &slider->rect);

    // Draw the current value of the slider
    int handle_x = slider->rect.x + (int)((slider->value - slider->min) / (slider->max - slider->min) * slider->rect.w);
    SDL_Rect handle = { handle_x - 5, slider->rect.y - 5, 10, slider->rect.h + 10 };
    SDL_SetRenderDrawColor(renderer, 100, 100, 255, 255);
    SDL_RenderFillRect(renderer, &handle);
}


bool handleSliderEvent(SDL_Event *event, Slider *slider, bool *sliderActive, Slider **activeSlider) {
    int x = event->button.x;
    int y = event->button.y;

    switch (event->type) {
        case SDL_MOUSEBUTTONDOWN:
            if (x >= slider->rect.x && x <= slider->rect.x + slider->rect.w &&
                y >= slider->rect.y && y <= slider->rect.y + slider->rect.h) {
                *sliderActive = true;
                *activeSlider = slider;
                slider->value = slider->min + (double)(x - slider->rect.x) / slider->rect.w * (slider->max - slider->min);
                return true;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            if (*sliderActive && *activeSlider == slider) {
                *sliderActive = false;
                *activeSlider = NULL;
                return true;
            }
            break;

        case SDL_MOUSEMOTION:
            if (*sliderActive && *activeSlider == slider) {
                slider->value = slider->min + (double)(x - slider->rect.x) / slider->rect.w * (slider->max - slider->min);
                return true;
            }
            break;
    }
    return false;
}


void renderText(SDL_Renderer *renderer, TTF_Font *font, const char *text, int x, int y, SDL_Color color) {
    SDL_Surface *surface = TTF_RenderText_Blended(font, text, color);
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);
    SDL_Rect dstrect = { x, y, surface->w, surface->h };
    SDL_RenderCopy(renderer, texture, NULL, &dstrect);
    SDL_FreeSurface(surface);
    SDL_DestroyTexture(texture);
}



int main()
{
    if (SDL_Init (SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL initialization failed : %s \n", SDL_GetError());
        return 1;
    }
    if (TTF_Init() == -1) {
    fprintf(stderr, "SDL_ttf initialization failed: %s\n", TTF_GetError());
    SDL_Quit();
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
    loadTexture(renderer, &texture,"fish.png");

    SDL_Texture *pred_texture;
    loadTexture(renderer, &pred_texture,"predateur.png");


    // Load font
    TTF_Font *font = TTF_OpenFont("ArialMT.ttf", 24);
    if (font == NULL) {
        fprintf(stderr, "Failed to load font: %s\n", TTF_GetError());
        TTF_Quit();
        SDL_Quit();
        return 1;
    }

    // Sliders initialization
    Slider speedSlider = { { SLIDER_X, SLIDER_Y, SLIDER_WIDTH, SLIDER_HEIGHT }, S, 1.0, 80.0 };   //Slider pour la vitesse
    Slider alphaSlider = { { SLIDER_X, SLIDER_Y - 70, SLIDER_WIDTH, SLIDER_HEIGHT }, ALPHA, 0, 5 };   // Slider pour alpha (angle mort)
    Slider rayon_repulSlider = { { SLIDER_X, SLIDER_Y - 140, SLIDER_WIDTH, SLIDER_HEIGHT }, S, 0.0, 20.0 };   //Slider pour le rayon de la zone de répulsion
    Slider rayon_alignSlider = { { SLIDER_X, SLIDER_Y - 210, SLIDER_WIDTH, SLIDER_HEIGHT }, S, RAYON_REPULSION, 60.0 };   //Slider pour le rayon de la zone d'alignement
    Slider rayon_attracSlider = { { SLIDER_X, SLIDER_Y - 280, SLIDER_WIDTH, SLIDER_HEIGHT }, S, RAYON_ALIGN, 100.0 };   //Slider pour le rayon de la zone d'attraction

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

        poissons[i].type = 0;
    }

    //Création du prédateur
    struct poisson predateur;
    predateur.x = (double)(rand() % WINDOW_WIDTH);
    predateur.y = (double)(rand() % WINDOW_HEIGHT);
    predateur.v.i = 1;
    predateur.v.j = 1;
    predateur.type = 1;

    SDL_Event event;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }

            // Handle slider events
            if (handleSliderEvent(&event, &speedSlider, &sliderActive, &activeSlider)) {
                S = speedSlider.value;
            }

            if (handleSliderEvent(&event, &alphaSlider, &sliderActive, &activeSlider)){
                ALPHA = alphaSlider.value;
            }

            if (handleSliderEvent(&event, &rayon_repulSlider, &sliderActive, &activeSlider)){
                RAYON_REPULSION = rayon_repulSlider.value;
            }

            if (handleSliderEvent(&event, &rayon_alignSlider, &sliderActive, &activeSlider)){
                RAYON_ALIGN = rayon_alignSlider.value;
            }

            if (handleSliderEvent(&event, &rayon_attracSlider, &sliderActive, &activeSlider)){
                RAYON_ATTRAC = rayon_attracSlider.value;
            }
        }

        // Simulation du mouvement des poissons : on met a jour le tableau des poissons
        simulation(poissons, 0.1,&predateur);


        // Rendu des positions mises à jour
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        for (int i = 0; i < NB_POISSONS; i++) {
            render(renderer, &texture, &pred_texture, poissons + i);
        }

        render(renderer, &texture, &pred_texture, &predateur);

        // Draw the sliders
        SDL_Color textColor = {0, 0, 0, 255};  // Black color
        drawSlider(renderer, &speedSlider);
        renderText(renderer, font, "Vitesse S", SLIDER_X, SLIDER_Y - 30, textColor);

        drawSlider(renderer, &alphaSlider);
        renderText(renderer, font, "alpha", SLIDER_X, SLIDER_Y - 100, textColor);

        drawSlider(renderer, &rayon_repulSlider);
        renderText(renderer, font, "rayon zone de repulsion", SLIDER_X, SLIDER_Y - 170, textColor);

        drawSlider(renderer, &rayon_alignSlider);
        renderText(renderer, font, "rayon zone d'alignement", SLIDER_X, SLIDER_Y - 240, textColor);

        drawSlider(renderer, &rayon_attracSlider);
        renderText(renderer, font, "rayon zone d'attraction", SLIDER_X, SLIDER_Y - 310, textColor);

        SDL_RenderPresent(renderer);
    }
    
    // Delay to control the frame rate
    // SDL_Delay(30);

    TTF_CloseFont(font);
    TTF_Quit();
    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    free(poissons);
    return 0;
}