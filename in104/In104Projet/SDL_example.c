#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <stdio.h>
#include <stdlib.h>

#define IMAGE_WIDTH 1200
#define IMAGE_HEIGHT 800

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
    SDL_Surface *surface = IMG_Load( "fish_bitmap.jpg ") ;
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


int main ( ) {
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


    SDL_Event event ;
    int quit = 0 ;
    while ( ! quit ) {
        while ( SDL_PollEvent(&event ) != 0) {
            if ( event.type == SDL_QUIT) {
                quit = 1 ;
            }
        }



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