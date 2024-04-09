#include <stdio.h>
#include <math.h>

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

struct vecteur distance_priv_tau(struct poisson p, struct poisson* zor, unsigned int nr, struct poisson* zoo, unsigned int no, struct poisson* zoa, unsigned int na){
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

int main(){
    printf("poisson");
}