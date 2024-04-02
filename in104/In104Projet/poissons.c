#include <stdio.h>
#include <math.h>

struct vecteur{
    double i;
    double j;
    double k;
};

struct poisson{
    double x;
    double y;
    double z;
    struct vecteur v;
};

double norm2(struct vecteur v){
    double n = sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
    return n;
}

struct vecteur r(struct poisson pi, struct poisson pj){
    struct vecteur rij = {pj.x-pi.x, pj.y-pi.y, pj.z-pi.z};
    struct vecteur *r = &rij;
    double n = norm2(rij);
    r->i = r->i/n;
    r->j = r->j/n;
    r->k = r->k/n;
    return rij;
}
//Fonction qui détermine la direction privilégiée d_i

struct vecteur distance_priv_tau(struct poisson p, struct poisson* zor, int nr, struct poisson* zoo, int no, struct poisson* zoa, int na){
    struct vecteur d_i;
    struct vecteur d_r;
    struct vecteur d_o;
    struct vecteur d_a;

}

int main(){
    printf("poisson");
}