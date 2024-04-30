# In104Projet
# in104
# Adèle Migozzi
# Noémie Lafine

Le but est de simuler les mouvements d'un banc de poissons. Pour cela, on crée un tableau renseignant les vecteurs vitesses et les positions pour chaque poisson. A chaque intervalle de temps tau, on met à jour ce tableau et on transmet l'information à un logiciel graphique représentant le banc de poisson. Pour tau assez petit, on voit le banc évoluer presque en continu. 
Pour mettre à jour le tableau, on s'intéresse à chaque poisson chacun leur tour. Soit un poisson p quelconque qu'on étudie. On va récupérer quels sont ses voisins (en fonction de leur proximité au poisson). On fait alors évoluer les paramètres de chaque poisson voisin v. Il y a trois zones différentes autour du poisson p, en fonction de celle où se trouve le poisson v, il devra adopter une atitude unique. 
