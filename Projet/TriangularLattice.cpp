#include "TriangularLattice.h"

#include <iostream>
#include <random>
#include <stdexcept>
#include <algorithm> //min,max
#include <fstream>

std::random_device rd;
std::mt19937 gen(rd());

//__________________________________________________________Classe Lattice__________________________________________________________\\


//La Classe Lattice est une classe nous permettant de créer un réseau triangulaire sur lequel siègent
//des hexagones réguliers qui peuvent se déplacer et s'orienter selon 6 directions (k\pi/3).
//
//Un site vide est repéré par un 0
//Un site occupé est repéré par un entier k compris entre 1 et 6 repérant son orientation (k\pi/3)
//par rapport à un axe horizontal orienté vers la droite
//
//Les données sont stockées dans un entier signé sur 8 bits, cela permet de réduire l'occupation de la mémoire


//_____________________________Méthodes

void Lattice::afficher() const{
   for(int i=0; i<nx; i++){
    if(i%2==1){
        std::cout << "  ";
    }    
    for(int j=0; j<ny; j++){
        std::cout << (int) data[i*ny+j] << "   ";
    };
    std::cout << "\n";
    
   };
}

//Crée un hexagone de côté a pixels
void hexagon(sf::RenderWindow& window, int8_t phi0, float a, float pos_x, float pos_y){

    sf::ConvexShape triangle;
    triangle.setPointCount(3);
    triangle.setPoint(0, sf::Vector2f(0, 0));
    triangle.setPoint(1, sf::Vector2f(a*cos(M_PI/6), a*sin(M_PI/6)));
    triangle.setPoint(2, sf::Vector2f(a*cos(M_PI/6), -a*sin(M_PI/6)));

    triangle.setPosition(pos_x, pos_y);
    triangle.setRotation((-phi0%6)*60);

    for(int i=0;i<6;i++){
        
        triangle.rotate(60);
        
        sf::Color triangleColor = sf::Color(30 * i, 50 * i, 70 * i);
        triangle.setFillColor(triangleColor);
        
        window.draw(triangle);
    }
}

//Affiche le Lattice dans une window SFML, les sites occupés contiennent des hexagones définis ci-dessus
void Lattice::affiche_SFML(sf::RenderWindow& window,const float a) const{
    float taille=a*sqrt(3)/2;

    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            int8_t valeur=data[i*ny+j];
            if(valeur!=0){
                hexagon(window,valeur,a,taille*(2*j+1+i%2),a*(1.5*i+1));
            }
        }        
    }
}

//Permet d'accéder au site (x,y) du Lattice
//Les conditions aux limites périodiques de Born-Von Kármán sont implémentées ici
Site Lattice::site_xy(const int x,const int y) const{
    int xm=x%nx;
    int ym=y%ny;
    if (xm < 0 && ym < 0){
        xm+=nx; ym+=ny;
        return Site(ny*xm+ym,xm,ym);
    }
    else if (xm < 0){
        xm+=nx;
        return Site(ny*xm+ym,xm,ym);
    }
    else if (ym < 0){
        ym+=ny;
        return Site(ny*xm+ym,xm,ym);
    }
    else{
        return Site(ny*xm+ym,xm,ym);
    }
}

//Permet de tirer un site aléatoirement
Site Lattice::site_aleatoire() const{
    
    std::uniform_int_distribution<int> distribution(0, nx*ny);
    
    int _index=distribution(gen);
    int y=_index%ny;
    if(y < 0){
        return Site(_index,_index/ny,y+ny);
    }
    else{
        return Site(_index,_index/ny,y);
    }
}

//Permet d'accéder au site situé à la k-ème position dans le pointeur
//double* data,  k = ny*x+y
Site Lattice::site_index(const int index) const{
    int y=index%ny;
    if(y < 0){
        return Site(index,index/ny,y+ny);
    }
    else{
        return Site(index,index/ny,y);
    }
}

//Renvoie une liste des premiers voisins de contact du Site s (définis dans la partie Renormalisation du rapport)
//Tourne dans le sens trigo, en partant de la droite
std::array<Site,6> Lattice::voisins(const Site s) const{
    int x=s._x, y=s._y;
    if(x%2 ==0){
        return {this->site_xy(x,y+1), this->site_xy(x-1,y), this->site_xy(x-1,y-1), this->site_xy(x,y-1), this->site_xy(x+1,y-1), this->site_xy(x+1,y)};
    }
    else{
        return {this->site_xy(x,y+1), this->site_xy(x-1,y+1), this->site_xy(x-1,y), this->site_xy(x,y-1), this->site_xy(x+1,y), this->site_xy(x+1,y+1)};
    }
}

//Renvoie une liste des seconds voisins de contact du Site s (définis dans la partie Renormalisation du rapport)
//Tourne dans le sens trigo, en partant de la droite
std::array<Site,6> Lattice::second_voisins(const Site s) const{
    int x=s._x, y=s._y;
    if(x%2 ==0){
        return {this->site_xy(x,y+2), this->site_xy(x-2,y+1), this->site_xy(x-2,y-1), this->site_xy(x,y-2), this->site_xy(x+2,y-1), this->site_xy(x+2,y+1)};
    }
    else{
        return {this->site_xy(x,y+2), this->site_xy(x-2,y+1), this->site_xy(x-2,y-1), this->site_xy(x,y-2), this->site_xy(x+2,y-1), this->site_xy(x+2,y+1)};
    }
}

//Compte le nombre de particules présentes sur le Lattice
//Tous les sites ne sont pas forcément occupé
int Lattice::Particule_Count() const{
    int count =0;

    for(int k=0;k<nx*ny;k++){
        if(data[k]!=0){
            count+=1;
        }
    }
    return count;
}

//Renvoie une Matrix Particules de taille N*2 contenant les coordonées des N particules
//On accède aux coordonnées x et y de la particule i par Particules(i, 0) et Particules(i, 1) respectivement
//La connaissance de la positions des particules est primordial car nous travaillons à basse densité (N_parts = 0.1 * N_Sites)
//Cela nous permettra dans l'algorithme de Metropolis  de tirer aléatoirement un site occupé
Matrix Lattice::Pos_Particules(const Lattice &L) const{
    int Nbr_particules= L.Particule_Count();    
    Matrix Particules = Matrix(Nbr_particules,2);

    int compteur = 0;
    for(int i=0;i<L.nx;i++){
        for(int j=0;j<L.ny;j++){
            if((int) L[L.site_xy(i,j)] != 0){
                Particules(compteur,0)=i;
                Particules(compteur,1)=j;
                compteur++;
            }
        }
    }

    return Particules;
}

//Cette méthode nous permet de renormaliser notre Lattice (Voir rapport partie Renormalisation)
Lattice Lattice::Renormalisation() const{
    Lattice R = Lattice(nx/2,ny/2);

    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            if((i%2==0) && (j%2==0)){
                if(i%4==2){
                    R[R.site_xy(i/2,j/2)] = (*this)[this->site_xy(i,j+1)];
                }
                else{
                    R[R.site_xy(i/2,j/2)] = (*this)[this->site_xy(i,j)];
                }
            }
        }
    }
    return R;
}

//_____________________________Opérateurs

int8_t Lattice::operator[] (Site s) const{
    return data[s._index];
}

int8_t& Lattice::operator[] (Site s){
    return data[s._index];
}

int ConComp::operator[] (SiteC s) const{
    return data[s._index];
}

int& ConComp::operator[] (SiteC s){
    return data[s._index];
}
//_____________________________Constructeurs

Lattice::Lattice(int nx_, int ny_) : nx(nx_), ny(ny_), data(nullptr){
    data= new int8_t[nx*ny];
    
    if(ny%2==1){
        throw std::runtime_error("ny has to be even to enforce BVK periodic conditions");
    }
    for(int i=0;i<nx*ny;i++){
        data[i]=0;
    }
}

Lattice::~Lattice(){
    if (data != nullptr){
        delete[] data;
    }
    data = nullptr;
}


//__________________________________________________________Classe ConComp__________________________________________________________\\


//La Classe ConComp est une classe nous permettant de repérer et numéroter les différentes composantes connexes formées par les
//particules hexagonales qui siègent sur un Lattice L
//
//Un site vide est repéré par un 0
//Un site occupé est repéré par un entier k > 1 qui renseigne sur l'appartenance de la particule à la composante connexe numéro k
//
//Isoler les composantes connexes va se montrer crucial dans l'algorithme générationnel afin de caractériser la phase qui règne 
//sur le Lattice
//
//Cette fois, les données sont stockées dans un int car le nombre de composantes connexes peut rapidement dépasser 127


//_____________________________Méthodes

//Echelle de couleur qui simule la colormap Viridis de Python car elle est jolie
//Le code ne viens pas de nous
// t dans [0, 1]
sf::Color viridisColor(float t) {
    // Viridis color map in RGB
    const std::vector<std::vector<unsigned char>> viridisColors = {
        {253, 231, 37}, {246, 230, 32}, {236, 229, 27}, {226, 228, 24}, {218, 227, 25},
        {208, 225, 28}, {197, 224, 33}, {189, 223, 38}, {178, 221, 45}, {168, 219, 52},
        {157, 217, 59}, {149, 216, 64}, {139, 214, 70}, {129, 211, 77}, {122, 209, 81},
        {112, 207, 87}, {103, 204, 92}, {94, 201, 98}, {88, 199, 101}, {80, 196, 106},
        {72, 193, 110}, {66, 190, 113}, {59, 187, 117}, {53, 183, 121}, {47, 180, 124},
        {44, 177, 126}, {39, 173, 129}, {36, 170, 131}, {34, 167, 133}, {32, 163, 134},
        {31, 160, 136}, {30, 156, 137}, {31, 153, 138}, {31, 149, 139}, {32, 146, 140},
        {33, 143, 141}, {34, 139, 141}, {36, 135, 142}, {37, 132, 142}, {38, 129, 142},
        {40, 125, 142}, {41, 121, 142}, {42, 118, 142}, {44, 114, 142}, {46, 111, 142},
        {47, 108, 142}, {49, 104, 142}, {50, 100, 142}, {52, 96, 141}, {54, 93, 141},
        {56, 89, 140}, {58, 84, 140}, {59, 81, 139}, {61, 77, 138}, {63, 72, 137},
        {65, 68, 135}, {66, 64, 134}, {68, 59, 132}, {69, 55, 129}, {70, 51, 127},
        {71, 46, 124}, {72, 41, 121}, {72, 36, 117}, {72, 32, 113}, {72, 27, 109},
        {72, 22, 104}, {72, 17, 100}, {71, 11, 94}, {70, 5, 89}, {68, 1, 84}
    };

    // Ensure t is in the [0, 1] range
    t = 1 - std::max(0.0f, std::min(1.0f, t));

    // Scale t to the range [0, viridisColors.size() - 1]
    float index = t * (viridisColors.size() - 1);

    // Get the lower and upper color indices
    int lowerIndex = static_cast<int>(std::floor(index));
    int upperIndex = std::min(lowerIndex + 1, static_cast<int>(viridisColors.size()) - 1);

    // Calculate the interpolation factor
    float factor = index - lowerIndex;

    // Interpolate between the lower and upper colors
    sf::Color lowerColor(viridisColors[lowerIndex][0], viridisColors[lowerIndex][1], viridisColors[lowerIndex][2]);
    sf::Color upperColor(viridisColors[upperIndex][0], viridisColors[upperIndex][1], viridisColors[upperIndex][2]);

    sf::Color interpolatedColor(
        static_cast<unsigned char>((1.0f - factor) * lowerColor.r + factor * upperColor.r),
        static_cast<unsigned char>((1.0f - factor) * lowerColor.g + factor * upperColor.g),
        static_cast<unsigned char>((1.0f - factor) * lowerColor.b + factor * upperColor.b)
    );

    return interpolatedColor;
}

//Crée un hexagone de côté a pixels avec, en son centre le numéro de la composante connexe dont la particule fait partie
void hexagonOneColor(sf::RenderWindow& window, float viridis, float a, float pos_x, float pos_y, int ConCompNumber){

    sf::ConvexShape triangle;
    triangle.setPointCount(3);
    triangle.setPoint(0, sf::Vector2f(0, 0));
    triangle.setPoint(1, sf::Vector2f(a*cos(M_PI/6), a*sin(M_PI/6)));
    triangle.setPoint(2, sf::Vector2f(a*cos(M_PI/6), -a*sin(M_PI/6)));

    triangle.setPosition(pos_x, pos_y);

    for(int i=0;i<6;i++){
        
        triangle.rotate(60);
        
        sf::Color triangleColor = viridisColor(viridis);
        triangle.setFillColor(triangleColor);
        
        window.draw(triangle);
    }

    sf::Text text;
    sf::Font font;
    
    if (!font.loadFromFile("Arial.ttf")) {
        throw "The font didn't load !, check for a file Arial.ttf";
        return;
    }

    font.loadFromFile("Arial.ttf");
    text.setFont(font);

    text.setString(std::to_string(ConCompNumber));
    text.setCharacterSize(3*a/4);
    text.setFillColor(sf::Color::Black);

    sf::FloatRect textBounds = text.getLocalBounds();
    text.setPosition(pos_x - textBounds.width / 2.0f, pos_y - textBounds.height / 1.0f);

    window.draw(text);
}

//Affiche les composantes connexes dans une window SFML
void ConComp::Show_Connected_Components(const float a) const{

    float taille=a*sqrt(3)/2;

    sf::RenderWindow window(sf::VideoMode(1500, 900), "Connected Components of The Lattice");

    window.clear(sf::Color(0,0,0));
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            int valeur=data[i*ny+j];
            if(valeur!=0){
                hexagonOneColor(window, (float)valeur /NbrCC , a, taille*(2*j+1+i%2), a*(1.5*i+1), (int) valeur);
            }
        }        
    }
    window.display();

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
    }
}

SiteC ConComp::site_xy(int x, int y) const{
    int xm=x%nx;
    int ym=y%ny;
    if (xm < 0 && ym < 0){
        xm+=nx; ym+=ny;
        return SiteC(ny*xm+ym,xm,ym);
    }
    else if (xm < 0){
        xm+=nx;
        return SiteC(ny*xm+ym,xm,ym);
    }
    else if (ym < 0){
        ym+=ny;
        return SiteC(ny*xm+ym,xm,ym);
    }
    else{
        return SiteC(ny*xm+ym,xm,ym);
    }
}

SiteC ConComp::site_aleatoire() const{
    
    std::uniform_int_distribution<int> distribution(0, nx*ny);
    
    int _index=distribution(gen);
    int y=_index%ny;
    if(y < 0){
        return SiteC(_index,_index/ny,y+ny);
    }
    else{
        return SiteC(_index,_index/ny,y);
    }
}

//index=ny*x+y
SiteC ConComp::site_index(int index) const{
    int y=index%ny;
    if(y < 0){
        return SiteC(index,index/ny,y+ny);
    }
    else{
        return SiteC(index,index/ny,y);
    }
}

//Renvoie une liste des premiers voisins de contact du SiteC s
//Tourne dans le sens trigo, en partant de la droite
std::array<SiteC,6> ConComp::voisins(const SiteC s) const{
    int x=s._x, y=s._y;
    if(x%2 ==0){
        return {this->site_xy(x,y+1), this->site_xy(x-1,y), this->site_xy(x-1,y-1), this->site_xy(x,y-1), this->site_xy(x+1,y-1), this->site_xy(x+1,y)};
    }
    else{
        return {this->site_xy(x,y+1), this->site_xy(x-1,y+1), this->site_xy(x-1,y), this->site_xy(x,y-1), this->site_xy(x+1,y), this->site_xy(x+1,y+1)};
    }
}

//Donne la longueur de la bordure extérieure d'une composante connexe
//Notre algorithme est inspiré de l'algorithme de Moore (Voir partie correspondante dans le rapport)
int ConComp::OuterBorderLength(const int ConCompNumber) const{
    if(ConCompNumber > NbrCC){
        throw std::invalid_argument("There are not that many connected components");
    }

    
    int k = nx*ny - 1;
    for(k = nx*ny - 1; data[k] != ConCompNumber; k--){}
    int Length = 0;

    SiteC StartingSite = (*this).site_index(k);
    std::array<SiteC,6> Voisins = (*this).voisins(StartingSite);
    int j=0;
    for(int i = 0; i<6 && (*this)[Voisins[i]] == 0; i++){
        j++;
    }
    if(j == 6){
        return 6;
    }
    
    int ExitFace = 5;
    int FirstExitFace = 5;
    int EntryFace = 5;
    while((*this)[Voisins[(ExitFace-1)%6]] == 0){
        ExitFace--;
        FirstExitFace--;
        EntryFace--;
    }
    while((*this)[Voisins[EntryFace%6]] == 0){
        EntryFace++;
    }
    Length+= EntryFace - ExitFace;

    SiteC s = Voisins[EntryFace%6];
    ExitFace = (EntryFace + 4)%6;

    while ( !((s._index == StartingSite._index) && (ExitFace == FirstExitFace))){
        Voisins = (*this).voisins(s);
        
        EntryFace = ExitFace;
        
        while((*this)[Voisins[EntryFace%6]] == 0){
            EntryFace++;
        }
        Length+= EntryFace - ExitFace;
        s = Voisins[EntryFace%6];
        ExitFace = (EntryFace + 4)%6;
    }
    
    return Length;

}

//Renvoie une matrice M où M(i, 0) renvoie la taille de la composante connexe i+1
Matrix ConComp::SizeConComps() const{
    Matrix Size = Matrix(NbrCC, 1);
    //Size(i-1,0) -> Taille de la composante connexe i

    for(int k=0; k < nx*ny; k++){
        if(!(data[k] == 0 || data[k] == -1)){
            Size(data[k]-1,0)++;
        }
    }

    return Size;
}

//Check si la composante connexe est sur le bord de gauche
bool ConComp::IsOnLeftEdge(const int ConCompNumber) const{
    for(int x=0; x<nx; x++){
        if ((*this)[(*this).site_xy(x,0)] == ConCompNumber){
            return true;
        }
    }
    return false;
}

//Check si la composante connexe est sur le bord supérieur
bool ConComp::IsOnTopEdge(const int ConCompNumber) const{
    for(int y=0; y<ny; y++){
        if ((*this)[(*this).site_xy(0,y)] == ConCompNumber){
            return true;
        }
    }
    return false;
}

//Isole la composante connexe sur un nouveau réseau dédié
ConComp ConComp::isolateConComp(const int ConCompNumber) const{
    if(ConCompNumber > NbrCC){
        throw std::invalid_argument("There are not that many connected components");
    }
    
    int x=0;
    int y=0;
    int k=0;

    //Si la composante connexe se trouve sur l'un de ces bords, on shift notre réseau de la moitié de sa longueur
    //afin de ramener la composante connexe au centre du réseau
    //On suppose que l'on a pas de composante connexe qui fait plus de la moitiée de la taille du réseau, sinon la décaler de la moitié la ramène sur le côté :( 
    if ((*this).IsOnLeftEdge(ConCompNumber)){
        y = y - ny/2;
    }
    if ((*this).IsOnTopEdge(ConCompNumber)){
        x = x - nx/2;
    }

    while(data[k]!=ConCompNumber){
        k++;
    }

    //Initialisation des positions extrémales de la composante pour créer le nouveau réseau à la bonne taille
    SiteC s = (*this).site_index(k);
    int xmin = s._x - x, xmax = s._x - x, ymin = s._y - y, ymax = s._y - y;

    for(int xp=0; xp < nx; xp++){
        for(int yp=0; yp < ny; yp++){
           if((*this)[(*this).site_xy(x + xp, y + yp)] == ConCompNumber){
                if(xp < xmin){
                    xmin = xp;
                }
                if(xp > xmax){
                    xmax = xp;
                }
                if(yp < ymin){
                    ymin = yp;
                }
                if(yp > ymax){
                    ymax = yp;
                }
            }
        }
    }

    int shift = (xmin + 1)%2; //Nécessaire pour conserver la forme de la structure !
    ConComp Isolated = ConComp(xmax - xmin + 3 + shift, ymax - ymin + 3 + (ymax - ymin + 1)%2);//On rajoute de la marge sur les côtés comme ça on a de la place
    Isolated.NbrCC = ConCompNumber;

    //Recopiage de la forme dans le nouveau réseau isolé
    for(int xp = 0; xp < xmax - xmin + 1; xp++){
       for(int yp = 0; yp < ymax - ymin + 1; yp++){
          if ((*this)[(*this).site_xy(xmin + x + xp, ymin + yp + y)] == ConCompNumber){
                Isolated[Isolated.site_xy(shift + 1 + xp, 1+yp)] = (*this)[(*this).site_xy(xmin + x + xp, ymin + yp + y)];
            }
        }
    }

    return Isolated;
}

//Remplace les cases vides situées à l'extérieur d'une composante connexe isolée par des -1
void dfs_complementary(const SiteC s, const ConComp& Isolated, ConComp& Complementary, ConComp& visited){ 
    if (visited[s] == 1 || Isolated[s] != 0){
        return;
    }

    visited[s] = 1;
    Complementary[s] = -1;

    for(SiteC voisin : Isolated.voisins(s)){
        dfs_complementary(voisin, Isolated, Complementary, visited);
    }
}

//Algorithme basé sur Depth-first search algorithm afin de trouver de manière récursive les composantes connexes
//Voir rapport
void dfs(const SiteC s, ConComp& Complementary, ConComp& visited, const int CurrentLabel){ 
    if (visited[s] == 1 || Complementary[s] == -1){
        return;
    }

    visited[s] = 1;
    Complementary[s] = CurrentLabel;

    for(SiteC voisin : Complementary.voisins(s)){
        dfs(voisin, Complementary, visited, CurrentLabel);
    }
}

//Renvoie le complémentaire d'une composante connexe isolée, un site vide est représenté par un -1
//Les trous, autrefois représenté par des zéros sont désormais repérés avec le numéro de leur composante connexe
ConComp ConComp::Complementary() const{
    ConComp visited = ConComp(nx, ny);
    ConComp Complementary = ConComp(nx, ny);

    for(int k=0; k < nx*ny; k++){
        if(data[k] != 0){
            Complementary[Complementary.site_index(k)] = -1;
        }
    }
    dfs_complementary((*this).site_index(0), *this, Complementary, visited);

    int CurrentLabel = 1;

    for(int index = 0; index < nx*ny; index++){
        if (visited[visited.site_index(index)] == 0 && Complementary[Complementary.site_index(index)] != -1){
            dfs(Complementary.site_index(index), Complementary, visited, CurrentLabel);
            CurrentLabel++;
        }
    }

    Complementary.NbrCC = CurrentLabel;

    return Complementary;
}

//Renvoie la taille de tous les trous d'une composante connexe isolée
int ConComp::SizeOfHoles() const{
    Matrix SizeHoles = (*this).Complementary().SizeConComps();
    int Size = 0;

    for(int ConCompNumber = 0; ConCompNumber < SizeHoles.nx; ConCompNumber++){
        Size+= SizeHoles(ConCompNumber, 0);
    }

    return Size;
}

//Renvoie les caractéristiques de chaque composante connexe
//Size, SizeHoles, Volume, Porosity, Surface to volume ratio, Sphericity
//Parameters(CCN-1,i) -> param i de la CCN
Matrix ConComp::ClustersParameters() const{
    Matrix Size = (*this).SizeConComps();//Size(CC-1,0)
    Matrix Parameters = Matrix(NbrCC - 1, 6);

    for(int CCN = 1; CCN < NbrCC; CCN++){
        ConComp Isolated = (*this).isolateConComp(CCN);
        int OuterSurface = Isolated.OuterBorderLength(Isolated.NbrCC);

        Parameters(CCN - 1, 0) = Size(CCN -1, 0);
        Parameters(CCN - 1, 1) = Isolated.SizeOfHoles();
        Parameters(CCN - 1, 2) = Parameters(CCN - 1, 0) + Parameters(CCN - 1, 1);
        Parameters(CCN - 1, 3) = Parameters(CCN - 1, 1) / ((double) Parameters(CCN - 1, 2));
        Parameters(CCN - 1, 4) = OuterSurface / ((double) 2*Parameters(CCN - 1, 2));    
        Parameters(CCN - 1, 5) = sqrt(Parameters(CCN - 1, 2) * (3 * sqrt(3)) * (2*M_PI)) / OuterSurface;
    }

    return Parameters;
}


//_____________________________Constructeurs

//Algorithme basé sur Depth-first search algorithm afin de trouver de manière récursive les composantes connexes
//Voir rapport
void dfs(const Site s, const Lattice& L, int* data, Lattice& visited, const int CurrentLabel){ 
    if (visited[s] == 1 || L[s] == 0){
        return;
    }

    visited[s] = 1;
    data[s._index] = CurrentLabel;

    for(Site voisin : L.voisins(s)){
        dfs(voisin, L, data, visited, CurrentLabel);
    }
}

//Assigne à chaque particule située sur un Lattice sa composante connexe
ConComp::ConComp(Lattice& L) : nx(L.nx), ny(L.ny), data(nullptr){
    data = new int[nx*ny];
    Lattice visited = Lattice(nx,ny);

    //Pour retirer les cas chiants quand on calcule les outer surface des composantes connexes
    //On va pas retirer beaucoup de particules anyway ~ 0.1 * ny << NbrParts
    for(int y=0; y<ny; y++){
        L[L.site_xy(nx-1, y)] = 0;
    }

    int CurrentLabel = 1;

    for(int index = 0; index < nx*ny; index++){
        if (visited[visited.site_index(index)] == 0 && L[L.site_index(index)] != 0){
            dfs(L.site_index(index), L, data, visited, CurrentLabel);
            CurrentLabel++;
        }
        else if(L[L.site_index(index)] == 0){
            data[index] = 0;
        }
    }

    NbrCC = CurrentLabel;
}

ConComp::ConComp(int nx_, int ny_) : nx(nx_), ny(ny_), data(nullptr){
    data= new int[nx*ny];
    
    if(ny%2==1){
        throw std::runtime_error("ny has to be even to enforce BVK periodic conditions");
    }
    for(int i=0;i<nx*ny;i++){
        data[i]=0;
    }
}

ConComp::~ConComp(){
    if (data != nullptr){
        delete[] data;
    }
    data = nullptr;
}