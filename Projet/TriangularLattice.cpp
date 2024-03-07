#include "TriangularLattice.h"

#include <iostream>
#include <random>
#include <stdexcept>
#include <algorithm> //min,max
#include <fstream>

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

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

void Lattice::affiche_SFML(sf::RenderWindow& window, float a) const{
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

void Lattice::write(std::string Name) const{
    std::ofstream fich(Name);

    for(int k=0;k<nx*ny;k++){
        if(k > 0 && k%ny == 0){
            fich << "\n";
        }
        fich << (int) data[k] << " ";
    }
}

Site Lattice::site_xy(int x, int y) const{
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

//index=ny*x+y
Site Lattice::site_index(int index) const{
    int y=index%ny;
    if(y < 0){
        return Site(index,index/ny,y+ny);
    }
    else{
        return Site(index,index/ny,y);
    }
}

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

int Lattice::Particule_Count() const{
    int count =0;

    for(int k=0;k<nx*ny;k++){
        if(data[k]!=0){
            count+=1;
        }
    }
    return count;
}

//Renvoie une matrix N*2 contenant les coordonées des N particules
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


//_____________________________Constructeurs

Lattice::Lattice(int nx_, int ny_) : nx(nx_), ny(ny_), data(nullptr){
    data= new int8_t[nx*ny];
    
    if(ny%2==1){
        throw std::runtime_error("ny doit être pair pour respecter les conditions de BVK");
    }
    for(int i=0;i<nx*ny;i++){
        data[i]=0;
    }
}

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

ConComp::ConComp(const Lattice& L) : nx(L.nx), ny(L.ny), data(nullptr){
    data = new int[nx*ny];
    Lattice visited = Lattice(nx,ny);
    Lattice Connected = Lattice(nx,ny);

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

//Destructeur
Lattice::~Lattice(){
    if (data != nullptr){
        delete[] data;
    }
    data = nullptr;
}

ConComp::~ConComp(){
    if (data != nullptr){
        delete[] data;
    }
    data = nullptr;
}
//_____________________________Privée

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

void hexagonOneColor(sf::RenderWindow& window, float viridis, float a, float pos_x, float pos_y, int number){

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
    
        // Draw the number in the middle of the hexagon
    sf::Text text;
    sf::Font font;

    // Load a font file (Arial) - replace with the path to your font file
    if (!font.loadFromFile("Arial.ttf")) {
        // Handle font loading error
        return;
    }

    text.setFont(font);
    text.setString(std::to_string(number));
    text.setCharacterSize(3*a/4);
    text.setFillColor(sf::Color::Black);

    // Center the text in the middle of the hexagon
    sf::FloatRect textBounds = text.getLocalBounds();
    text.setPosition(pos_x - textBounds.width / 2.0f, pos_y - textBounds.height / 1.0f);

    window.draw(text);
}

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

void ConComp::write(const std::string Name) const{
    std::ofstream fich(Name);

    for(int k=0;k<nx*ny;k++){
        if(k > 0 && k%ny == 0){
            fich << "\n";
        }
        fich << data[k] << " ";
    }
}