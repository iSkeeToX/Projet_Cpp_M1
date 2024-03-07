# arguments pour l'étape de compilation : tous les avertissements
CPPFLAGS := -std=c++2b -Wall -Wextra -Ofast  #ou 23 #2b

# arguments pour l'édition de liens, bibliothèques externes
LDFLAGS := -lm -lsfml-graphics -lsfml-system -lsfml-window -flto

# nom du compilateur
CPP := g++

# cible par défaut
all: monprogramme

# cible monprogramme : compilation des .cpp (dépendances .o), édition des liens, puis → exécution

monprogramme: main.o Matrix.o TriangularLattice.o IsingModel.o
	$(CPP) -o prog $^ $(LDFLAGS)
	./prog

# compilation de main.cpp
main.o: main.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)

# compilation de Matrix.cpp
Matrix.o: Matrix.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)

# compilation de TriangularLattice.cpp
TriangularLattice.o: TriangularLattice.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)

# compilation de IsingModel.cpp
IsingModel.o: IsingModel.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)
