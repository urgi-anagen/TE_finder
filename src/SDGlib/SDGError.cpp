
#include<SDGError.h>
#include<unistd.h>
#include <iostream>
#include<stdlib.h>

void SDGError(int code)
{ 
  std::cerr << "Erreur #" << code << " lors de l'execution du programme\n";
  exit(code);
}
