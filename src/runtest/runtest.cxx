#include <iostream>

int main( int argc, char **argv )
{
   std::cout << "test 1" << std::endl;

   double zstep = 5.0;

   for(int i = 0; i< 100; i++ ) {
      double xv = 100.0;
      double yv = 200.0;
      double zv = 100.0 + double(i)*zstep;
      std::cout << "(" << xv << ", " << yv<< ", " << zv<< ")" << std::endl;
   }
   return 0;
}

