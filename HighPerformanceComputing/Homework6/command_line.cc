#include <iostream>  // for <<, >>, cout, cin, etc.
#include <cstdlib>   // for atoi, atol, atof, etc.
#include <cstdio>   // for EOF
#include <unistd.h>  // needed for getopt()


// -----------------------------------------------------------
// PARSELINE - parse the argument line

int parseline(int argc, char **argv, bool &a_set, int &count)
{
   const char *optlist = "an:";
   int opt;
   int err = 0;

   // assign some default values, either here or in the caller
   a_set = false;
   count = 0;

   // parse the argument line
   while((opt = getopt(argc, argv, optlist)) != EOF)  {
      switch(opt)  {
      case 'a':
         a_set = true;
         break;
      case 'n':
         count = atoi(optarg);
         break;
      default:  // getopt prints message
         err = 1;
         break;
      }
   }
   return(err);
}
// -----------------------------------------------------------

int main(int argc, char **argv)
{
   int count;
   bool have_a;
   using std::cout;
   using std::endl;

   if(parseline(argc, argv, have_a, count))  // error in arguments
      exit(1);   
   if(have_a)
      cout << "-a is specified" << endl;

   cout << "count is: " << count << endl;

   if(optind < argc)     // remaining argument is a file name
      cout << "file name argument is " << argv[optind] << endl;
   else
      cout << "no file name specified" << endl;
   return(0);
}
