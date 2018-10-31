// Skeleton C++ program for the quartile utility.

#include <iostream>   // needed for cin, cout
#include <fstream>    // needed for ifstream
#include <algorithm>  // needed for sort
#include <vector>     // needed for vector
#include <unistd.h>   // for getopt

using namespace std;

// -----------------------------------------------------------
//  PROCESS_FILE - read the data, sort it, and compute the quartiles
//  You may add other input arguments as necessary to adjust the behavior.

int process_file(istream& in)
{
   int retcode = EXIT_SUCCESS;  // unless an error occurs during processing
   vector<double> v;
   double value;

   while(in >> value)
     v.push_back(value);  // add to the end of the existing vector

   //  Sort the vector if necessary and compute the requested statistics.
   return(retcode);
}
// -----------------------------------------------------------
int parseline(int argc, char **argv)  // other arguments as needed
{
   // your code here
   return(0);
}
// -----------------------------------------------------------

int main(int argc, char **argv)
{
   //  RETCODE lets shell scripts check whether the program has completed
   //  successfully.
   int retcode;

   if(parseline(argc, argv))  // error in arguments
      retcode = EXIT_FAILURE;
   if(optind < argc) {  // read from the file named on the command line
      ifstream infile(argv[optind]);
      if(infile.good())  // file successfully opened
         retcode = process_file(infile);
      else
         retcode = EXIT_FAILURE;
   } else  // read from stdin
      retcode = process_file(cin);

   return(retcode);
}
