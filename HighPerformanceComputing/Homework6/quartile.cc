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
void extremes(vector<double> v,double &min,double &max)
{
  min=v[0];
  max=v[0];
  for(int i=1;i<v.size();i++)
  {
    if(v[i]<min)
      min=v[i];
    if(v[i]>max)
      max=v[i];
  }
}

double median(vector<double> v,int begin,int end)
{
  int N=end-begin+1;
  double median;
  if(N % 2 == 0)
    median=(v[begin+N/2]+v[begin+N/2-1])/2;
  else
    median=v[begin+N/2];
  return(median);
}

void quartiles(vector<double> v,double &Q1,double &Q2, double &Q3)
{
  int index,N;
  N=v.size();
  if(N % 2 == 0)
  {
    Q1=median(v,0,N/2-1);
    Q2=median(v,0,N-1);
    Q3=median(v,N/2,N-1);
  }
  else
  {
    Q1=median(v,0,N/2-1);
    Q2=median(v,0,N-1);
    Q3=median(v,N/2+1,N-1);
  }
}

int process_file(istream& in)
{
   int retcode = EXIT_SUCCESS;  // unless an error occurs during processing
   vector<double> v;
   double value,vmin,vmax,Q1,Q2,Q3;

   while(in >> value){
     v.push_back(value);  // add to the end of the existing vector
   }
   sort(v.begin(),v.end());
   extremes(v,vmin,vmax);
   quartiles(v,Q1,Q2,Q3);
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
