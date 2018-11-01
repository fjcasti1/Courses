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
  if (N==1)
  {
    Q1=v[0];
    Q2=v[0];
    Q3=v[0];
  }
  else
  {
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
}

int process_file(istream& in, bool e_opt)
{
   int retcode = EXIT_SUCCESS;  // unless an error occurs during processing
   vector<double> v;
   double value,vmin,vmax,Q1,Q2,Q3;

   while(in >> value){
     v.push_back(value);  // add to the end of the existing vector
   }
   if(v.size()!=0)
   {
     sort(v.begin(),v.end());
     extremes(v,vmin,vmax);
     if (!e_opt)
       quartiles(v,Q1,Q2,Q3);

//     cout << "Min :" << vmin << endl;
//     cout << "Q1 : " << Q1 << endl;
//     cout << "Q2 : " << Q2 << endl;
//     cout << "Q3 : " << Q3 << endl;
//     cout << "Max :" << vmax << endl;
     cout << vmin << endl;
     if (!e_opt)
     {
       cout << Q1 << endl;
       cout << Q2 << endl;
       cout << Q3 << endl;
     }
     cout << vmax << endl;
   }
   //  Sort the vector if necessary and compute the requested statistics.
   return(retcode);
}
// -----------------------------------------------------------
int parseline(int argc, char **argv, bool &e_opt, bool &h_opt)  // other arguments as needed
{
  const char *optlist = "eh:";
  int opt;
  int err = 0;

  e_opt = false;
  h_opt = false;

  while((opt = getopt(argc, argv, optlist)) != EOF)
  {
    switch(opt)
    {
      case 'e':
        e_opt = true;
        break;
      case 'h':
        h_opt = true;
        break;
      default:
        err = 1;
        break;
    }
  }
  return(err);
}
// -----------------------------------------------------------

int main(int argc, char **argv)
{
   //  RETCODE lets shell scripts check whether the program has completed
   //  successfully.
   int retcode;
   bool e_opt, h_opt;

   if(parseline(argc, argv,e_opt,h_opt))  // error in arguments
      retcode = EXIT_FAILURE;
   if(h_opt)
   {
     //print help
     cout << "HELP MESSAGE" << endl;
   }
   else
   {
     if(optind < argc)   // read from the file named on the command line
     {
        ifstream infile(argv[optind]);
        if(infile.good())  // file successfully opened
           retcode = process_file(infile,e_opt);
        else
           retcode = EXIT_FAILURE;
     } 
     else  // read from stdin
        retcode = process_file(cin,e_opt);
   }

   return(retcode);
}
