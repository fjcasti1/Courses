#!/usr/bin/awk -f

# We assume the first and second field of the data file are the name and
# last name, respectively.

/^#/{print $1,"\t\t",$3,"\tBMI"}  # Display desired labels.
!/^#/{
  split($5,ht,/[-]/)
  inches = ht[1]*12+ht[2]
  print $1,"\t\t",inches
  }
