#!/usr/bin/awk -f

# We assume the first and second field of the data file are the name and
# last name, respectively.
# We assume the height is given in a format H-h where H is in feet and h
# in inches.

/^#/{print $1,"\t\t",$4}  # Display desired labels.
!/^#/{
  split($5,ht,/[-]/)
  inches = ht[1]*12+ht[2]
  print $1,$2,"\t\t",inches
  }
