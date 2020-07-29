#!/usr/bin/awk -f

# We assume the first and second field of the data file are the name and
# last name, respectively.

/^#/{print $1,"\t\t",$3,"\tBMI"}  # Display desired labels.
!/^#/{
  split($5,ht,/[-]/)
  inches = ht[1]*12+ht[2]
  h = inches*0.0254
  w = $6/2.2
  bmi = w/h^2
  print $1,$2,"\t\t",$4,"\t",bmi
  }
