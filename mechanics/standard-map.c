/* -*-C-*-

$Id: copyright.c,v 1.5 2005/09/25 01:28:09 cph Exp $

Copyright 2005 Massachusetts Institute of Technology

This file is part of MIT/GNU Scheme.

MIT/GNU Scheme is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

MIT/GNU Scheme is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with MIT/GNU Scheme; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

*/
#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846
#define TWOPI 6.283185307179586476925287

main(argc, argv)
     int argc;
     char **argv;
{
  double x, y, K, yp;
  int i, n;
  double pv();

  if(argc != 5) {
    printf("args: K x y n\n");
    exit(-1);
  }
  sscanf(argv[1], "%lf", &K);
  sscanf(argv[2], "%lf", &x);
  sscanf(argv[3], "%lf", &y);
  sscanf(argv[4], "%d", &n);

  for(i=0; i<n; i++) {
    yp = pv(y + K*sin(x));
    x = pv(x + yp);
    y=yp;
    printf("%.16lf %.16lf\n", x, y);
  }
}

double pv(x)
     double x;
{
  double y;

  y = x - TWOPI*floor(x/TWOPI);
  return(y);
}
