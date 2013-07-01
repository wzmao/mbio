#include <stdio.h>

void qs(double a[],int l, int h)
{
/*
It is a quick sorting function
it could qs the array a from l to h.

example:
qs.argtypes=[POINTER(c_double),c_int,c_int]
qs.restype=c_void_p
a=range(100000)
from random import *
shuffle(a)
b=(c_double * len(a))()
for i in range(len(a)):
    b[i]=a[i]
qs(b,0,len(b)-1)
print list(b)==range(100000)
*/

  if (l >= h) 
  {
    return;
  }
  int i, j;
  double key;
  i = l;
  j = h;
  key = a[i];
  while (i < j) 
  {
    while (i < j && a[j] >= key) 
    {
      j--;
    }
    while (i < j && a[i] < key) 
    {
      i++;
    }
    if (i < j) 
    {
      double f = a[j];
      a[j] = a[i];
      a[i] = f;
    }
  }
  // printf("%d %d %d %d\n",l,h,i,j);
  if (l < j ) 
  {
    qs(a, l, j );
  }
  if (j+1 < h) 
  {
    qs(a, j+1 , h);
  }
  return;
}
