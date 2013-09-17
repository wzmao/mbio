#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define m(i,j) m[i*l+j]
#define om(i,j) om[i*l+j]

double* calcOMES(char m[],int n,int l)
{
/*
This function is used to calculate the OMES matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.

example:
*/

  int i,j,k,k1,k2;
  char reslist[]={'-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
      'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
  int resnum=sizeof(reslist)/sizeof(char);
  double p[l][resnum];
  double *om;
  om = malloc( l*l * sizeof(double) );
  for (i=0;i<l;i++)
    for (j=0;j<resnum;j++)
      p[i][j]=0;
  for (i=0;i<l;i++)
  {
    double count=0;
    for (j=0;j<n;j++)
      for (k=0;k<resnum;k++)
        if (m(j,i)==reslist[k])
        {
          p[i][k]++;
          count++;
          break;
        }
    for (k=0;k<resnum;k++)
      p[i][k]=p[i][k]/n;
  }
  for (i=0;i<l;i++)
  {
    om(i,i)=0;
    for (j=i+1;j<l;j++)
    {
      double add=0;
      for (k1=0;k1<resnum;k1++)
        for (k2=0;k2<resnum;k2++)
        {
          if ((p[i][k1]!=0.0)&&(p[j][k2]!=0.0))
          {
            double count=0;
            for (k=0;k<n;k++)
            {
              if ((m(k,i)==reslist[k1])&&(m(k,j)==reslist[k2]))
                count++;
            }
            if (count!=0.0)
            {
              count=count/n;
              add+=(count-p[i][k1]*p[j][k2])*(count-p[i][k1]*p[j][k2])/(p[i][k1]*p[j][k2]);
            }
          }
        }
      om(i,j)=add;
      om(j,i)=add;
    }
  }
  return om;
}


double* calcOMES_line(char m[],int n,int l,int pos)
{
/*
This function is used to calculate the OMES matrix for only one line.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
pos is the return column.
*/
  int i,j,k,k1,k2;
  char reslist[]={'-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
      'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
  int resnum=sizeof(reslist)/sizeof(char);
  double p[l][resnum];
  double *om;
  om = malloc( l * sizeof(double) );
  for (i=0;i<l;i++)
    for (j=0;j<resnum;j++)
      p[i][j]=0;
  for (i=0;i<l;i++)
  {
    double count=0;
    for (j=0;j<n;j++)
      for (k=0;k<resnum;k++)
        if (m(j,i)==reslist[k])
        {
          p[i][k]++;
          count++;
          break;
        }
    for (k=0;k<resnum;k++)
      p[i][k]=p[i][k]/n;
  }
  for (i=pos;i<pos+1;i++)
  {
    om[i]=0;
    for (j=0;j<l;j++)
    {
      double add=0;
      for (k1=0;k1<resnum;k1++)
        for (k2=0;k2<resnum;k2++)
        {
          if ((p[i][k1]!=0.0)&&(p[j][k2]!=0.0))
          {
            double count=0;
            for (k=0;k<n;k++)
            {
              if ((m(k,i)==reslist[k1])&&(m(k,j)==reslist[k2]))
                count++;
            }
            if (count!=0.0)
            {
              count=count/n;
              add+=(count-p[i][k1]*p[j][k2])*(count-p[i][k1]*p[j][k2])/(p[i][k1]*p[j][k2]);
            }
          }
        }
      om[j]=add;
    }
  }
  return om;
}

int* calcOMES_P(char m[],int n,int l,int pos,int times, int cutoff)
{
/*
This function is used to calculate P value for the OMES matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
pos is the position to shuffle.
times is the shuffle times.
cutoff is the lower OMES cutoff for shuffle.
*/
  srand((int)time(0));
  int i, j,k,time;
  char s[n*l];
  for (i=0;i<n*l;i++)
  {
    s[i]=m[i];
  }
  int *p;
  p=malloc( l * sizeof(int) );
  double *saveom;
  saveom=calcOMES_line(m,n,l,pos);
  for (i=0;i<l;i++)
  {
    p[i]=0;
  }
  for (time=0;time<times;time++)
  {
    //shuffle
		char ch= ' ';
		for (i=0;i<n-1;i++)
		{
			j = i + rand()%(n-1-i);
			ch=s[j*l+pos];
			s[j*l+pos]=s[i*l+pos];
			s[i*l+pos]=ch;
		}
    //Calcate
    double *newom;
    newom=calcOMES_line(s,n,l,pos);
    for (i=0;i<l;i++)
    {
      if (newom[i]>saveom[i])
      {
        p[i]++;
      }
    }
    free(newom);
  }
  free(saveom);
  return p;
}
