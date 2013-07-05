#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define m(i,j) m[i*l+j]
#define mi(i,j) mi[i*l+j]

double* calcMI(char m[],int n,int l)
{
/*
This function is used to calculate the MI matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.

example:
c_CalcMI = M.calcMI
c_CalcMI.argtypes = [ct.POINTER(ct.c_char), ct.c_int, ct.c_int]
c_CalcMI.restype = ct.POINTER(ct.c_double)
allsequence = ''.join(sequences)
m = (ct.c_char * len(allsequence))()
for i in range(len(allsequence)):
    m[i] = allsequence[i]
l = len(sequences[0])
result = c_CalcMI(m, len(sequences), l)
mi = []
for i in range(l**2):
    if i % l == 0:
        mi.append([])
    mi[-1].append(result[i])
*/

  int i,j,k,k1,k2;
  char reslist[]={'-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
      'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
  int resnum=sizeof(reslist)/sizeof(char);
  double p[l][resnum];
  double *mi;
  mi = malloc( l*l * sizeof(double) );
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
    mi(i,i)=0;
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
              add+=count*log(count/p[i][k1]/p[j][k2]);
            }
          }
        }
      mi(i,j)=add;
      mi(j,i)=add;
    }
  }
  return mi;
}


double* calcMI_line(char m[],int n,int l,int pos)
{
/*
This function is used to calculate the MI matrix for only one line.
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
  double *mi;
  mi = malloc( l * sizeof(double) );
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
    mi[i]=0;
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
              add+=count*log(count/p[i][k1]/p[j][k2]);
            }
          }
        }
      mi[j]=add;
    }
  }
  return mi;
}

int* calcMI_P(char m[],int n,int l,int pos,int times, int cutoff)
{
/*
This function is used to calculate P value for the MI matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
pos is the position to shuffle.
times is the shuffle times.
cutoff is the lower MI cutoff for shuffle.
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
  double *savemi;
  savemi=calcMI_line(m,n,l,pos);
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
    //Calculate
    double *newmi;
    newmi=calcMI_line(s,n,l,pos);
    for (i=0;i<l;i++)
    {
      if (newmi[i]>savemi[i])
      {
        p[i]++;
      }
    }
    free(newmi);
  }
  free(savemi);
  return p;
}

double* calcMIp(char m[],int n,int l)
{
/*
This function is used to calculate the MIp matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
*/
  double *mip=calcMI(m,n,l);
  int i,j,k;
  double mean[l],allmean=0;
  k=0;
  for (i=0;i<l;i++)
  {
    mean[i]=0;
    for (j=0;j<l;j++)
    {
      mean[i]+=mip[k];
      k++;
    }
    mean[i]=mean[i]/(l-1);
    allmean+=mean[i];
  }
  allmean=allmean/l;
  for (i=0;i<l;i++)
  {
    for (j=0;j<l;j++)
    {
      mip[i*l+j]-=((mean[i]*mean[j])/allmean);
    }
  }
  return mip;
}

double* calcMIp_line(char m[],int n,int l, int pos)
{
/*
This function is used to calculate the MIp matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
*/
  double *mi=calcMI(m,n,l);
  int i,j,k;
  double mean[l],allmean=0;
  k=0;
  for (i=0;i<l;i++)
  {
    mean[i]=0;
    for (j=0;j<l;j++)
    {
      mean[i]+=mi[k];
      k++;
    }
    mean[i]=mean[i]/(l-1);
    allmean+=mean[i];
  }
  allmean=allmean/l;
  double *mip;
  mip=malloc( l * sizeof(double) );
  for (j=0;j<l;j++)
  {
    mip[j]=mi[pos*l+j]-((mean[pos]*mean[j])/allmean);
  }
  return mip;
}

int* calcMIp_P(char m[],int n,int l,int pos,int times, int cutoff,double savemi[])
{
/*
This function is used to calculate P value for the MIp matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
pos is the position to shuffle.
times is the shuffle times.
cutoff is the lower MI cutoff for shuffle.
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
  for (i=0;i<l;i++)
  {
    p[i]=0;
  }
  double savemi1[l];
  double linesum[l],allsum=0.0,linesumtemp[l],allsumtemp=0.0;
  for (i=0;i<l;i++)
  {
    linesum[i]=0.0;
    linesumtemp[i]=0.0;
    for (j=0;j<l;j++)
    {
      linesum[i]+=savemi[i*l+j];
      allsum+=savemi[i*l+j];
    }
  }
  double *savemip;
  savemip = malloc(l* sizeof(double));
  for (i=0;i<l;i++)
  {
    savemip[i]=savemi[pos*l+i]-(linesum[pos]*linesum[i]/(l-1)*l/allsum);
    savemi1[i]=savemi[pos*l+i];
  }
  double tempmip=0.0;
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
    //Calculate
    double *newmi;
    newmi=calcMI_line(s,n,l,pos);
    for (i=0;i<l;i++)
    {
      linesumtemp[i]=linesum[i]-savemi1[i]+newmi[i];
      allsumtemp+=-savemi1[i]+newmi[i];
    }
    linesumtemp[pos]+=allsumtemp+savemi1[pos]-newmi[pos];
    allsumtemp+=allsum;
    for (i=0;i<l;i++)
    {
      tempmip=newmi[i]-(linesumtemp[pos]*linesumtemp[i]/(l-1)*l/allsumtemp);
      if (tempmip>savemip[i])
      {
        p[i]++;
      }
    }
    free(newmi);
  }
  free(savemip);
  return p;
}