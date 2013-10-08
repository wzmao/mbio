#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define m(i,j) m[i*l+j]
#define mi(i,j) mi[i*l+j]
#define om(i,j) om[i*l+j]
#define sca(i,j) sca[i*l+j]

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
    allsumtemp=0;
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
              add+=n*(count-p[i][k1]*p[j][k2])*(count-p[i][k1]*p[j][k2])/(p[i][k1]*p[j][k2]);
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
              add+=n*(count-p[i][k1]*p[j][k2])*(count-p[i][k1]*p[j][k2])/(p[i][k1]*p[j][k2]);
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

double* calcSCA(char m[],int n,int l)
{
/*
This function is used to calculate the OMES matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.

example:
*/

  int i,j,k,k1,k2;
  char reslist[]={ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
      'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
  double q[]={.073, .025, .050, .061, .042, .072, .023, .053, .064, .089,
         .023, .043, .052, .040, .052, .073, .056, .063, .013, .033};
  int resnum=sizeof(reslist)/sizeof(char);
  double *p;
  p= malloc(l*resnum*sizeof(double));
  double *sca;
  sca = malloc( l*l * sizeof(double) );
  for (i=0;i<l;i++)
    for (j=0;j<resnum;j++)
      p[i*resnum+j]=0;
  for (i=0;i<l;i++)
  {
    double count=0;
    for (j=0;j<n;j++)
      for (k=0;k<resnum;k++)
        if (m(j,i)==reslist[k])
        {
          p[i*resnum+k]++;
          count++;
          break;
        }
    double sumall=0.0,temp=0.0;
    sumall=0;
    for (k=0;k<resnum;k++)
    {
      p[i*resnum+k]=p[i*resnum+k]*1.0/n;
      if (p[i*resnum+k]>=1.0-1e-10 || p[i*resnum+k]<=1e-10)
      {
        p[i*resnum+k]=0.0;
        temp=0;
      }
      else
      {
        temp=fabs(log((1-q[k])/q[k]*p[i*resnum+k]/(1-p[i*resnum+k])))*p[i*resnum+k];
        p[i*resnum+k]=temp*fabs(log((1-q[k])/q[k]*p[i*resnum+k]/(1-p[i*resnum+k])));
      }
        sumall=sumall+temp*temp;
    }
    sumall=sqrt(sumall);
    if (sumall<1e-10)
      for(k=0;k<resnum;k++)
        p[i*resnum+k]=0;
    else
      for (k=0;k<resnum;k++)
      {
        p[i*resnum+k]=p[i*resnum+k]/sumall;
      }
  }
  double *x;
  x=malloc(n*l*sizeof(double));
  double *sumx;
  sumx=malloc(l*sizeof(double));
  for (j=0;j<l;j++)
    sumx[j]=0.0;
  for (i=0;i<n;i++)
  {
    for (j=0;j<l;j++)
    {
      for (k=0;k<resnum;k++)
      {
        if (reslist[k]==m(i,j))
        {
          x[i*l+j]=p[j*resnum+k];
          sumx[j]+=x[i*l+j];
        }
      }
    }
  }
  for (i=0;i<l;i++)
  {
    for (j=i;j<l;j++)
    {
      double add=0;
      for (k=0;k<n;k++)
      {
        add+=x[k*l+i]*x[k*l+j];
      }
      add=(add/n-sumx[i]*sumx[j]/n/n);
      if (add<0)
        add=-add;
      sca(i,j)=add;
      sca(j,i)=add;
    }
  }
  free(p);
  free(x);
  free(sumx);
  return sca;
}


double* calcMy(char m[],int n,int l)
{
/*
This function is used to calculate the My matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
*/

  int i,j,k,k1,k2,k3,k4;
  char reslist[]={ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
      'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
  int resnum=sizeof(reslist)/sizeof(char);
  double *simi;
  simi=malloc(n*n*sizeof(double));
  for (i=0;i<n;i++)
  {
    simi[i*n+i]=0.0;
    for (j=i+1;j<n;j++)
    {
      double allcount=0,same=0;
      for (k=0;k<l;k++)
      {
        if (m(i,k)!='-' || m(j,k)!='-')
        {
          allcount++;
          if (m(i,k)==m(j,k))
            same++;
        }
      }
      simi[i*n+j]=same/allcount;
      simi[j*n+i]=simi[i*n+j];
    }
  }
  int *edge,edgenow=0;
  edge=malloc(5*n*2*sizeof(int));
  for (i=0;i<5*n*2;i++)
  {
    edge[i]=-1;
  }
  for (i=0;i<n;i++)
  {
    int has=0,now=0,maxn[5],max5=-1;
    while (has<5)
    {
      if (now!=i)
      {
        maxn[has]=now;
        has++;
      }
      now++;
    }
    for(j=0;j<5;j++)
      if (max5==-1||simi[i*n+maxn[max5]]>simi[i*n+maxn[j]])
        max5=j;
    for (j=now;j<n;j++)
    {
      if (simi[i*n+j]>simi[i*n+maxn[max5]])
      {
        maxn[max5]=j;
        max5=-1;
        for(k=0;k<5;k++)
          if (max5==-1||simi[i*n+maxn[max5]]>simi[i*n+maxn[k]])
            max5=k;
      }
    }
    for (j=0;j<5;j++)
    {
      edge[edgenow*2+0]=i;
      edge[edgenow*2+1]=maxn[j];
      edgenow++;
    }
  }
  free(simi);
  double *my;
  int n1=0,n2=0;
  my=malloc(l*l*sizeof(double));
  for (i=0;i<l*l;i++)
    my[i]=0.0;
  for (i=0;i<l;i++)
  {
    for (j=i+1;j<l;j++)
    {
      double *prop,*pro1,*pro2,allnumber=0,add=0;
      prop=malloc(20*20*20*20*sizeof(double));
      pro1=malloc(20*20*sizeof(double));
      pro2=malloc(20*20*sizeof(double));
      int count1=0,count2=0,count3=0,count4=0;
      for (k=0;k<160000;k++)
        prop[k]=0.0;
      for (k=0;k<400;k++)
      {
        pro1[k]=0;
        pro2[k]=0;
      }

      for (k=0;k<5*n;k++)
      {
        n1=edge[k*2];
        n2=edge[k*2+1];
        if (m(n1,i)!='-' && m(n1,j)!='-' && m(n2,i)!='-' && m(n2,j)!='-')
        {
          count1=0;
          count2=0;
          count3=0;
          count4=0;
          for (k1 = 0; k1 < 20; k1++)
          {
            if (reslist[k1]==m(n1,i))
            {
              count1=k1;
              break;
            }
          }
          for (k1 = 0; k1 < 20; k1++)
          {
            if (reslist[k1]==m(n2,i))
            {
              count2=k1;
              break;
            }
          }
          for (k1 = 0; k1 < 20; k1++)
          {
            if (reslist[k1]==m(n1,j))
            {
              count3=k1;
              break;
            }
          }
          for (k1 = 0; k1 < 20; k1++)
          {
            if (reslist[k1]==m(n2,j))
            {
              count4=k1;
              break;
            }
          }
          prop[(count1*20+count3)*400+(count2*20+count4)]+=1;
          prop[(count2*20+count4)*400+(count1*20+count3)]+=1;
          pro1[count1*20+count2]+=1;
          pro1[count2*20+count1]+=1;
          pro2[count3*20+count4]+=1;
          pro2[count4*20+count3]+=1;
          allnumber+=2;
        }
      }
      if (allnumber!=0)
      {
        for (k=0;k<160000;k++)
          prop[k]=prop[k]/allnumber;
        for (k=0;k<400;k++)
        {
          pro1[k]=pro1[k]/allnumber;
          pro2[k]=pro2[k]/allnumber;
        }
      }
      allnumber=0.0;
      for (k1=0;k1<20;k1++)
        for (k2=0;k2<20;k2++)
        {
          if (pro1[k1*20+k2]==0)
          {
            continue;
          }
          else
          {
            for (k3=0;k3<20;k3++)
              for (k4=0;k4<20;k4++)
              {
                if (pro2[k3*20+k4]==0||prop[(k1*20+k3)*400+(k2*20+k4)]==0)
                {
                  continue;
                }
                else
                {
                  allnumber+=prop[(k1*20+k3)*400+(k2*20+k4)]*log(prop[(k1*20+k3)*400+(k2*20+k4)]/pro1[k1*20+k2]/pro2[k3*20+k4]);                  
                }
              }
          }
        }
      my[i*l+j]=allnumber;
      my[j*l+i]=allnumber;
      free(prop);
      free(pro1);
      free(pro2);
    } 
  }
  free(edge);
  return my;
}

double* calcMyp(char m[],int n,int l)
{
/*
This function is used to calculate the Myp matrix.
m is the fastas sequences which has been concanated to one array
n is the number of sequences and l is the length.
len(m) must eaqul to l*n.
*/

  double *myp=calcMy(m,n,l);
  int i,j,k;
  double mean[l],allmean=0;
  k=0;
  for (i=0;i<l;i++)
  {
    mean[i]=0;
    for (j=0;j<l;j++)
    {
      mean[i]+=myp[k];
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
      myp[i*l+j]-=((mean[i]*mean[j])/allmean);
    }
  }
  return myp;
}