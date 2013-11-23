#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "../Sequence/correlation.c"
#include "../Application/sort.c"

typedef int bool;

#define OUTPUT
#define TRUE 1
#define FALSE 0
#define True 1
#define False 0
#define true 1
#define false 0
#define seqnum
#define lennum
#define cutoff 
#define times 

char reslist[]={'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'};
int resnum=sizeof(reslist);

int main (int argc, char *argv[])
{
  MPI_Init(NULL,NULL);
  int id,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  srand((int)time(0));
  if ((id==0)&&  OUTPUT)
  {
    printf("* There are all %d process.\n",size);
  }
  int i,j,k;
  char lines[seqnum][lennum];
  for (i=0;i<seqnum;i++)
  {
    lines[i][0]='\0';
  }
  char ch;
  FILE *f;
  if (f=fopen("file.fasta","r"))
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if ((id==0)&&  OUTPUT)
      printf("* Reading FASTA\n");
    MPI_Barrier(MPI_COMM_WORLD);
    i=0;
    j=0;
    ch=fgetc(f);
    while (!feof(f))
    {
      if (ch=='\n')
      {
        lines[i][j]='\0';
        i++;
        j=0;
      }
      else
      {
        if ((ch!='\n')&&(ch<=127))
        {
          lines[i][j]=ch;
          j++;
        }
      }
      ch=fgetc(f);
    }
    fclose(f); 
  }
  else
  {
    printf("* Procdss %d couldn't get FASTA.\n",id);
    return 0;
  } 
  lines[i][j]='\0';
  int number=i+1,length=j;
  MPI_Barrier(MPI_COMM_WORLD);
  if (OUTPUT)
    printf("* Processe %d get %d sequences with length of %d\n",id,number,length);
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0;i<number;i++)
  {
    if (strlen(lines[i])!=length)
    {
      printf ("* There are something wrong with No.%d sequence length.\n",i);
      printf("%s\n",lines[i]);
    }
  }
  double mi[length][length];
  int p[length][length];

  char temp[number][2];
  MPI_Barrier(MPI_COMM_WORLD);
  if ((id==0)&&  OUTPUT)
    printf("* Calculating the frequency.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  double freq[number][resnum];
  for (i=0;i<length;i++)
  {
    long fre[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},all=0;
    for (j=0;j<number;j++)
    {
      for (k=0;k<resnum;k++)
      {
        if (lines[j][i]==reslist[k])
        {
          fre[k]++;
          all++;
          break;
        }
      }
    }
    for (k=0;k<resnum;k++)
    {
      freq[i][k]=1.0*fre[k]/all;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if ((id==0)&&  OUTPUT)
    printf("* Frequency calculation finished.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  if ((id==0)&&  OUTPUT)
    printf("* Calculating the MIp value.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  char llll[length*number];
  k=0;
  for (i=0;i<number;i++)
    for (j=0;j<length;j++)
    {
      llll[k]=lines[i][j];
      k++;
    }
  double *mitemp=calcMIp(llll,number,length);
  k=0;
  for (i=0;i<length;i++)
    for (j=0;j<length;j++)
    {
      mi[i][j]=mitemp[k];
      k++;
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if ((id==0)&&  OUTPUT)
    printf("* MIp calculation finished.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  if ((id==0)&&  OUTPUT)
    printf("* Saving MIp result.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  if (id==0)
  {
    if (f=fopen("MIpsave.save","wb"))
    {
     fwrite(&mi[0][0],sizeof(double),length*length,f);
    }
    fclose(f);
  }
  if ((id==0)&&  OUTPUT)
    printf("* MIp Saving finished.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  double save=0;
  long all=length*length,done=0,cccc=0;
  if (id==0)
  { 
    double sortmi[length*length];
    k=0;
    for (i=0;i<length;i++)
    {
      for (j=0;j<length;j++)
      {
        sortmi[k]=mi[i][j];
        k++;
      }
    }
    qs(sortmi,0,k-1);
    int rank=k*cutoff-1;
    save=sortmi[rank];
  }
  MPI_Bcast(&save,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (OUTPUT==1)  
    printf("* Processe %d calculating the P for each positions.\n",id);
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0;i<length;i++)
    for (j=0;j<length;j++)
    {
      p[i][j]=-1;
    }
  double *misave=calcMI(llll,number,length);
  for (i=id;i<length;i+=size)
  {
    int j,k;
    p[i][i]=-100;
    done++;
    if (OUTPUT==1)  
      {printf("\r* %ld positions have been done, %ld/%ld=%.4lf%% from %d",done,done,all,100.0*done/all,id);fflush(stdout);}
    int *tempp;
    tempp=calcMIp_P(llll,number,length,i,times,save,misave);
    j=0;
    while (j<length)
    {
      p[i][j]=tempp[j];
      j++;
    }
    done+=2*(length-i);
    free(tempp);
  }
  free(misave);
  MPI_Barrier(MPI_COMM_WORLD);
  if (id!=0)
    for (i=id;i<length;i+=size)
      MPI_Send(&p[i][0],length,MPI_INT,0,i,MPI_COMM_WORLD);
  else
  {
    MPI_Status status;
    for (i=0;i<length;i++)
      if (i%size!=0)
        MPI_Recv(&p[i][0],length,MPI_INT,(i%size),i,MPI_COMM_WORLD,&status);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if ((id==0)&&(OUTPUT==1))
    printf("\n* All P value calculations have been done.\n");
  if ((id==0)&&(OUTPUT==1))
    printf("* Saving P value.\n");
  if (id==0)
  {
    if (f=fopen("Psave.save","wb"))
     fwrite(&p[0][0],sizeof(int),length*length,f);
    fclose(f);
  }
  if ((id==0)&&(OUTPUT==1))
    printf("* P value saved.\n");
  MPI_Finalize();
  return 0;
}
