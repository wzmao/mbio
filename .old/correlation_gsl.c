#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#define m(i,j) m[i*l+j]


double* calcMys(char m[],int n,int l)
{
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
  double *mys;
  int n1=0,n2=0;
  mys=malloc(l*l*sizeof(double));
  for (i=0;i<l*l;i++)
    mys[i]=0.0;
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
          allnumber+=2;
        }
      }

      gsl_matrix_view m = gsl_matrix_view_array (prop, 400, 400);
      gsl_vector *eval = gsl_vector_alloc (400);
      gsl_matrix *evec = gsl_matrix_alloc (400, 400);
      gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (400);
      gsl_eigen_symmv (&m.matrix, eval, evec, w);
      gsl_eigen_symmv_free (w);
      gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
      gsl_matrix *mid = gsl_matrix_alloc (400, 400);
      gsl_matrix *right = gsl_matrix_alloc (400, 400);
      for (k1 = 0; k1 < 400; k1++)
      {
        double ser=gsl_vector_get(eval,k1);
        for (k2=0;k2<400;k2++)
        {
          gsl_matrix_set(mid,k1,k2,0);
          gsl_matrix_set(right,k1,k2,gsl_matrix_get(evec,k2,k1));
        }
        gsl_matrix_set(mid,k1,k1,ser/(ser+1));
      }
      gsl_matrix *temp= gsl_matrix_alloc (400, 400);
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, evec, mid,
                  0.0, temp);
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, right,
                  0.0, evec);
      gsl_matrix_free(mid);
      gsl_matrix_free(right);
      gsl_vector_free (eval);
      for (k1=0;k1<400;k1++)
        for(k2=0;k2<400;k2++)
        {
          prop[k1*400+k2]=gsl_matrix_get(evec,k1,k2);
          if (prop[k1*400+k2]<0)
            prop[k1*400+k2]=-prop[k1*400+k2];
        }
      gsl_matrix_free (evec);
      gsl_matrix_free (temp);
      for (k1=0;k1<20;k1++)
        for(k2=0;k2<20;k2++)
          for(k3=0;k3<20;k3++)
            for(k4=0;k4<20;k4++)
            {
              pro1[k1*20+k2]+=prop[(k1*20+k3)*400+(k2*20+k4)];
              pro2[k3*20+k4]+=prop[(k1*20+k3)*400+(k2*20+k4)];
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
      mys[i*l+j]=allnumber;
      mys[j*l+i]=allnumber;
      free(prop);
      free(pro1);
      free(pro2);
    } 
  }
  free(edge);
  return mys;
}