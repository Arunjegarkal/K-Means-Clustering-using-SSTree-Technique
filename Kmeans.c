#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/**********************************************************************************************************
Function to find second cluster's centroid
**********************************************************************************************************/

void findCentroid_2(int *data1, int **data2, int dim,int size,int current_cluster_size);

/**********************************************************************************************************
Function to find i'th cluster's centroid
**********************************************************************************************************/

void findCentroid_i(int *data1, int **data2, int dim,int size,int current_cluster_size);

/**********************************************************************************************************
Function to find mean cluster's centroid
**********************************************************************************************************/

void findMidpoint(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int ** temp_centroid);

/**********************************************************************************************************
Function to create cluster using Kmeans algorithm
**********************************************************************************************************/

void KMEANS(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *searchData,int *clusterstart);

/**********************************************************************************************************
Function to assign data points to cluster
**********************************************************************************************************/

void ClusterAssign(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *temp_array);

/**********************************************************************************************************
Function to search a datapoint in a cluster
**********************************************************************************************************/

void Search(int dim,int ndata,int *data,int *searchData,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *clusterstart,int **cluster_assign_level1,int **cluster_centroid_level1,int *clustersize_level1,int *clusterstart_level1,double *cluster_raidus_level1,int temp_array[][2]);

/**********************************************************************************************************
Function to clear Cluster size, cluster assign and cluster radious before calculating new centroids
**********************************************************************************************************/

void initializeZero(int ndata,int k,int *clustersize,int **cluster_assign,double *cluster_raidus);

/**********************************************************************************************************
Sort the data array according to cluster assign
**********************************************************************************************************/

void sort(int dim,int ndata,int *data,int k,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *clusterstart);

void KMEANS(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *searchData,int *clusterstart)
{
    int i=0,y,j;
    int x=(rand()%ndata);
    y=x-(x%dim);
    int **temp_centroid;
    int *temp_array=(int*)malloc(ndata * sizeof(int));
    temp_centroid=(int **)malloc(k * sizeof(int *));
    for (i=0; i<k; i++)
         *(temp_centroid+i) = (int *)malloc(dim * sizeof(int));
    //Find cluster1 centroids
	for(i=0;i<dim;i++)
    {
        cluster_centroid[0][i]=y+i;
    }
	//Find centroid 2
    findCentroid_2(data,cluster_centroid,dim,ndata,1);
    //Find centroids of i clusters
	for(i=2;i<k;i++)
        findCentroid_i(data,cluster_centroid,dim,ndata,i);

    //ClusterAssign function call to assign all the data points to assign to recently calculated clusters
	ClusterAssign(dim,ndata,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_array);
	//find the mid point of all the clusters
	findMidpoint(dim,ndata,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_centroid);
    int flag=0;
	//Loop untill centroids of K clusters wont change
    do
    {
        flag=0;
        for(i=0;i<k;i++)
        {
            for(j=0;j<dim;j++)
            {
                if(cluster_centroid[i][j]!=temp_centroid[i][j])
                    flag=1;
            }
        }
        for(i=0;i<k;i++)
        {
            for(j=0;j<dim;j++)
            {
                cluster_centroid[i][j]=temp_centroid[i][j];
            }
        }

        initializeZero(ndata,k,clustersize,cluster_assign,cluster_raidus);
        ClusterAssign(dim,ndata,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_array);
        findMidpoint(dim,ndata,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_centroid);

    }while(flag==1);
    sort(dim,ndata,data,k,cluster_centroid,cluster_assign,clustersize,clusterstart);

}

void sort(int dim,int ndata,int *data,int k,int ** cluster_centroid, int **cluster_assign,int *clustersize,int * clusterstart)
{
    int i,j,l,count=0;
    int *temp_array = (int*)malloc(ndata * sizeof(int));
    /*repeat for K clusters*/
    for(i=0;i<k;i++)
    {
        clusterstart[i]=count;
        /*read the data points assigned to the ith cluster and move it to temp */
        for(j=0;j<clustersize[i];j++)
        {
            for(l=0;l<dim;l++)
            {
                temp_array[count]=data[cluster_assign[i][j]+l];
                count++;
            }
        }
    }
    /*Move the data from temp to data*/
    for(i=0;i<ndata;i++)
        data[i]=temp_array[i];
}


void Search(int dim,int ndata,int *data,int *searchData,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *clusterstart,int **cluster_assign_level1,int **cluster_centroid_level1,int *clustersize_level1,int *clusterstart_level1,double *cluster_raidus_level1,int temp_array[][2])
{
    int i,ii,j,jj,m,val1,val2,search_cluster[k],nearest_node,visited_nodes=0,search_cluster_size=0,startPos;
    double temp,dist,mindist=9999999999,search_cluster_dist[k];
    /*Find the cluster at level 1 in which search node fall in based on radius*/
    for(i=0;i<k;i++)
    {
        temp=0.0;
        for(j=0;j<dim;j++)
        {
            val1=cluster_centroid[i][j];
            val2=searchData[j];
            temp+=((val1-val2)*(val1-val2));
        }
        dist=sqrt(temp);
        //Store the cluster number if the search data point falls with in the radious of cluster
		if ((cluster_raidus[i]-dist)>0)
        {
            mindist=dist;
            search_cluster_dist[search_cluster_size]=mindist;
            search_cluster[search_cluster_size]=i;
            search_cluster_size++;
        }
    }
    mindist=99999999999;
    for(i=0;i<search_cluster_size;i++)
    {
        /*Find the cluster at level 2 in which search node fall in based on radius*/
        for(j=temp_array[search_cluster[i]][1];j<(temp_array[search_cluster[i]][1]+temp_array[search_cluster[i]][0]);j++)
        {
            /*Skip the search in cluster if cluster size is zero*/
            if(clustersize_level1[j]>0)
            {
                temp=0.0;
                for(ii=0;ii<dim;ii++)
                {
                    val1=cluster_centroid_level1[j][ii];
                    val2=searchData[ii];
                    temp+=((val1-val2)*(val1-val2));
                }
                dist=sqrt(temp);
                /*Perform exhaustive search in level2 sub cluster if the search node falls with in the radious*/
                if ((cluster_raidus_level1[j]-dist)>0)
                 {
                     for(jj=clusterstart_level1[j];jj<((clustersize_level1[j]*dim)+clusterstart_level1[j]);jj=jj+dim)
                     {
                        temp=0.0;
                        /*compute the distance between search node and the data node*/
                        for(m=0;m<dim;m++)
                        {
                            val1=data[jj+m];
                            val2=searchData[m];
                            temp+=(val1-val2)*(val1-val2);
                        }
                        dist=sqrt(temp);
                        if (mindist>dist)
                        {
                            mindist=dist;
                            nearest_node=jj;
                        }
                        visited_nodes++;
                     }
                 }
            }
        }
    }
    printf("\nSearch node is ");
    for(i=0;i<dim;i++)
    {
        printf("[%d] ",searchData[i]);
    }
    printf("\nNearest node is ");
    for(i=0;i<dim;i++)
    {
        printf("[%d] ",data[nearest_node+i]);
    }
    printf("Number of Nodes visited %d \nDistance is %f",visited_nodes,mindist);

}

void findMidpoint(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int ** temp_centroid)
{
    int i,j,l,val=0;

    for(i=0;i<k;i++)
    {
        for(j=0;j<dim;j++)
        {
            temp_centroid[i][j]=0;
        }
    }
	/*Repeat for K clusters*/
    for(i=0;i<k;i++)
    {
        for(l=0;l<dim;l++)
        {
            for(j=0;j<clustersize[i];j++)
            {
				/*Calculate Sum of all data points in a cluster*/
                val=val+data[cluster_assign[i][j]+l];
            }
            if(val!=0)
                temp_centroid[i][l]=(val/clustersize[i]);
            else
                temp_centroid[i][l]=0;
            val=0;
        }

    }
}
void ClusterAssign(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *temp_array)
{

    int i=0,j,m,val1,val2;
    double temp,dist;
     for(i=0;i<k;i++)
    {
        for(j=0;j<dim;j++)
        {
			cluster_raidus[i]=0;
        }
    }


    for(j=0;j<ndata;j=j+dim)
    {
        double min_dist=99999999999;
        for(m=0;m<k;m++)
        {
            /*Loop to calculate the distance*/
            temp=0;
            for(i=0;i<dim;i++)
            {
                val1=data[j+i];
                val2=cluster_centroid[m][i];
                temp+=(val1-val2)*(val1-val2);
            }
            dist=sqrt(temp);
			/*store the distance if the calculated distance is less than previously stored distance*/
            if(min_dist>dist)
            {
                min_dist=dist;
                temp_array[j]=m;
				/*Store the radius*/
                if(cluster_raidus[m]<min_dist)
                {
                    cluster_raidus[m]=min_dist;
                }
            }
        }
    }
    for(j=0;j<ndata;j=j+dim)
    {
        cluster_assign[temp_array[j]][clustersize[temp_array[j]]]=j;
        clustersize[temp_array[j]]=clustersize[temp_array[j]]+1;
    }
}
void initializeZero(int ndata,int k,int *clustersize,int **cluster_assign,double *cluster_raidus)
{
    int i,j;
    for(i=0;i<k;i++)
    {
        clustersize[i]=0;
        cluster_raidus[i]=0;
    }
    for(i=0;i<k;i++)
    {
        for(j=0;j<ndata;j++)
            cluster_assign[i][j]=0;
    }

}
/*Select the data point as the centroid for 2nd cluster which is far from cluster 1 centroid*/
void findCentroid_2(int *data1, int **data2, int dim,int size,int current_cluster_size)
{
    int i,j,*max_dist_pt;
    double dist,mindist=0,tmp=0.0;
    for(j=0;j<size;j=j+dim)
    {
        tmp=0.0;
        for(i=0;i<dim;i++)
        {
            int val=data1[i+j];
            int val2=data1[data2[0][i]];
            tmp+=((val-val2)*(val-val2));
        }
        dist=sqrt(tmp);
        if(mindist<dist)
        {
            mindist=dist;
            int m;
            for(m=0;m<dim;m++)
            {
                data2[current_cluster_size][m]=j+m;
            }
        }

    }
}

void findCentroid_i(int *data1, int **data2, int dim,int size,int current_cluster_size)
{
    int i,ii,j,*max_dist_pt;
    double dist,mindist=0,tmp=0.0;
    for(j=0;j<size;j=j+dim)
    {

        for(ii=0;ii<current_cluster_size;ii++)
        {
            tmp=0.0;
            for(i=0;i<dim;i++)
            {
                int val1=data1[j+i];
                int val2=data1[data2[ii][i]];
                tmp+=((val1-val2)*(val1-val2));
            }
            dist=sqrt(tmp);
            if(mindist<dist)
            {
                mindist=dist;
                int m,n,flag=1;
                for(n=0;n<current_cluster_size;n++)
                {
                    for(m=0;m<dim;m++)
                    {
                        if(data2[n][m]!=j+m)
                        {
                            flag=0;
                            m=dim;
                        }
                        else
                            flag=1;
                    }
                }
                if(flag==0)
                {
                    for(m=0;m<dim;m++)
                    {
                        data2[current_cluster_size][m]=j+m;
                    }
                }
            }
        }
    }
}

int main()
{
    int ndata=10000,i,dim=5,*data,k,j,l,size,n=2,level1_clustercount=0;
    int **cluster_assign,**cluster_centroid,*searchData,*clustersize,*clusterstart;
    double *cluster_raidus;
    int *dist;
    int *ptr_data,*ptr_search,temp=0;

    printf("Enter the Number of Data points\n");
    scanf("%d",&ndata);
    printf("Enter the Number of Dimention\n");
    scanf("%d",&dim);
    size=ndata*dim;
    k=(int)ceil(cbrt(ndata));
    //Dynamic memory allocation
    data = (int*)malloc(size * sizeof(int));
    clusterstart = (int*)malloc(k * sizeof(int));
    cluster_assign = (int**)malloc(k * sizeof(int));
    for (i=0; i<k; i++)
    {
        *(cluster_assign+i) = (int *)malloc(size * sizeof(int));
    }
    cluster_raidus =(double*)malloc(k * sizeof(double));
    searchData = (int*)malloc(dim * sizeof(int));
    cluster_centroid=(int **)malloc(k * sizeof(int *));
    for (i=0; i<k; i++)
         *(cluster_centroid+i) = (int *)malloc(dim * sizeof(int));
    clustersize = (int*)malloc(k * sizeof(int));
    for(i=0;i<k;i++)
    {
        clustersize[i]=0;
    }
    //Generating random dataset
    for(i=0;i<size;i++)
    {
        data[i]=((rand()%1000)+(rand() / (double)RAND_MAX));
    }
    //Initializing all the radius to 0
    for(i=0;i<k;i++)
    {
        cluster_raidus[i]=0;
    }
    for(i=0;i<dim;i++)
        searchData[i]=((rand()%1000)+(rand() / (double)RAND_MAX));

    KMEANS(dim,size,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,searchData,clusterstart);
    for(i=0;i<k;i++)
    {
        //printf("clustersize %d %d %d clusterstart %d \n",i,clustersize[i],(int)ceil(cbrt(clustersize[i])),clusterstart[i]);
        level1_clustercount=level1_clustercount + (int)ceil(cbrt(clustersize[i]));
        if(temp<clustersize[i])
            temp=clustersize[i];
    }
    //Data structures for level 2
    int array[level1_clustercount][2],*temp_data,count=0,**temp_array;
    int **cluster_assign_level1,**cluster_centroid_level1,*clustersize_level1,*clusterstart_level1;
    double *cluster_raidus_level1;
    //Dynamic memory allocation of level 2 data structures
    clusterstart_level1 = (int*)malloc(level1_clustercount * sizeof(int));
    cluster_assign_level1 = (int**)malloc(level1_clustercount * sizeof(int));
    for (i=0; i<level1_clustercount; i++)
    {
        *(cluster_assign_level1+i) = (int *)malloc(temp*dim * sizeof(int));
    }
    cluster_raidus_level1 =(double*)malloc(level1_clustercount * sizeof(double));
    cluster_centroid_level1=(int **)malloc(level1_clustercount * sizeof(int *));
    for (i=0; i<level1_clustercount; i++)
         *(cluster_centroid_level1+i) = (int *)malloc(dim * sizeof(int));
    clustersize_level1 = (int*)malloc(level1_clustercount * sizeof(int));

    int **temp_cluster_centroid_level1,*temp_clustersize_level1,*temp_clusterstart_level1,temp_start=0;
    double *temp_cluster_raidus_level1;
    temp_clustersize_level1 = (int*)malloc(level1_clustercount * sizeof(int));
    temp_cluster_centroid_level1=(int **)malloc(level1_clustercount * sizeof(int *));
    for (i=0; i<level1_clustercount; i++)
         *(temp_cluster_centroid_level1+i) = (int *)malloc(dim * sizeof(int));
    temp_cluster_raidus_level1 =(double*)malloc(level1_clustercount * sizeof(double));
    temp_clusterstart_level1 = (int*)malloc(level1_clustercount * sizeof(int));
    //initialize to zero
    for(i=0;i<level1_clustercount;i++)
    {
        clustersize_level1[i]=0;
        temp_clustersize_level1[i]=0;
        cluster_raidus_level1[i]=0;
        temp_cluster_raidus_level1[i]=0;
    }
    /*Find the total number of cluster in level 2 based on number data points in level 1th cluster */
    for(i=0;i<k;i++)
    {
        array[i][0]=(int)ceil(cbrt(clustersize[i]));
    }
    for(i=0;i<level1_clustercount;i++)
        for(j=0;j<temp;j++)
            cluster_assign_level1[i][j]=0;
    /*Sub Cluster level 1 clusters*/
    for(i=0;i<k;i++)
    {
        /*Provided level 1 cluster size is greater than Zero*/
        if(clustersize[i]!=0)
        {
            temp_data = (int*)malloc(clustersize[i]*dim * sizeof(int));
            for(j=0;j<(clustersize[i]*dim);j++)
            {
                temp_data[j]=data[clusterstart[i]+j];
            }size=clustersize[i]*dim;
            int m=array[i][0];
            KMEANS(dim,size,temp_data,m,temp_cluster_raidus_level1,temp_cluster_centroid_level1,cluster_assign_level1,temp_clustersize_level1,searchData,temp_clusterstart_level1);
            for(j=0;j<(clustersize[i]*dim);j++)
            {
                data[clusterstart[i]+j]=temp_data[j];
            }
            array[i][1]=count;
            for(j=0;j<array[i][0];j++)
            {
                clustersize_level1[count]=temp_clustersize_level1[j];
                cluster_raidus_level1[count]=temp_cluster_raidus_level1[j];
                clusterstart_level1[count]=temp_clusterstart_level1[j]+clusterstart[i];
                for(l=0;l<dim;l++)
                    cluster_centroid_level1[count][l]=temp_cluster_centroid_level1[j][l];

                count++;
            }


        }

    }
    count=0;
    for(i=0;i<k;i++)
    {
        printf(" Level 1 cluster %d clustersize %d \n",i,clustersize[i]);
        int n;
        if(i==(k-1))
            n=level1_clustercount;
        else
            n=array[i+1][1];
        for(j=array[i][1];j<n;j++)
        {
            printf("        Level 2 Sub Cluster size of %d \n",clustersize_level1[j]);

        }
    }

    Search(dim,ndata,data,searchData,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,clusterstart,cluster_assign_level1,cluster_centroid_level1,clustersize_level1,clusterstart_level1,cluster_raidus_level1,array);
    return 0;
}

