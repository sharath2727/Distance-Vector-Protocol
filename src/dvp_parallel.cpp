#include<stdlib.h>
#include<stdio.h>
#include<cassert>
#include<string>
#include<fstream>
#include<iostream>
#include<limits.h>
#include<vector>
#include<sstream>
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[]){
   if(argc != 2){
      printf("Usage: Please pass a file name with graph information.\n");
   }
   char* filename = argv[1];
   
   ifstream graphdata(filename);
   
   std::string line;
   int num_nodes;
   int num_edges;
   
   graphdata.is_open();
   getline(graphdata, line);
   num_nodes = atoi(line.c_str());
   int routing_tables[num_nodes][num_nodes];
   int next_hop_table[num_nodes][num_nodes];
   getline(graphdata, line);
   num_edges = atoi(line.c_str());

   for(int j=0;j<num_nodes;j++){
     for(int k=0;k<num_nodes;k++){
       routing_tables[j][k] = INT_MAX/2;
       next_hop_table[j][k] = -1;
     }
   }
  //The weightes are placed into the routing tables
  //We have a hot table to calculate the next hop and it is initialised here.
   while(getline(graphdata, line)){
     if(line.length() == 0){
       continue;
     }
     std::stringstream edge_ss(line);
     std::string edge_token;
     
     getline(edge_ss, edge_token, ' ');
     int left = atoi(edge_token.c_str());
     getline(edge_ss, edge_token, ' ');
     int right = atoi(edge_token.c_str());
     getline(edge_ss, edge_token, ' ');
     int weight = atoi(edge_token.c_str());
     
     routing_tables[left][right] = weight;
     next_hop_table[left][right] = right;
     
     routing_tables[left][left] = 0;
     next_hop_table[left][left] = 0; 
   }
   graphdata.close();
   
   int updates_available = 1;
   int rank, size;
   
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
  //Timing analysis- start
   double t1=0.0,t2=0.0;
   MPI_Barrier(MPI_COMM_WORLD);
   t1=MPI_Wtime();

   int remainder = num_nodes % size;
   int blocksize = num_nodes/size;

   int num_rows = rank < remainder ? blocksize+1 : blocksize;

   int* counts = (int *)malloc(sizeof(int)*size);
   int* displs = (int *)malloc(sizeof(int)*size);

   //Here we handle the condition when the number of nodes is not divisible by the number of processors-Load Balancing
   for(int i=0;i<size;i++){
     if(i<remainder){
       counts[i] = (blocksize+1)*num_nodes;
       displs[i] = i*(blocksize+1)*num_nodes;
     }else{
       counts[i] = blocksize*num_nodes;
       displs[i] = remainder*(blocksize+1)*num_nodes + (i-remainder)*blocksize*num_nodes;
     }
   }

   int start_row = -1;
   if(rank < remainder){
     start_row = rank*(blocksize+1);
   }else{
     start_row = remainder*(blocksize+1) + (rank-remainder)*blocksize;
   }
   /* DVP Starts here */
   while(updates_available){
     updates_available = 0;
     for(int i=0;i<num_rows;i++){
       int node_num = start_row + i;
       for(int k=0;k<num_nodes;k++){
         int already_existing_value = routing_tables[node_num][k];
         for(int l=0;l<num_nodes;l++){
           int node_to_L = routing_tables[node_num][l];
           int l_to_k = routing_tables[l][k];
           int total_reroute_val = node_to_L + l_to_k;
           if(already_existing_value > total_reroute_val){
             routing_tables[node_num][k] = total_reroute_val;
             next_hop_table[node_num][k] = l;
             updates_available = 1;
           }
         }
       }
     }
     MPI_Allgatherv(routing_tables[start_row], num_rows*num_nodes, MPI_INT, routing_tables, counts, displs, MPI_INT, MPI_COMM_WORLD); 
     MPI_Allgatherv(next_hop_table[start_row], num_rows*num_nodes, MPI_INT, next_hop_table, counts, displs, MPI_INT, MPI_COMM_WORLD); 
     MPI_Allreduce(&updates_available, &updates_available, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   }
    
   MPI_Barrier(MPI_COMM_WORLD);
   t2=MPI_Wtime();
   //timing ends here
   double mytime=t2-t1;
   double max,min,avgtime;
   MPI_Reduce(&mytime,&avgtime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&mytime,&max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&mytime,&min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
 
   
   if(rank==0){
     printf("Average time: %f \t Maximum time: %f \t Minimum time: %f\n",(avgtime/size),max,min);
   }

   //Routing tables printed in this section

  /* if(rank == 0){
     printf("\nRouting table of node %d\n", rank);
     printf("------------------------------\n");
     for(int j=0;j<num_nodes;j++){
       for(int k=0;k<num_nodes;k++){
         int routing_val = routing_tables[j][k];
         if(routing_val != INT_MAX/2){
           printf("%d/%d\t", routing_tables[j][k], next_hop_table[j][k]);
         }else{
           printf("-/-\t");
         }
       }
       printf("\n");
     }
   }
   */
   MPI_Finalize();
   return 0;
}
