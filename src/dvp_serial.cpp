#include<stdlib.h>
#include<stdio.h>
#include<cassert>
#include<map>
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

const char* SPACE = " ";

ifstream graphdata(filename);

std::string line;
std::map<int, int**> routing_tables;
int num_nodes;
int num_edges;

//  Initialsation of routing tables
if(graphdata.is_open()){
   getline(graphdata, line);
   num_nodes = atoi(line.c_str());
   getline(graphdata, line);
   num_edges = atoi(line.c_str());
   cout << "Num of nodes is:" << num_nodes << "Num of Edges is: " << num_edges << "\n";
   for(int node = 0; node < num_nodes; node++){
      int **routing_table = (int**)malloc(sizeof(int*)*num_nodes);
      for(int i=0;i<num_nodes;i++){
         routing_table[i] = (int *)malloc(sizeof(int)*num_nodes);
         for(int j=0;j<num_nodes;j++){
            routing_table[i][j] = INT_MAX/2;
         }
      }
      routing_tables[node] = routing_table;
   }
//Weights are obtained from the input file and placed into the routing tables
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

      int **left_node_routing_table = routing_tables[left];
      left_node_routing_table[left][right] = weight;

      left_node_routing_table[left][left] = 0; 
   }
   graphdata.close();
}
double t1=0.0,t2=0.0;
t1=MPI_Wtime();
int updates_available = 1;

//DVP starts here
while(updates_available){
   updates_available = 0;
   for(int i=0;i<num_nodes;i++){
      int** routing_table_of_node = routing_tables[i];
      for(int j=0;j<num_nodes;j++){
         for(int k=0;k<num_nodes;k++){
            if(j==i){
                int already_existing_value = routing_table_of_node[j][k];
                for(int l=0;l<num_nodes;l++){
                   int node_to_L = routing_table_of_node[i][l];
                   int** l_to_k_table = routing_tables[l];
                   int l_to_k = l_to_k_table[l][k];
                   int total_reroute_val = node_to_L + l_to_k;
                   if(already_existing_value > total_reroute_val){
                      routing_table_of_node[j][k] = total_reroute_val;
                      updates_available = 1;
                   }
               }
           }else{
               int** routing_table_of_j = routing_tables[j];
               routing_table_of_node[j][k] = routing_table_of_j[j][k];
           }
         }
     } 
   }
}
t2=MPI_Wtime();
printf("Avg my_bcast time=%f",t2-t1);
printf("\n");

//Printing the routing tables

/*   int** routing_table = routing_tables[0];
   printf("Routing table of node 0\n");
   printf("------------------------------\n");
   for(int j=0;j<num_nodes;j++){
      for(int k=0;k<num_nodes;k++){
         int routing_val = routing_table[j][k];
         if(routing_val != INT_MAX/2){
            printf("%d\t", routing_table[j][k]);
         }else{
            printf("-\t");
         }
   }
      printf("\n");
}*/
}
