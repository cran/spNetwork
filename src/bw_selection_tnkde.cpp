#include "spNetwork.h"
#include "base_kernel_funtions.h"
#include "matrices_functions.h"

#include <cmath>


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTIONS TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR SIMPLE TNKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The worker function to calculate simple TNKDE likelihood cv
//' @name ess_kernel_loo_tnkde
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
//' @keywords internal
arma::cube ess_kernel_loo_tnkde(fptros kernel_func, arma::sp_imat &edge_mat,
                                IntegerVector &events,
                                IntegerVector &events_wid,
                                NumericVector &time_events,
                                List &neighbour_list,
                                int v, int wid, double v_time,
                                arma::vec &bws_net, arma::vec &bws_time,
                                NumericVector &line_weights, int max_depth){

  //step0 : generate the queue
  int depth = 0;
  std::queue <List> data_holder;
  arma::cube kvalues(bws_net.n_elem, bws_time.n_elem, events.length());
  kvalues.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0
  );
  data_holder.push(cas1);

  double max_bw_net = arma::max(bws_net);
  double bw_net, bw_time;

  // avant de lancer les iteration, nous devons vérifier s'il est necessaire
  // de rajouter la densite a des evenement qui se situeraient sur le meme point
  // de depart !
  std::vector<int> index = get_all_indeces_int(events,v);
  if(index.size() > 1){
    for(int ii = 0; ii < bws_net.n_elem ; ii++){
      bw_net = bws_net(ii);
      double kernel_net =  kernel_func(0,bw_net);
      for(int j = 0 ; j < bws_time.n_elem; j ++){
        bw_time = bws_time(j);
        // NOTE : we are doing the bw scaling here
        for (int zz : index){
          if(time_events[zz] != v_time){
            double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
            kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
          }
        }
      }
    }
  }

  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.front();
    data_holder.pop();
    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    int prev_node = cas["prev_node"];

    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          // find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double d2 = line_weights[edge_id-1] + d;
          std::vector<int> index = get_all_indeces_int(events,v2);
          if(index.size() >0 ){
            // il semble que v2 soit un noeud pour lequel au moins un evenement est present
            for(int ii = 0; ii < bws_net.n_elem ; ii++){
              bw_net = bws_net(ii);
              double kernel_net =  kernel_func(d2,bw_net);
              for(int j = 0 ; j < bws_time.n_elem; j ++){
                bw_time = bws_time(j);
                // NOTE : we are doing the bw scaling here
                for (int zz : index){
                  double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
                  kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
                }

              }
            }
          }
          //evaluating for the next move
          if (d2< max_bw_net and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return kvalues;
}


//' @title The worker function to calculate simple TNKDE likelihood cv (adaptive case)
//' @name ess_kernel_loo_tnkde_adpt
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::mat with the network bandwidths to consider
//' @param bws_time an arma::mat with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
//' @keywords internal
 arma::cube ess_kernel_loo_tnkde_adpt(fptros kernel_func, arma::sp_imat &edge_mat,
                                 IntegerVector &events,
                                 IntegerVector &events_wid,
                                 NumericVector &time_events,
                                 List &neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::mat &bws_net, arma::mat &bws_time,
                                 NumericVector &line_weights, int max_depth){

   //step0 : generate the queue
   // Rcout << "starting the new function !\n";
   int depth = 0;
   std::queue <List> data_holder;
   arma::cube kvalues(bws_net.n_rows, bws_net.n_cols, events.size());
   kvalues.fill(0.0);
   //step1 : generate the first case
   List cas1 = List::create(Named("d")=0.0,
                            Named("v") = v,
                            Named("prev_node") = -999,
                            Named("depth") = 0
   );
   data_holder.push(cas1);

    // NOTE : POSSIBLE ENHANCEMENT HERE, I COULD RESET THE max_bw_net
    // DURING THE ITERATIONS TO REDUCE A BIT THE CALCULATION TIME. WE WOULD
    // AVOID SOME ITERATION THIS WAY.
   double max_bw_net = bws_net.max();
   double bw_net, bw_time;

   // avant de lancer les iterations, il est important de verifier que
   // certains evenement ne se trouvent pas sur la meme vertex d'origine.
   // si c'est le cas, nous devons leur assigner la valeur de densite correspondante
   std::vector<int> index = get_all_indeces_int(events,v);
   if(index.size() > 1){
     // Rcout << "See me ! I must add some density at the origin ! \n";
     for (int zz : index){
       if(time_events[zz] != v_time){
         // arma::mat slice_bw_net = bws_net.slice(zz);
         // arma::mat slice_bw_time = bws_time.slice(zz);
         //Rcout << "step8\n";
         for(int j = 0 ; j < bws_net.n_rows; j ++){
           for(int jj = 0 ; jj < bws_time.n_cols; jj ++){
             // NOTE : we are doing the bw scaling here
             double bw_net = bws_net(j,jj);
             double bw_time = bws_time(j,jj);
             double kernel_net =  kernel_func(0.0, bw_net);
             double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
             kvalues(j,jj,zz) = kvalues(j,jj,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
           }
         }
       }
     }
   }



   //lancement des iterations
   while(data_holder.empty()==FALSE){
     //unpacking (imagine some loop unrolling here with a function to deal with.)
     List cas = data_holder.front();
     data_holder.pop();
     int v = cas["v"];
     int depth = cas["depth"];
     double d = cas["d"];
     int prev_node = cas["prev_node"];
     //step1 : find all the neighbours
     IntegerVector neighbours = neighbour_list[v-1];

     //step2 : iterate over the neighbours
     int cnt_n = neighbours.length();
     int new_depth;
     if(cnt_n>2){
       new_depth = depth+1;
     }else{
       new_depth = depth;
     }
     //if we have only one neighbour, we must stop
     if(cnt_n>1 or prev_node<=0){
       for(int i=0; i < cnt_n; ++i){
         int v2 = neighbours[i];
         // comment retrouver le potentiel wid du noeud ?
         //on ne veut pas revenir en arriere !
         if(v2!=prev_node){
           //Rcout << "step4\n";
           // find the edge between the two nodes
           int edge_id = edge_mat(v,v2);
           double d2 = line_weights[edge_id-1] + d;
           //Rcout << "step5\n";
           std::vector<int> index = get_all_indeces_int(events,v2);
           //Rcout << "step6\n";
           if(index.size() >0 ){
             // Rcout << "step7\n";
             // il semble que v2 soit un noeud pour lequel au moins un evenement est present
             // NB : nrow pour les BW network et ncol pour le bw temporelles
             for (int zz : index){
               //arma::mat slice_bw_net = bws_net.slice(zz);
               //arma::mat slice_bw_time = bws_time.slice(zz);
               //Rcout << "step8\n";
               for(int j = 0 ; j < bws_net.n_rows; j ++){
                 for(int jj = 0 ; jj < bws_time.n_cols; jj ++){
                   // NOTE : we are doing the bw scaling here
                   double bw_net = bws_net(j,jj);
                   double bw_time = bws_time(j,jj);
                   double kernel_net =  kernel_func(d2, bw_net);
                   double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
                   //Rcout << "step9\n";
                   kvalues(j,jj,zz) = kvalues(j,jj,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
                 }
               }
             }

           }
           //evaluating for the next move
           if (d2< max_bw_net and new_depth<max_depth){
             List new_cas = List::create(Named("d")=d2,
                                         Named("v") = v2,
                                         Named("prev_node") = v,
                                         Named("depth") = new_depth
             );
             data_holder.push(new_cas);
           }
         }
       }
     }
   }
   // Rcout << "OUT OF FUNCTION\n";
   return kvalues;
 }



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTION TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR DISCONTINUOUS NKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The worker function to calculate discontinuous TNKDE likelihood cv
//' @name esd_kernel_loo_tnkde
//' @description The worker function to calculate discontinuous TNKDE likelihood cv (INTERNAL)
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
arma::cube esd_kernel_loo_tnkde(fptros kernel_func, arma::sp_imat &edge_mat,
                                 IntegerVector &events,
                                 IntegerVector &events_wid,
                                 NumericVector &time_events,
                                 List &neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::vec &bws_net, arma::vec &bws_time,
                                 NumericVector &line_weights, int max_depth){

  //step0 : generate the queue
  int depth = 0;
  std::queue <List> data_holder;
  arma::cube kvalues(bws_net.n_elem, bws_time.n_elem, events.length());
  kvalues.fill(0.0);
  //step1 : generate the first case
  List cas1 = List::create(Named("d")=0.0,
                           Named("v") = v,
                           Named("prev_node") = -999,
                           Named("depth") = 0,
                           Named("alpha") = 1
  );
  data_holder.push(cas1);

  double max_bw_net = arma::max(bws_net);
  double bw_net, bw_time;

  // avant de lancer les iteration, nous devons vérifier s'il est necessaire
  // de rajouter la densite a des evenement qui se situeraient sur le meme point
  // de depart !
  std::vector<int> index = get_all_indeces_int(events,v);
  if(index.size() > 1){
    for(int ii = 0; ii < bws_net.n_elem ; ii++){
      bw_net = bws_net(ii);
      double kernel_net =  kernel_func(0,bw_net);
      for(int j = 0 ; j < bws_time.n_elem; j ++){
        bw_time = bws_time(j);
        // NOTE : we are doing the bw scaling here
        for (int zz : index){
          if(time_events[zz] != v_time){
            double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
            kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
          }
        }
      }
    }
  }


  //lancement des iterations
  while(data_holder.empty()==FALSE){
    //unpacking (imagine some loop unrolling here with a function to deal with.)
    List cas = data_holder.front();
    data_holder.pop();

    int v = cas["v"];
    int depth = cas["depth"];
    double d = cas["d"];
    double alpha = cas["alpha"];
    int prev_node = cas["prev_node"];


    //step1 : find all the neighbours
    IntegerVector neighbours = neighbour_list[v-1];

    //step2 : iterate over the neighbours
    int cnt_n = neighbours.length();
    int new_depth;
    if(cnt_n>2){
      new_depth = depth+1;
    }else{
      new_depth = depth;
    }

    //step 3 prepare the new alpha value
    double new_alpha;
    if((prev_node < 0)  && (cnt_n > 2)){
      new_alpha = 2.0/(cnt_n);
    }else if((prev_node < 0) && (cnt_n == 1)){
      new_alpha = 1.0;
    }else{
      new_alpha = alpha * (1.0/(cnt_n-1.0));
    }

    //if we have only one neighbour, we must stop
    if(cnt_n>1 or prev_node<=0){
      for(int i=0; i < cnt_n; ++i){
        int v2 = neighbours[i];
        //on ne veut pas revenir en arriere !
        if(v2!=prev_node){
          //find the edge between the two nodes
          int edge_id = edge_mat(v,v2);
          double d2 = line_weights[edge_id-1] + d;

          //est ce que v2 est un evenement pour lequel on doit garder la valeur

          std::vector<int> index = get_all_indeces_int(events,v2);
          if(index.size() >0 ){

            for(int ii = 0; ii < bws_net.n_elem ; ii++){
              bw_net = bws_net(ii);
              double kernel_net =  kernel_func(d2,bw_net) * new_alpha;
              for(int j = 0 ; j < bws_time.n_elem; j ++){
                bw_time = bws_time(j);
                for (int zz : index){
                  double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
                  kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
                }
              }
            }
          }

          //evaluating for the next move
          if (d2< max_bw_net and new_depth<max_depth){
            List new_cas = List::create(Named("d")=d2,
                                        Named("v") = v2,
                                        Named("prev_node") = v,
                                        Named("depth") = new_depth,
                                        Named("alpha")=new_alpha
            );
            data_holder.push(new_cas);
          }
        }
      }
    }
  }

  return kvalues;
}


//' @title The worker function to calculate discontinuous TNKDE likelihood cv (adaptive case)
//' @name esd_kernel_loo_tnkde_adpt
//' @description The worker function to calculate discontinuous TNKDE likelihood cv (INTERNAL)
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::mat with the network bandwidths to consider
//' @param bws_time an arma::mat with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
//' @keywords internal
 arma::cube esd_kernel_loo_tnkde_adpt(fptros kernel_func, arma::sp_imat &edge_mat,
                                 IntegerVector &events,
                                 IntegerVector &events_wid,
                                 NumericVector &time_events,
                                 List& neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::mat &bws_net,
                                 arma::mat &bws_time,
                                 NumericVector &line_weights, int max_depth){

   //step0 : generate the queue
   int depth = 0;
   std::queue <List> data_holder;
   arma::cube kvalues(bws_net.n_rows, bws_time.n_cols, events.length());
   kvalues.fill(0.0);
   //step1 : generate the first case
   List cas1 = List::create(Named("d")=0.0,
                            Named("v") = v,
                            Named("prev_node") = -999,
                            Named("depth") = 0,
                            Named("alpha") = 1
   );
   data_holder.push(cas1);

   double max_bw_net = bws_net.max();
   double bw_net, bw_time;

   // avant de lancer les iterations, il est important de verifier que
   // certains evenement ne se trouvent pas sur la meme vertex d'origine.
   // si c'est le cas, nous devons leur assigner la valeur de densite correspondante
   std::vector<int> index = get_all_indeces_int(events,v);
   if(index.size() > 1){
     // Rcout << "See me ! I must add some density at the origin ! \n";
     for (int zz : index){
       if(time_events[zz] != v_time){
         //arma::mat slice_bw_net = bws_net.slice(zz);
         //arma::mat slice_bw_time = bws_time.slice(zz);
         //Rcout << "step8\n";
         for(int j = 0 ; j < bws_net.n_rows; j ++){
           for(int jj = 0 ; jj < bws_time.n_cols; jj ++){
             // NOTE : we are doing the bw scaling here
             double bw_net = bws_net(j,jj);
             double bw_time = bws_time(j,jj);
             double kernel_net =  kernel_func(0.0, bw_net);
             double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
             kvalues(j,jj,zz) = kvalues(j,jj,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
           }
         }
       }
     }
   }



   //lancement des iterations
   while(data_holder.empty()==FALSE){
     //unpacking (imagine some loop unrolling here with a function to deal with.)
     List cas = data_holder.front();
     data_holder.pop();
     int v = cas["v"];
     int depth = cas["depth"];
     double d = cas["d"];
     double alpha = cas["alpha"];
     int prev_node = cas["prev_node"];


     //step1 : find all the neighbours
     IntegerVector neighbours = neighbour_list[v-1];

     //step2 : iterate over the neighbours
     int cnt_n = neighbours.length();
     int new_depth;
     if(cnt_n>2){
       new_depth = depth+1;
     }else{
       new_depth = depth;
     }

     //step 3 prepare the new alpha value
     double new_alpha;
     if((prev_node < 0)  && (cnt_n > 2)){
       new_alpha = 2.0/(cnt_n);
     }else if((prev_node < 0) && (cnt_n == 1)){
       new_alpha = 1.0;
     }else{
       new_alpha = alpha * (1.0/(cnt_n-1.0));
     }

     //if we have only one neighbour, we must stop
     if(cnt_n>1 or prev_node<=0){
       for(int i=0; i < cnt_n; ++i){
         int v2 = neighbours[i];
         //on ne veut pas revenir en arriere !
         if(v2!=prev_node){
           //find the edge between the two nodes
           int edge_id = edge_mat(v,v2);
           double d2 = line_weights[edge_id-1] + d;

           //est ce que v2 est un evenement pour lequel on doit garder la valeur

           std::vector<int> index = get_all_indeces_int(events,v2);
           if(index.size() >0 ){
             for (int zz : index){
               //arma::mat slice_bw_net = bws_net.slice(zz);
               //arma::mat slice_bw_time = bws_time.slice(zz);
               for(int j = 0 ; j < bws_net.n_rows; j ++){
                 for(int jj = 0 ; jj < bws_time.n_cols; jj ++){
                   // NOTE : we are doing the bw scaling here
                   double bw_net = bws_net(j,jj);
                   double bw_time = bws_time(j,jj);
                   double kernel_net =  kernel_func(d2,bw_net) * new_alpha;
                   double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
                   kvalues(j,jj,zz) = kvalues(j,jj,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
                 }
               }
             }
             // for(int ii = 0; ii < bws_net.n_elem ; ii++){
             //   bw_net = bws_net(ii);
             //   double kernel_net =  kernel_func(d2,bw_net) * new_alpha;
             //   for(int j = 0 ; j < bws_time.n_elem; j ++){
             //     bw_time = bws_time(j);
             //     for (int zz : index){
             //       double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
             //       kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
             //     }
             //   }
             // }
           }

           //evaluating for the next move
           if (d2< max_bw_net and new_depth<max_depth){
             List new_cas = List::create(Named("d")=d2,
                                         Named("v") = v2,
                                         Named("prev_node") = v,
                                         Named("depth") = new_depth,
                                         Named("alpha")=new_alpha
             );
             data_holder.push(new_cas);
           }
         }
       }
     }
   }

   return kvalues;
 }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE FUNCTION TO CALCULATE BW SELECTION CV LOO TEMPORAL FOR CONTINUOUS NKDE
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//' @title The worker function to calculate continuous TNKDE likelihood cv
//' @name esc_kernel_loo_tnkde
//' @description The worker function to calculate continuous TNKDE likelihood cv (INTERNAL)
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
arma::cube esc_kernel_loo_tnkde(fptros kernel_func, arma::sp_imat &edge_mat,
                                 IntegerVector &events,
                                 IntegerVector &events_wid,
                                 NumericVector &time_events,
                                 List &neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::vec &bws_net, arma::vec &bws_time,
                                 NumericVector &line_weights, int max_depth){
  //step0 : generate the queue
  int depth = 0;
  struct acase{
    int v;
    int prev_node;
    int depth;
    double d;
    double alpha;
  };

  // prepare the data cube and matrix
  arma::cube kvalues(bws_net.n_elem, bws_time.n_elem, events.length());
  arma::mat net_k_values(events.length(),bws_net.n_elem);
  arma::mat time_k_values(events.length(),bws_time.n_elem);
  kvalues.fill(0.0);
  net_k_values.fill(0.0);
  time_k_values.fill(0.0);

  double max_bw_net = arma::max(bws_net);
  double bw_net, bw_time;

  //queue<acase> data_holder;
  std::vector<acase> data_holder;

  // let us prepare the first cases
  IntegerVector v_neighbours = neighbour_list[v-1];
  double alpha = 2.0/v_neighbours.length();

  acase el = {v,-999,0,0.0,alpha};
  data_holder.push_back(el);

  // avant de lancer les iteration, nous devons vérifier s'il est necessaire
  // de rajouter la densite a des evenement qui se situeraient sur le meme point
  // de depart !
  std::vector<int> index = get_all_indeces_int(events,v);
  if(index.size() > 1){
    for(int ii = 0; ii < bws_net.n_elem ; ii++){
      bw_net = bws_net(ii);
      double kernel_net =  kernel_func(0,bw_net);
      for(int j = 0 ; j < bws_time.n_elem; j ++){
        bw_time = bws_time(j);
        // NOTE : we are doing the bw scaling here
        for (int zz : index){
          if(time_events[zz] != v_time){
            double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
            kvalues(ii,j,zz) = kvalues(ii,j,zz) + ((kernel_net * kernel_time) * (1.0/(bw_net*bw_time)));
          }
        }
      }
    }
  }



  //lancement des iterations
  while(data_holder.empty()==FALSE){

    //unpacking (imagine some loop unrolling here with a function to deal with.)
    acase cas = data_holder.back();
    data_holder.pop_back();
    // int v = cas.v;
    // int prev_node = cas.prev_node;
    // int depth = cas.depth;
    // double d = cas.d;
    // double alpha = cas.alpha;

    //Rcout << "v_old is : "<< prev_node <<"\n";
    //Rcout << "v is : "<< v <<"\n";
    //Rcout << "d is : "<< d<<"\n";
    //Rcout << "alpha is : "<< alpha <<"\n";
    //Rcout << "\n\n";

    // we will update the densities on v
    // but only if v is a vertex on wich I can find an event

    std::vector<int> index = get_all_indeces_int(events,cas.v);

    // NOTE : WE SHOULD FIRST CALCULATE THE COMPLETE NKDE DENSITY CAUSED BY AN EVENT
    // AND ONLY AFTER MULTIPLYING IT WITH THE TIME DENSITY
    if(index.size() >0 ){
      for(int ii = 0; ii < bws_net.n_elem ; ii++){
        bw_net = bws_net(ii);
        double kernel_net =  kernel_func(cas.d,bw_net) * cas.alpha;
        // NOTE : we are not doing the bw scaling here but later
        for (int zz : index){
          net_k_values(zz,ii) = net_k_values(zz,ii) + kernel_net;
        }
      }
    }

    // we can now prepare the next steps
    if(max_bw_net >= cas.d){
      IntegerVector v_neighbours = neighbour_list[cas.v-1];
      int n = v_neighbours.length();
      int new_depth;

      //updating depth
      if(n>2){
        new_depth = cas.depth+1;
      }else{
        new_depth = cas.depth;
      }

      if(n>1){
        // we must continue only if we are not too deep
        if(new_depth <= max_depth){
          // and iterate over the neighbours
          for(int j = 0; j < n ; j++){
            // ---- first, the simple case, we are not going back ----
            int v2 = v_neighbours[j];
            int l2 = edge_mat(cas.v,v2);
            double d2 = cas.d + line_weights[l2-1];
            // first case, we must back fire
            if(v2 == cas.prev_node){
              if(n>2){
                double p2 = (2.0-n)/n;
                double new_alpha = cas.alpha * p2;
                acase new_case = {v2,cas.v,new_depth,d2,new_alpha};
                data_holder.push_back(new_case);
                //data_holder.push_back((struct acase){v2,cas.v,new_depth,d2,new_alpha});
              }
            }else{
              double new_alpha = cas.alpha * (2.0/n);
              acase new_case = {v2,cas.v,new_depth,d2,new_alpha};
              data_holder.push_back(new_case);
              //data_holder.push_back((struct acase){v2,cas.v,new_depth,d2,new_alpha});
            }
          }
        }
      }
    }
  }

  //Rcout << "HERE ARE THE CALCULATED NETWORK DENSITIES\n";
  //Rcout << net_k_values <<"\n\n";

  // LET US CALCULATE THE TEMPORAL DENSITIES
  for(int e = 0; e < events.length() ; e++){
    double d = std::abs(v_time - time_events(e));
    for(int t = 0; t < bws_time.n_elem; t++){
      time_k_values(e,t) = kernel_func(d,bws_time(t));
    }
  }
  //Rcout << "HERE ARE THE CALCULATED TIME DENSITIES\n";
  //Rcout << time_k_values <<"\n\n";


  // BEFORE THE SCALING HERE, I MUST APPLY THE TEMPORAL DENSITIES TOO

  for(int e = 0; e < events.length() ; e++){
    for(int n = 0; n < bws_net.n_elem; n++){
      for(int t = 0; t < bws_time.n_elem; t++){
        kvalues(n,t,e) = net_k_values(e,n) * time_k_values(e,t);
      }
    }
  }

  //Rcout << "HERE ARE THE UNSCALED DENSITIES\n";
  //Rcout << kvalues <<"\n\n";

  // and now we can apply the scaling
  arma::mat scale_mat(bws_net.n_elem, bws_time.n_elem);

  for(int i = 0; i < bws_time.n_elem; i++){
    scale_mat.col(i) = 1.0/(bws_net * bws_time(i));
  }

  for(int i = 0; i < events.length(); i++){
    kvalues.slice(i) = kvalues.slice(i) % scale_mat;
  }


  return kvalues;
}

//' @title The worker function to calculate continuous TNKDE likelihood cv (adaptive case)
//' @name esc_kernel_loo_tnkde_adpt
//' @description The worker function to calculate continuous TNKDE likelihood cv (INTERNAL)
//' @param kernel_func a cpp pointer function (selected with the kernel name)
//' @param edge_mat matrix, to find the id of each edge given two neighbours.
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param time_events a NumericVector indicating the timestamp of each event
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param v the actual node to consider (int)
//' @param v_time the time of v (double)
//' @param bws_net an arma::mat with the network bandwidths to consider
//' @param bws_time an arma::mat with the time bandwidths to consider
//' @param line_weights a vector with the length of the edges
//' @param depth the actual recursion depth
//' @param max_depth the maximum recursion depth
//' @return a cube with the impact of the event v on each other event for
//' each pair of bandwidths (cube(bws_net, bws_time, events))
//' @keywords internal
 arma::cube esc_kernel_loo_tnkde_adpt(fptros kernel_func, arma::sp_imat &edge_mat,
                                 IntegerVector &events,
                                 IntegerVector &events_wid,
                                 NumericVector &time_events,
                                 List &neighbour_list,
                                 int v, int wid, double v_time,
                                 arma::mat &bws_net,
                                 arma::mat &bws_time,
                                 NumericVector &line_weights, int max_depth){
   //step0 : generate the queue
   int depth = 0;
   struct acase{
     int v;
     int prev_node;
     int depth;
     double d;
     double alpha;
   };

   // prepare the data cube and matrix
   int n_bw_net = bws_net.n_rows;
   int n_bw_time = bws_time.n_cols;
   arma::cube kvalues(bws_net.n_rows, bws_time.n_cols, events.length());
   arma::cube net_k_values(bws_net.n_rows, bws_time.n_cols, events.length());
   arma::cube time_k_values(bws_net.n_rows, bws_time.n_cols, events.length());
   kvalues.fill(0.0);
   net_k_values.fill(0.0);
   time_k_values.fill(0.0);

   double max_bw_net = bws_net.max();
   double bw_net, bw_time;

   //queue<acase> data_holder;
   std::vector<acase> data_holder;

   // let us prepare the first cases
   IntegerVector v_neighbours = neighbour_list[v-1];
   double alpha = 2.0/v_neighbours.length();

   acase el = {v,-999,0,0.0,alpha};
   data_holder.push_back(el);

   // avant de lancer les iterations, il est important de verifier que
   // certains evenement ne se trouvent pas sur la meme vertex d'origine.
   // si c'est le cas, nous devons leur assigner la valeur de densite correspondante
   std::vector<int> index = get_all_indeces_int(events,v);
   if(index.size() > 1){
     // Rcout << "See me ! I must add some density at the origin ! \n";
     for (int zz : index){
       if(time_events[zz] != v_time){
         //arma::mat slice_bw_net = bws_net.slice(zz);
         //arma::mat slice_bw_time = bws_time.slice(zz);
         //Rcout << "step8\n";
         for(int j = 0 ; j < bws_net.n_rows; j ++){
           for(int jj = 0 ; jj < bws_time.n_cols; jj ++){
             // NOTE : we are NOT doing the bw scaling here
             double bw_net = bws_net(j,jj);
             double bw_time = bws_time(j,jj);
             double kernel_net =  kernel_func(0.0, bw_net);
             double kernel_time = kernel_func(std::abs(v_time-time_events[zz]), bw_time);
             kvalues(j,jj,zz) = kvalues(j,jj,zz) + ((kernel_net * kernel_time));
           }
         }
       }
     }
   }


   //lancement des iterations
   while(data_holder.empty()==FALSE){

     //unpacking (imagine some loop unrolling here with a function to deal with.)
     acase cas = data_holder.back();
     data_holder.pop_back();
     int v = cas.v;
     int prev_node = cas.prev_node;
     int depth = cas.depth;
     double d = cas.d;
     double alpha = cas.alpha;

     // Rcout << "v_old is : "<< prev_node <<"\n";
     // Rcout << "v is : "<< v <<"\n";
     // Rcout << "d is : "<< d<<"\n";
     // Rcout << "alpha is : "<< alpha <<"\n";
     // Rcout << "\n\n";

     // we will update the densities on v
     // but only if v is a vertex on wich I can find an event

     std::vector<int> index = get_all_indeces_int(events,v);

     // NOTE : WE SHOULD FIRST CALCULATE THE COMPLETE NKDE DENSITY CAUSED BY AN EVENT
     // AND ONLY AFTER MULTIPLYING IT WITH THE TIME DENSITY
     if(index.size() >0 ){
       for(int j = 0; j < bws_net.n_rows ; j++){
         for(int jj = 0 ; jj < bws_time.n_cols ; jj ++){
           bw_net = bws_net(j,jj);
           // NOTE : we are not doing the bw scaling here but later
           for (int zz : index){
             double kernel_net =  kernel_func(d,bw_net) * alpha;
             net_k_values(j,jj,zz) = net_k_values(j,jj,zz) + kernel_net;
           }
         }
       }
     }

     // we can now prepare the next steps
     if(max_bw_net >= d){
       IntegerVector v_neighbours = neighbour_list[v-1];
       int n = v_neighbours.length();
       int new_depth;

       //updating depth
       if(n>2){
         new_depth = depth+1;
       }else{
         new_depth = depth;
       }

       if(n>1){
         // we must continue only if we are not too deep
         if(new_depth <= max_depth){
           // and iterate over the neighbours
           for(int j = 0; j < n ; j++){
             // ---- first, the simple case, we are not going back ----
             int v2 = v_neighbours[j];
             int l2 = edge_mat(v,v2);
             double d2 = d + line_weights[l2-1];
             // first case, we must back fire
             if(v2 == prev_node){
               if(n>2){
                 double p2 = (2.0-n)/n;
                 double new_alpha = alpha * p2;
                 acase new_case = {v2,v,new_depth,d2,new_alpha};
                 //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
                 //data_holder.push(new_case);
                 data_holder.push_back(new_case);
               }
             }else{
               double new_alpha = alpha * (2.0/n);
               acase new_case = {v2,v,new_depth,d2,new_alpha};
               //Rcout << "  adding this cas : d="<<d2<<", v2="<<v2<<", v3="<<v3<<", l2="<<l2<<" alpha="<<new_alpha<<"\n";
               data_holder.push_back(new_case);
             }
           }
         }
       }
     }
   }

   // Rcout << "HERE ARE THE CALCULATED NETWORK DENSITIES\n";
   // Rcout << net_k_values <<"\n\n";

   // Rcout << "HERE IS THE MATRIX OF TIME BANDWIDTHS\n";
   // Rcout << bws_time;
   // Rcout << "\n\n";

   // Rcout << "HERE IS THE MATRIX OF NETWORK BANDWIDTHS\n";
   // Rcout << bws_net;
   // Rcout << "\n\n";

   int n_net = bws_net.n_rows;
   // Rcout << "here is bws_net.n_rows "<<bws_net.n_rows << "\n";
   // LET US CALCULATE THE TEMPORAL DENSITIES
   for(int e = 0; e < events.length() ; e++){
     // Rcout << "this is e : "<< e <<"\n";
     double d = std::abs(v_time - time_events(e));
     for(int n = 0 ; n < bws_net.n_rows ; n++){
       // Rcout << "this is n : "<< n <<"\n";
       for(int t = 0; t < bws_time.n_cols ; t++){
         // Rcout << "this is t : "<< t <<"\n";
         double bw_time = bws_time(n,t);
         // Rcout << "here is bw_time" << bw_time << "\n";
         // Rcout << "here is old time_keval" << time_k_values(n,t,e) << "\n";
         time_k_values(n,t,e) = kernel_func(d,bw_time);
       }
     }
   }
   // Rcout << "HERE ARE THE CALCULATED TIME DENSITIES\n";
   // Rcout << time_k_values <<"\n\n";


   // BEFORE THE SCALING HERE, I MUST APPLY THE TEMPORAL DENSITIES TOO

   kvalues = time_k_values % net_k_values;

   // Rcout << "HERE ARE THE UNSCALED DENSITIES\n";
   // Rcout << kvalues <<"\n\n";

   // and now we can apply the scaling
   arma::mat scale_mat = arma::pow((bws_net % bws_time),-1);
   for(int i = 0; i < events.length(); i ++){
     kvalues.slice(i) = kvalues.slice(i) % scale_mat;
   }


   return kvalues;
 }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// THE EXPOSED FUNCTIONS
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//' @title The exposed function to calculate TNKDE likelihood cv
//' @name tnkde_get_loo_values
//' @description The exposed function to calculate TNKDE likelihood cv (INTERNAL)
//' @param method a string, one of "simple", "continuous", "discontinuous"
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param sel_events a Numeric vector indicating the selected events (id of nodes)
//' @param sel_events_wid a Numeric Vector indicating the unique if of the selected events
//' @param sel_events_time a Numeric Vector indicating the time of the selected events
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param events_wid a NumericVector indicating the unique id of all the events
//' @param events_time a NumericVector indicating the timestamp of each event
//' @param weights a cube with the weights associated with each event for each
//' bws_net and bws_time.
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param kernel_name a string with the name of the kernel to use
//' @param line_list a DataFrame describing the lines
//' @param max_depth the maximum recursion depth
//' @param min_tol a double indicating by how much 0 in density values must be replaced
//' @return a matrix with the CV score for each pair of bandiwdths
//' @export
//' @examples
//' # no example provided, this is an internal function
// [[Rcpp::export]]
arma::cube tnkde_get_loo_values(std::string method, List &neighbour_list,
                               IntegerVector &sel_events,
                               IntegerVector &sel_events_wid,
                               NumericVector &sel_events_time,
                               IntegerVector &events,
                               IntegerVector &events_wid,
                               NumericVector &events_time,
                               arma::cube &weights,
                               arma::vec &bws_net,
                               arma::vec &bws_time,
                               std::string kernel_name,
                               DataFrame &line_list, int max_depth, double min_tol){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0 (bw_net 0 et bw_time 0)
  // NOTE WE calculate the values only for the events in sel_events
  arma::cube base_k(bws_net.n_elem, bws_time.n_elem, sel_events.length());

  //calculer la matrice des lignes
  //arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
  arma::cube k;
  //step2 : iterer sur chaque event de la zone d'etude
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    int wid = events_wid[i];
    double v_time = events_time[i];
    // launching recursion
    // here we got the the influences of the vertex y on each other selected event in quadra
    if(method == "simple"){
      k = ess_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               sel_events_time, neighbour_list,
                               y, wid, v_time,
                               bws_net, bws_time,
                               line_weights, max_depth);

    }else if (method == "discontinuous"){
      k = esd_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);
    }else{
      k = esc_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);
    }
    // NOTE : the scaling by bws is applied above
    // and we must now remove its influence on itself
    arma::mat w = weights.slice(i);
    for(int ii = 0; ii < sel_events.length(); ii++){
      k.slice(ii) = k.slice(ii) % w;
    }

    // if y was a selected event, its own weight must be set to 0
    int index = get_first_index_int(sel_events_wid,wid);
    if(index >= 0){
      k.slice(index).fill(0);
    }

    // and summing the values at each iteration (% is the element wise product)
    base_k += k;
  };

  return base_k;
  // // and calculate the final values
  // arma::mat result;
  //
  // arma::uvec neg_elems = arma::find(base_k <= 0);
  // base_k.elem(neg_elems).fill(min_tol);
  //
  //
  // result = arma::sum(arma::log(base_k),2);
  //
  //
  // return result;
}


//' @title The exposed function to calculate TNKDE likelihood cv
//' @name tnkde_get_loo_values2
//' @description The exposed function to calculate TNKDE likelihood cv (INTERNAL) when an adaptive bandwidth is used
//' @param method a string, one of "simple", "continuous", "discontinuous"
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param sel_events a Numeric vector indicating the selected events (id of nodes)
//' @param sel_events_wid a Numeric Vector indicating the unique if of the selected events
//' @param sel_events_time a Numeric Vector indicating the time of the selected events
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param events_wid a NumericVector indicating the unique id of all the events
//' @param events_time a NumericVector indicating the timestamp of each event
//' @param weights a cube with the weights associated with each event for each
//' bws_net and bws_time.
//' @param bws_net an arma::cube of three dimensions with the network bandwidths calculated for each observation for each global time and network bandwidths
//' @param bws_time an arma::cube of three dimensions with the time bandwidths calculated for each observation for each global time and network bandwidths
//' @param kernel_name a string with the name of the kernel to use
//' @param line_list a DataFrame describing the lines
//' @param max_depth the maximum recursion depth
//' @param min_tol a double indicating by how much 0 in density values must be replaced
//' @return a matrix with the CV score for each pair of global bandiwdths
//' @export
//' @examples
//' # no example provided, this is an internal function
// [[Rcpp::export]]
 arma::cube tnkde_get_loo_values2(std::string method, List &neighbour_list,
                                 IntegerVector &sel_events,
                                 IntegerVector &sel_events_wid,
                                 NumericVector &sel_events_time,
                                 IntegerVector &events,
                                 IntegerVector &events_wid,
                                 NumericVector &events_time,
                                 arma::cube &weights,
                                 arma::cube &bws_net,
                                 arma::cube &bws_time,
                                 std::string kernel_name,
                                 DataFrame &line_list, int max_depth, double min_tol){

   //selecting the kernel function
   //Rcout << "print stage 1 ... \n";
   fptros kernel_func = select_kernelos(kernel_name);

   //step0 extract the columns of the dataframe
   NumericVector line_weights = line_list["weight"];

   //step 1 : mettre toutes les valeurs a 0 (bw_net 0 et bw_time 0)
   // NOTE WE calculate the values only for the events in sel_events
   arma::cube base_k(bws_net.n_rows, bws_time.n_cols, sel_events.length());

   //calculer la matrice des lignes
   //arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
   arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
   arma::cube k;
   // Rcout << "Here is the neighbours wids : \n";
   // Rcout << events_wid << "\n\n";
   // Rcout << "Here is the neighbours time : \n";
   // Rcout << events_time << "\n\n";

   // Rcout << "iterating to calculate the densities from each event \n";
   //step2 : iterer sur chaque event de la zone d'etude
   int cnt_e = events.length()-1;
   for(int i=0; i <= cnt_e; ++i){
     // Rcout << "    "<< i<<"\n";
     //preparer les differentes valeurs de departs pour l'event y
     int y = events[i];
     int wid = events_wid[i];
     double v_time = events_time[i];
     // launching recursion
     // here we got the the influences of the vertex y on each other selected event in quadra

     if(method == "simple"){
       k = ess_kernel_loo_tnkde_adpt(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net.slice(i), bws_time.slice(i),
                                line_weights, max_depth);

     }else if (method == "discontinuous"){
       // NOTE TO MYSELF HERE : If we are considering the densities caused by
       // the event y, then we only need the bandwidths of the event y
       k = esd_kernel_loo_tnkde_adpt(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net.slice(i), bws_time.slice(i),
                                line_weights, max_depth);
     }else{
       // NOTE TO MYSELF HERE : I MUST IMPLEMENT A VERSION OF THIS FUNCTION
       // TO WORK WITH ADAPTIVE BANDWIDTH (SEE THE ONE PREPARED ABOVE FOR
       // THE SIMPLE METHOD)
       k = esc_kernel_loo_tnkde_adpt(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net.slice(i), bws_time.slice(i),
                                line_weights, max_depth);
     }
     // NOTE : the scaling by bws is applied above
     arma::mat w = weights.slice(i);

     for(int ii = 0; ii < sel_events.length(); ii++){
       k.slice(ii) = k.slice(ii) % w;
     }

     // if y was a selected event, its own weight must be set to 0
     int index = get_first_index_int(sel_events_wid,wid);
     if(index >= 0){
       k.slice(index).fill(0);
     }

     // and summing the values at each iteration (% is the element wise product)
     base_k += k;
   };

   return base_k;

 }



//// DEVELOPMENT /////



//' @title The exposed function to calculate adaptive bandwidth with space-time
//' interaction for TNKDE (INTERNAL)
//' @name adaptive_bw_tnkde_cpp
//' @param method a string, one of "simple", "continuous", "discontinuous"
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param sel_events a Numeric vector indicating the selected events (id of nodes)
//' @param sel_events_wid a Numeric Vector indicating the unique if of the selected events
//' @param sel_events_time a Numeric Vector indicating the time of the selected events
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param events_wid a NumericVector indicating the unique id of all the events
//' @param events_time a NumericVector indicating the timestamp of each event
//' @param weights a cube with the weights associated with each event for each
//' bws_net and bws_time.
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param kernel_name a string with the name of the kernel to use
//' @param line_list a DataFrame describing the lines
//' @param max_depth the maximum recursion depth
//' @param min_tol a double indicating by how much 0 in density values must be replaced
//' @return a vector witht the estimated density at each event location
//' @export
//' @examples
//' # no example provided, this is an internal function
// [[Rcpp::export]]
arma::rowvec adaptive_bw_tnkde_cpp(std::string method,
                                List &neighbour_list,
                                IntegerVector &sel_events,
                                IntegerVector &sel_events_wid,
                                NumericVector &sel_events_time,
                                IntegerVector &events,
                                IntegerVector &events_wid,
                                NumericVector &events_time,
                                arma::vec &weights,
                                arma::vec &bws_net,
                                arma::vec &bws_time,
                                std::string kernel_name,
                                DataFrame &line_list,
                                int max_depth,
                                double min_tol){

  //selecting the kernel function
  fptros kernel_func = select_kernelos(kernel_name);

  //step0 extract the columns of the dataframe
  NumericVector line_weights = line_list["weight"];

  //step 1 : mettre toutes les valeurs a 0 (bw_net 0 et bw_time 0)
  // NOTE WE calculate the values only for the events in sel_events
  arma::cube base_k(bws_time.n_elem, bws_net.n_elem, sel_events.length());

  //calculer la matrice des lignes
  //arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
  arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
  arma::cube k;
  //step2 : iterer sur chaque event de la zone d'etude
  int cnt_e = events.length()-1;
  for(int i=0; i <= cnt_e; ++i){
    //preparer les differentes valeurs de departs pour l'event y
    int y = events[i];
    int wid = events_wid[i];
    double v_time = events_time[i];
    // launching recursion
    // here we got the the influences of the vertex y on each other selected event in quadra
    if(method == "simple"){
      k = ess_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               sel_events_time, neighbour_list,
                               y, wid, v_time,
                               bws_net, bws_time,
                               line_weights, max_depth);

    }else if (method == "discontinuous"){
      k = esd_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               sel_events_time, neighbour_list,
                               y, wid, v_time,
                               bws_net, bws_time,
                               line_weights, max_depth);
    }else{
      k = esc_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                               sel_events_time, neighbour_list,
                               y, wid, v_time,
                               bws_net, bws_time,
                               line_weights, max_depth);
    }
    // NOTE : the scaling by bws is applied above
    float w = weights(i);
    for(int ii = 0; ii < sel_events.length(); ii++){
      k.slice(ii) = k.slice(ii) * w;
    }

    // and summing the values at each iteration (% is the element wise product)
    base_k += k;
  };
  // and calculate the final values

  arma::uvec neg_elems = arma::find(base_k <= 0);
  base_k.elem(neg_elems).fill(min_tol);

  arma::rowvec result = base_k(arma::span(0), arma::span(0), arma::span::all);

  return result;
}


//' @title The exposed function to calculate adaptive bandwidth with space-time
//' interaction for TNKDE (INTERNAL)
//' @name adaptive_bw_tnkde_cpp2
//' @param method a string, one of "simple", "continuous", "discontinuous"
//' @param neighbour_list a List, giving for each node an IntegerVector with
//' its neighbours
//' @param sel_events a Numeric vector indicating the selected events (id of nodes)
//' @param sel_events_wid a Numeric Vector indicating the unique if of the selected events
//' @param sel_events_time a Numeric Vector indicating the time of the selected events
//' @param events a NumericVector indicating the nodes in the graph being events
//' @param events_wid a NumericVector indicating the unique id of all the events
//' @param events_time a NumericVector indicating the timestamp of each event
//' @param weights a cube with the weights associated with each event for each
//' bws_net and bws_time.
//' @param bws_net an arma::vec with the network bandwidths to consider
//' @param bws_time an arma::vec with the time bandwidths to consider
//' @param kernel_name a string with the name of the kernel to use
//' @param line_list a DataFrame describing the lines
//' @param max_depth the maximum recursion depth
//' @param min_tol a double indicating by how much 0 in density values must be replaced
//' @return a vector with the estimated density at each event location
//' @export
//' @examples
//' # no example provided, this is an internal function
// [[Rcpp::export]]
 arma::cube adaptive_bw_tnkde_cpp2(std::string method,
                                    List &neighbour_list,
                                    IntegerVector &sel_events,
                                    IntegerVector &sel_events_wid,
                                    NumericVector &sel_events_time,
                                    IntegerVector &events,
                                    IntegerVector &events_wid,
                                    NumericVector &events_time,
                                    arma::vec &weights,
                                    arma::vec &bws_net,
                                    arma::vec &bws_time,
                                    std::string kernel_name,
                                    DataFrame &line_list,
                                    int max_depth,
                                    double min_tol){
   // Rcout << "step0\n";
   //selecting the kernel function
   fptros kernel_func = select_kernelos(kernel_name);
   // Rcout << "step1\n";
   //step0 extract the columns of the dataframe
   NumericVector line_weights = line_list["weight"];
   // Rcout << "step2\n";
   //step 1 : mettre toutes les valeurs a 0 (bw_net 0 et bw_time 0)
   // NOTE WE calculate the values only for the events in sel_events
   arma::cube base_k(bws_net.n_elem, bws_time.n_elem, sel_events.length());
   // Rcout << "step3\n";
   //calculer la matrice des lignes
   //arma::sp_mat edge_mat = make_matrix_sparse(line_list,neighbour_list);
   arma::sp_imat edge_mat = make_imatrix_sparse(line_list, neighbour_list);
   arma::cube k;
   // Rcout << "step3\n";
   //step2 : iterer sur chaque event de la zone d'etude
   int cnt_e = events.length()-1;
   for(int i=0; i <= cnt_e; ++i){
     // Rcout << "Iterating on " <<  i<<"\n";
     //preparer les differentes valeurs de departs pour l'event y
     int y = events[i];
     int wid = events_wid[i];
     double v_time = events_time[i];
     // launching recursion
     // here we got the the influences of the vertex y on each other selected event in quadra
     if(method == "simple"){
       k = ess_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);

     }else if (method == "discontinuous"){
       k = esd_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);
     }else{
       k = esc_kernel_loo_tnkde(kernel_func, edge_mat, sel_events,sel_events_wid,
                                sel_events_time, neighbour_list,
                                y, wid, v_time,
                                bws_net, bws_time,
                                line_weights, max_depth);
     }
     // Rcout << "Here is k : \n";
     // Rcout << k <<"\n\n";
     // NOTE : the scaling by bws is applied above
     float w = weights(i);
     for(int ii = 0; ii < sel_events.length(); ii++){
       k.slice(ii) = k.slice(ii) * w;
     }
     // Rcout << "weight applied !\n";
     // and summing the values at each iteration (% is the element wise product)
     base_k += k;
   };
   // and calculate the final values
   // Rcout << "end of iterations !\n";
   arma::uvec neg_elems = arma::find(base_k <= 0);
   base_k.elem(neg_elems).fill(min_tol);

   return base_k;
 }




