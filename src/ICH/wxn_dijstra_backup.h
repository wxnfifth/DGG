#ifndef _WXN_GRAPH_H_
#define _WXN_GRAPH_H_

#include<iostream>
#include<vector>
#include<string>
#include <queue>
#include <algorithm>
#include <list>
#include <map>
#include <fstream>
#include <assert.h>

#include "RichModel.h"

template<typename T>  class SparseGraph{
private:
  struct HeadOfSVG {
    int begin_vertex_index;
    int end_vertex_index;
    int num_of_vertex;
    HeadOfSVG(int _begin_vertex_index , 
      int _end_vertex_index , 
      int _num_of_vertex)
      :
    begin_vertex_index(_begin_vertex_index) ,
      end_vertex_index(_end_vertex_index) , 
      num_of_vertex(_num_of_vertex){}
    HeadOfSVG(){}
    void print(){
      printf("head of svg : num_of_vertex %d\n" , num_of_vertex);
    }
  };

  struct BodyHeadOfSVG{
    int source_index;
    int neighbor_num;
    BodyHeadOfSVG(int _source_index , int _neighbor_num)
    {
      source_index = _source_index;
      neighbor_num = _neighbor_num;
    }
    BodyHeadOfSVG(){}
    void print(){
      printf("source_index %d\n" , source_index);
      printf("neighbor_num %d\n" , neighbor_num);
    }
  };

  struct BodyPartOfSVG{
    int dest_index;
    float dest_dis;
    BodyPartOfSVG(){}
    BodyPartOfSVG(int _dest_index , float _dest_dis):
      dest_index(_dest_index),
      dest_dis(_dest_dis),
    {
      dest_index = _dest_index;
      dest_dis = _dest_dis;
    }
  };

  struct BodyPartOfSVGWithAngle{
    int dest_index;
    float dest_dis;
    float angle;
    int begin_pos;
    int end_pos;
    BodyPartOfSVGWithAngle(){}
    BodyPartOfSVGWithAngle(int _dest_index , float _dest_dis, float _angle, int _begin_pos, int _end_pos):
      dest_index(_dest_index),
      dest_dis(_dest_dis),
      angle(_angle),
      begin_pos(_begin_pos),
      end_pos(_end_pos)
    {
    }
    bool operator<(const BodyPartOfSVGWithAngle& o) const{
      return angle < o.angle;
    }
  };

  struct BodyPartOfSVGWithK : BodyPartOfSVG{
    int rank_k;
    BodyPartOfSVGWithK(){}
    BodyPartOfSVGWithK(int _dest_index , float _dest_dis , int _rank_k):
      BodyPartOfSVG(_dest_index , _dest_dis),
      rank_k(_rank_k){}
    bool operator<(const BodyPartOfSVGWithK& other)const{
      return rank_k < other.rank_k;
    }
  };


protected:
  struct GraphEdge {
    int neigh;
    T dis;
    GraphEdge(){}
    GraphEdge(int _neigh, T _dis) :neigh(_neigh),dis(_dis){}
    //bool operator<(const GraphEdge& other)const{
    //  return dis > other.dis;
    //}
  };


  int node_number_;

  std::vector<GraphEdge> graph_edge;
  std::vector<int> graph_edge_degree;


  std::vector<std::vector<int>> graph_neighbor;
  std::vector<std::vector<int>> graph_pos_in_neighbor;
  std::vector<std::vector<T>> graph_neighbor_dis;
  std::vector<std::vector<T>> graph_neighbor_angle;
  std::vector<std::vector<pair<unsigned short,unsigned short>>> graph_neighbor_begin_end_pos;
  std::vector<std::map<int,int>> graph_neighbor_map;


public:
  //double theta;
  CRichModel* model_ptr;

public:
  SparseGraph(){
    node_number_ = -1;
  }

  void initialize(int _node_number) {
    node_number_ = _node_number;
    graph_neighbor.reserve(_node_number);
    graph_neighbor.resize(_node_number);
    graph_pos_in_neighbor.reserve(_node_number);
    graph_pos_in_neighbor.resize(_node_number);
    graph_neighbor_dis.reserve(_node_number);
    graph_neighbor_dis.resize(_node_number);
    graph_neighbor_angle.reserve(_node_number);
    graph_neighbor_angle.resize(_node_number);
    graph_neighbor_begin_end_pos.reserve(_node_number);
    graph_neighbor_begin_end_pos.resize(_node_number);
    graph_neighbor_map.reserve(_node_number);
    graph_neighbor_map.resize(_node_number);



  }

  void allocate_for_neighbor(int u , int number_of_neighbor) {
    graph_neighbor[u].reserve(number_of_neighbor);
    graph_pos_in_neighbor[u].reserve(number_of_neighbor);
    graph_neighbor_dis[u].reserve(number_of_neighbor);
    graph_neighbor_angle[u].reserve(number_of_neighbor);
    graph_neighbor_begin_end_pos[u].reserve(number_of_neighbor);
  }
  void addedge(int u , int v , T w) {
    //u , v is the two edge
    // w is the distance
    assert(u < node_number_ && v < node_number_);
    graph_neighbor[u].push_back(v);
    graph_neighbor_dis[u].push_back(w);
  }

  void addedge(int u , int v , T w, T angle, int begin_pos , int end_pos) {
    //u , v is the two edge
    // w is the distance
    assert(u < node_number_ && v < node_number_);
    graph_neighbor[u].push_back(v);
    graph_neighbor_dis[u].push_back(w);
    graph_neighbor_map[u][v] = graph_neighbor_angle[u].size();
    graph_neighbor_angle[u].push_back(angle);
    graph_neighbor_begin_end_pos[u].push_back(make_pair(begin_pos,end_pos));
  }

  int NodeNum() {
    return node_number_;
  }


  int read_svg_file_binary_new(const std::string& svg_file_name) {
    read_svg_file_binary(svg_file_name);
    //    std::vector<Edge> graph_edge;
    //std::vector<int> graph_edge_pos;
    int total_edge_num = 0;
    for (auto &neigh:graph_neighbor) {
      total_edge_num += neigh.size();
    }
    graph_edge.reserve(total_edge_num+10);
    graph_edge.resize(total_edge_num+10);
    graph_edge_degree.resize(node_number_+1);
    graph_edge_degree[0] = 0;
    for (int i = 0; i < node_number_; ++i) {
      graph_edge_degree[i+1] = graph_edge_degree[i] + graph_neighbor[i].size();
      for (int j = graph_edge_degree[i]; j < graph_edge_degree[i+1]; ++j) {
        int pos = j-graph_edge_degree[i];
        graph_edge[j] = GraphEdge(graph_neighbor[i][pos], graph_neighbor_dis[i][pos]);
      }
    }
  }



  int read_svg_file_binary(const std::string& svg_file_name) {

    std::ifstream input_file (svg_file_name, std::ios::in | std::ios::binary);
    HeadOfSVG head_of_svg;
    input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
    head_of_svg.print();
    initialize(head_of_svg.num_of_vertex);

    double average_neighbor_number(0.0);
    double max_radius_in_file(0.0);
    double average_radius(0.0);

    for (int i = 0; i < head_of_svg.num_of_vertex;++i) {
      BodyHeadOfSVG body_head;
      input_file.read( (char*)&body_head , sizeof(body_head));
      average_neighbor_number += (double)body_head.neighbor_num;
      std::vector<BodyPartOfSVG> body_parts;
      for(int j = 0; j < body_head.neighbor_num;++j){
        BodyPartOfSVG body_part;
        input_file.read((char*)&body_part , sizeof(body_part));
        body_parts.push_back(body_part);
      }
      allocate_for_neighbor(body_head.source_index , body_parts.size());
      for (int j = 0; j < body_parts.size();++j) {
        addedge(body_head.source_index ,
          body_parts[j].dest_index ,
          body_parts[j].dest_dis);
      }
      if( i > 0 && i % (head_of_svg.num_of_vertex / 10 ) == 0 ){
        std::cerr << "read " << i * 100 / head_of_svg.num_of_vertex  << " percent \n";
      }
    }

    average_neighbor_number /= head_of_svg.num_of_vertex;
    printf("average neigh %lf\n" , average_neighbor_number);
    std::cerr << "reading done..\n";
    input_file.close();

    return 0;
  }

  int read_svg_file_with_angle(const std::string& svg_file_name) {

    std::ifstream input_file (svg_file_name, std::ios::in | std::ios::binary);
    HeadOfSVG head_of_svg;
    input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
    head_of_svg.print();
    initialize(head_of_svg.num_of_vertex);

    double average_neighbor_number(0.0);
    double max_radius_in_file(0.0);
    double average_radius(0.0);

    for (int i = 0; i < head_of_svg.num_of_vertex;++i) {
      BodyHeadOfSVG body_head;
      input_file.read( (char*)&body_head , sizeof(body_head));
      average_neighbor_number += (double)body_head.neighbor_num;
      std::vector<BodyPartOfSVGWithAngle> body_parts;
      for(int j = 0; j < body_head.neighbor_num;++j){
        BodyPartOfSVGWithAngle body_part;
        input_file.read((char*)&body_part , sizeof(body_part));
        body_parts.push_back(body_part);
      }
      allocate_for_neighbor(body_head.source_index , body_parts.size());
      for (int j = 0; j < body_parts.size();++j) {
        addedge(body_head.source_index ,
          body_parts[j].dest_index ,
          body_parts[j].dest_dis,
          body_parts[j].angle,
          body_parts[j].begin_pos,
          body_parts[j].end_pos
          );
      }
      if( i > 0 && i % (head_of_svg.num_of_vertex / 10 ) == 0 ){
        std::cerr << "read " << i * 100 / head_of_svg.num_of_vertex  << " percent \n";
      }
    }

    average_neighbor_number /= head_of_svg.num_of_vertex;
    printf("average neigh %lf\n" , average_neighbor_number);
    std::cerr << "reading done..\n";
    input_file.close();

    for (int i = 0; i < graph_pos_in_neighbor.size(); ++i) {
      graph_pos_in_neighbor[i].resize(graph_neighbor[i].size());
      for (int j = 0; j < graph_neighbor[i].size(); ++j) {
        int neigh = graph_neighbor[i][j];
        int pos = graph_neighbor_map[neigh][i];
        graph_pos_in_neighbor[i][j] = pos;
      }
    }

    return 0;
  }


  virtual void findShortestDistance(int source)=0;
  //virtual int getSource(int v)=0;
  virtual T distanceToSource(int index)=0;
  virtual void getPath(const int& dest , std::vector<int>& path_nodes)=0;

};


template <typename T> class LC_FIM:public SparseGraph<T> {
  enum LabelType { SeedPoint, ActivePoint, FarPoint, DeadPoint };

  std::vector<T> dis_;
  std::vector<LabelType> label_;

  T distanceToSource(int index) {
    return dis_[index];
  }

  void getPath(const int& dest , std::vector<int>& path_nodes){

      printf("error in line 248!\n");

  }


  T upwind(int u) {
    T new_dis = FLT_MAX;
    for (int i = 0; i < graph_neighbor[u].size(); ++i) {
      int v = graph_neighbor[u][i];
      T d = dis_[v] + graph_neighbor_dis[u][i];
      if (new_dis > d) {
        new_dis = d;
      }
    }
    return new_dis;
  }

  void findShortestDistance(int source)
  {
    dis_.resize(node_number_);
    fill(dis_.begin(), dis_.end(), FLT_MAX);
    label_.resize(node_number_);
    fill(label_.begin(), label_.end(), FarPoint);

    label_[source] = SeedPoint;
    dis_[source] = 0;
    vector<bool> label_set(node_number_,false);
    label_set[source] = true;

    list<int> active_list;
    for (int i = 0; i < graph_neighbor[source].size(); ++i) {
      active_list.push_back(graph_neighbor[source][i]);
      label_[graph_neighbor[source][i]] = ActivePoint;
    }
    int times = 0;

    while( !active_list.empty()) {
      //printf("active list size %d\n" , active_list.size());
      auto itr = active_list.begin();
      while (itr != active_list.end()) {
        int tmp_index = *itr;
        T old_dis = dis_[tmp_index];
        T new_dis = upwind(tmp_index);
        if (fabs(old_dis - new_dis) < 1e-5) {
          if (old_dis > new_dis) {
            dis_[tmp_index] = new_dis;
          }
          for (int i = 0; i < graph_neighbor[tmp_index].size();++i) {
            int v = graph_neighbor[tmp_index][i];

            if (label_[v] == FarPoint) {
              T old_dis_v = dis_[v];
              T new_dis_v = upwind(v);
              if (true) {  
                if (old_dis_v > new_dis_v) {
                  dis_[v] = new_dis_v;
                }
                if (label_[v] != ActivePoint) {
                  active_list.insert(itr, v);
                  label_[v] = ActivePoint;
                }
              }
            }
          }
          itr = active_list.erase(itr);
          label_[tmp_index] = DeadPoint;
        } else {
          if (old_dis > new_dis) {
            dis_[tmp_index] = new_dis;
          }
        }
      }

    }
     
  }

  
};

template <typename T> class LC_LLL:public SparseGraph<T>{

private:

  std::vector<T> dis;
  std::vector<bool> visited;
  std::vector<int> fathers;

public:
  LC_LLL(){}

  void getPath(const int& dest , std::vector<int>& path_nodes) {
    path_nodes.clear();
    int u = dest;
    while( u != fathers[u] ){
      path_nodes.push_back(u);
      u = fathers[u];
    }
    path_nodes.push_back(u);
    std::reverse(path_nodes.begin() , path_nodes.end());
  }

  
  void findShortestDistance(int source)
  {
    fathers.resize(node_number_);
    fill(fathers.begin(),fathers.end(),-1);
    dis.resize(node_number_);
    fill(dis.begin(), dis.end(), FLT_MAX);
    visited.resize(node_number_);
    fill(visited.begin(), visited.end(), false);
    std::deque<int> que;
    double dis_sum = 0.0;
    dis[source] = 0;
    que.push_back(source);
    fathers[source] = source;
    visited[source] = true;

    int processed_ele_ment =0;
    int total_neigh = 0;
    while (!que.empty()) {
      int u = -1;
      processed_ele_ment++;
      while (true) {

        if (dis[que.front()] > dis_sum / que.size()) {
          que.push_back(que.front());
          que.pop_front();
        } else {
          u = que.front();
          que.pop_front();
          dis_sum -= dis[u];
          visited[u] = false;
          break;
        }
      }
      for (int i = 0; i < graph_neighbor[u].size(); ++i) {
        if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
          total_neigh++;
          int v = graph_neighbor[u][i];
          T w = graph_neighbor_dis[u][i];
          double old_v = dis[v];
          dis[v] = dis[u] + w;
          fathers[v] = u;
          if (visited[v] == false) {
            que.push_back(v);
            visited[v] = true;
            dis_sum += dis[v];
          } else {
            dis_sum -= old_v;
            dis_sum += dis[v];
          }
        }
      }
    }
        printf("********* processed_element %d total_vert %d percenter %.2lf%%\n" , processed_ele_ment, node_number_, processed_ele_ment / (double) node_number_ * 100.0); 
        printf("********* total_neigh %d total_vert %d percenter %.2lf%%\n" , total_neigh, node_number_, total_neigh / (double) node_number_ * 100.0); 


  }

  inline int getSource(int v){
    if( fathers[v] == v ){
      return v;
    }else{
      fathers[v] = getSource(fathers[v]);
      return fathers[v];
    }
  }

  inline T distanceToSource(int index) {
    if(index < 0 || index >= node_number_ ){
      std::cerr << "wrong index " << index << "\n";
      return 0;
    }
    return dis[index];
  }

};

template <typename T> class LC_HY:public SparseGraph<T>{

private:

  std::vector<T> dis;
  std::vector<bool> visited;
  std::vector<int> fathers;

public:
  LC_HY(){}

  void getPath(const int& dest , std::vector<int>& path_nodes) {
    path_nodes.clear();
    int u = dest;
    while( u != fathers[u] ){
      path_nodes.push_back(u);
      u = fathers[u];
    }
    path_nodes.push_back(u);
    std::reverse(path_nodes.begin() , path_nodes.end());
  }
  
  void findShortestDistance(int source)
  {
    printf("???\n");
    vector<int> father_in_neighbor_pos(node_number_);

    fathers.resize(node_number_);
    fill(fathers.begin(),fathers.end(),-1);
    dis.resize(node_number_);
    fill(dis.begin(), dis.end(), FLT_MAX);
    visited.resize(node_number_);
    fill(visited.begin(), visited.end(), false);
    std::deque<int> que;
    double dis_sum = 0.0;
    dis[source] = 0;
    que.push_back(source);
    fathers[source] = source;
    visited[source] = true;

    int processed_ele_ment = 0;
    int total_neigh = 0;
    while (!que.empty()) {
      processed_ele_ment++;
      int u = -1;
      while (true) {

        if (dis[que.front()] > dis_sum / que.size()) {
          que.push_back(que.front());
          que.pop_front();
        } else {
          u = que.front();
          que.pop_front();
          dis_sum -= dis[u];
          visited[u] = false;
          break;
        }
      }


      int father_u = fathers[u];
      if (father_u == u) {
        for (int i = 0; i < graph_neighbor[u].size(); ++i) {
          if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
            int v = graph_neighbor[u][i];
            T w = graph_neighbor_dis[u][i];
            double old_v = dis[v];
            dis[v] = dis[u] + w;
            fathers[v] = u;
            father_in_neighbor_pos[v] = graph_pos_in_neighbor[u][i];
            if (visited[v] == false) {
              que.push_back(v);
              visited[v] = true;
              dis_sum += dis[v];
            } else {
              dis_sum -= old_v;
              dis_sum += dis[v];
            }
          }
        }
      } else {
        //bool is_convex_vert = model_ptr->IsConvexVert(u);
        ////if (is_convex_vert) {
        ////  continue;
        ////}
        //T angle_sum = model_ptr->AngleSum(u);
        //T father_angle = graph_neighbor_angle_map[u][father_u];
        //int father_pos = graph_neighbor_map[u][father_u];
       int father_pos = father_in_neighbor_pos[u];
        // printf("father_pos %d pos_in_neigh %d\n" , father_pos, father_in_neighbor_pos[u]); 
        //T father_angle = graph_neighbor_angle[u][father_pos];

        //based on father_angle as 0 
        //T start_angle = M_PI - theta;
        //T end_angle = angle_sum - (M_PI - theta);
        //printf("******** theta %lf du\n" , theta / M_PI * 180.0);
        int begin_pos = graph_neighbor_begin_end_pos[u][father_pos].first;
        int end_pos = graph_neighbor_begin_end_pos[u][father_pos].second;

        if (begin_pos == -1) continue;
        if (begin_pos <= end_pos) {
          for (int i = begin_pos; i <= end_pos; ++i) {
            if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
              total_neigh++;
              int v = graph_neighbor[u][i];
              T w = graph_neighbor_dis[u][i];
              double old_v = dis[v];
              dis[v] = dis[u] + w;
              fathers[v] = u;
              father_in_neighbor_pos[v] = graph_pos_in_neighbor[u][i];
              if (visited[v] == false) {
                que.push_back(v);
                visited[v] = true;
                dis_sum += dis[v];
              } else {
                dis_sum -= old_v;
                dis_sum += dis[v];
              }
            }
          }
        } else {
          for (int i = begin_pos; i < graph_neighbor[u].size(); ++i) {
            if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
              total_neigh++;
              int v = graph_neighbor[u][i];
              T w = graph_neighbor_dis[u][i];
              double old_v = dis[v];
              dis[v] = dis[u] + w;
              fathers[v] = u;
              father_in_neighbor_pos[v] = graph_pos_in_neighbor[u][i];
              if (visited[v] == false) {
                que.push_back(v);
                visited[v] = true;
                dis_sum += dis[v];
              } else {
                dis_sum -= old_v;
                dis_sum += dis[v];
              }
            }
          }
          for (int i = 0; i <= end_pos; ++i) {
            if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
              total_neigh++;
              int v = graph_neighbor[u][i];
              T w = graph_neighbor_dis[u][i];
              double old_v = dis[v];
              dis[v] = dis[u] + w;
              fathers[v] = u;
              father_in_neighbor_pos[v] = graph_pos_in_neighbor[u][i];
              if (visited[v] == false) {
                que.push_back(v);
                visited[v] = true;
                dis_sum += dis[v];
              } else {
                dis_sum -= old_v;
                dis_sum += dis[v];
              }
            }
          }
        }


      }
    }
    printf("********* processed_element %d total_vert %d percenter %.2lf%%\n" , processed_ele_ment, node_number_, processed_ele_ment / (double) node_number_ * 100.0); 
    printf("********* total_neigh %d total_vert %d percente %.2lf%%\n" , total_neigh, node_number_, total_neigh / (double) node_number_ * 100.0);
  }

  //  void findShortestDistance_backup(int source)
  //{
  //  printf("???\n");
  //  fathers.resize(node_number_);
  //  fill(fathers.begin(),fathers.end(),-1);
  //  dis.resize(node_number_);
  //  fill(dis.begin(), dis.end(), FLT_MAX);
  //  visited.resize(node_number_);
  //  fill(visited.begin(), visited.end(), false);
  //  std::deque<int> que;
  //  double dis_sum = 0.0;
  //  dis[source] = 0;
  //  que.push_back(source);
  //  fathers[source] = source;
  //  visited[source] = true;
  //  int processed_ele_ment = 0;
  //  int total_neigh = 0;
  //  while (!que.empty()) {
  //    processed_ele_ment++;
  //    int u = -1;
  //    while (true) {
  //      if (dis[que.front()] > dis_sum / que.size()) {
  //        que.push_back(que.front());
  //        que.pop_front();
  //      } else {
  //        u = que.front();
  //        que.pop_front();
  //        dis_sum -= dis[u];
  //        visited[u] = false;
  //        break;
  //      }
  //    }
  //    int father_u = fathers[u];
  //    if (father_u == u) {
  //      for (int i = 0; i < graph_neighbor[u].size(); ++i) {
  //        if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
  //          int v = graph_neighbor[u][i];
  //          T w = graph_neighbor_dis[u][i];
  //          double old_v = dis[v];
  //          dis[v] = dis[u] + w;
  //          fathers[v] = u;
  //          if (visited[v] == false) {
  //            que.push_back(v);
  //            visited[v] = true;
  //            dis_sum += dis[v];
  //          } else {
  //            dis_sum -= old_v;
  //            dis_sum += dis[v];
  //          }
  //        }
  //      }
  //    } else {
  //      //bool is_convex_vert = model_ptr->IsConvexVert(u);
  //      ////if (is_convex_vert) {
  //      ////  continue;
  //      ////}
  //      T angle_sum = model_ptr->AngleSum(u);
  //      T father_angle = graph_neighbor_angle_map[u][father_u];
  //      //based on father_angle as 0
  //      T start_angle = M_PI - theta;
  //      T end_angle = angle_sum - (M_PI - theta);
  //      //printf("******** theta %lf du\n" , theta / M_PI * 180.0);
  //      for (int i = 0; i < graph_neighbor[u].size(); ++i) {
  //        double angle = graph_neighbor_angle[u][i];
  //        angle -= father_angle;
  //        if (angle < 0) {
  //          angle += angle_sum;
  //        }
  //        if (angle < start_angle || angle > end_angle) {
  //          continue;
  //        }
  //        
  //        if (dis[graph_neighbor[u][i]] > dis[u] + graph_neighbor_dis[u][i]) {
  //        total_neigh++;
  //        int v = graph_neighbor[u][i];
  //          T w = graph_neighbor_dis[u][i];
  //          double old_v = dis[v];
  //          dis[v] = dis[u] + w;
  //          fathers[v] = u;
  //          if (visited[v] == false) {
  //            que.push_back(v);
  //            visited[v] = true;
  //            dis_sum += dis[v];
  //          } else {
  //            dis_sum -= old_v;
  //            dis_sum += dis[v];
  //          }
  //        }
  //      }
  //    }
  //  }
  //  printf("********* processed_element %d total_vert %d percenter %.2lf%%\n" , processed_ele_ment, node_number_, processed_ele_ment / (double) node_number_ * 100.0); 
  //  printf("********* total_neigh %d total_vert %d percente %.2lf%%\n" , total_neigh, node_number_, total_neigh / (double) node_number_ * 100.0);
  //}


  inline int getSource(int v){
    if( fathers[v] == v ){
      return v;
    }else{
      fathers[v] = getSource(fathers[v]);
      return fathers[v];
    }
  }

  inline T distanceToSource(int index) {
    if(index < 0 || index >= node_number_ ){
      std::cerr << "wrong index " << index << "\n";
      return 0;
    }
    return dis[index];
  }

};

template <typename T> class  Dijstra_vector:public SparseGraph<T> {
public:
private:
    struct QueueNode{
        T dis;
        int node_index;
        QueueNode(){}
        QueueNode(int _node_index,T _dis){
            dis = _dis;
            node_index = _node_index;
        }
        bool operator<(const QueueNode& other)const{
            return dis > other.dis;
        }
    };
    std::vector<T> dis;
    std::vector<bool>visited;
    std::vector<int> fathers;

public:
  Dijstra_vector(){}

  void getPath(const int& dest , std::vector<int>& path_nodes){
    path_nodes.clear();
    int u = dest;
    while( u != fathers[u] ){
      path_nodes.push_back(u);
      u = fathers[u];
    }
    path_nodes.push_back(u);
    std::reverse(path_nodes.begin() , path_nodes.end());
  }

	void findShortestDistance(int source)
	{
		fathers.resize(node_number_);
		fill(fathers.begin(),fathers.end(),-1);
		dis.resize(node_number_);
    fill(dis.begin(),dis.end(),FLT_MAX);
		visited.resize(node_number_);
		fill(visited.begin(),visited.end(),false);
    priority_queue<QueueNode> que;
		dis[source] = 0;
    que.push(QueueNode(source,0));
    fathers[source] = source;
		int max_que_size = 0;

    while(!que.empty()){
			QueueNode u = que.top();
			max_que_size = max((int)que.size(), max_que_size);
			que.pop();
			if(visited[u.node_index]) continue;
			visited[u.node_index] = true;
		
        //for(int i = 0; i < graph_neighbor[u.node_index].size();++i){
      GraphEdge* edge_ptr = &graph_edge[graph_edge_degree[u.node_index]];
      GraphEdge* edge_ptr_end = &graph_edge[graph_edge_degree[u.node_index+1]];
      while (edge_ptr != edge_ptr_end) {
          int v = edge_ptr->neigh;//  graph_neighbor[u.node_index][i];
          //T w = graph_neighbor_dis[u.node_index][i];
          T w = edge_ptr->dis;
          if( !visited[v] && dis[v] > dis[u.node_index] + w) {
            dis[v] = dis[u.node_index] + w;
            que.push(QueueNode(v,dis[v]));
            fathers[v] = u.node_index;
          }
          edge_ptr++;
        }

		}
		//printf("max queue size %d\n" , max_que_size );
	}

	int getSource(int v){
		if( fathers[v] == v ){
			return v;
		}else{
			fathers[v] = getSource(fathers[v]);
			return fathers[v];
		}
	}



  inline T distanceToSource(int index){
    if(index < 0 || index >= node_number_ ){
      cerr << "wrong index " << index << "\n";
      return 0;
    }
    return dis[index];
  }


};


#endif