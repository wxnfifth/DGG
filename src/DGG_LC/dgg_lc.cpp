
#include <string>
#include <ctime>
#include <random>
#include <windows.h>
#include <psapi.h>
#include <stdio.h>   
#include <tchar.h>
#include "wxn\wxnTime.h"
#include "wxn\wxn_dijstra.h"
#include "svg_definition.h"

using namespace std;




std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

int main(int argc, char** argv)
{
  bool test_performance = false;
  string obj_file_name;
  string svg_file_name;

  if (argc < 4) {
    printf("parameter insufficient !\n usage SVG_LC.exe [obj_file_name] [svg_file_name] [dij  or lll or fim or hy]\n");
    return 1;
  }
  obj_file_name = argv[1];

  CRichModel rich_model(obj_file_name);
  rich_model.Preprocess();

  svg_file_name = argv[2];
  string method = argv[3];
  SparseGraph<float>* s_graph = NULL;

  if (method == "dij") {
    s_graph = new Dijstra_vector<float>();
    printf("size before %d\n" , sizeof(*s_graph));
    s_graph->read_svg_file_with_angle((string)svg_file_name);
    printf("size after %d\n" , sizeof(*s_graph));
  } else if (method == "lll") {
    s_graph = new LC_LLL<float>();
    s_graph->read_svg_file_with_angle((string)svg_file_name);//load the precomputed infomation 
  } else if(method == "fim") {
    s_graph = new LC_FIM<float>();//temp modify
    s_graph->read_svg_file_binary((string)svg_file_name);//load the precomputed infomation 
  }  else if(method == "hy") {
    s_graph = new LC_HY<float>();
    s_graph->read_svg_file_with_angle((string)svg_file_name);
    //double vg_eps;
    //if (argc < 5) {
    //  printf("xx.exe xx.obj xx.binary hy 0.01\n");
    //  exit(1);
    //}
    //sscanf(argv[4],"%lf",&vg_eps);
    //double theta = asin(sqrt(vg_eps));
    //s_graph->theta = theta * 5;
    //printf("******** theta %lf du\n" , s_graph->theta  / M_PI * 180.0);
    s_graph->model_ptr = &rich_model;
  }
  else{
    printf("invalid choice\n"); 
    return 1;
  }

  vector<pair<double,double>> erros_list;//<dis,error>
  std::mt19937 rng;
  std::uniform_int_distribution<int> uint_dist(0,s_graph->NodeNum()-1);
  double average_time=0;
  int iteration_times = 10;
  double max_dis = 0;
  for (int itr = 0; itr < iteration_times; ++itr) 
  {

    int source_vert = uint_dist(rng);//index of source vertex
    //double t0 = clock();
    ElapasedTime t;
    s_graph->findShortestDistance(source_vert);  //find the geodesic distance from a single source to all vertex of a mesh
    //printf("t0 %lf\n" , (clock() - t0)/1000.0);
    average_time += t.getTime();
    t.printTime("time");


    char buf[256];
    sprintf(buf, "c:\\util\\XW_geodesic_to_file.exe %s %d" , obj_file_name.c_str(), source_vert);
    system(buf);
    string exact_file_name = obj_file_name.substr(0, obj_file_name.length()-4) + "_exact_distance.txt";
    //here is an example to output the distance field to a txt file
    FILE* correct_dis_file = fopen(exact_file_name.c_str(), "r");
    if (correct_dis_file == NULL) {
      printf("file not found\n");
      return 0;
    }
    vector<double> correct_dis(s_graph->NodeNum());
    for (int i = 0; i < s_graph->NodeNum(); ++i) {
      fscanf(correct_dis_file, "%lf", &correct_dis[i]);
    }
    fclose(correct_dis_file);

    max_dis = std::max(max_dis,*std::max_element(correct_dis.begin(), correct_dis.end()));


    FILE* out_file = fopen("svg_dis.txt","w");
    for (int i = 0; i < s_graph->NodeNum(); ++i) {
      double geodesic_distance = s_graph->distanceToSource(i);//get distance from source vert to dest vert
      fprintf(out_file, "%lf\n" , geodesic_distance);
    }
    fclose(out_file);



    //double max_error=0;
    //double ave_error=0;
    int cnt_error_dis = 0;
    for (int i = 0; i < correct_dis.size(); ++i) {
      if (fabs(correct_dis[i]) < 1e-6 ) continue;
      double error = fabs(s_graph->distanceToSource(i) - correct_dis[i])/correct_dis[i];
      //max_error = std::max(error,max_error);
      //ave_error += error;
      //if (error < 1 && error > 0.9) {
      //  printf("************source %d dest num %d\n" , source_vert, i);
      //  vector<int> path_nodes;
      //  s_graph->getPath(i,path_nodes);
      //  printf("path\n");
      //  for (auto& v:path_nodes){
      //    printf("** %d " , v);
      //  }
      //  printf("\n");
      //}
      if (error > 1) {
        cnt_error_dis ++;
        continue; 
      }
      erros_list.push_back(make_pair(correct_dis[i],error));
    }
    printf("________error dis percent %.2lf\n" , (double)cnt_error_dis / correct_dis.size());
  }

  //ave_error /= correct_dis.size();
  ////printf("max_error %lf\n" , max_error);
  //printf("ave_error %lf\n" , ave_error);


  double longest_dis = max_dis;
  int interval_num = 100;
  vector<double> average_error_sep(interval_num+1, 0.0);
  vector<double> average_error_cnt(interval_num+1, 0);
  double all_average_error = 0;

  for (int i = 0; i < erros_list.size(); ++i) {
    double dis = erros_list[i].first;
    double error = erros_list[i].second;
    double percent =  dis / longest_dis;
    int pos = percent * interval_num;
    average_error_sep[pos] += error;
    average_error_cnt[pos]++;
    all_average_error += erros_list[i].second;
  }

  for (int i = 0; i < average_error_sep.size(); ++i) {
    if (average_error_cnt[i] != 0) {
      printf("dis percent %d to %d: average error %.10lf\n" , i , i + 1 , average_error_sep[i] / average_error_cnt[i]);
    }else{
       printf("dis percent %d to %d: average error non\n" , i , i + 1 );
    }
  }
  printf("total_average_running_time %lf\n" , average_time / iteration_times);
  printf("total_average_error %.10lf\n" , all_average_error / erros_list.size());

  if (method == "hy") {
    int cnt = 0;
    
    auto& lst = split(svg_file_name, '_');
    string eps_str = lst[lst.size()-3].substr(2);
    printf("eps %s\n",  eps_str.c_str());
    double eps = atof(eps_str.c_str());
    for (auto& e:erros_list) {
      if (e.second > eps) {
        cnt++;
      }
    }
    printf("percenter %lf\n" , (double)cnt / erros_list.size());
  }

   HANDLE current_process = GetCurrentProcess();
   PROCESS_MEMORY_COUNTERS pmc;

   if ( GetProcessMemoryInfo( current_process, &pmc, sizeof(pmc)) )
   {
     printf("PeakWorkingSetSize: %d M\n" , pmc.PeakWorkingSetSize / 1024 / 1024);
   }

  delete s_graph;

 /* system("pause");
 */ return 0;
}