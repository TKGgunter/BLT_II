#include <stdio.h>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>

struct Tree{
		int max_depth;
		int n_nodes;
		std::vector<float> threshold;	
		std::vector<float> value0;
		std::vector<float> value1;
		std::vector<int>   child_right;
		std::vector<int>   child_left;
		std::vector<int> 	feature;
};


struct Forest{
		int n_trees;
		std::vector<Tree> trees;	
		std::vector<std::string> features;
};

float score_forest(std::vector<float> arr, Forest* f){
		 	float rt = 0.0;	
			for (int i = 0 ; i < f->n_trees; i++){

					Tree* tree = &(f->trees[i]);
					int cursor = 0;

					while (true){
							if ( cursor < 0 ) break;
							if ( tree->child_left[cursor] < 0 ){

									//Does is this really how scores are handled?
									//https://github.com/scikit-learn/scikit-learn/blob/7389dbac82d362f296dc2746f10e43ffa1615660/sklearn/tree/_tree.pyx
									rt +=  (float) tree->value0[cursor] / (float) (  tree->value1[cursor] + tree->value0[cursor] );
									break;
							}


							int c_feature = tree->feature[cursor];
							if ( arr[c_feature] > tree->threshold[cursor])	{
									cursor = tree->child_right[cursor];
							} else{
									cursor = tree->child_left[cursor];
							}
					}
			}
			return rt / (float) f->n_trees;
}




void print_tree(Tree* tree){
		printf("Tree: \n");
		printf("max depth: %d \n", tree->max_depth);
		printf("number of nodes: %d \n", tree->n_nodes);
	
		if ( tree->threshold.size() != tree->feature.size()){
				printf("This tree has a problem the length of the threshold and value0 arrays are not equal.\n");
				printf("threshold %d\n", (int) tree->threshold.size());
				printf("value0 %d\n", (int) tree->value0.size());
				printf("value1 %d\n", (int) tree->value1.size());
				printf("child right %d\n", (int) tree->child_right.size());
				printf("child left %d\n", (int) tree->child_left.size());
				printf("feature %d\n", (int) tree->feature.size());
		}
		for (int i = 0; i != 5; i++ ){
				if (i == 0 )printf("Threshold: ");
				printf("%f ", tree->threshold[i]);
		} printf("\n");

		for (int i = 0; i != 5; i++ ){
				if (i == 0 )printf("Value0: ");
				printf("%f ", tree->value0[i]);
		} printf("\n");

		for (int i = 0; i != 5; i++ ){
				if (i == 0 )printf("Value1: ");
				printf("%f ", tree->value1[i]);
		} printf("\n");

		for (int i = 0; i != 5; i++ ){
				if (i == 0 )printf("child right: ");
				printf("%d ", tree->child_right[i]);
		} printf("\n");

		for (int i = 0; i != 5; i++ ){
				if (i == 0 )printf("child left: ");
				printf("%d ", tree->child_left[i]);
		} printf("\n");

		for (int i = 0; i != 5; i++ ){
				if (i == 0 )printf("feature: ");
				printf("%d ", tree->feature[i]);
		} printf("\n");
		
}

void print_forest(Forest* forest){
		
		printf("Forest\nNumber of trees: %d == %d\n", forest->n_trees, (int)forest->trees.size());
		printf("Features: ");
		for(auto f = forest->features.begin(); f != forest->features.end(); f++){
				printf( " %s ", f->c_str());
		}
		printf( "\n");

		int i = 0;
		for(auto t = forest->trees.begin(); t != forest->trees.end(); t++){
				printf( "Tree number: %d\n", i);
				print_tree(&(*t));
				i++;
				if (i == 1) break;
		}
			
}




void load_tree( std::string filename, Tree* tree){
		std::vector<std::string> members_tree  = {"max_depth", "n_nodes", "child_left", "child_right", "threshold", "value0", "value1", "feature"};
		std::fstream file;
		file.open(filename.c_str(), std::fstream::in);

		if(!file.good()){
				printf("Tree reader could not open file: %s!\n", filename.c_str());
				return;
		}

		
		std::string line, text;
		while(!file.eof()) {
				std::getline(file, line);
				text.append(line+"\n");
				line.clear();
		}
		file.close();
		
		tree->max_depth = 0;
		tree->n_nodes = 0;
		tree->threshold.clear();
		tree->value0.clear();
		tree->value1.clear();
		tree->child_left.clear();
		tree->child_right.clear();
		tree->feature.clear();

		std::string temp_buffer;
		int cursor = -1;
		
		for(auto it = text.begin(); it != text.end(); it++){
				if (*it == ':' && temp_buffer.length() > 0) { 
						int k = 0;
						for(auto kt = members_tree.begin(); kt != members_tree.end(); kt++){
								if (strcmp(temp_buffer.c_str(), kt->c_str()) == 0){
										temp_buffer.clear();
										cursor = k;
								}
								k++; 
						}
						continue;
				}
				if ((*it == ' ' || *it == '\n') && cursor != -1 && temp_buffer.length() > 0) { 

						if(members_tree[cursor] == "max_depth"){
								tree->max_depth = (int)atoi(temp_buffer.c_str());
						}
						else if(members_tree[cursor] == "n_nodes"){
								tree->n_nodes = (int)atoi(temp_buffer.c_str());
						}
						else if(members_tree[cursor] == "threshold")	{
								tree->threshold.push_back((float)atof(temp_buffer.c_str()));
						}
						else if(members_tree[cursor] == "value0")	{
								tree->value0.push_back((float)atof(temp_buffer.c_str()));
						}
						else if(members_tree[cursor] == "value1")	{
								tree->value1.push_back((float)atof(temp_buffer.c_str()));
						}
						else if(members_tree[cursor] == "child_left")	{
								tree->child_left.push_back((int)atoi(temp_buffer.c_str()));
						}
						else if(members_tree[cursor] == "child_right")	{
								tree->child_right.push_back((int)atoi(temp_buffer.c_str()));
						}
						else if(members_tree[cursor] == "feature")	{
								tree->feature.push_back((int)atoi(temp_buffer.c_str()));
						}
						temp_buffer.clear();

						if (*it == '\n') { cursor = -1; }
						continue;
				}
				if(*it == ' ' || *it == '\n')	{ continue;}
				temp_buffer += *it;
		}
}

void load_forest( std::string path, Forest* forest){
		std::vector<std::string> members_forest= {"number_trees", "trees", "features"};

		std::fstream file;
		std::string filename = path + "/rf.txt";
		file.open(filename.c_str(), std::fstream::in);

		if(!file.good()){
				printf("Forest reader could not open file: %s!\n", filename.c_str());
				return;
		}

		
		std::string line, text;
		while(!file.eof()) {
				std::getline(file, line);
				text.append(line+"\n");
				line.clear();
		}
		file.close();

		printf("%s\n", text.c_str());
		
		forest->n_trees = 0;
		forest->trees.clear();
		forest->features.clear();

		std::string temp_buffer;
		int cursor = -1;
		for(auto it = text.begin(); it != text.end(); it++){
				if (*it == ':' && temp_buffer.length() > 0) { 


						int k = 0;
						for(auto kt = members_forest.begin(); kt != members_forest.end(); kt++){
								if (strcmp(temp_buffer.c_str(), kt->c_str()) == 0){
										temp_buffer.clear();
										cursor = k;
										break;
								}
								k++; 
						}
						continue;
				}
				if ((*it == ' ' || *it == '\n') && cursor != -1 && temp_buffer.length() > 0) { 

						if(members_forest[cursor] == "number_trees"){

								printf("len %s\n", temp_buffer.c_str());
								forest->n_trees = (int)atoi(temp_buffer.c_str());
						}
						else if(members_forest[cursor] == "features"){
								forest->features.push_back( temp_buffer.c_str() );
						}
						temp_buffer.clear();
						continue;
						if (*it == '\n') { cursor = -1; }
				}
				if(*it == ' ' || *it == '\n')	{ continue;}
				
				temp_buffer += *it;
		}

		for(int i = 0; i < forest->n_trees; i++){
				Tree tree;
				filename = path + "/tree" + std::to_string(i) + ".txt";
				load_tree(filename, &tree);
				forest->trees.push_back(tree);
		}
}




