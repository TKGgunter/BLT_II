#ifndef CSVREADER_HH
#define CSVREADER_HH


#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>

struct CSV{
		std::vector<std::string> column_names;	
		std::vector<std::vector<float>> datas;
		int nfeatures;	
		int nentries;	
};

std::vector<float>* get_column(CSV* csv, std::string column){
		int i = 0;
		for( auto it = csv->column_names.begin(); it != csv->column_names.end(); it++){
				if ( *it == column ) return &(csv->datas[i]);
				i++;
		}
		printf("get column failed. %s .. we need a Option enum", column.c_str());
		exit(0);
		return &(csv->datas[0]);
}

void print_csv(CSV* csv){

		printf("CSV\n");
		printf("number of features %d %d\n", (int)csv->column_names.size(), csv->nfeatures);
		for(int i= 0; i< csv->nfeatures; i++){
				printf("~%s~%d\t", csv->column_names[i].c_str(), (int)csv->datas[i].size());
		}

		printf("\n");

		for(int i= 0; i < csv->nentries; i++){
				for(int j= 0; j < csv->nfeatures; j++){

						printf("%f\t", csv->datas[j][i]);
				}
				printf("\t %d\n", i);
		}
}

void read_csv( std::string filename, CSV* csv){
		std::fstream file;
		file.open(filename.c_str(), std::fstream::in);

		if(!file.good()){
				printf("CSV reader could not open file: %s!", filename.c_str());
				return;
		}

		
		std::string line, csv_text;
		while(!file.eof()) {
				std::getline(file, line);
				csv_text.append(line+"\n");
				line.clear();
		}
		file.close();

		//NOTE
		//This is a test to be sure that we have copied the entire file into memory
		//TODO: Remove after we test Jan 9, 2018
		//printf("this is a test\n%s", csv_text.c_str());

		bool first_line = true;
		std::vector<char> delimiters = {'\t', ' '};

		int cursor = 0;
		csv->column_names.clear();
		csv->datas.clear();
		csv->nfeatures = 0;
		csv->nentries  = 0;
		std::string temp_buffer;
		for (auto it = csv_text.begin(); it != csv_text.end(); it++){

				//@Robustness this could be moved within the delimiter loop 
				//not sure if that makes things easier to read or not
				if (*it == '\n'){ 
						if (first_line) {
								//I'm scared that they is a copy of pointers not a true clone of data
					  		csv->column_names.push_back(temp_buffer);
								csv->nfeatures++;
								for(int i = 0; i < csv->nfeatures; i++){
										csv->datas.push_back(std::vector<float>());
								}
								first_line = false; 
						} else{
								if (temp_buffer.length() > 0) csv->datas[cursor].push_back((float)atof(temp_buffer.c_str()));
						}
						cursor = 0;
						temp_buffer.clear();
						continue;
				}
				bool skip_char = false;
				for(auto _it = delimiters.begin(); _it != delimiters.end(); _it++){

						if( *_it == *it  && temp_buffer.length() > 0) { 
								if (first_line) {
										csv->column_names.push_back(temp_buffer);
										csv->nfeatures++;

								} else {
										csv->datas[cursor].push_back((float)atof(temp_buffer.c_str()));
										cursor++;
								}
								temp_buffer.clear();	
								break;
						} 
						if (*_it == *it && (int)temp_buffer.length() == 0){
							skip_char = true;
						}
				}
				//CLEANUP
				//May not need this
				if (skip_char){ 
						continue;
				}
				if ( ' ' == *it &&  0 == (int)temp_buffer.length()) {
						continue;
				}
				temp_buffer += *it;
				//printf("%s %d\n", temp_buffer.c_str(), (int)temp_buffer.length());
		} 
		csv->nentries = (int)csv->datas[0].size() ;
	
}





#endif






