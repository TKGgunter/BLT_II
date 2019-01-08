#include <string>
#include <vector>
#include <stdio.h>

struct CSV{
		std::vector<std::string> column_names;	
		std::vector<std::vector<float>> datas;
		int nfeatures;	
		int nentries;	
};
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
				csv_text.append(line);
		}
		file.close();

		bool first_line = true;
		std::vector<char> delimiters = {'\t', ' '};
		char newline_chars = '\n';

		int cursor = 0;
		csv->column_names.clear();
		csv->datas.clear();
		csv->nfeatures = 0;
		csv->nentries  = 0;
		std::string temp_buffer;
		for (auto it = csv_text.begin(); it != csv_text.end(); it++){

				//@Robustness this could be moved within the delimiter loop 
				//not sure if that makes things easier to read or not
				if *it == newline_chars{ 
						if first_line {
								//I'm scared that they is a copy of pointers not a true clone of that data
					  		csv->column_names.push_back(temp_buffer);
								csv->nfeatures++;
								for(int i = 0; i < csv->nfeatures; i++){
										csv->datas.push_back(std::vector());
								}
								first_line = false; 
						} else {
								csv->nentries++;
						}
						cursor = 0;
						temp_buffer.clear();
						continue;
				}
				for(auto _it = delimiters.begin(); _it != delimeter.end(); _it++){

						if( *_it == *it ) { 
								if first_line {
										//I'm scared that they is a copy of pointers not a true clone of that data
										csv->column_names.push_back(temp_buffer);
										csv->nfeatures++;

								} else {
										csv->datas[cursor].push_back((float)atof(temp_buffer));
								}
								cursor++;
								temp_buffer.clear();	
								break;	
						}
				}
				temp_buffer += *it;
		} 
	
};












