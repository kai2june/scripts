// /// It seemed to contain segmentation fault now.
// #include <iostream>
// #include <string>
// #include <vector>
// #include <fstream>
// #include <sstream>
// #include <algorithm>
// #include <map>
// #include <stdexcept>
// #include <set>
// #include <utility>

// std::map<std::string, std::string> txps_genes_map_;

// void split(std::string line, std::vector<std::string>& ll, char delimiter=' ')
// {
//     std::istringstream ss(line);
//     std::string elem;
//     if (delimiter == ' ')
//         while( ss >> elem )
//             ll.emplace_back(elem);
//     else
//     {
//         while ( std::getline(ss, elem, delimiter) )
//             ll.emplace_back(elem);
//     }
// }

// void generateGeneGff3(std::string txpGff3FileName, std::string outputGeneGff3FileName)
// {
//     std::vector<std::vector<std::string>> gene_records;
//     std::ifstream ifs(txpGff3FileName);
//     std::string line, last_found, found;
//     while(std::getline(ifs, line))
//     {
//         if (line[0] == '#')
//             continue;
        
//         std::string txpName;
//         if ( line.find("FBgn") != std::string::npos )
//         {
//             std::vector<std::string> words;
//             std::vector<std::string> tags;
//             split(line, words);
//             split(words.back(), tags, ';');
//             for(auto iter=tags.begin(); iter!=tags.end(); ++iter)
//             {
//                 if ( iter->find("geneID=") != std::string::npos )
//                     found = iter->substr(iter->find("geneID=") + 7);
//                 else if ( iter->find("ID=", 0, 3) != std::string::npos )
//                     txpName = iter->substr(iter->find("ID=", 0, 3) + 3);
//             }
//             txps_genes_map_[txpName] = found;

//             words[2] = "exon";
//             words.back() = "Parent=" + found;

//             if (last_found == found)
//             {
//                 gene_records.back()[3] = 
//                     stoi(gene_records.back()[3]) < stoi(words[3]) ? gene_records.back()[3] : words[3];
//                 gene_records.back()[4] = 
//                     stoi(gene_records.back()[4]) > stoi(words[4]) ? gene_records.back()[4] : words[4];
//             }
//             else
//             {
//                 gene_records.emplace_back(words);
//             }
//             last_found = found;
//         }
//         else
//         {
//             found.clear();
//         }
//     }

//     /// @brief sort for the sake of unsorted gff3 input file 
//     std::sort(gene_records.begin(), gene_records.end(), 
//         [](const std::vector<std::string>& lhs, const std::vector<std::string>& rhs)
//         {
//             return lhs.back() < rhs.back();
//         }
//     );
//     std::vector<std::vector<std::string>> unique_gene_records;
//     last_found.clear();
//     found.clear();
//     for(auto iter=gene_records.begin(); iter!=gene_records.end(); ++iter)
//     {
//         found = iter->back();
//         if (last_found == found)
//         {
//             unique_gene_records.back()[3] = 
//                 stoi(unique_gene_records.back()[3]) < stoi((*iter)[3]) ? unique_gene_records.back()[3] : (*iter)[3];
//             unique_gene_records.back()[4] = 
//                 stoi(unique_gene_records.back()[4]) > stoi((*iter)[4]) ? unique_gene_records.back()[4] : (*iter)[4];
//         }
//         else
//         {
//             unique_gene_records.emplace_back(*iter);
//         }
//         last_found = found;
//     }
    
//     std::sort(unique_gene_records.begin(), unique_gene_records.end(), 
//         [](const std::vector<std::string>& lhs, const std::vector<std::string>& rhs)
//         {
//             if (lhs[0] == rhs[0])
//                 return stoi(lhs[3]) < stoi(rhs[3]);
//             return lhs[0] < rhs[0];
//         }
//     );

//     std::ofstream ofs(outputGeneGff3FileName);
//     for(auto iter=unique_gene_records.begin(); iter!=unique_gene_records.end(); ++iter)
//     {
//         std::vector<std::string> txp_line = (*iter);
//         txp_line[2] = "transcript";
//         std::string id = txp_line.back().substr(txp_line.back().find("Parent=") + 7);
//         txp_line.back() = "ID=" + id + ";geneID=" + id;

//         std::string s1, s2;
//         for(size_t i=0; i<txp_line.size(); ++i)
//         {
//             s1 = s1 + txp_line[i];
//             s2 = s2 + (*iter)[i];
//             if ( i != txp_line.size()-1 )
//             {
//                 s1 = s1 + "\t";
//                 s2 = s2 + "\t";
//             }
//             else
//             {
//                 s1 = s1 + "\n";
//                 s2 = s2 + "\n";
//             }
//         }
//         ofs << s1 << s2;
//     }
// }

// int main(int argc, char* argv[])
// {
//     generateGeneGff3(argv[1], argv[2]);
    
//     return 0;
// }