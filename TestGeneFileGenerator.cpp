#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <set>
#include <utility>

void split(std::string line, std::vector<std::string>& ll, char delimiter=' ')
{
    std::istringstream ss(line);
    std::string elem;
    if (delimiter == ' ')
        while( ss >> elem )
            ll.emplace_back(elem);
    else
    {
        while ( std::getline(ss, elem, delimiter) )
            ll.emplace_back(elem);
    }
}

struct GeneInfo
{
    std::string chromosomeName;
    std::string geneName;
    int32_t start;
    int32_t end;
    char strand;
};

void generateGeneTxpFasta(std::string genomeFileName,
                            std::string txpFastaFileName, 
                            std::string outputGeneGff3FileName, 
                            std::string outputGeneTxpFastaFileName)
{
    /// @brief parse genome.fa positive strand
    std::ifstream ifs_genome(genomeFileName);
    std::map<std::string, std::string> genome;
    std::string current_chromosome, line;
    while (std::getline(ifs_genome, line))
    {
        if (line[0] == '>')
        {
            std::istringstream ss(line);
            ss >> current_chromosome;
            current_chromosome = current_chromosome.substr(1);
        }
        else
            genome[current_chromosome] += line;
    }
    ifs_genome.close();
    /// @brief reverse it to get negative strand of genome
    std::map<std::string, std::string> genome_negative_strand;
    for(auto iter=genome.begin(); iter!=genome.end(); ++iter)
    {
        std::cerr << "Generating complement strand of chromosome " << iter->first << std::endl;
        genome_negative_strand[iter->first] = iter->second;
        for(size_t i=0, sz = genome_negative_strand[iter->first].size(); i<sz; ++i)
        {
            char ch = genome_negative_strand[iter->first][i];
            if(ch == 'A' || ch == 'a')
                genome_negative_strand[iter->first][i] += ('T' - 'A');
            else if (ch == 'C' || ch == 'c')
                genome_negative_strand[iter->first][i] += ('G' - 'C');
            else if (ch == 'G' || ch == 'g')
                genome_negative_strand[iter->first][i] += ('C' - 'G');
            else if (ch == 'T' || ch == 't')
                genome_negative_strand[iter->first][i] += ('A' - 'T');
        }
    }


    std::vector<GeneInfo> gene_info_;
    std::ifstream ifs_geneMap(outputGeneGff3FileName);
    while (std::getline(ifs_geneMap, line))
    {
        if(line[0] == '#')
            continue;

        std::vector<std::string> ll;
        split(line, ll);
        if(ll[2] != "transcript")
            continue;
        if( ll[8].find("ID=", 0, 3) == std::string::npos )
            throw std::runtime_error("gff3 record should contain ID!!!");

        gene_info_.emplace_back(GeneInfo{});
        gene_info_.back().chromosomeName = ll[0];
        gene_info_.back().start = stoi(ll[3]);
        gene_info_.back().end = stoi(ll[4]);
        gene_info_.back().strand = ll[6][0];
        std::vector<std::string> ids;
        split(ll[8], ids, ';');
        for(auto elem : ids)
            if (elem.find("ID=", 0, 3) != std::string::npos)
            {
                gene_info_.back().geneName = elem.substr(3);
                break;
            }
    }
    ifs_geneMap.close();
    // for (auto elem : gene_info_)
    //     std::cerr << "chrname=" << elem.chromosomeName << " start=" << elem.start << " end=" << elem.end << " genename=" << elem.geneName << std::endl;

    std::ofstream ofs_fasta(outputGeneTxpFastaFileName, std::ios::out);
    std::ifstream ifs_txp(txpFastaFileName);
    while(std::getline(ifs_txp, line))
    {
        ofs_fasta << line << "\n";
    }
    for(const auto& elem : gene_info_)
    {
        ofs_fasta << ">" << elem.geneName;
        std::string s;
        if (elem.strand == '+')
            s = genome[elem.chromosomeName].substr(elem.start-1, elem.end-elem.start+1);
        else if (elem.strand == '-')
        {
            s = genome_negative_strand[elem.chromosomeName].substr(elem.start-1, elem.end-elem.start+1);
            std::reverse(s.begin(), s.end());
        }
        else
            std::cerr << ">>>>>>>>>>>>Error: strand should be either '+' or '-', not '" << elem.strand << "'" << std::endl;
        for(int32_t i(0); i<s.size(); ++i)
        {
            if(i%70==0)
                ofs_fasta << "\n";
            ofs_fasta << s[i];
        }
        ofs_fasta << "\n";
    }
    ofs_fasta.close();
    ifs_txp.close();
}


int main()
{
    std::string genomeFileName("/home/0309meeting/0413/txp_vs_genetxp/txp/Drosophila_melanogaster.BDGP6.dna.chromosome.all.fa");
    std::string txpFastaFileName("fake.fasta");
    std::string outputGeneGff3FileName("gene.gff3");
    std::string outputGeneTxpFastaFileName("jimmy.gene.fasta");
    generateGeneTxpFasta(genomeFileName, txpFastaFileName, outputGeneGff3FileName, outputGeneTxpFastaFileName);

    return 0;
}